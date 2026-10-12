/*
 * mapsim - simulate exhaustive paired-end reads from a reference and score
 *          per-base mappability from the resulting alignments.
 *
 * Subcommands:
 *   simulate  write interleaved FASTQ for every fragment start in a set of regions
 *   evaluate  read SAM (stdin), record per-base best MAPQ of correctly placed reads
 *   merge     combine chunk outputs into a best-MAPQ bedgraph + unmappable BED
 *
 * Read names encode their origin as  <contig>:<frag_start0>:<strand>
 * where strand is F (normal), OT or OB (bisulfite, directional library).
 *
 * Fragment at 0-based start f, length F, read length L:
 *   left  window = ref[f, f+L)        right window = ref[f+F-L, f+F)
 *   normal: R1 = left (fwd)             R2 = rc(right)
 *   OT    : R1 = C2T(left) (fwd)        R2 = rc(C2T(right))
 *   OB    : R1 = rc(G2A(right))         R2 = G2A(left) (fwd)
 * Pairs where both reads are entirely N are skipped (they cannot align).
 * Single-end (--layout single): F = L and only R1 is written, so each
 * position yields one read (normal) or an OT and an OB read (bisulfite).
 *
 * A read is "correctly placed" if it is primary, on the origin contig and
 * strand, and its unclipped start is within --tolerance bp of its origin.
 * Each correctly placed read credits the reference bases it aligns to with
 * its MAPQ; a base's score is the max over reads. A base is unmappable if
 * its score is below --min-mapq. Reference N bases are treated like any
 * other base: an N covered by a good read is mappable.
 *
 * In bisulfite mode OT and OB reads are tracked separately; the combined score
 * is the max of the two, and merge also reports bases only one strand covers.
 *
 * Chunk files (evaluate output) are gzip and hold one or more records:
 *   "MAPSIM2\t<contig>\t<start0>\t<n>\t<F|OT|OB>\n" followed by n bytes, where
 *   byte i is best MAPQ+1 at position start0+i from reads of that strand type
 *   (0 = no correctly placed read covered it).
 */
#include <ctype.h>
#include <errno.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>

#define VERSION "1.4.0"

static void die(const char *fmt, const char *a) {
    fprintf(stderr, "[mapsim] error: ");
    fprintf(stderr, fmt, a ? a : "");
    fprintf(stderr, "\n");
    exit(1);
}

static void *xmalloc(size_t n) {
    void *p = malloc(n ? n : 1);
    if (!p) die("out of memory%s", NULL);
    return p;
}

/* ------------------------------------------------------------------ */
/* FASTA index (.fai) and region loading                               */
/* ------------------------------------------------------------------ */

typedef struct {
    char *name;
    int64_t len, offset, linebases, linewidth;
} fai_rec_t;

typedef struct {
    fai_rec_t *r;
    int n;
} fai_t;

static fai_t fai_load(const char *ref) {
    char path[4096];
    snprintf(path, sizeof path, "%s.fai", ref);
    FILE *fp = fopen(path, "r");
    if (!fp) die("cannot open %s (reference must be uncompressed and samtools-faidx'd)", path);
    fai_t f = {NULL, 0};
    int cap = 0;
    char line[4096];
    while (fgets(line, sizeof line, fp)) {
        char name[2048];
        long long len, off, lb, lw;
        if (sscanf(line, "%2047s %lld %lld %lld %lld", name, &len, &off, &lb, &lw) != 5) continue;
        if (f.n == cap) {
            cap = cap ? cap * 2 : 256;
            f.r = realloc(f.r, cap * sizeof *f.r);
            if (!f.r) die("out of memory%s", NULL);
        }
        f.r[f.n].name = strdup(name);
        f.r[f.n].len = len;
        f.r[f.n].offset = off;
        f.r[f.n].linebases = lb;
        f.r[f.n].linewidth = lw;
        f.n++;
    }
    fclose(fp);
    return f;
}

static int fai_idx(const fai_t *f, const char *name) {
    for (int i = 0; i < f->n; i++)
        if (strcmp(f->r[i].name, name) == 0) return i;
    return -1;
}

/* Load ref[beg, end) uppercased; anything not ACGT becomes N. */
static char *load_seq(const char *fa_path, const fai_rec_t *r, int64_t beg, int64_t end) {
    if (beg < 0) beg = 0;
    if (end > r->len) end = r->len;
    int64_t n = end > beg ? end - beg : 0;
    char *s = xmalloc(n + 1);
    s[n] = 0;
    if (n == 0) return s;
    FILE *fp = fopen(fa_path, "rb");
    if (!fp) die("cannot open %s", fa_path);
    static char buf[1 << 20];
    setvbuf(fp, buf, _IOFBF, sizeof buf);
    int64_t off = r->offset + (beg / r->linebases) * r->linewidth + beg % r->linebases;
    if (fseeko(fp, off, SEEK_SET) != 0) die("seek failed in %s", fa_path);
    int64_t got = 0;
    int c;
    while (got < n && (c = getc(fp)) != EOF) {
        if (c == '\n' || c == '\r') continue;
        if (c == '>') break;
        c = toupper(c);
        s[got++] = (c == 'A' || c == 'C' || c == 'G' || c == 'T') ? (char)c : 'N';
    }
    fclose(fp);
    if (got != n) die("short read from %s (gzipped, or .fai out of date?)", fa_path);
    return s;
}

/* ------------------------------------------------------------------ */
/* Region lists: --regions FILE [--chunk ID] or --contig/--start/--end */
/* FILE lines: "contig start end" or "chunk_id contig start end"       */
/* ------------------------------------------------------------------ */

typedef struct {
    char *contig;
    int64_t start, end;
} region_t;

typedef struct {
    const char *ref, *regions, *chunk, *contig, *out, *stats, *mode;
    int64_t start, end;
    int L, F, tol, minq, se;
} opts_t;

static void parse_opts(int argc, char **argv, opts_t *o, const char *cmd) {
    memset(o, 0, sizeof *o);
    o->start = o->end = -1;
    o->L = 150; o->F = 500; o->tol = 5; o->minq = 20; o->mode = "normal";
    for (int i = 1; i < argc; i++) {
        const char *a = argv[i];
        if (i + 1 >= argc) die("missing value for %s", a);
        const char *v = argv[++i];
        if (!strcmp(a, "--ref")) o->ref = v;
        else if (!strcmp(a, "--regions")) o->regions = v;
        else if (!strcmp(a, "--chunk")) o->chunk = v;
        else if (!strcmp(a, "--contig")) o->contig = v;
        else if (!strcmp(a, "--start")) o->start = atoll(v);
        else if (!strcmp(a, "--end")) o->end = atoll(v);
        else if (!strcmp(a, "--read-len")) o->L = atoi(v);
        else if (!strcmp(a, "--frag-len")) o->F = atoi(v);
        else if (!strcmp(a, "--tolerance")) o->tol = atoi(v);
        else if (!strcmp(a, "--min-mapq")) o->minq = atoi(v);
        else if (!strcmp(a, "--mode")) o->mode = v;
        else if (!strcmp(a, "--layout")) {
            if (!strcmp(v, "single")) o->se = 1;
            else if (strcmp(v, "paired")) die("--layout must be paired or single, not %s", v);
        }
        else if (!strcmp(a, "-o") || !strcmp(a, "--out")) o->out = v;
        else if (!strcmp(a, "--stats")) o->stats = v;
        else { fprintf(stderr, "[mapsim %s] ", cmd); die("unknown option %s", a); }
    }
    if (o->se) o->F = o->L; /* single-end: the "fragment" is just the read */
    if (o->L <= 0 || o->F < o->L) die("need 0 < read-len <= frag-len%s", NULL);
    if (o->minq < 0 || o->minq > 254) die("min-mapq must be 0..254%s", NULL);
}

static int load_regions(const opts_t *o, region_t **out) {
    region_t *r = NULL;
    int n = 0, cap = 0;
    if (o->regions) {
        FILE *fp = fopen(o->regions, "r");
        if (!fp) die("cannot open %s", o->regions);
        char line[4096];
        while (fgets(line, sizeof line, fp)) {
            char a[2048], b[2048];
            long long s, e;
            const char *name;
            if (sscanf(line, "%2047s %2047s %lld %lld", a, b, &s, &e) == 4) {
                if (o->chunk && strcmp(a, o->chunk)) continue;
                name = b;
            } else if (sscanf(line, "%2047s %lld %lld", a, &s, &e) == 3) {
                name = a;
            } else continue;
            if (n == cap) {
                cap = cap ? 2 * cap : 64;
                r = realloc(r, cap * sizeof *r);
                if (!r) die("out of memory%s", NULL);
            }
            r[n].contig = strdup(name);
            r[n].start = s;
            r[n].end = e;
            n++;
        }
        fclose(fp);
    } else if (o->contig && o->start >= 0 && o->end >= 0) {
        r = xmalloc(sizeof *r);
        r[0].contig = strdup(o->contig);
        r[0].start = o->start;
        r[0].end = o->end;
        n = 1;
    }
    if (n == 0) die("no regions given (use --regions FILE [--chunk ID] or --contig/--start/--end)%s", NULL);
    *out = r;
    return n;
}

/* ------------------------------------------------------------------ */
/* simulate                                                            */
/* ------------------------------------------------------------------ */

static char comp_tab[256];

static void init_comp(void) {
    for (int i = 0; i < 256; i++) comp_tab[i] = 'N';
    comp_tab['A'] = 'T'; comp_tab['C'] = 'G'; comp_tab['G'] = 'C'; comp_tab['T'] = 'A';
}

static void revcomp(char *dst, const char *src, int n) {
    for (int i = 0; i < n; i++) dst[i] = comp_tab[(unsigned char)src[n - 1 - i]];
}

static void convert(char *s, int n, char from, char to) {
    for (int i = 0; i < n; i++)
        if (s[i] == from) s[i] = to;
}

static void write_fq(FILE *out, const char *name, const char *seq, const char *qual, int L) {
    fputc('@', out);
    fputs(name, out);
    fputc('\n', out);
    fwrite(seq, 1, L, out);
    fputs("\n+\n", out);
    fwrite(qual, 1, L, out);
    fputc('\n', out);
}

static int cmd_simulate(int argc, char **argv) {
    opts_t o;
    parse_opts(argc, argv, &o, "simulate");
    if (!o.ref)
        die("usage: mapsim simulate --ref ref.fa (--regions chunks.tsv [--chunk ID] | --contig C --start S --end E)\n"
            "         [--mode normal|bisulfite] [--layout paired|single] [--read-len 150] [--frag-len 500]\n"
            "         > reads.fq (interleaved if paired)%s", NULL);
    int bs = !strcmp(o.mode, "bisulfite");
    if (!bs && strcmp(o.mode, "normal")) die("unknown mode %s", o.mode);
    const int L = o.L, F = o.F, se = o.se;
    fai_t fai = fai_load(o.ref);
    region_t *reg;
    int nreg = load_regions(&o, &reg);
    init_comp();

    static char obuf[1 << 22];
    setvbuf(stdout, obuf, _IOFBF, sizeof obuf);
    char *qual = xmalloc(L + 1), *a = xmalloc(L + 1), *b = xmalloc(L + 1), *t = xmalloc(L + 1);
    memset(qual, 'I', L);
    char name[4096];

    for (int ri = 0; ri < nreg; ri++) {
        const region_t *g = &reg[ri];
        int ci = fai_idx(&fai, g->contig);
        if (ci < 0) die("contig %s not in reference", g->contig);
        const fai_rec_t *r = &fai.r[ci];
        int64_t start = g->start, end = g->end;
        if (end > r->len - F + 1) end = r->len - F + 1; /* fragment must fit */
        int64_t npairs = 0, nskip = 0;
        if (start < end) {
            int64_t rbeg = start, rend = end - 1 + F, n = rend - rbeg;
            char *seq = load_seq(o.ref, r, rbeg, rend);
            int64_t *acgt = xmalloc((n + 1) * sizeof *acgt);
            acgt[0] = 0;
            for (int64_t i = 0; i < n; i++) acgt[i + 1] = acgt[i] + (seq[i] != 'N');

            for (int64_t f = start; f < end; f++) {
                int64_t lo = f - rbeg, ro = f + F - L - rbeg;
                if (acgt[lo + L] == acgt[lo] && acgt[ro + L] == acgt[ro]) { nskip++; continue; }
                const char *left = seq + lo, *right = seq + ro;
                if (!bs) {
                    snprintf(name, sizeof name, "%s:%lld:F", g->contig, (long long)f);
                    revcomp(b, right, L);
                    write_fq(stdout, name, left, qual, L);
                    if (!se) write_fq(stdout, name, b, qual, L);
                    npairs++;
                } else {
                    snprintf(name, sizeof name, "%s:%lld:OT", g->contig, (long long)f);
                    memcpy(a, left, L); convert(a, L, 'C', 'T');
                    memcpy(t, right, L); convert(t, L, 'C', 'T'); revcomp(b, t, L);
                    write_fq(stdout, name, a, qual, L);
                    if (!se) write_fq(stdout, name, b, qual, L);

                    snprintf(name, sizeof name, "%s:%lld:OB", g->contig, (long long)f);
                    memcpy(t, right, L); convert(t, L, 'G', 'A'); revcomp(a, t, L);
                    memcpy(b, left, L); convert(b, L, 'G', 'A');
                    write_fq(stdout, name, a, qual, L);
                    if (!se) write_fq(stdout, name, b, qual, L);
                    npairs += 2;
                }
            }
            free(acgt);
            free(seq);
        }
        fprintf(stderr, "[mapsim simulate] %s:%lld-%lld mode=%s layout=%s %s=%lld skipped_allN=%lld\n", g->contig,
                (long long)g->start, (long long)g->end, o.mode, se ? "single" : "paired", se ? "reads" : "pairs", (long long)npairs, (long long)nskip);
    }
    if (fflush(stdout) != 0) die("write failed: %s", strerror(errno));
    return 0;
}

/* ------------------------------------------------------------------ */
/* evaluate                                                            */
/* ------------------------------------------------------------------ */

/* strand types: F (normal), OT, OB (bisulfite) */
#define NSTRAND 3
static const char *strand_label[NSTRAND] = {"F", "OT", "OB"};

static int strand_idx(const char *s) {
    for (int i = 0; i < NSTRAND; i++)
        if (!strcmp(s, strand_label[i])) return i;
    return -1;
}

typedef struct {
    const char *contig;
    size_t clen;
    int64_t fbeg, fend; /* fragment starts [fbeg, fend) */
    int64_t abeg, n;    /* credited positions [abeg, abeg+n) */
    uint8_t *best[NSTRAND];
} eval_reg_t;

/* parse leading integer, advance pointer */
static long parse_int(const char **p) {
    long v = 0;
    while (**p >= '0' && **p <= '9') v = v * 10 + (*(*p)++ - '0');
    return v;
}

static int cmd_evaluate(int argc, char **argv) {
    opts_t o;
    parse_opts(argc, argv, &o, "evaluate");
    if (!o.ref || !o.out)
        die("usage: mapsim evaluate --ref ref.fa (--regions chunks.tsv [--chunk ID] | --contig C --start S --end E)\n"
            "         -o chunk.gz [--mode normal|bisulfite] [--layout paired|single] [--stats FILE]\n"
            "         [--read-len 150] [--frag-len 500]\n"
            "         [--tolerance 5] [--min-mapq 20 (stats only)]  < aln.sam%s", NULL);
    const int L = o.L, F = o.F, tol = o.tol, minq = o.minq;
    int bs = !strcmp(o.mode, "bisulfite");
    if (!bs && strcmp(o.mode, "normal")) die("unknown mode %s", o.mode);
    const int s_lo = bs ? 1 : 0, s_hi = bs ? 2 : 0; /* strand indices used by this mode */
    fai_t fai = fai_load(o.ref);
    region_t *reg;
    int nreg = load_regions(&o, &reg);

    eval_reg_t *er = xmalloc(nreg * sizeof *er);
    for (int i = 0; i < nreg; i++) {
        int ci = fai_idx(&fai, reg[i].contig);
        if (ci < 0) die("contig %s not in reference", reg[i].contig);
        int64_t len = fai.r[ci].len;
        er[i].contig = reg[i].contig;
        er[i].clen = strlen(reg[i].contig);
        er[i].fbeg = reg[i].start;
        er[i].fend = reg[i].end;
        er[i].abeg = reg[i].start;
        int64_t aend = reg[i].end - 1 + F;
        if (aend > len) aend = len;
        if (reg[i].end > len) aend = len;
        er[i].n = aend > er[i].abeg ? aend - er[i].abeg : 0;
        memset(er[i].best, 0, sizeof er[i].best);
        for (int st = s_lo; st <= s_hi; st++) {
            er[i].best[st] = calloc(er[i].n ? er[i].n : 1, 1);
            if (!er[i].best[st]) die("out of memory%s", NULL);
        }
    }

    int64_t nprim = 0, nunmapped = 0, nwrong = 0, ncorrect = 0, ngood = 0;
    char *line = NULL;
    size_t cap = 0;
    ssize_t len;
    static char ibuf[1 << 22];
    setvbuf(stdin, ibuf, _IOFBF, sizeof ibuf);
    int last = 0;

    while ((len = getline(&line, &cap, stdin)) > 0) {
        if (line[0] == '@') continue;
        /* QNAME FLAG RNAME POS MAPQ CIGAR */
        char *fld[6];
        char *p = line;
        int k = 0;
        for (; k < 6; k++) {
            fld[k] = p;
            char *tab = strchr(p, '\t');
            if (!tab) break;
            *tab = 0;
            p = tab + 1;
        }
        if (k < 6) die("malformed SAM line%s", NULL);
        int flag = atoi(fld[1]);
        if (flag & 0x900) continue; /* secondary / supplementary */
        nprim++;
        if (flag & 0x4) { nunmapped++; continue; }

        /* QNAME = contig:f:strand (contig may contain ':') */
        char *qn = fld[0];
        char *c2 = strrchr(qn, ':');
        if (!c2) die("unexpected read name %s", qn);
        *c2 = 0;
        char *c1 = strrchr(qn, ':');
        if (!c1) die("unexpected read name %s", qn);
        *c1 = 0;
        const char *strand = c2 + 1;
        int sti = strand_idx(strand);
        if (sti < s_lo || sti > s_hi) die("read strand type %s does not match --mode", strand);
        int64_t f = atoll(c1 + 1);
        size_t qlen = (size_t)(c1 - qn);

        /* find this read's region (usually same as last) */
        eval_reg_t *g = &er[last];
        if (!(g->clen == qlen && f >= g->fbeg && f < g->fend && !memcmp(g->contig, qn, qlen))) {
            int found = -1;
            for (int i = 0; i < nreg; i++)
                if (er[i].clen == qlen && f >= er[i].fbeg && f < er[i].fend && !memcmp(er[i].contig, qn, qlen)) {
                    found = i;
                    break;
                }
            if (found < 0) die("read %s does not belong to this chunk", qn);
            last = found;
            g = &er[last];
        }

        int is_r1 = o.se || (flag & 0x40) != 0; /* single-end reads behave like R1 */
        int is_rev = (flag & 0x10) != 0;
        int64_t exp_start;
        int exp_rev;
        if (sti == 2) { /* OB */
            exp_start = is_r1 ? f + F - L : f;
            exp_rev = is_r1;
        } else { /* F or OT */
            exp_start = is_r1 ? f : f + F - L;
            exp_rev = !is_r1;
        }
        if (strcmp(fld[2], g->contig) != 0 || is_rev != exp_rev) { nwrong++; continue; }

        int64_t pos0 = atoll(fld[3]) - 1;
        const char *cg = fld[5];
        int64_t lead = 0;
        for (const char *q = cg; *q;) {
            long v = parse_int(&q);
            char op = *q++;
            if (op == 'S' || op == 'H') lead += v;
            else break;
        }
        int64_t d = pos0 - lead - exp_start;
        if (d < -tol || d > tol) { nwrong++; continue; }

        int mapq = atoi(fld[4]);
        if (mapq > 254) mapq = 254;
        ncorrect++;
        if (mapq >= minq) ngood++;
        uint8_t v = (uint8_t)(mapq + 1);

        uint8_t *best = g->best[sti];
        int64_t rp = pos0;
        for (const char *q = cg; *q && *q != '*';) {
            long ln = parse_int(&q);
            char op = *q++;
            if (op == 'M' || op == '=' || op == 'X' || op == 'D') {
                for (long j = 0; j < ln; j++, rp++) {
                    int64_t i = rp - g->abeg;
                    if (i >= 0 && i < g->n && best[i] < v) best[i] = v;
                }
            } else if (op == 'N') {
                rp += ln;
            }
        }
    }
    free(line);
    if (ferror(stdin)) die("error reading SAM from stdin%s", NULL);

    char tmp[4096];
    snprintf(tmp, sizeof tmp, "%s.tmp", o.out);
    gzFile gz = gzopen(tmp, "wb1");
    if (!gz) die("cannot write %s", tmp);
    for (int i = 0; i < nreg; i++)
        for (int st = s_lo; st <= s_hi; st++) {
            gzprintf(gz, "MAPSIM2\t%s\t%lld\t%lld\t%s\n", er[i].contig, (long long)er[i].abeg,
                     (long long)er[i].n, strand_label[st]);
            for (int64_t off = 0; off < er[i].n;) {
                unsigned w = (unsigned)((er[i].n - off) > (1 << 30) ? (1 << 30) : (er[i].n - off));
                if (gzwrite(gz, er[i].best[st] + off, w) != (int)w) die("gzwrite failed for %s", tmp);
                off += w;
            }
        }
    if (gzclose(gz) != Z_OK) die("gzclose failed for %s", tmp);
    if (rename(tmp, o.out) != 0) die("rename to %s failed", o.out);

    FILE *sf = o.stats ? fopen(o.stats, "w") : stderr;
    if (!sf) die("cannot write %s", o.stats);
    fprintf(sf, "chunk\tregions\tprimary_reads\tunmapped\tmisplaced\tcorrect\tcorrect_mapq_ge_%d\n", minq);
    fprintf(sf, "%s\t%d\t%lld\t%lld\t%lld\t%lld\t%lld\n", o.chunk ? o.chunk : "-", nreg, (long long)nprim,
            (long long)nunmapped, (long long)nwrong, (long long)ncorrect, (long long)ngood);
    if (o.stats) fclose(sf);
    if (nprim == 0) fprintf(stderr, "[mapsim evaluate] warning: no alignment records seen\n");
    return 0;
}

/* ------------------------------------------------------------------ */
/* merge                                                               */
/* ------------------------------------------------------------------ */

static int cmd_merge(int argc, char **argv) {
    opts_t o;
    parse_opts(argc, argv, &o, "merge");
    const char *list = o.regions, *prefix = o.out;
    if (!o.ref || !list || !prefix)
        die("usage: mapsim merge --ref ref.fa --regions chunk_files.txt --out PREFIX [--min-mapq 20]\n"
            "  chunk_files.txt: one evaluate output path per line\n"
            "  writes PREFIX.unmappable.bed, PREFIX.bestmapq.bedgraph.gz, PREFIX.N.bed and, for\n"
            "  bisulfite chunks, PREFIX.single_strand.bed%s", NULL);
    const int minq = o.minq;
    fai_t fai = fai_load(o.ref);

    /* best[strand][contig], allocated on first use */
    uint8_t **best[NSTRAND];
    for (int st = 0; st < NSTRAND; st++) best[st] = calloc(fai.n, sizeof(uint8_t *));
    int *seen = calloc(fai.n, sizeof *seen);
    int strand_seen[NSTRAND] = {0};

    FILE *lf = fopen(list, "r");
    if (!lf) die("cannot open %s", list);
    char path[4096];
    int nch = 0;
    uint8_t *buf = NULL;
    int64_t bufcap = 0;
    while (fgets(path, sizeof path, lf)) {
        path[strcspn(path, "\r\n")] = 0;
        if (!path[0]) continue;
        gzFile gz = gzopen(path, "rb");
        if (!gz) die("cannot open chunk %s", path);
        gzbuffer(gz, 1 << 20);
        char hdr[4096];
        int nrec = 0;
        while (gzgets(gz, hdr, sizeof hdr)) {
            char magic[16], name[2048], lab[16];
            long long s, n;
            if (sscanf(hdr, "%15s %2047s %lld %lld %15s", magic, name, &s, &n, lab) != 5 || strcmp(magic, "MAPSIM2"))
                die("bad record header in %s (chunk files from an older mapsim must be regenerated)", path);
            int st = strand_idx(lab);
            if (st < 0) die("bad strand label in %s", path);
            int ci = fai_idx(&fai, name);
            if (ci < 0) die("chunk contig not in reference (%s)", name);
            if (s < 0 || s + n > fai.r[ci].len) die("chunk record out of contig bounds in %s", path);
            if (n > bufcap) { bufcap = n; free(buf); buf = xmalloc(bufcap); }
            for (int64_t got = 0; got < n;) {
                unsigned w = (unsigned)((n - got) > (1 << 30) ? (1 << 30) : (n - got));
                int rr = gzread(gz, buf + got, w);
                if (rr <= 0) die("truncated chunk %s", path);
                got += rr;
            }
            if (!best[st][ci]) {
                best[st][ci] = calloc(fai.r[ci].len ? fai.r[ci].len : 1, 1);
                if (!best[st][ci]) die("out of memory%s", NULL);
            }
            strand_seen[st] = 1;
            uint8_t *dst = best[st][ci] + s;
            for (int64_t i = 0; i < n; i++)
                if (dst[i] < buf[i]) dst[i] = buf[i];
            seen[ci] = 1;
            nrec++;
        }
        int err;
        gzerror(gz, &err);
        if (err != Z_OK && err != Z_STREAM_END) die("gzip error reading %s", path);
        gzclose(gz);
        if (nrec == 0) die("empty chunk file %s", path);
        nch++;
    }
    fclose(lf);
    free(buf);
    const int bs = strand_seen[1] || strand_seen[2];
    if (bs && strand_seen[0]) die("chunk list mixes normal and bisulfite chunks%s", NULL);

    snprintf(path, sizeof path, "%s.unmappable.bed", prefix);
    FILE *bed = fopen(path, "w");
    snprintf(path, sizeof path, "%s.N.bed", prefix);
    FILE *nbed = fopen(path, "w");
    snprintf(path, sizeof path, "%s.bestmapq.bedgraph.gz", prefix);
    gzFile bg = gzopen(path, "wb6");
    FILE *sbed = NULL;
    if (bs) {
        snprintf(path, sizeof path, "%s.single_strand.bed", prefix);
        sbed = fopen(path, "w");
    }
    if (!bed || !bg || !nbed || (bs && !sbed)) die("cannot write outputs with prefix %s", prefix);

    const uint8_t thr = (uint8_t)(minq + 1);
    int64_t tot_len = 0, tot_unmap = 0, tot_n = 0, tot_n_unmap = 0, tot_nomap = 0, tot_ot = 0, tot_ob = 0;
    for (int ci = 0; ci < fai.n; ci++) {
        const fai_rec_t *r = &fai.r[ci];
        const int64_t clen = r->len;
        /* combined score = max over strand types */
        uint8_t *b = calloc(clen ? clen : 1, 1);
        if (!b) die("out of memory%s", NULL);
        for (int st = 0; st < NSTRAND; st++) {
            const uint8_t *x = best[st][ci];
            if (!x) continue;
            for (int64_t i = 0; i < clen; i++)
                if (b[i] < x[i]) b[i] = x[i];
        }
        if (!seen[ci]) fprintf(stderr, "[mapsim merge] warning: no chunk data for %s (all unmappable)\n", r->name);
        char *seq = load_seq(o.ref, r, 0, clen);
        int64_t i = 0;
        while (i < clen) { /* unmappable = best MAPQ < minq; N bases are not special */
            if (b[i] < thr) {
                int64_t j = i;
                while (j < clen && b[j] < thr) {
                    if (seq[j] == 'N') tot_n_unmap++;
                    j++;
                }
                fprintf(bed, "%s\t%lld\t%lld\n", r->name, (long long)i, (long long)j);
                tot_unmap += j - i;
                i = j;
            } else i++;
        }
        i = 0;
        while (i < clen) {
            if (seq[i] == 'N') {
                int64_t j = i;
                while (j < clen && seq[j] == 'N') j++;
                fprintf(nbed, "%s\t%lld\t%lld\n", r->name, (long long)i, (long long)j);
                tot_n += j - i;
                i = j;
            } else i++;
        }
        if (bs) { /* bases good on exactly one of OT/OB; label = the strand that maps */
            const uint8_t *ot = best[1][ci], *ob = best[2][ci];
            i = 0;
            while (i < clen) {
                int a = ot && ot[i] >= thr, c = ob && ob[i] >= thr;
                if (a == c) { i++; continue; }
                int64_t j = i;
                while (j < clen && (ot && ot[j] >= thr) == a && (ob && ob[j] >= thr) == c) j++;
                fprintf(sbed, "%s\t%lld\t%lld\t%s\n", r->name, (long long)i, (long long)j, a ? "OT_only" : "OB_only");
                if (a) tot_ot += j - i; else tot_ob += j - i;
                i = j;
            }
        }
        i = 0;
        while (i < clen) { /* best MAPQ, -1 = no correctly placed read */
            int64_t j = i;
            while (j < clen && b[j] == b[i]) j++;
            if (b[i] == 0) tot_nomap += j - i;
            gzprintf(bg, "%s\t%lld\t%lld\t%d\n", r->name, (long long)i, (long long)j, (int)b[i] - 1);
            i = j;
        }
        tot_len += clen;
        free(seq);
        free(b);
        for (int st = 0; st < NSTRAND; st++) free(best[st][ci]);
    }
    if (fclose(bed) || fclose(nbed) || gzclose(bg) != Z_OK || (sbed && fclose(sbed)))
        die("error closing outputs%s", NULL);
    fprintf(stderr,
            "[mapsim merge] chunk_files=%d contigs=%d total_bp=%lld\n"
            "[mapsim merge] unmappable_bp (best MAPQ < %d) = %lld (%.3f%%); no correctly placed read at all = %lld\n"
            "[mapsim merge] reference N bp = %lld, of which unmappable = %lld\n",
            nch, fai.n, (long long)tot_len, minq, (long long)tot_unmap,
            tot_len ? 100.0 * tot_unmap / tot_len : 0.0, (long long)tot_nomap, (long long)tot_n,
            (long long)tot_n_unmap);
    if (bs)
        fprintf(stderr, "[mapsim merge] single-strand bp: OT_only = %lld, OB_only = %lld\n", (long long)tot_ot,
                (long long)tot_ob);
    return 0;
}

int main(int argc, char **argv) {
    if (argc < 2 || !strcmp(argv[1], "-h") || !strcmp(argv[1], "--help")) {
        fprintf(stderr, "mapsim %s\nusage: mapsim <simulate|evaluate|merge> [options]\n", VERSION);
        return argc < 2;
    }
    if (!strcmp(argv[1], "simulate")) return cmd_simulate(argc - 1, argv + 1);
    if (!strcmp(argv[1], "evaluate")) return cmd_evaluate(argc - 1, argv + 1);
    if (!strcmp(argv[1], "merge")) return cmd_merge(argc - 1, argv + 1);
    if (!strcmp(argv[1], "--version") || !strcmp(argv[1], "version")) { printf("%s\n", VERSION); return 0; }
    fprintf(stderr, "unknown subcommand %s\n", argv[1]);
    return 1;
}
