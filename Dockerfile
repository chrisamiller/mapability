# Exhaustive paired-end mappability (normal + bisulfite) with bwa-mem3.
#
# Build for the cluster (x86_64):
#   docker buildx build --platform linux/amd64 -t mapability:1.4.0 .
FROM debian:bookworm-slim AS build
RUN apt-get update && apt-get install -y --no-install-recommends gcc libc6-dev make zlib1g-dev \
    && rm -rf /var/lib/apt/lists/*
COPY src /src
RUN make -C /src CFLAGS="-O3 -Wall -Wextra -Wno-unused-parameter" LDFLAGS="-static" mapsim

FROM mambaorg/micromamba:2.3.0-debian12-slim
USER root
RUN micromamba install -y -n base -c conda-forge -c bioconda \
        bwa-mem3=0.14.0 samtools=1.24 htslib=1.24 bedtools=2.31.1 python=3.12 curl \
    && micromamba clean -a -y \
    && apt-get update && apt-get install -y --no-install-recommends procps \
    && rm -rf /var/lib/apt/lists/*
COPY --from=build /src/mapsim /usr/local/bin/mapsim
COPY scripts/*.sh scripts/*.py /usr/local/bin/
RUN chmod 755 /usr/local/bin/*.sh /usr/local/bin/*.py /usr/local/bin/mapsim \
    && for t in bwa-mem3 samtools bgzip tabix bedtools python3; do ln -s /opt/conda/bin/$t /usr/local/bin/$t; done
# Work without the micromamba entrypoint or the image PATH (e.g. LSF docker passes the host PATH);
# the scripts also prepend these directories themselves
ENV PATH=/opt/conda/bin:/usr/local/bin:$PATH
ENTRYPOINT []
CMD ["/bin/bash"]
