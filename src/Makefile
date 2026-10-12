CC ?= gcc
CFLAGS ?= -O3 -Wall -Wextra -Wno-unused-parameter
LDLIBS = -lz

mapsim: mapsim.c
	$(CC) $(CFLAGS) -o $@ $< $(LDFLAGS) $(LDLIBS)

clean:
	rm -f mapsim
