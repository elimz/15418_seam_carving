CC=gcc
CFLAGS=-g -Wall -DDEBUG=$(DEBUG)

CFILES = seam.c cycletimer.c
HFILES = seam.h cycletimer.h
LDFLAGS= -lm

all:seam

seam: $(CFILES) $(HFILES)
	$(CC) $(CFLAGS) -o seam $(CFILES) $(LDFLAGS)
clean:
	rm seam
	rm tower_out.ppm
