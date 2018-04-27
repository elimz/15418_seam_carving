CC=gcc
CFLAGS=-g -Wall -DDEBUG=$(DEBUG)

CFILES = seam.c cycletimer.c sort.c
HFILES = seam.h cycletimer.h sort.h
LDFLAGS= -lm

all:seam

seam: $(CFILES) $(HFILES)
	$(CC) $(CFLAGS) -o seam $(CFILES) $(LDFLAGS)
clean:
	rm seam
	rm output.ppm
	rm output_seam.ppm
