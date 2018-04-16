CC=gcc
CFLAGS=-g -Wall -DDEBUG=$(DEBUG)

CFILES = seam.c
HFILES = seam.h

all:seam

seam: $(CFILES) $(HFILES)
	$(CC) $(CFLAGS) -o seam $(CFILES)
clean:
	rm seam
