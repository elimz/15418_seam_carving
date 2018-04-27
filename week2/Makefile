CC=gcc
CFLAGS=-g -Wall -DDEBUG=$(DEBUG)

CFILES = seam.c cycletimer.c
HFILES = seam.h cycletimer.h
LDFLAGS= -lm
OMP=-fopenmp -DOMP

all:seam seam-omp

seam: $(CFILES) $(HFILES)
	$(CC) $(CFLAGS) -o seam $(CFILES) $(LDFLAGS)
seam-omp: $(CFILES) $(HFILES)
	$(CC) $(CFLAGS) $(OMP) -o seam-omp $(CFILES) $(LDFLAGS)
clean:
	rm seam
	rm seam-omp
