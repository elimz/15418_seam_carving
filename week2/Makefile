CC=gcc
DEBUG=0
CFLAGS=-g -O3 -Wall -DDEBUG=$(DEBUG)

CFILES = seam.c cycletimer.c sort.c
HFILES = seam.h cycletimer.h sort.h
LDFLAGS= -lm
OMP=-fopenmp -DOMP

# all: seam seam-omp
all: seam seam-omp

seam: $(CFILES) $(HFILES)
	$(CC) $(CFLAGS) -o seam $(CFILES) $(LDFLAGS)

seam-omp: $(CFILES) $(HFILES)
	$(CC) $(CFLAGS) $(OMP) -o seam-omp $(CFILES) $(LDFLAGS)

clean:
	rm seam
	rm seam-omp
