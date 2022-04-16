# -*- Mode: Makefile; -*-

CCPLUS=g++
MPICC = mpic++
CFLAGS= -O3 -Wall -DRK2 -DAUSM -DTIMING
BINS=./bin/euler
SOURCES= main.cpp Field.cpp Mesh.cpp Cell.cpp Flux.cpp Solver.cpp timeAdvance.cpp boundary.cpp debug.cpp
all: $(BINS)

$(BINS): $(SOURCES)
	$(CCPLUS) $(CFLAGS) $^ -o $@

clean:
	 rm -f $(BINS) ./bin/result.dat