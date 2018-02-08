#PMETISROOT=./parmetis
PMETISROOT=./pmetis
INCLUDE = -I. -I$(PMETISROOT)/include -I$(PMETISROOT)/metis/GKlib -I$(PMETISROOT)/metis/include -I$(PMETISROOT)/programs/.
LIBS = libaztec.a $(PMETISROOT)/lib/libparmetis.a $(PMETISROOT)/lib/libmetis.a -lm

BINDIR = ../bin
CFLAGS = -DLINUX -std=c99 -DNDEBUG -DNDEBUG2 -DHAVE_EXECINFO_H -DHAVE_GETLINE -O3
FFLAGS = -extend-source -O2 -mcmodel=medium -shared-intel

FC = mpif77
CC = mpicc
.SUFFIXES: .o .f .c

OBJS = module.o pcs.o partmesh.o mylib.o partsub0.o partcoor.o solvsub.o qsort.o

all: pcs

pcs: $(OBJS)
	$(FC) $(FFLAGS) -o $@ $^ $(LIBS)
	cp $@ $(BINDIR)

.f.o:
	$(FC) $(FFLAGS) -c $<
.c.o:
	$(CC) $(CFLAGS) $(INCLUDE) -c $<

.PHONY: clean
clean:
	rm -rf *.o *.mod

realclean:
	rm -rf pcs *.o *.mod
