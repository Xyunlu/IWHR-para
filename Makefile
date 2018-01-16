PMETISROOT=/md3220i/home/supery/working/pmetis
INCLUDE = -I. -I$(PMETISROOT)/include -I$(PMETISROOT)/metis/GKlib -I$(PMETISROOT)/metis/include -I$(PMETISROOT)/programs/.
LIBS = libaztec.a $(PMETISROOT)/lib/libparmetis.a $(PMETISROOT)/lib/libmetis.a -lm
#LIBS = /opt/itplibs/libfepg.a /opt/itplibs/libblas.a /opt/itplibs/libmetis.a /opt/itplibs/libfepgsolv.a /opt/itplibs/superlu_linux.a /opt/itplibs/blas_linux.a /md3220i/home/xiaojie/super/pmetis/lib/libparmetis.a /md3220i/home/xiaojie/super/pmetis/lib/libmetis.a -lm

BINDIR = ../bin
CFLAGS = -DLINUX -std=c99 -DNDEBUG -DNDEBUG2 -DHAVE_EXECINFO_H -DHAVE_GETLINE -O3 -shared-intel
FFLAGS = -extend-source -O2 -mcmodel=medium

FC = mpif77
CC = mpicc
.SUFFIXES: .o .f .c

OBJS = module.o pcs.o partmesh.o mylib.o partsub0.o partcoor.o solvsub.o

all: pcs

pcs: $(OBJS)
	$(FC) -o $@ $^ $(LIBS)
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
