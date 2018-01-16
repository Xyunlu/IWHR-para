#TOP LEVEL MAKEFILE FOR AZTEC LIBRARY

#====================================================================
# ------------------------
# | CVS File Information |
# ------------------------
#
# $RCSfile: Makefile_dist_lib,v $
#
# $Author: tuminaro $
#
# $Date: 2000/06/02 16:49:21 $
#
# $Revision: 1.28 $
#
# $Name:  $
#====================================================================
include ../dir.inc
include ../Make.com

#
# Machine dependent timing routine
#
TIME_GENERIC    = md_timer_generic.c
TIME_SUN4       = md_timer_sun.c
TIME_SOLARIS    = md_timer_sun.c
TIME_SGI        = md_timer_sun.c
TIME_SGIM4      = md_timer_sun.c
TIME_SGI10K     = md_timer_sun.c
TIME_DEC        = md_timer_sun.c
TIME_I860       = md_timer_intel.c
TIME_HP         = md_timer_generic.c
TIME_SUNMOS     = md_timer_intel.c
TIME_NCUBE      = md_timer_ncube.c
TIME_SP2        = md_timer_generic.c
TIME_T3E        = md_timer_mpi.c
TIME_LINUX      = md_timer_sun.c
TIME_TFLOP      = md_timer_intel.c

#
# Communication package
#
COMM_SERIAL      = md_wrap_scalar_c.c
COMM_I860        = md_wrap_intel_c.c
COMM_MPI         = md_wrap_mpi_c.c
COMM_NCUBE       = md_wrap_ncube_c.c
COMM_SUNMOS      = md_wrap_puma_c.c

#
# Blas routines NOT supplied by vendor
#
BLAS_GENERIC = blas1.f  blas2.f blas3.f xerbla.f lsame.f
BLAS_SUN4    = blas1.f  blas2.f blas3.f xerbla.f lsame.f
BLAS_SOLARIS = blas1.f  blas2.f blas3.f xerbla.f lsame.f
BLAS_SGI     = blas1.f  blas2.f blas3.f xerbla.f lsame.f
BLAS_SGIM4   = blas1.f  blas2.f blas3.f xerbla.f lsame.f
BLAS_SGI10K  = blas1.f  blas2.f blas3.f xerbla.f lsame.f
BLAS_DEC     = blas1.f  blas2.f blas3.f xerbla.f lsame.f
BLAS_I860    =
BLAS_HP      = blas1.f  blas2.f blas3.f xerbla.f lsame.f
BLAS_SUNMOS  =
BLAS_NCUBE   =          blas2.f blas3.f xerbla.f lsame.f
BLAS_SP2     = blas1.f  blas2.f blas3.f xerbla.f lsame.f
BLAS_T3E     =
BLAS_LINUX   = blas1.f  blas2.f blas3.f xerbla.f lsame.f
BLAS_TFLOP   = 

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#
# No need to change the rest of this file when adding
# a new machine to the makefile
#
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

CC              = $(CC_LINUX)
FC              = $(FC_LINUX)
AR              = $(AR_LINUX)
BLAS            = $(BLAS_LINUX)
TIMER           = $(TIME_LINUX)
COMMUNICATION   = $(COMM_MPI)
RNLIB           = $(RNLIB_LINUX)
CFORT           = $(CFORT_LINUX)
CFLAGS2         = $(CFLAGS_LINUX)
FFLAGS          = $(FFLAGS_LINUX)
CFLAGS          = $(CFLAGS2) $(CFORT) -DAZ_MPI $(MPI_INCLUDE_DIR)

# Aztec source code

AZ_SRC = az_bilu.c              az_cg.c		        az_cgs.c           \
         az_cgstab.c            az_check.c              az_comm.c          \
         az_converge.c          az_dd_overlap.c         az_dgemv2.c        \
         az_dgemv3.c            az_domain_decomp.c      az_fortran_wrap.c  \
         az_scaling.c           az_flop_cnt.c           az_gmres.c         \
         az_gmresr.c            az_ilu_util.c           az_ilut.c          \
         az_interface.c         az_lu_y12.c             az_matrix_util.c   \
         az_matvec_mult.c       az_old_matvec_mult.c    az_pad_utils.c     \
         az_poly.c              az_precond.c            az_qmrcgs.c        \
         az_rilu.c              az_solve.c              az_sort.c          \
         az_subdomain_solver.c  az_tools.c              az_util.c          \
         az_icc.c 		az_fix_pt.c


Y12M =   y12m.f			y12mae.f		y12maf.f \
	 y12mbe.f 		y12mbf.f 		y12mce.f \
	 y12mcf.f		y12mde.f		y12mdf.f \
	 y12mfe.f		y12mge.f		y12mhe.f \
	 y12cck.f


LA_SRC = la_dgetrf.f		la_dgetri.f		la_ilaenv.f \
	 la_dgetf2.f		la_dlaswp.f		la_dtrtri.f \
	 la_dtrti2.f		la_dgetrs.f		la_dgeco.f \
	 la_dgedi.f		la_dgefa.f		la_dpotf2.f \
	 la_dpotrf.f		la_dlamch.f		la_dlaic1.f \
	 la_ieeeck.f

SRC_F = $(LA_SRC) $(Y12M)  $(BLAS) az_reorder.f
SRC_C = $(AZ_SRC) $(TIMER) $(COMMUNICATION)

OBJ = 	$(SRC_C:.c=.o) $(SRC_F:.f=.o)

all: aztec


#####################################
# Define implicit compilation rules #
#####################################

#%.o:%.f
#	$(FC) $(FFLAGS) -c $<
#%.o:%.c
#	$(CC) $(CFLAGS) -c $<

# solver library

aztec: $(OBJ)
	@echo "Building library $(TARGET)"
	@rm -f libaztec.a
	$(AR) ruv libaztec.a $(OBJ)
	$(RNLIB) libaztec.a

clean:
	@echo "cleaning ..."
	/bin/rm -f *.o

realclean:
	@echo "real cleaning ..."
	/bin/rm -f *.o *.a 
