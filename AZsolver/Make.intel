#====================================================================
COMM    = SERIAL
MACHINE = GENERIC
MPI_INCLUDE_DIR = -I${MpiDir}/intel64/include -I${SrcDir}/include
#MPI_INCLUDE_DIR = -I${SrcDir}/include
MPI_LIB         = -L${MpiDir} -lmpi
#
# ANSI C compiler
#
CC_GENERIC      = cc
CC_SUN4         = acc
CC_SOLARIS      = cc
CC_SGI          = cc
CC_SGIM4        = cc
CC_SGI10K	= cc
CC_DEC          = cc
CC_I860         = icc
CC_HP           = cc
CC_SUNMOS	= sicc
CC_NCUBE        = ncc
CC_SP2          = mpcc
CC_T3E          = cc
CC_LINUX        = ${MPICC}
CC_TFLOP        = pgcc -cougar

#
# Fortran compiler
#
FC_GENERIC      = f77
FC_SUN4         = f77
FC_SOLARIS      = f77
FC_SGI          = f77
FC_DEC          = f77
FC_SGIM4        = f77
FC_SGI10K  	= f77
FC_I860         = if77
FC_HP           = f77
FC_SUNMOS	= sif77
FC_NCUBE        = ncc
FC_SP2          = mpxlf
FC_T3E          = f90
FC_LINUX        = ${MPIFC}
FC_TFLOP        = pgf77 -cougar

#
# Archive program
#
AR_GENERIC      = ar
AR_SUN4 	= ar
AR_SOLARIS    = ar
AR_SGI  	= ar
AR_DEC		= ar
AR_SGIM4 	= ar
AR_SGI10K 	= ar
AR_I860 	= ar860
AR_HP           = ar
AR_SUNMOS       = ar860
AR_NCUBE	= nar
AR_SP2          = ar
AR_T3E          = ar
AR_LINUX        = ar
AR_TFLOP        = xar

#
# Ranlib program
#
RNLIB_GENERIC   = touch
RNLIB_SUN4      = ranlib
RNLIB_SOLARIS   = ranlib
RNLIB_SGI       = touch
RNLIB_SGIM4     = touch
RNLIB_SGI10K    = touch
RNLIB_DEC       = touch
RNLIB_I860      = touch
RNLIB_HP        = touch
RNLIB_SUNMOS    = touch
RNLIB_NCUBE     = touch
RNLIB_SP2       = touch
RNLIB_T3E       = ranlib
RNLIB_LINUX     = ranlib
RNLIB_TFLOP     = xranlib

#
# Machine dependent fortran/C interface
#
CFORT_GENERIC     = -Dappend_
CFORT_SUN4        = -Dappend_
CFORT_SOLARIS     = -Dappend_
CFORT_SGI         = -Dappend_
CFORT_SGIM4       = -Dappend_
CFORT_SGI10K      = -Dappend_
CFORT_DEC         = -Dappend_
CFORT_I860        = -Dappend_
CFORT_HP          = -Dmatched
CFORT_SUNMOS      = -Dappend_
CFORT_NCUBE       = -Dcaps
CFORT_SP2         = -Dmatched
CFORT_T3E         = -Dcaps
CFORT_LINUX       = -Dappend_
CFORT_TFLOP       = -Dappend_

#
# Compilation flags
#
CFLAGS_GENERIC = -O
CFLAGS_SUN4   = -O2 -vc -Xc
CFLAGS_SOLARIS= -O -vc
CFLAGS_SGI    = -O2 -n32
CFLAGS_SGIM4  = -O  -n32
CFLAGS_SGI10K = -O  -64 -r10000 # 64 bit mips processors {R10000}
CFLAGS_DEC    = -O
CFLAGS_I860   = -O4
CFLAGS_HP     = -O4
CFLAGS_SUNMOS = -O4
CFLAGS_NCUBE  = -O
CFLAGS_SP2    = -O2
CFLAGS_T3E    = -O3 -DT3E
CFLAGS_LINUX  = -O2 -mcmodel=medium -shared-intel
CFLAGS_TFLOP  = -O3

FFLAGS_GENERIC = -O
FFLAGS_SUN4   = -O2
FFLAGS_SOLARIS= -O2
FFLAGS_SGI    = ${CFLAGS_SGI}
FFLAGS_SGIM4  = ${CFLAGS_SGIM4}
FFLAGS_SGI10K = ${CFLAGS_SGI10K}
FFLAGS_DEC    = ${CFLAGS_DEC}
FFLAGS_I860   = ${CFLAGS_I860}
FFLAGS_HP     = -O
FFLAGS_SUNMOS = ${CFLAGS_SUNMOS}
FFLAGS_NCUBE  = ${CFLAGS_NCUBE}
FFLAGS_SP2    = ${CFLAGS_SP2}
FFLAGS_T3E    = -O3 -dp
FFLAGS_LINUX  = -132 -O2 -mcmodel=medium -shared-intel
FFLAGS_TFLOP  = -O3

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#
# No need to change the rest of this file when adding
# a new machine to the makefile
#
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
