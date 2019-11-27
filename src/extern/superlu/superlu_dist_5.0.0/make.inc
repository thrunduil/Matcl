############################################################################
#
#  Program:         SuperLU_DIST
#
#  Module:          make.inc
#
#  Purpose:         Top-level Definitions
#
#  Creation date:   March 1, 2016	version 5.0.0
#
#  Modified:	    
#		    
#
############################################################################
#
#  The name of the libraries to be created/linked to
#
SuperLUroot	= /home/xiaoye/Dropbox/Codes/SuperLU/SuperLU_DIST_5.0.0/build
DSUPERLULIB   	= $(SuperLUroot)/SRC/libsuperlu_dist.so

# BLASDEF 	= -DUSE_VENDOR_BLAS

LIBS		= $(DSUPERLULIB) /usr/lib/libf77blas.so /usr/lib/libatlas.so /home/xiaoye/lib/dynamic/parmetis-4.0.3/build/Linux-x86_64/libparmetis/libparmetis.so /home/xiaoye/lib/dynamic/parmetis-4.0.3/build/Linux-x86_64/libmetis/libmetis.so

#
#  The archiver and the flag(s) to use when building archive (library)
#  If your system has no ranlib, set RANLIB = echo.
#
ARCH         = /usr/bin/ar
ARCHFLAGS    = cr
RANLIB       = /usr/bin/ranlib

CC           = /home/xiaoye/mpich-install/bin/mpicc
CFLAGS 	     = -O3 -DNDEBUG   -DUSE_VENDOR_BLAS -DDEBUGlevel=0 -DPRNTlevel=0 -std=c99 -g 
# CFLAGS       += -D
# CFLAGS       +=  
NOOPTS       = -O0
FORTRAN	     = /usr/bin/gfortran

LOADER       = $(CC)
LOADOPTS     = -Wl,-rpath,/home/xiaoye/Dropbox/Codes/SuperLU/SuperLU_DIST_5.0.0/lib 
