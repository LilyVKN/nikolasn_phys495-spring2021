FC 		:= gfortran
MPIFC 	:= mpifortran

SHELL 	:= /bin/sh
FFLAGS 	:= -O2

AR 		:= ar
ARFLAGS := cr
RANLIB 	:= ranlib

BUILDDIR 	:= $(TOPSRCDIR)/build
OBJDIR 		:= $(BUILDDIR)/obj
LIBDIR		:= $(BUILDDIR)/lib

FMMLIB 		:= $(LIBDIR)/libfmm.a
SPHARMLIB	:= $(LIBDIR)/libspharm.a