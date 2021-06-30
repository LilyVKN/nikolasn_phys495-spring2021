TOPSRCDIR := ../..
include $(TOPSRCDIR)/config.mk

.SUFFIXES: .f90 .o
$(OBJDIR)/%.o: %.f90
	$(FC) $(FFLAGS) -c -o $@ $<

DIRECTSRC := $(wildcard *.f90)
DIRECTOBJ := $(patsubst %.f90, $(OBJDIR)/%.o, $(DIRECTSRC))

default: all

.PHONY: all
all: $(DIRECTOBJ)