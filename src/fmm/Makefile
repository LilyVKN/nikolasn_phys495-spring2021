TOPSRCDIR := ../..
include $(TOPSRCDIR)/config.mk

.SUFFIXES: .f90 .o
$(OBJDIR)/%.o: %.f90
	$(FC) $(FFLAGS) -c -o $@ $<

FMMSRC := $(wildcard *.f90)
FMMOBJ := $(patsubst %.f90, $(OBJDIR)/%.o, $(FMMSRC))

default: all

.PHONY: all
all: $(FMMOBJ)