TOPSRCDIR := ../..
include $(TOPSRCDIR)/config.mk

.SUFFIXES: .f .f90 .o
$(OBJDIR)/%.o: %.f90
	$(FC) $(FFLAGS) -c -o $@ $^

.PHONY: all
all: $(SPHARMLIB)

# retrieve all object files ===================================================
ALLSRC = $(wildcard *.f90)
ALLOBJ = $(patsubst %.f90, $(OBJDIR)/%.o, $(ALLSRC))

$(SPHARMLIB): $(ALLOBJ)
	$(AR) $(ARFLAGS) $@ $^
	$(RANLIB) $@

# clean =======================================================================
.PHONY: clean cleanobj
clean: cleanobj

cleanobj:
	$(RM) -f *.o 