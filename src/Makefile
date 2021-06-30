TOPSRCDIR := ..
include $(TOPSRCDIR)/config.mk

default: all

.PHONY: all
all: $(FMMLIB)

.PHONY: fmm naive
fmm:
	$(MAKE) -C fmm

naive:
	$(MAKE) -C naive

$(FMMLIB): fmm naive
	$(AR) $(ARFLAGS) $@ $(OBJDIR)/*.o
	$(RANLIB) $@

# initialization targets ======================================================
.PHONY: init
init:
	mkdir -p $(OBJDIR) $(LIBDIR)

# cleaning targets ============================================================
.PHONY: clean cleanobj
clean:
	$(MAKE) -C fmm clean
	$(MAKE) -C naive clean

cleanobj:
	$(RM) *.o *.mod