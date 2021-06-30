TOPSRCDIR := ..
include $(TOPSRCDIR)/config.mk
PROJECTNAME := spharm-library

.PHONY: all
all: spharm

.PHONY: spharm
spharm:
	$(MAKE) -C src

# initialization targets ======================================================
.PHONY: init
init:
	mkdir -p $(OBJDIR) $(LIBDIR)

# cleaning targets ============================================================
.PHONY: clean cleanobj
clean: cleanobj

cleanobj: 
	$(MAKE) -C src cleanobj