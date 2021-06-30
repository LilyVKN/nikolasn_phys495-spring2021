TOPSRCDIR := .
PROJECTNAME := fmm-library
include $(TOPSRCDIR)/config.mk

default: fmmlib

.PHONY: all
all: init spharmlib fmmlib

.PHONY: fmmlib spharmlib
fmmlib: init
	$(MAKE) -C src

spharmlib: init
	$(MAKE) -C spharm

# initialization targets ======================================================
.PHONY: init
init:
	mkdir -p $(BUILDDIR)
	$(file > $(BUILDDIR)/.gitignore,*)

# cleaning targets ============================================================
.PHONY: clean
clean:
	$(MAKE) -C src clean