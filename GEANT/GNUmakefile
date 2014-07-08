# --------------------------------------------------------------
# Makefile for building G4 Root aplication
# --------------------------------------------------------------

name := fcalor
G4TARGET := $(name)
G4EXLIB := true

ifndef G4INSTALL
  G4INSTALL = ../../..
endif

.PHONY: all
all: lib bin

export ROOTSYS=/sharesoft/cmssw/slc5_amd64_gcc434/lcg/root/5.27.06b-cms25

include $(G4INSTALL)/config/binmake.gmk

ifdef HEPMC_DIR
  INCFLAGS  += -I$(HEPMC_DIR)/include
  LDLIBS    += -L$(HEPMC_DIR)/lib -lHepMC
endif

ifdef ROOTSYS
  CPPFLAGS += -I$(shell $(ROOTSYS)/bin/root-config --incdir)
  LDLIBS   += $(shell $(ROOTSYS)/bin/root-config --libs)
endif

CPPFLAGS += -g

ifeq ($(G4SYSTEM),Linux-g++)
  LDLIBS   += -Wl,-rpath,$(shell $(ROOTSYS)/bin/root-config --libdir)
endif

visclean:
	rm -f g4*.prim g4*.eps g4*.wrl
	rm -f .DAWN_*

