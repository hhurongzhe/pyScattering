SHELL := /bin/sh

PYTHON ?= python3
F2PY := $(PYTHON) -m numpy.f2py

MODULE := force_module
WRAPPER := wrapforce.f90
INTERFACE := wrapforce.pyf

FORTRAN_SOURCES := \
	idahopot.f \
	av18.f \
	cdbonn.f \
	lo450.f lo500.f lo550.f \
	nlo450.f nlo500.f nlo550.f \
	n2lo450.f n2lo500.f n2lo550.f \
	n3lo450new.f n3lo500new.f n3lo550new.f \
	n4lo450.f n4lo500.f n4lo550.f \
	nnlo_opt.f nnlo_sat.f \
	n3lo.f ichiral.f90 smschiral_RS.f90

.PHONY: all clean

all:
	$(F2PY) $(WRAPPER) -h $(INTERFACE) -m$(MODULE) --overwrite-signature
	$(F2PY) -c -m$(MODULE) $(INTERFACE) $(WRAPPER) $(FORTRAN_SOURCES)

clean:
	rm -f *.o *.mod $(INTERFACE) $(MODULE)*.so $(MODULE)*.dylib $(MODULE)*.pyd
