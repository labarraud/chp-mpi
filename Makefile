FC = mpif90
FCFLAGS=-Ofast

SRCDIR=Src
LIBDIR=Object
MODDIR=Module

PROG=chppar

SOURCES  := $(wildcard $(SRCDIR)/*.f90)
OBJECTS  := $(SOURCES:$(SRCDIR)/%.f90=$(LIBDIR)/%.o)

all:$(PROG)

$(PROG):$(OBJECTS)
	$(FC) $(FCFLAGS) $(OBJECTS) -o $(PROG) 

$(OBJECTS): $(LIBDIR)/%.o : $(SRCDIR)/%.f90
	$(FC) $(FCFLAGS) -c $< -o $@ -J$(MODDIR)
clean:
	rm -f $(LIBDIR)/*.o $(MODDIR)/*.mod

cleanall:
	rm -rf $(PROG) $(LIBDIR)/*.o $(MODDIR)/*.mod Result/*.dat Stat/*.dat *.gnu *.dat *.eps Doc/*

doc:
	doxygen DoxygenConfig

$(LIBDIR)/main.o : $(LIBDIR)/mathconst.o  $(LIBDIR)/gradiantconj.o $(LIBDIR)/boundarycondition.o $(LIBDIR)/modcharge.o $(LIBDIR)/write_read.o
