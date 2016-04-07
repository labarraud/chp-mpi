FC = gfortran
FCFLAGS=-g -fcheck=all -Warray-bounds

PROG=chp
SOURCES= precision.o  boundarycondition.o GradiantConj.o main.o

all:$(PROG)

$(PROG):$(SOURCES)
	$(FC) $(FCFLAGS) -o $@ $^

%.o:%.f90
	$(FC) $(FCFLAGS) -c $<

clean:
	rm -f *.o *.mod

cleanall:
	rm -f $(PROG) *.o *.mod
