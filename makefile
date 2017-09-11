FC = $(shell which gfortran)
CPP = $(shell which g++)

FF = -O4 -ffixed-line-length-130
CF = -O4

FOBJ = main.o read_file.o ejecta.o orbital_parameter.o classification.o omega.o

F90OBJ = time_elapsed.o

COBJ = roche_potential.o

LIB = -L/usr/lib -lstdc++ #-lc

%.o:%.f common_sph_var.h
	$(FC) -c $(FF) $<
%.o:%.cpp 
	$(CPP) -c $(CF) $<
%.o:%.f90
	$(FC) $(LIB) -c $(CF) $<

EXEC = $(shell basename $(shell pwd))

OBJ = $(FOBJ) $(COBJ) $(F90OBJ)

$(EXEC):$(OBJ)
	$(FC) -o $@ $(OBJ) $(LIB)

clean:
	$(shell which rm) *.o $(EXEC)
