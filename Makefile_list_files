# variables
FC=gfortran
CFLAGS=-c -g -Og -Wall

# linking
a.out: poltrian.o
	$(FC) poltrian.o list_files.o
	clean

# compiling
poltrian.o: poltrian.f90
	$(FC) $(CFLAGS) poltrian.f90 list_files.f90

# cleanup
clean:
	rm *.o a.out

# run
run:
	make