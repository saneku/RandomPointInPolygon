# variables
FC=gfortran
CFLAGS=-c -g -Og -Wall

# linking
a.out: poltrian.o
	$(FC) poltrian.o polygon_triangulate_test.o

# compiling
poltrian.o: poltrian.f90
	$(FC) $(CFLAGS) poltrian.f90 polygon_triangulate_test.f90

# cleanup
clean:
	rm *.o a.out

# run
run:
	make
	./a.out