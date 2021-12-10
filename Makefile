compiler = mpifort
compiler_flags = -O3 -fbackslash

libraries = -L/usr/lib -llapack -lblas -lsymspg
modules = spglib_f08.f90 cpu.f90 timer.f90 import.f90 route.f90 potential.f90 transform.f90 bulk.f90 matrix.f90 spectral.f90 density.f90 export.f90
target = main.f90
target_test = test.f90
output = a.out

.PHONEY : run debug clean test

default :
	$(compiler) $(compiler_flags) $(modules) $(target) -o $(output) $(libraries)
	
debug : run
	mpirun -np 4 ./a.out
	make clean
	
clean :
	rm *.o *.mod
