FC = gfortran		# standard GNU FORTRAN compiler
MPIFC = mpifort		# MPI FORTRAN compiler

default: parallel

# parallel - main program =====================================================
parallel: part_init_default.o force_grav_submat.o
	$(MPIFC) -o parallel.out part_init_default.o force_grav_submat.o parallel.f90

part_init_default.o: part_init_default.f90
	$(FC) -c part_init_default.f90

force_grav_submat.o: force_grav_submat.f90
	$(FC) -c force_grav_submat.f90

parallel_clean:
	$(RM) parallel.out

# clean =======================================================================
clean:  parallel_clean
	$(RM) *.out *.o *~