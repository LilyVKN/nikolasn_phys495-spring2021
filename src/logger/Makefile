FC = gfortran		# standard GNU FORTRAN compiler
MPIFC = mpifort		# MPI FORTRAN compiler

default: log_test

# log test ====================================================================
log_test: log_test.f90 sim_log.o sim_params.o
	$(MPIFC) -o log_test.out ../naive/part_init_default.o \
		../naive/force_grav_submat.o sim_params.o sim_log.o log_test.f90

sim_params.o: ../naive/sim_params.f90
	$(FC) -c ../naive/sim_params.f90

# test ========================================================================
test: sim_log.o test.f90
	$(FC) -o test.out sim_log.o test.f90

sim_log.o: sim_log.f90
	$(MPIFC) -c sim_log.f90

# clean =======================================================================
clean:
	$(RM) *.out *.o *~