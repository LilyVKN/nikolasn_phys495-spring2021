FC = gfortran		# standard GNU FORTRAN compiler

default: test

# parallel - main program =====================================================
test: sph_density.o sph_du_dt.o sph_dv_dt.o sph_eos.o sph_kernel.o sph_rho.o \
	test_sph.f90
	$(FC) -o test_sph.out sph_kernel.o sph_density.o sph_du_dt.o sph_dv_dt.o \
		sph_eos.o sph_rho.o ../naive/part_init_default.o \
		../naive/force_grav_submat.f90 ../naive/sim_params.o test_sph.f90

sph_density.o: sph_density.f90 sph_kernel.o
	$(FC) -c sph_density.f90

sph_du_dt.o: sph_du_dt.f90
	$(FC) -c sph_du_dt.f90

sph_dv_dt.o: sph_dv_dt.f90
	$(FC) -c sph_dv_dt.f90

sph_eos.o: sph_eos.f90
	$(FC) -c sph_eos.f90

sph_kernel.o: sph_kernel.f90
	$(FC) -c sph_kernel.f90

sph_rho.o: sph_rho.f90
	$(FC) -c sph_rho.f90

# clean =======================================================================
clean:
	$(RM) *.mod *.out *.o *~

clear_cache:
	$(RM) -r ./SPH_fortran_* ./tmp