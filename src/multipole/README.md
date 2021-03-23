# Multipole Procedure Library
**Nikolas Nguyen** (nikolas@usc.edu)

This collection of procedures allows for the calculation of spherical harmonics
and multipole moments with the intention to be used for the fast multipole
method and other N-body algorithms.

## Building the Code
Since these are mainly utility subroutines, only compile whatever is needed.
This can be done using `$ gfortran -c file.f90`. The dependencies between the
subroutines are as follows:

> `dyml_full.f90: dpmlcos_full.f90`  
> `dyml.f90: dpmlcos.f90`  
> `multipole.f90: yml_full.f90 pmlcos_full.f90`  
> `yml_full.f90: pmlcos_full.f90`
> `yml.f90: pmlcos.f90`

## Running the Code
Each file contains a separate subroutine like many common FORTRAN libraries, so
simply linking the compiled files properly allows for calling the subroutines in
further code. Descriptions of arguments and behavior are provided at the heading
comment of each file.

## Known Issues
Because all the code will first be tested on 32-bit computers, `REAL(4)` values
are used to do all the calculations. Therefore, large errors accumulate for
high l-values of the multipole moments, ~5% in some cases. There are plans to
increase the floating-point precision in future.

## References
[1] Wiggins, R., & Saito, M. (1971). Evaluation of Computational
        Algorithms for the Associated Legendre Polynomiasl By Interval Analysis.
        Bulletin of the Seismological Society of America, 61(2), pp.375-381.  
[2] Bosch, W. (2000) On the Computation of Derivatives of Legendre
        Functions. Physics and Chemistry of the Earth (A), 25(9-11), pp.655-659.  
[3] Pfalzner, S., & Gibbon, P. (2005). The Fast Multipole Method. In
        Many-body tree methods in physics (pp. 126-148). Cambridge: Cambridge
        University Press.