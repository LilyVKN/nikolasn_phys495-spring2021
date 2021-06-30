# Spherical Harmonics Subprogram Library
**Nikolas Nguyen** (nikolas@usc.edu)

This collection of procedures allows for the calculation of spherical harmonics
and multipole moments with the intention to be used for the fast multipole
method and other N-body algorithms.

For a more comprehensive overview of the procedures and the details of the
calculation methods, see the [`DESIGNDOC`](./DESIGNDOC.md)

## Building the Code
Both CMake and GNU Make compiling is supported. All compiling commands should be
accessed from the top directory of the project (i.e. the folder containing
`./spharm/` and `./src/`) as the subdirectory CMake and Make configurations
depend on shared settings. 

To run CMake create a build directory and run `cmake` on the top source
directory.

```
mkdir build
cd build
cmake ..
```

For standalone compiling with GNU Make, run `make spharmlib` from the top
directory. The Makefile will automatically create a new directory in the top
level (./build) for compiling the library and output the final library to
`./build/lib/libspharm.a`. The Makefile settings are accessed in the file 
[`config.mk`](../config.mk) in the top directory. Details on each option
including the path to the build directory are described in that file.

```
make spharmlib
```

Alternatively, the library can be compiled with the rest of the library using
the Makefile `all` target. This will also compile the sample programs.

```
make all
```

Since these are mainly utility subroutines, you can compile a la carte as well.
This can be done using `$ gfortran -c file.f90`. The dependencies between the
subroutines are as follows:

> `dyml_full.f90: dpmlcos_full.f90`  
> `dyml.f90: dpmlcos.f90`  
> `multipole.f90: yml_full.f90 pmlcos_full.f90`  
> `yml_full.f90: pmlcos_full.f90`
> `yml.f90: pmlcos.f90`

The [`DESIGNDOC`](./DESIGNDOC.md) contains detailed descriptions on each
procedure that should make such dependencies clear. Additionally, the source
files' documentation list each of their dependencies.

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
[3] Pfalzner, S., & Gibbon, P. (2005). Many-body tree methods in physics. 
        Cambridge: Cambridge University Press.