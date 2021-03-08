# Naive Parallel N-Body Simulation

**Nikolas Nguyen** (nikolas@usc.edu)

## Building the Code
The software is currently hard-coded with the simulation parameters held in the
FORTRAN module file [`sim_params.f90`](./sim_params.f90). In order to operate
properly, the particle count, `PCOUNT`, is optimized to be divisible by the 
number of worker nodes and the number of subdivistions you run the program with
(i.e. if `PCOUNT` is set to `19` and there are `3` workers with `2`
subdivisions,`6` subdivisions total, the particle count will be pared to `18`). 

The timesteps represent the number of frames at the standard 24 FPS which are
further broken into substeps for refining the simulation which is given by
`SUBSTEPS_PER_FRAME`. It is not recommended to change the `DOMAIN` array since,
currently, the video output is limited to a 4 x 4 x 4 grid.

For a more comprehensive overview of the simulation parameters and the details
of the calculation methods, see the [`DESIGNDOC`](./DESIGNDOC.md)

A Makefile is included in the directory compatible with Linux (see 
[Dependencies](#dependencies)). To build, simply run:
    
    $ make

## Running the Code
Since this program utilizes MPI for parallelization, it must be run using
multiple process. Keep in mind that the particle count will be rounded down to
be divisible by the number of nodes. To run the default, use:

    $ mpirun -n 4 ./parallel.out

The code will output the total current energy for each frame.

## Retrieving the Output
When you run the program, the output of the simulation is stored in a 
timestamped folder `./SPH_fortran_dd-mm-yy_HHMM` as `*.fcache` files. The format
of the file is one particle position per line ordered x, y, z delimited only by
whitespace. This can then be read by the utility program
[`ft2py.py`](../../util/ft2py.py). For detailed usage, use the built-in help
display, `$ python3 ./ft2py.py -h`. As an example usage:

    $ python3 ./ft2py.py -f ./folder/SPH_fortran_%04d.fcache -o ./output.mp4

will compile all the frames of the `.fcache` files into the video file 
`./output.mp4` using the python `matplotlib` 3D libraries.

For convenience *all* cache files including the simulation output can be deleted
using `$ make clear_cache`. Be sure to save any required data before using this
option.

## Dependencies <a name=dependencies> </a>
This software requires an MPI library installed on Linux computer(s). This code
was built and tested using [MPICH 3.3.2](https://www.mpich.org). The Makefile
uses the `gfortran` compiler and the MPI-specific compiler `mpifort` provided by
the MPICH distribution. 

## Known Issues
The currently implemented file-lock system of preventing concurrent access to
data occasionally fails especially for small data sets. The lock file fails to
be deleted in FORTRAN's call to `CLOSE(18,status='DELETE')`. This will sometimes
cause further fatal errors in the next substep for nodes trying to read the
cached submatrix data while it's being written. For larger data sets, the nodes
will rarely overlap in submatrix calculations and so fatal error is especially
less likely. This behavior has only been observed in the Visual Studio Code
terminal; however there is no reason to believe that it is exclusive to that
environment.