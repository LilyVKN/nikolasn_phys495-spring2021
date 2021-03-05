# Naive Parallel N-Body Simulation

**Nikolas Nguyen** (nikolas@usc.edu)

## Building the Code
The software is currently hard-coded with the simulation parameters held in the
FORTRAN module file [`SIM_PARAMS.f90`](./SIM_PARAMS.f90). In order to operate
properly, the particle count, `PCOUNT`, **MUST** be divisible by the number of
worker nodes you intend to run the program with (i.e. running `mpirun` with
`-n 4` will create 3 worker nodes and one root). The other parameters are
described in the module file. The timesteps represent the number of frames at
the standard 24 FPS which are further broken into substeps for refining the 
simulation which given by `SUBSTEPS_PER_FRAME`. It is not recommended to change
the `DOMAIN` array since, currently, the video output is limited to a 4 x 4 x 4
grid.

A Makefile is included in the directory compatible with Linux (see 
[Dependencies](#dependencies)). To build, simply run:
    
    $ make

## Running the Code
Since this program utilizes MPI for parallelization, it should be run using
multiple process. Since the simulation parameters must be set such that the
particle count is divisible by the number of worker nodes, the program must be
run with the appropriate number of processes (i.e. # of worker nodes + the root)
The default setup uses 27 particles and 3 worker nodes. To run the default, use:

    $ mpirun -n 4 ./parallel

## Retrieving the Output
When you run the program, the output of the simulation is stored in a 
timestamped folder `./SPH_fortran_dd-mm-yy_HHMM` as `*.fcache` files. The format
of the file is one particle position per line ordered x, y, z delimited only by
whitespace. This can then be read by the utility program
[`ft2py.py`](../../util/ft2py.py). For detailed usage, use the built-in help
display, `-h`. For example:

    $ python3 ./ft2py.py -f ./folder/SPH_fortran_%04d.fcache -o ./output.mp4

will compile all the frames of `.fcache` files into the video file 
`./output.mp4` using the python `matplotlib` 3D libraries.

## Dependencies <a name=dependencies> </a>
This software requires an MPI library installed on Linux computer(s). This code
was built and tested using [MPICH 3.3.2](https://www.mpich.org). The Makefile
uses the GNU FORTRAN compiler, `gfortran` and the MPI-specific compiler 
`mpifort` provided by the MPICH distribution. 