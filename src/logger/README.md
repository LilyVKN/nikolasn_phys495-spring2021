# Simulation Logging System

**Nikolas Nguyen** (nikolas@usc.edu)

## Structure of the Code
All the data types and functions are contained in [`sim_log.f90`](./sim_log.f90)
in a single FORTRAN `MODULE`. There are plans to compartmentalize them into a
more logical organization in the near future. Included is a logger testing file
which is dependent on the code in the `feature-naive-parallel`. For convenience,
a Makefile is included. Simply run `$ make` and run with the same MPI options.
(See [Naive Parallel N-Body Simulation](../naive/README.md) for more details)

## Basic Design Overview
Currently, the `simulation_logger` acts as a convenient wrapper object for 
organizing a formatted log format. Log lines are encapsulated by the 
`simlog_element` abstract data type for which derived types must implement a 
`print()` and `update()` function. The abstract data type is further extended
into organizing data structures such as fixed-sized collections which are not
intended to be directly handled by the user. The whole logging system likely has
more overhead than may be desired for heavy simulations.

The data types labeled by `simlog_entry*` are all meant to track values by
pointer. For example, if there a user wishes to track the energy, the variable
which stores the energy should be provided to an instance of `simlog_entry_real`
or `simlog_entry_scinum`. Because the entries use pointers to passively
track the values, the variables meant to be tracked by the entries should be
declared with the `TARGET` keyword.

Additional data types are included for formatting including the `simlog_report`
which creates a timestamped list of entries. The `simulation_logger` has by
default a start report run before the simulation, per-frame reports which can
provide vital data such as energy values, intermediate reports which print every
certain number of frames, and a completion report to be run at the end. A timer
which can be lapped and paused can be used to track either a continuous run or
a specific part of the calcuations.

A per-node data tracker has been added which can appropriately synchronize data
and update the root-side entries; however the organizing structure is not
currently updating and printing the correct entries. Fixes are incoming.

Full design documentation is still incoming as all the functionalities are
appropriately implemented.

## Using the Logger
To use an instance of the `simulation_logger` data type, first call the `init()`
function with rank of the node, the *total* number of processes, and the total
frame count. (The rank tracking was intended for organizing per-node data;
however this may be removed in future version. Right now it also prevents
printing from non-root nodes). 

All the desired entries should then be initialized, organized, and added to the
appropriate parts of the logger (e.g. start, mid, or end reports, frame
entries). Only the per-node data should be initialized in all nodes, everything
else can left to the root.

Right before starting the simulation loop, call `start()` which runs the start
report and begins the global simulation timer. Each frame `lap_frame()` will
tell the logger to increment the frame counter and to track the last frame's
interval. Additionally if the memeber `log_frames` is flagged as true (default)
a frame report will also be printed.

At the end of the simulation loop, call `finalize()` which will stop the timer
and indicate the end of the simulation in the log. All the necessary tracked
values can then be updated before calling `run_end_report()` which updates all
entries and prints the completion report.

## Dependencies <a name=dependencies> </a>
This software requires an MPI library installed on Linux computer(s). This code
was built and tested using [MPICH 3.3.2](https://www.mpich.org). The Makefile
uses the `gfortran` compiler and the MPI-specific compiler `mpifort` provided by
the MPICH distribution. 

## Known Issues
The organization of the data types and functions need heavy reworking.

The `simlog_node_tracker` data type does not appropriately call updates. So
while the data and entries can be synchronized, they aren't called properly by
the root node.