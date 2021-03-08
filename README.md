# Parallel N-Body Simulation
**PHYS 495 - Spring 2021** : Nikolas Nguyen (nikolasn@usc.edu)

A software library for my PHYS 495 thesis project, Spring 2021, intended to
provide tools for astrophysical N-body simulations and visualization with
parallel capabilities 

## Contents
**Source**
1. [Naive Parallel N-Body Simulation](./src/naive/README.md) (4 March 2021)
    - MPI-parallel implementation of naive N-body gravitational SPH in FORTRAN

**Utilities**
1. [FORTRAN Cache to Python](./util/ft2py.py)  
    - Program for taking simulation output stored in a frame-ordered series of
    `.fcache` files to `.mp4` using `matplotlib` 3D plotting libraries

**Samples**
1. [Python Quad-Tree](./samples/ConceptVisualization/Quad-Tree_26-January-2021.py)
(26 January 2021)
    - Visualizes the concept of the quad-tree algorithm
2. [Python Oct-Tree](./samples/ConceptVisualization/Oct-Tree_26-January-2021.py)
(26 January 2021)
    - Visualizes the concept of the oct-tree algorithm
3. [Python 2D SPH](./samples/ConceptVisualization/SPH_2D_21-January-2021.py)
(21 January 2021)
    - Visualizes 2D gravitation SPH in python
3. [Python 3D SPH](./samples/ConceptVisualization/SPH_3D_23-January-2021.py)
(21 January 2021)
    - Visualizes 3D gravitation SPH in python