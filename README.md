# nikolasn_phys495-spring2021
A software library for my PHYS 495 thesis project, Spring 2021, intended to
provide tools for astrophysical N-body simulations and visualization with
parallel capabilities

Please check the `develop` and `feature-*` branches for in-progress changes to
code! Some of the features are described below
- `feature-fmm-parallel`: This branch contains code for the Fast Multipole 
Method (FMM) which is intended to allow O(N) calculations on electrostatic
N-body problems
- `feature-sph-parallel`: This branch implements Smooth Particle Hydrodynamic
calculations for particles to handle internal pressures
- `feature-naive-parallel`: This is the first feature branch which used the
direct O(N^2) calculation approach to N-body simulations
