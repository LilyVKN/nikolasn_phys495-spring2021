# Fast Multipole Method N-Body Simulation - Design Documentation
**PHYS 495** | Nikolas Nguyen | nikolasvknguyen@outlook.com

The fast multipole method (FMM)was developed by Greengard and Rokhlin for
quicker, accurate solutions to electromagnetic numberical simulations. Under
ideal conditions, the method approaches O(N) complexity. It is, however,
computationally intensive to calculate the multipole expansion to suitable
levels of error; therefore, FMM should only used for sufficiently large
simulations.

*Note: This document is best viewed using a Markdown reader; however the
document was written to be as easily read without it as possible. The figures
used throughout are supplementary and can be found in the github repository.*

### Table of Contents
1. [Overview of FMM](#overview)
2. [Multipole Expansions](#multipole)
3. [Local Expansions](#local)
4. [Translations](#translations)
3. [MPI Implementation](#mpi_implement)

## Overview of FMM <a name=overview></a>
FMM is based on multiple levels of subdivision of the simulation domain. This
tree of subdivision is traversed downwards (from finest level to root), then
upwards (from root to finest level), before calculating the forces on each
particle and updating the positions.

## Multipole Expansions <a name=multipole></a>
The multipole expansions are derived from the electic potential in 3-dimensional
spherical coordinates. The potential can be described by an infinite series of
spherical harmonics which may be separated into a directionally-dependent
component and a distance-dependent component. When evaluated at positions far
more distant from the multipole center than the radius that encapsulates the
particles of the multipole, the multipole expansion converges.

## Local Expansions <a name=local></a>
The local expansions are derived from the same electrical potential as the
multipole expansion but for values within the multipole. In order for the
expansion to converge, the particles contributing to the multipole expansion
must be sufficiently far from the center relative to the positions where the
potential is being evaluated.

## Translations <a name=translations></a>
Since potentials follow the laws of superposition, each element of two multipole
expansions may be added together assuming they are measured in the same
coordinate space. Additionally, the expansions can be evaluated relative to a
new center by translating the coordinate space in which it was calculated.

## MPI Implementation <a name=mpi_implement></a>
As noted by [insert reference!] there is a bottleneck between the downward and
upward pass when number of nodes exceeds the number of cells to be calculated.
This means much of the computation in contributing daughter cells to the lower
level parents is done almost sequentially as is propogating the local expansions
upward.