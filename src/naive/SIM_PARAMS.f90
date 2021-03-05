! this module contains all the shared simulation parameters
! TODO: add configuration file support, so I don't have to recompile each time
! TODO: add command line arguments?
MODULE SIM_PARAMS
    IMPLICIT NONE
    ! operation variables

    ! domain variables
    REAL, DIMENSION(3), PARAMETER :: DOMAIN = (/ 4.0, 4.0, 4.0 /)
    INTEGER, PARAMETER :: TIMESTEPS = 240           ! assumes 24 fps
    INTEGER, PARAMETER :: SUBSTEPS_PER_FRAME = 20

    ! simulation variables
    ! PCOUNT must be divisible by the number of nodes running on it
    ! TODO: make the method of partitioning PCOUNT more automatic
    INTEGER, PARAMETER :: PCOUNT = 81           ! number of particles
    REAL, PARAMETER :: MASS = 0.1               ! mass per particle
    REAL, PARAMETER :: PI = 4.0 / ATAN(1.0)     ! retrieve pi from computer
    REAL, PARAMETER :: G = 0.01                 ! gravitational constant
    REAL, PARAMETER :: SOFTENING = 0.04         ! softening parameter
END MODULE