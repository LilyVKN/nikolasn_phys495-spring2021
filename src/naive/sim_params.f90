! this module contains all the shared simulation parameters and the functions
! to be used to update them for the current simulation run
!
! TODO: add configuration file support, so I don't have to recompile each time
! TODO: add command line arguments?
MODULE sim_params
    IMPLICIT NONE
    ! operation variables

    ! domain variables
    REAL, DIMENSION(3), PARAMETER :: DOMAIN = (/ 4.0, 4.0, 4.0 /)
    INTEGER, PARAMETER :: TIMESTEPS = 240           ! assumes 24 fps
    INTEGER, PARAMETER :: SUBSTEPS_PER_FRAME = 20

    ! simulation variables
    INTEGER :: PCOUNT = 81                      ! number of particles
    REAL, PARAMETER :: MASS = 0.1               ! mass per particle
    REAL, PARAMETER :: PI = 4.0 / ATAN(1.0)     ! retrieve pi from computer
    REAL, PARAMETER :: G = 0.01                 ! gravitational constant
    REAL, PARAMETER :: SOFTENING = 0.04         ! softening parameter

CONTAINS

    ! this function is meant to perform the (sometimes necessary) optimizations
    ! of the simulation parameters based on the number of processes. currently
    ! it pares the values to prevent wasted process divisions - i.e. such that
    ! the particle count is divisible by the number of processes
    !
    !   arguments (in) ========================================================
    !     num_procs : INTEGER
    !       the number of worker processes contributing to the calculations
    SUBROUTINE optimize_params(num_procs)
        INTEGER, INTENT(IN) :: num_procs

        ! make sure the input number of processes is valid
        IF(num_procs <= 0) RETURN

        ! pare the extra values that would not be evenly divided among workers 
        PCOUNT = PCOUNT - MODULO(PCOUNT, num_procs)

    END SUBROUTINE optimize_params
END MODULE SIM_PARAMS