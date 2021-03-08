! this module contains all the shared simulation parameters and the functions
! to be used to update them for the current simulation run
!
! TODO: add configuration file support, so I don't have to recompile each time
! TODO: add command line arguments?
MODULE sim_params
    IMPLICIT NONE
    ! operation variables
    INTEGER, PARAMETER :: SUBDIV = 2            ! number of partitions per node

    ! domain variables
    REAL, DIMENSION(3), PARAMETER :: DOMAIN = (/ 4.0, 4.0, 4.0 /)
    INTEGER, PARAMETER :: TIMESTEPS = 120           ! assumes 24 fps
    INTEGER, PARAMETER :: SUBSTEPS_PER_FRAME = 20

    ! simulation variables
    INTEGER :: PCOUNT = 486                     ! number of particles
    REAL, PARAMETER :: MASS = 0.1               ! mass per particle
    REAL, PARAMETER :: PI = 4.0 / ATAN(1.0)     ! retrieve pi from computer
    REAL, PARAMETER :: G = 0.01                 ! gravitational constant
    REAL, PARAMETER :: SOFTENING = 0.04         ! softening parameter

CONTAINS

    ! this function is meant to perform the (sometimes necessary) optimizations
    ! of the simulation parameters based on the number of processes. currently
    ! it pares the values to prevent wasted process divisions - i.e. such that
    ! the particle count is divisible by the number of processes and the number
    ! of subdivided partitions per node
    !
    !   arguments (in) ========================================================
    !     num_procs : INTEGER
    !       the number of worker processes contributing to the calculations
    SUBROUTINE optimize_params(num_procs)
        INTEGER, INTENT(IN) :: num_procs

        ! make sure the input number of processes is valid
        IF(num_procs <= 0) RETURN

        ! pare the extra values that would not be evenly divided among workers 
        PCOUNT = PCOUNT - MODULO(PCOUNT, num_procs * SUBDIV)

    END SUBROUTINE optimize_params

    ! this function is meant to perform optimization for node memory usage.
    ! using the given memory allowed per node, it calculates (roughly) the
    ! necessary subdivisions of the matrix calculations. there is no padding, so
    ! the memory overhead should be considered in advance of this function
    !
    !   arguments (in) ========================================================
    !     num_procs : INTEGER
    !       the number of worker processes contributing to the calculations
    !     mem_per_node : INTEGER
    !       the desired maximum amount of memory (in bytes) to be used for
    !       calculating in each worker node
    !
    !   notes =================================================================
    !       this function overrides changes by optimize_params(), so only one of
    !       the two should be used for one program
    SUBROUTINE optimize_memory(num_procs,mem_per_node)
        INTEGER, INTENT(IN) :: num_procs
        INTEGER, INTENT(IN) :: mem_per_node

        INTEGER :: isize

        ! make sure the inputs are valid
        IF (num_procs <= 0) RETURN
        IF (mem_per_node < 24) RETURN   ! at least large enough for one particle

        ! calculate the number of data entries that can fit in a node given the
        ! base size per data entry is 12 (3 REALs for a 3D vector)
        isize = mem_per_node / 12

        

    END SUBROUTINE optimize_memory
END MODULE sim_params