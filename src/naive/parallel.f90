! Parallel N-Body Simulation - 4 March 2021 ===================================
!   This program simulates an N-body gravitational simulation by creating a
!   pairwise force matrix (i.e. f_ij represents the force between particles i
!   and j). The matrix is partitioned into square submatrices to parallelize
!   among MPI processes. This version of the program is naive, because it runs
!   pairwise calculations, it takes O(N^2) time. The simulation frames are
!   stored in a timestamped folder in custom .fcache files.
!
!   The MPI root collates the partitioned calculations from each node. As it
!   stands, each node is given a partition of the particles to calculate the
!   forces for. Therefore, they are given a row of submatrices to calculate and
!   they are each responsible for performing the leap-frog integration of its
!   respective particle dynamics.
!   
!   Dependencies --------------------------------------------------------------
!    - MPI library - This code was tested on MPICH
!   
!   References ----------------------------------------------------------------
!   .. [1] Pfalzner, S., & Gibbon, P. (2005). Basic Principles of the 
!       Hierarchical Tree Method. In Many-body tree methods in physics 
!       (pp. 9-18). Cambridge: Cambridge University Press.
!
! TODO: Create submatrix caching system to prevent redundancies
! TODO: Implement different parallelization
! TODO: Add progress readouts - e.g. system energy
PROGRAM main
    ! load the shared simulation parameters
    use SIM_PARAMS
    ! load the MPI FORTRAN header to access the MPI library
    include 'mpif.h'
    
    ! TODO: remove implicit variables
    
    ! particle simulation data
    REAL, DIMENSION(PCOUNT,3) :: pos, vel
    REAL, DIMENSION(:,:), ALLOCATABLE :: data_new       ! temporary data array
    REAL, DIMENSION(:,:,:), ALLOCATABLE :: force_mat

    ! file management variables
    INTEGER, DIMENSION(8) :: date_values
    CHARACTER(len=15) :: timestamp
    CHARACTER(len=27) :: dir_name
    CHARACTER(len=80) :: filename

    ! MPI utility variables
    INTEGER :: status(MPI_STATUS_SIZE)

    ! initialize MPI
    call MPI_INIT(ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD,irank,ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD,inum_procs,ierr)

    ! make sure mpirun called multiple nodes to make use of the parallelization
    IF (inum_procs <= 2) THEN
        WRITE(*,*) "no single process support yet :<"
        GOTO 10
    END IF

    ! DETERMINE METHOD ========================================================
    !   The root should decide the appropriate parallelization method based on
    !   node and particle counts

    ! calculate the partition
    ipartition = PCOUNT / (inum_procs - 1)

    ! determine which process is the root - i.e. the main node
    IF (irank .eq. 0) THEN
        ! SETUP SIMULATION ====================================================

        ! create directory for caching
        ! TODO: Separate file management into module for caching
        call date_and_time(values=date_values)
        WRITE(timestamp,'(I2.2,"-",I2.2,"-",I4.4,"_",I4.4)') date_values(3), &
            date_values(2), date_values(1), &
            ((date_values(5) * 100) + date_values(6))
        dir_name = 'SPH_fortran_' // timestamp
        call execute_command_line('mkdir -p ' // dir_name)

        ! create a directory for temporary files
        call execute_command_line('mkdir -p /tmp')

        ! run initialization on the particles
        call part_init_default(pos,vel,PCOUNT,(DOMAIN(1)/2.0))

        ! allocate an array for storing temporary values for particle data
        ALLOCATE(data_new(PCOUNT,3))

        ! DELEGATE ============================================================
        !   might be necessary for complicated partitioning or job distribution
        !   perhaps user-determined delegation (e.g. some nodes do gravity
        !   calculations others do EM calculations)

        ! simulation loop (over all timesteps) using leap frog integration
        DO i = 1, TIMESTEPS
            ! substep loop for refining the simulation between frames
            DO j = 1, SUBSTEPS_PER_FRAME

                ! SYNC ========================================================

                ! Global Data -------------------------------------------------
                
                ! signal the timestep
                call MPI_BCAST((i*SUBSTEPS_PER_FRAME)+j,1,MPI_INTEGER,0, &
                    MPI_COMM_WORLD,ierr)

                ! send collated position data
                call MPI_BCAST(pos,PCOUNT*3,MPI_REAL,0,MPI_COMM_WORLD,ierr)
                ! send collated velocity data
                call MPI_BCAST(vel,PCOUNT*3,MPI_REAL,0,MPI_COMM_WORLD,ierr)

                ! Node-Specific Data ------------------------------------------
                !   there is no node-specific data to be distributed here

                ! CALCULATE ===================================================
                !   usually only other nodes will be doing granular work

                ! MERGE =======================================================

                ! retrieve each node's partition of position data
                data_new = 0.0      ! clear the temporary data array
                DO inode = 1, inum_procs-1
                    call MPI_RECV(data_new(((inode-1)*ipartition)+1: & 
                        (inode*ipartition),:),ipartition*3,MPI_REAL,inode, &
                        MPI_ANY_TAG,MPI_COMM_WORLD,status,ierr)
                END DO
                pos = data_new      ! update the positions

                ! retrieve each node's partition of velocity data
                data_new = 0.0      ! clear the temporary data array
                DO inode = 1, inum_procs-1
                    call MPI_RECV(data_new(((inode-1)*ipartition)+1: & 
                        (inode*ipartition),:),ipartition*3,MPI_REAL,inode, &
                        MPI_ANY_TAG,MPI_COMM_WORLD,status,ierr)
                END DO
                vel = data_new      ! update the velocities

                ! GLOBAL CALCULATIONS/MANAGEMENT ==============================
                !   none - the nodes do the integration

            END DO

            ! write the data to frame cache for rendering
            ! TODO: move to caching module
            WRITE(filename,'(A,"/SPH_fortran_",I4.4,".fcache")') dir_name, i
            OPEN(file=filename,unit=16)
            DO j = 1, PCOUNT
                WRITE(16,*) pos(j,:)
            END DO
            CLOSE(16)

        END DO

        ! send a quit message (-1) to all nodes
        i = -1
        call MPI_BCAST(i,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

        ! safely deallocate the temporary data array
        DEALLOCATE(data_new)

    ELSE
        ! allocate the force submatrix for calculation substeps
        ALLOCATE(force_mat(ipartition,ipartition,3))
        ! allocate a temporary array for intermediate particle calculations
        ALLOCATE(data_new(ipartition,3))

        ! calculate half-substep time period for updating dynamics
        rperiod = 1.0 / (48.0 * SUBSTEPS_PER_FRAME)
        
        ! keep performing calculations for root until the quit signal is sent
        DO WHILE(i >= 0)

            ! SYNC ============================================================
            
            ! Global ----------------------------------------------------------

            ! retrieve the timestep signal
            call MPI_BCAST(i,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
            IF (i < 0) EXIT     ! if the the signal is negative, that means quit

            ! retrieve collated particle data - position then velocity
            call MPI_BCAST(pos,PCOUNT*3,MPI_REAL,0,MPI_COMM_WORLD,ierr)
            call MPI_BCAST(vel,PCOUNT*3,MPI_REAL,0,MPI_COMM_WORLD,ierr)
            
            ! Node-Specific Data ------------------------------------------
            !   there is no node-specific data to be received here

            ! CALCULATE ===================================================
            !   use the leapfrog kick-drift-kick method

            ! calculate this node's force submatrix
            data_new = 0.0
            DO j = 1, inum_procs - 1
                call force_grav_submat(pos,irank,j,ipartition,force_mat)
                data_new = data_new + SUM(force_mat,1)
            END DO

            ! kick-drift-kick
            pos(((irank-1)*ipartition)+1:(irank*ipartition),:) = & 
                pos(((irank-1)*ipartition)+1:(irank*ipartition),:) + & 
                (vel(((irank-1)*ipartition)+1:(irank*ipartition),:) * rperiod)
            vel(((irank-1)*ipartition)+1:(irank*ipartition),:) = &
                vel(((irank-1)*ipartition)+1:(irank*ipartition),:) + & 
                (data_new * (2.0 * rperiod))
            pos(((irank-1)*ipartition)+1:(irank*ipartition),:) = &
                pos(((irank-1)*ipartition)+1:(irank*ipartition),:) + &
                (vel(((irank-1)*ipartition)+1:(irank*ipartition),:) * rperiod)

            ! MERGE =======================================================

            ! send this node's partition of particle data - pos then vel
            call MPI_SEND(pos(((irank-1)*ipartition)+1:(irank*ipartition),:), &
                ipartition*3,MPI_REAL,0,0,MPI_COMM_WORLD,ierr)
            call MPI_SEND(vel(((irank-1)*ipartition)+1:(irank*ipartition),:), &
                ipartition*3,MPI_REAL,0,0,MPI_COMM_WORLD,ierr)
        END DO

        ! securely deallocate the temporary data memory
        DEALLOCATE(data_new)
        DEALLOCATE(force_mat)
    END IF

    ! Safely exit the MPI program
10  call MPI_FINALIZE(ierr)
END PROGRAM main