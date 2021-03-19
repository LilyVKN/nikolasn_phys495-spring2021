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
! TODO: Implement different parallelization
! TODO: Add progress readouts - e.g. system energy, simulation timer, etc.
PROGRAM main
    ! load the shared simulation parameters
    use sim_params
    use sim_log
    ! load the MPI FORTRAN header to access the MPI library
    include 'mpif.h'
    
    ! TODO: remove implicit variables
    
    ! particle simulation data
    REAL, DIMENSION(:,:), ALLOCATABLE :: pos, vel
    REAL, DIMENSION(:,:), ALLOCATABLE :: data_new       ! temporary data array
    REAL, DIMENSION(:,:,:), ALLOCATABLE :: force_mat
    REAL, TARGET :: energy = 0, KE = 0, PE_grav = 0
    REAL, TARGET :: E_start = 0, E_err = 0.0

    ! file management variables
    INTEGER, DIMENSION(8) :: date_values
    CHARACTER(len=15) :: timestamp
    CHARACTER(len=27) :: dir_name
    CHARACTER(len=80) :: filename

    ! MPI utility variables
    INTEGER, TARGET :: irank, inum_procs
    INTEGER :: status(MPI_STATUS_SIZE)

    ! logging tools
    TYPE(simulation_logger) :: logger
    TYPE(simlog_entry_string) :: method_log
    TYPE(simlog_entry_integer) :: node_cnt_log
    TYPE(simlog_entry_scinum) :: E_log, KE_log, PE_grav_log
    TYPE(simlog_entry_percentage) :: E_err_log

    ! logging timers
    TYPE(simlog_timer) :: force_calc_timer, submat_calc_timer, integrate_timer
    TYPE(simlog_subsection) :: calc_speed_section
    TYPE(simlog_entry_string) :: root_label
    TYPE(simlog_node_tracker) :: calc_speed_log
    TYPE(simlog_entry_time), TARGET :: force_calc_log, submat_calc_log, &
        integrate_log
    TYPE(simlog_entry_count), TARGET :: cache_log
    REAL, TARGET :: avg_force_calc = 0, avg_submat_calc = 0, &
        avg_integrate_calc = 0
    INTEGER, TARGET :: icache_hits = 0

    ! initialize MPI
    call MPI_INIT(ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD,irank,ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD,inum_procs,ierr)

    ! make sure mpirun called multiple nodes to make use of the parallelization
    IF (inum_procs <= 2) THEN
        WRITE(*,*) "no single process support yet :<"
        GOTO 10
    END IF

    ! make sure to optimize the parameters before allocation and simulation
    call optimize_params(inum_procs - 1)
    ALLOCATE(pos(PCOUNT,3))
    ALLOCATE(vel(PCOUNT,3))

    ! PROGRAM SETUP ===========================================================
    !   The root should decide the appropriate parallelization method based on
    !   node and particle counts

    ! calculate the partition data
    idiv_count = (inum_procs - 1) * SUBDIV  ! total number of subdivisions
    ipart = PCOUNT / idiv_count             ! number of particles per submatrix
    inode_part = PCOUNT / (inum_procs - 1)  ! total number of particles per node

    ! setup the logger
    call logger%init(irank,inum_procs,TIMESTEPS)

    force_calc_log%value => avg_force_calc
    force_calc_log%label = "Force Calcs.:"
    force_calc_log%tab_count = 1

    submat_calc_log%value => avg_submat_calc
    submat_calc_log%label = "Submat. Calcs.:"
    submat_calc_log%tab_count = 2

    cache_log%value => icache_hits
    cache_log%label = "Avg. Cache Hit:"
    cache_log%total = idiv_count
    cache_log%tab_count = 2

    integrate_log%value => avg_integrate_calc
    integrate_log%label = "Integration:"
    integrate_log%tab_count = 1

    call calc_speed_log%init(irank,inum_procs)
    call calc_speed_log%add_entry(force_calc_log)
    call calc_speed_log%add_entry(submat_calc_log)
    call calc_speed_log%add_entry(cache_log)
    call calc_speed_log%add_entry(integrate_log)

    ! determine which process is the root - i.e. the main node
    IF (irank.EQ.0) THEN

        ! SIMULATION SETUP ====================================================

        ! setup the logger start report
        method_log%value = "Direct-Evaluation"
        method_log%label = "Method:"

        node_cnt_log%value => inum_procs
        node_cnt_log%label = "MPI Node Count:"
        node_cnt_log%use_newline = .TRUE.

        call logger%add_start_report_entry(method_log)
        call logger%add_start_report_entry(node_cnt_log)
        call report_sim_params(logger)

        ! setup the logger values to track
        E_log%value => energy
        E_log%label = "Energies    Total:"
        call logger%add_frame_entry(E_log)
        call logger%add_mid_report_entry(E_log)
    
        KE_log%value => KE
        KE_log%label = "KE:"
        KE_log%tab_count = 3
        call logger%add_frame_entry(KE_log)
        call logger%add_mid_report_entry(KE_log)
    
        PE_grav_log%value => PE_grav
        PE_grav_log%label = "PE_grav:"
        PE_grav_log%tab_count = 3
        PE_grav_log%use_newline = .TRUE.
        call logger%add_frame_entry(PE_grav_log)
        call logger%add_mid_report_entry(PE_grav_log)
    
        E_err_log%value => E_err
        E_err_log%label = "Energy Error:"
        E_err_log%use_newline = .TRUE.
        call logger%add_mid_report_entry(E_err_log)

        root_label%value = " "
        root_label%label = "Root:"
        force_calc_log%label = "Calculations:"
        integrate_log%label = "Global Calcs.:"
        integrate_log%use_newline = .TRUE.
        calc_speed_section%title = "Calculation Speeds"
        call calc_speed_section%add_entry(root_label)
        call calc_speed_section%add_entry(force_calc_log)
        call calc_speed_section%add_entry(integrate_log)
        call calc_speed_section%add_entry(calc_speed_log)
        call logger%add_end_report_section(calc_speed_section)
        
        ! create directory for caching
        ! TODO: Separate file management into module for caching
        call date_and_time(values=date_values)
        WRITE(timestamp,'(I2.2,"-",I2.2,"-",I4.4,"_",I4.4)') date_values(3), &
            date_values(2), date_values(1), &
            ((date_values(5) * 100) + date_values(6))
        dir_name = 'SPH_fortran_' // timestamp
        call execute_command_line('mkdir -p ' // dir_name)

        ! create a directory for temporary files
        call execute_command_line('mkdir -p ./tmp/fortran')
        call execute_command_line('rm -f ./tmp/fortran/*')

        ! run initialization on the particles
        call part_init_default(pos,vel,PCOUNT,(DOMAIN(1)/2.0))

        ! calculate the initial energy
        KE = 0.5 * MASS * SUM(vel**2)

        ! allocate an array for storing temporary values for particle data
        ALLOCATE(data_new(PCOUNT,3))

        ! DELEGATE ============================================================
        !   might be necessary for complicated partitioning or job distribution
        !   perhaps user-determined delegation (e.g. some nodes do gravity
        !   calculations others do EM calculations)

        ! simulation loop (over all timesteps) using leap frog integration
        call logger%start()
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
                
                call force_calc_timer%start()
                ! calculate the kinetic energy of the system before the update
                KE = 0.5 * MASS * SUM(vel**2)
                call force_calc_timer%lap(.TRUE.)

                ! MERGE =======================================================

                ! retrieve each node's partition of position data
                data_new = 0.0      ! clear the temporary data array
                DO inode = 1, inum_procs-1
                    call MPI_RECV(data_new(((inode-1)*inode_part)+1: & 
                        (inode*inode_part),:),inode_part*3,MPI_REAL,inode, &
                        MPI_ANY_TAG,MPI_COMM_WORLD,status,ierr)
                END DO
                pos = data_new      ! update the positions

                ! retrieve each node's partition of velocity data
                data_new = 0.0      ! clear the temporary data array
                DO inode = 1, inum_procs-1
                    call MPI_RECV(data_new(((inode-1)*inode_part)+1: & 
                        (inode*inode_part),:),inode_part*3,MPI_REAL,inode, &
                        MPI_ANY_TAG,MPI_COMM_WORLD,status,ierr)
                END DO
                vel = data_new      ! update the velocities

                PE_grav = 0.0
                energy = 0.0
                ! retrieve the sum of all the node's calculated potential energy
                call MPI_REDUCE(energy,PE_grav,1,MPI_REAL,MPI_SUM,0, &
                    MPI_COMM_WORLD,ierr)

                ! GLOBAL CALCULATIONS/MANAGEMENT ==============================
                
                call integrate_timer%start()

                ! add the potential energy to the total energy
                energy = PE_grav + KE
                IF (i.EQ.1) E_start = energy
                E_err = ((energy - E_start) / E_start) * 100.0

                ! delete all the temporary submatrix cache's for last substep
                call execute_command_line('rm -f ./tmp/fortran/*')

                call integrate_timer%lap(.TRUE.)

            END DO

!            WRITE(*,'("Current energy: ",D12.6)') energy

            ! write the data to frame cache for rendering
            ! TODO: move to caching module
            WRITE(filename,'(A,"/SPH_fortran_",I4.4,".fcache")') dir_name, i
            OPEN(file=filename,unit=16)
            DO j = 1, PCOUNT
                WRITE(16,*) pos(j,:)
            END DO
            CLOSE(16)

            call logger%lap_frame()

        END DO

        ! send a quit message (-1) to all nodes
        i = -1
        call MPI_BCAST(i,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

        avg_force_calc = force_calc_timer%get_avg()
        avg_integrate_calc = integrate_timer%get_avg()
        call calc_speed_log%sync()

        call logger%finalize()
        call logger%run_end_report()

        ! safely deallocate the temporary data array
        DEALLOCATE(data_new)

    ELSE
        ! SIMULATION SETUP ====================================================

        ! allocate the force submatrix for calculation substeps
        ALLOCATE(force_mat(ipart,ipart,3))
        ! allocate a temporary array for intermediate particle calculations
        ALLOCATE(data_new(inode_part,3))

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

            energy = 0.0

            ! calculate this node's force submatrix
            data_new = 0.0
            isubdiv = 1

            ! start by calculating the upper triangle to help with caching
            call force_calc_timer%start()
            DO k = (irank - 1) * SUBDIV + 1, irank * SUBDIV
                DO l = k, k + idiv_count
                    ! calculate (or retrieve from cache) this submatrix
                    call submat_calc_timer%start()
                    call force_grav_submat(pos,k,MODULO(l-1,idiv_count)+1, &
                        ipart,force_mat,PE_grav,ihit)
                    call submat_calc_timer%lap(.TRUE.)

                    icache_hits = icache_hits + ihit

                    ! add the force contributions of this submatrix to the total
                    data_new((isubdiv-1)*ipart+1:isubdiv*ipart,:) = &
                        data_new((isubdiv-1)*ipart+1:isubdiv*ipart,:) +&
                        SUM(force_mat,1)

                    ! update the energy
                    energy = energy + PE_grav
                END DO
                isubdiv = isubdiv + 1
            END DO
            call force_calc_timer%lap(.TRUE.)

            call integrate_timer%start()

            ! kick-drift-kick
            pos(((irank-1)*inode_part)+1:(irank*inode_part),:) = & 
                pos(((irank-1)*inode_part)+1:(irank*inode_part),:) + & 
                (vel(((irank-1)*inode_part)+1:(irank*inode_part),:) * rperiod)
            vel(((irank-1)*inode_part)+1:(irank*inode_part),:) = &
                vel(((irank-1)*inode_part)+1:(irank*inode_part),:) + & 
                (data_new * (2.0 * rperiod))
            pos(((irank-1)*inode_part)+1:(irank*inode_part),:) = &
                pos(((irank-1)*inode_part)+1:(irank*inode_part),:) + &
                (vel(((irank-1)*inode_part)+1:(irank*inode_part),:) * rperiod)

            call integrate_timer%lap(.TRUE.)

            ! MERGE =======================================================

            ! send this node's partition of particle data - pos then vel
            call MPI_SEND(pos(((irank-1)*inode_part)+1:(irank*inode_part),:), &
                inode_part*3,MPI_REAL,0,0,MPI_COMM_WORLD,ierr)
            call MPI_SEND(vel(((irank-1)*inode_part)+1:(irank*inode_part),:), &
                inode_part*3,MPI_REAL,0,0,MPI_COMM_WORLD,ierr)

            ! send the potential energy contributions this node has tracked
            call MPI_REDUCE(energy,PE_grav,1,MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD, &
                ierr)
        END DO

        avg_force_calc = force_calc_timer%get_avg()
        avg_submat_calc = submat_calc_timer%get_avg()
        icache_hits = icache_hits / (TIMESTEPS * SUBSTEPS_PER_FRAME)
        avg_integrate_calc = integrate_timer%get_avg()
        call calc_speed_log%update()
        call calc_speed_log%sync()

        ! securely deallocate the temporary data memory
        DEALLOCATE(data_new)
        DEALLOCATE(force_mat)
    END IF

    ! safely deallocate the particle data
    DEALLOCATE(pos)
    DEALLOCATE(vel)

    ! Safely exit the MPI program
10  call MPI_FINALIZE(ierr)
END PROGRAM main

SUBROUTINE report_sim_params(logger)
    use sim_params
    use sim_log
    TYPE(simulation_logger), INTENT(INOUT) :: logger
    
    TYPE(simlog_subsection) :: op_param, domain_vars, sim_vars
    TYPE(simlog_entry_memory) :: mem_log
    TYPE(simlog_entry_integer) :: subdiv_log, timestep_log, substep_log, &
        pcount_log
    TYPE(simlog_entry_real) :: domx_log, domy_log, domz_log, g_log, soft_log, &
        mass_log
    TYPE(simlog_entry_string) :: dom_label, init_log

    ! prepare the operation parameter report
    ALLOCATE(mem_log%value)
    mem_log%value = 0
    mem_log%label = "Node Memory:"

    ALLOCATE(subdiv_log%value)
    subdiv_log%value = SUBDIV
    subdiv_log%label = "Subdivisions:" 

    op_param%title = "Operation Parameters"
    call op_param%add_entry(mem_log)
    call op_param%add_entry(subdiv_log)
    call logger%add_start_report_section(op_param)

    ! prepare the domain variables report
    dom_label%label = "Domain Size:"
    dom_label%value = " "
    ALLOCATE(domx_log%value)
    domx_log%value = DOMAIN(1)
    ALLOCATE(domy_log%value)
    domy_log%value = DOMAIN(2)
    ALLOCATE(domz_log%value)
    domz_log%value = DOMAIN(3)
    domx_log%label = "x-axis"
    domy_log%label = "y-axis"
    domz_log%label = "z-axis"
    domx_log%tab_count = 1
    domy_log%tab_count = 1
    domz_log%tab_count = 1

    ALLOCATE(timestep_log%value)
    timestep_log%value = TIMESTEPS
    timestep_log%label = "Timesteps:"

    ALLOCATE(substep_log%value,source=SUBSTEPS_PER_FRAME)
    !substep_log%value = SUBSTEPS_PER_FRAME
    substep_log%label = "Substeps:"
    substep_log%tab_count = 1

    domain_vars%title = "Domain Variables"
    call domain_vars%add_entry(dom_label)
    call domain_vars%add_entry(domx_log)
    call domain_vars%add_entry(domy_log)
    call domain_vars%add_entry(domz_log)
    call domain_vars%add_entry(timestep_log)
    call domain_vars%add_entry(substep_log)
    call logger%add_start_report_section(domain_vars)

    ! setup simulation variables report
    ALLOCATE(g_log%value)
    g_log%value = G
    g_log%label = "Grav. Constant:"
    ALLOCATE(soft_log%value)
    soft_log%value = SOFTENING
    soft_log%label = "Softening Param:"
    soft_log%use_newline = .TRUE.

    ALLOCATE(pcount_log%value)
    pcount_log%value = PCOUNT
    pcount_log%label = "Particle Count:"
    ALLOCATE(mass_log%value)
    mass_log%value = MASS
    mass_log%label = "Particle Mass:"
    mass_log%use_newline = .TRUE.

    init_log%value = "DEFAULT"
    init_log%label = "Init. Method:"

    sim_vars%title = "Simulation Variables"
    call sim_vars%add_entry(g_log)
    call sim_vars%add_entry(soft_log)
    call sim_vars%add_entry(pcount_log)
    call sim_vars%add_entry(mass_log)
    call sim_vars%add_entry(init_log)
    call logger%add_start_report_section(sim_vars)

END SUBROUTINE report_sim_params