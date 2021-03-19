MODULE sim_log
    PRIVATE
    INTEGER, PARAMETER :: SIMLOG_MAX_ELEMENTS = 16

    ! make certain procedures public
    PUBLIC :: format_time, line_log, line_section, line_subsection

    TYPE, PUBLIC, ABSTRACT :: simlog_element
    CONTAINS
        PROCEDURE(element_print_func), DEFERRED :: print
        PROCEDURE(element_update_func), DEFERRED :: update
    END TYPE simlog_element

    ABSTRACT INTERFACE
        SUBROUTINE element_print_func(this)
            import simlog_element
            CLASS(simlog_element), INTENT(INOUT) :: this
        END SUBROUTINE element_print_func

        SUBROUTINE element_update_func(this)
            import simlog_element
            CLASS(simlog_element), INTENT(INOUT) :: this
        END SUBROUTINE element_update_func
    END INTERFACE

    TYPE, PUBLIC :: simlog_element_ptr
        CLASS(simlog_element), POINTER :: ptr
    END TYPE simlog_element_ptr

    TYPE, PUBLIC, EXTENDS(simlog_element) :: simlog_collection
        TYPE(simlog_element_ptr) :: elements(SIMLOG_MAX_ELEMENTS)
        INTEGER :: num_elements = 0
    CONTAINS
        PROCEDURE :: print => collection_print
        PROCEDURE :: update => collection_update

        PROCEDURE :: add_element => collection_add_element
    END TYPE simlog_collection

    TYPE, PUBLIC, EXTENDS(simlog_collection) :: simlog_report
        CHARACTER(len=58) :: title
    CONTAINS
        PROCEDURE :: print => report_print

        PROCEDURE :: add_entry => report_add_entry
        PROCEDURE :: add_subsection => report_add_subsection
    END TYPE simlog_report

    TYPE, PUBLIC, EXTENDS(simlog_collection) :: simlog_subsection
        CHARACTER(len=58) :: title
    CONTAINS
        PROCEDURE :: print => subsection_print
        PROCEDURE :: add_entry => subsection_add_entry
    END TYPE simlog_subsection

    TYPE, PUBLIC:: simlog_timer
        REAL :: start_time, last_time, last_interval, avg_time, total_time
        INTEGER :: ilaps
        LOGICAL, PRIVATE :: bpaused = .FALSE., binitialized = .FALSE.
    CONTAINS
        PROCEDURE :: init => timer_init
        PROCEDURE :: lap => timer_lap
        PROCEDURE :: pause => timer_pause
        PROCEDURE :: start => timer_start

        PROCEDURE :: get_avg => timer_get_avg
        PROCEDURE :: get_current => timer_get_current
        PROCEDURE :: get_last => timer_get_last
        PROCEDURE :: get_total => timer_get_total
    END TYPE simlog_timer

    TYPE, ABSTRACT, EXTENDS(simlog_element) :: simlog_entry
        CHARACTER(len=36) :: label = " "
        CHARACTER(len=18) :: content = " ", units = " "
        INTEGER :: tab_count = 0
        LOGICAL :: use_newline = .FALSE.
    CONTAINS
        PROCEDURE :: print => entry_print
    END TYPE simlog_entry

    TYPE, PUBLIC :: simlog_entry_ptr
        CLASS(simlog_entry), POINTER :: ptr
    END TYPE

    TYPE, PUBLIC, EXTENDS(simlog_entry) :: simlog_entry_integer
        INTEGER, POINTER :: value
    CONTAINS
        PROCEDURE :: update => entry_integer_update
    END TYPE simlog_entry_integer

    TYPE, PUBLIC, EXTENDS(simlog_entry_integer) :: simlog_entry_count
        INTEGER :: total
    CONTAINS
        PROCEDURE :: update => entry_count_update
    END TYPE simlog_entry_count

    TYPE, PUBLIC, EXTENDS(simlog_entry_integer) :: simlog_entry_memory
    CONTAINS
        PROCEDURE :: update => entry_memory_update
    END TYPE simlog_entry_memory

    TYPE, PUBLIC, EXTENDS(simlog_entry) :: simlog_entry_real
        REAL, POINTER :: value
    CONTAINS
        PROCEDURE :: update => entry_real_update
    END TYPE simlog_entry_real

    TYPE, PUBLIC, EXTENDS(simlog_entry_real) :: simlog_entry_percentage
    CONTAINS
        PROCEDURE :: update => entry_percentage_update
    END TYPE simlog_entry_percentage

    TYPE, PUBLIC, EXTENDS(simlog_entry_real) :: simlog_entry_scinum
    CONTAINS
        PROCEDURE :: update => entry_scinum_update
    END TYPE simlog_entry_scinum

    TYPE, PUBLIC, EXTENDS(simlog_entry_real) :: simlog_entry_time
        REAL :: time
    CONTAINS
        PROCEDURE :: update => entry_time_update
    END TYPE simlog_entry_time

    TYPE, PUBLIC, EXTENDS(simlog_entry) :: simlog_entry_string
        CHARACTER(len=18) :: value
    CONTAINS
        PROCEDURE :: update => entry_string_update
    END TYPE simlog_entry_string

    TYPE, PUBLIC, EXTENDS(simlog_entry) :: simlog_entry_mpi
        TYPE(simlog_entry_ptr) :: sync_entry
        INTEGER :: irank = 1, inum_procs = 1, isource = 1
        LOGICAL, PRIVATE :: binitialized = .FALSE.
    CONTAINS
        PROCEDURE :: init => entry_mpi_init
        PROCEDURE :: link => entry_mpi_link
        PROCEDURE :: sync => entry_mpi_sync
        PROCEDURE :: update => entry_mpi_update
    END TYPE simlog_entry_mpi

    TYPE, PUBLIC, EXTENDS(simlog_entry) :: simlog_node_tracker
        TYPE(simlog_entry_string), DIMENSION(:), ALLOCATABLE :: node_labels
        TYPE(simlog_entry_mpi), DIMENSION(:,:), ALLOCATABLE :: node_entries
        INTEGER :: irank = 1, inum_procs = 1, num_elements = 0
        LOGICAL, PRIVATE :: binitialized = .FALSE.
    CONTAINS
        PROCEDURE :: print => node_tracker_print
        PROCEDURE :: update => node_tracker_update

        PROCEDURE :: init => node_tracker_init
        PROCEDURE :: sync => node_tracker_sync

        PROCEDURE :: add_entry => node_tracker_add_entry
    END TYPE simlog_node_tracker

    TYPE, PUBLIC :: simulation_logger
        TYPE(simlog_timer) :: global_timer
        REAL :: avg_time, total_time

        INTEGER :: mid_report_trigger = 5
        LOGICAL :: log_frames = .TRUE.

        TYPE(simlog_collection), PRIVATE :: frame_report
        TYPE(simlog_report), PRIVATE :: start_report, mid_report, end_report

        INTEGER, PRIVATE :: irank, inum_procs, total_frames
        INTEGER, PRIVATE :: num_elements
    CONTAINS
        PROCEDURE :: init => logger_init
        PROCEDURE :: start => logger_start
        PROCEDURE :: lap_frame => logger_lap_frame
        PROCEDURE :: finalize => logger_finalize
        PROCEDURE :: run_end_report => logger_run_end_report

        PROCEDURE :: add_frame_entry => logger_add_frame_entry
        PROCEDURE :: add_start_report_entry => logger_add_start_entry
        PROCEDURE :: add_start_report_section => logger_add_start_section
        PROCEDURE :: add_mid_report_entry => logger_add_mid_entry
        PROCEDURE :: add_mid_report_section => logger_add_mid_section
        PROCEDURE :: add_end_report_entry => logger_add_end_entry
        PROCEDURE :: add_end_report_section => logger_add_end_section

        PROCEDURE :: print_frame_report => logger_print_frame_report
    END TYPE simulation_logger

CONTAINS

! collection procedures =======================================================

    SUBROUTINE collection_add_element(this,element)
        CLASS(simlog_collection), INTENT(INOUT) :: this
        CLASS(simlog_element), TARGET, INTENT(IN) :: element

        IF (this%num_elements >= SIMLOG_MAX_ELEMENTS) RETURN
        this%num_elements = this%num_elements + 1

        ALLOCATE(this%elements(this%num_elements)%ptr,source=element)
    END SUBROUTINE collection_add_element

    SUBROUTINE collection_print(this)
        CLASS(simlog_collection), INTENT(INOUT) :: this

        IF (this%num_elements <= 0) RETURN

        call this%update()
        DO i = 1, this%num_elements
            IF (ASSOCIATED(this%elements(i)%ptr)) THEN
                call this%elements(i)%ptr%print()
            ENDIF
        END DO
    END SUBROUTINE collection_print

    SUBROUTINE collection_update(this)
        CLASS(simlog_collection), INTENT(INOUT) :: this

        IF (this%num_elements <= 0) RETURN

        DO i = 1, this%num_elements
            IF (ASSOCIATED(this%elements(i)%ptr)) THEN
                call this%elements(i)%ptr%update()
            ENDIF
        END DO

    END SUBROUTINE collection_update

! entry procedures ============================================================

    SUBROUTINE entry_print(this)
        CLASS(simlog_entry), INTENT(INOUT) :: this

        CHARACTER(len=80) :: line

        line = " "
        WRITE(line,'(A,T63,A)') line, TRIM(this%units)
        WRITE(line,'(A,T44,A)') line, this%content

        SELECT CASE(this%tab_count)
        CASE (1)
            WRITE(line,'(A,T29,A)') line, TRIM(this%label)
        CASE (2)
            WRITE(line,'(A,T33,A)') line, TRIM(this%label)
        CASE (3)
            WRITE(line,'(A,T37,A)') line, TRIM(this%label)
        CASE DEFAULT
            WRITE(line,'(A,T25,A)') line, TRIM(this%label)
        END SELECT
        WRITE(*,*) line

        IF (this%use_newline) WRITE(*,*)

    END SUBROUTINE entry_print

    SUBROUTINE entry_count_update(this)
        CLASS(simlog_entry_count), INTENT(INOUT) :: this

        CHARACTER(len=15) :: count

        this%content = " "
        WRITE(this%content,'(I18)') this%value

        this%units = " "
        WRITE(count,'(I15)') this%total
        WRITE(this%units,'(A,T1,A,T3,A)') this%units, "/", ADJUSTL(count)
    END SUBROUTINE

    SUBROUTINE entry_integer_update(this)
        CLASS(simlog_entry_integer), INTENT(INOUT) :: this

        this%content = " "
        WRITE(this%content,'(I18)') this%value
    END SUBROUTINE entry_integer_update

    SUBROUTINE entry_memory_update(this)
        CLASS(simlog_entry_memory), INTENT(INOUT) :: this

        this%content = " "
        this%units = " "
        IF (this%value >= 10000000) THEN
            WRITE(this%content,'(I18)') (this%value / 1000000)
            this%units = "GB"
        ELSE IF (this%value >= 10000) THEN
            WRITE(this%content,'(I18)') (this%value / 1000)
            this%units = "MB"
        ELSE
            WRITE(this%content,'(I18)') this%value
            this%units = "KB"
        ENDIF
    END SUBROUTINE entry_memory_update

    SUBROUTINE entry_mpi_init(this,irank,inum_procs,isender)
        CLASS(simlog_entry_mpi), INTENT(INOUT) :: this
        INTEGER, INTENT(IN) :: irank, inum_procs, isender

        IF (irank < 0) RETURN
        IF (inum_procs < 0) RETURN
        IF (isender < 0) RETURN

        this%irank = irank
        this%inum_procs = inum_procs
        this%isource = isender
        this%binitialized = .TRUE.
    END SUBROUTINE entry_mpi_init

    SUBROUTINE entry_mpi_link(this,new_entry)
        CLASS(simlog_entry_mpi), INTENT(INOUT) :: this
        CLASS(simlog_entry), TARGET, INTENT(IN) :: new_entry

        IF (.NOT.this%binitialized) RETURN
        ALLOCATE(this%sync_entry%ptr,source=new_entry)

        IF (this%irank.EQ.0) THEN
            this%label = new_entry%label
            this%content = "null"
            this%units = new_entry%units
            this%tab_count = new_entry%tab_count
            this%use_newline = new_entry%use_newline
        ENDIF
    END SUBROUTINE entry_mpi_link

    SUBROUTINE entry_mpi_sync(this)
        include 'mpif.h'
        CLASS(simlog_entry_mpi), INTENT(INOUT) :: this

        CHARACTER(len=18) :: buffer

        IF (this%irank.EQ.0) THEN
            this%content = " "
            call MPI_RECV(buffer,18,MPI_CHARACTER,this%isource, &
                MPI_ANY_TAG,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
            WRITE(this%content,'(T1,A)') buffer

            this%units = " "
            call MPI_RECV(buffer,18,MPI_CHARACTER,this%isource, &
                MPI_ANY_TAG,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
            WRITE(this%units,'(T1,A)') buffer
        ELSE
            call MPI_SEND(this%sync_entry%ptr%content,18,MPI_CHARACTER,0,0, &
                MPI_COMM_WORLD,ierr)
            call MPI_SEND(this%sync_entry%ptr%units,18,MPI_CHARACTER,0,0, &
                MPI_COMM_WORLD,ierr)
        ENDIF
    END SUBROUTINE entry_mpi_sync

    SUBROUTINE entry_mpi_update(this)
        CLASS(simlog_entry_mpi), INTENT(INOUT) :: this

        CHARACTER(len=18) :: line

        IF (.NOT.this%binitialized) RETURN

        IF (this%irank.EQ.0) THEN
            this%content = ADJUSTR(this%content)
        ELSE
            IF (.NOT.ASSOCIATED(this%sync_entry%ptr)) RETURN
            call this%sync_entry%ptr%update()
        ENDIF
    END SUBROUTINE entry_mpi_update

    SUBROUTINE entry_percentage_update(this)
        CLASS(simlog_entry_percentage), INTENT(INOUT) :: this

        REAL :: percentage

        this%content = " "
        WRITE(this%content,'(A,T3,SP,F16.5)') this%content, this%value

        this%units = "%"
    END SUBROUTINE entry_percentage_update

    SUBROUTINE entry_real_update(this)
        CLASS(simlog_entry_real), INTENT(INOUT) :: this
        
        WRITE(this%content,'(F18.6)') this%value
        DO i = 18,1,-1
            IF (this%content(i:i) /= "0") THEN
                index = i
                IF (this%content(i:i) == ".") index = i + 1
                EXIT
            ENDIF
        END DO

        WRITE(this%content,*) this%content(1:index)
        this%content = ADJUSTR(this%content)
    END SUBROUTINE entry_real_update

    SUBROUTINE entry_scinum_update(this)
        CLASS(simlog_entry_scinum), INTENT(INOUT) :: this

        this%content = " "
        WRITE(this%content,'(A,T6,SP,ES13.6)') this%content, this%value
    END SUBROUTINE

    SUBROUTINE entry_string_update(this)
        CLASS(simlog_entry_string), INTENT(INOUT) :: this
        this%content = ADJUSTR(this%value)
    END SUBROUTINE entry_string_update

    SUBROUTINE entry_time_update(this)
        CLASS(simlog_entry_time), INTENT(INOUT) :: this

        CHARACTER(len=15) :: time_line

        call format_time(this%value,time_line)
        this%content = " "
        WRITE(this%content,'(A,T4,A)') this%content, time_line
    END SUBROUTINE entry_time_update

! logger procedures ===========================================================

    SUBROUTINE logger_add_end_entry(this,new_entry)
        CLASS(simulation_logger), INTENT(INOUT) :: this
        CLASS(simlog_entry), TARGET, INTENT(IN) :: new_entry
        call this%end_report%add_entry(new_entry)
    END SUBROUTINE logger_add_end_entry

    SUBROUTINE logger_add_end_section(this,new_section)
        CLASS(simulation_logger), INTENT(INOUT) :: this
        CLASS(simlog_subsection), INTENT(IN) :: new_section
        IF (this%irank == 0) call this%end_report%add_subsection(new_section)
    END SUBROUTINE logger_add_end_section

    SUBROUTINE logger_add_frame_entry(this,tracker)
        CLASS(simulation_logger), INTENT(INOUT) :: this
        CLASS(simlog_element), TARGET, INTENT(IN) :: tracker
        call this%frame_report%add_element(tracker)
    END SUBROUTINE logger_add_frame_entry

    SUBROUTINE logger_add_mid_entry(this,new_entry)
        CLASS(simulation_logger), INTENT(INOUT) :: this
        CLASS(simlog_entry), TARGET, INTENT(IN) :: new_entry
        call this%mid_report%add_entry(new_entry)
    END SUBROUTINE logger_add_mid_entry

    SUBROUTINE logger_add_mid_section(this,new_section)
        CLASS(simulation_logger), INTENT(INOUT) :: this
        CLASS(simlog_subsection), POINTER, INTENT(IN) :: new_section
        call this%mid_report%add_subsection(new_section)
    END SUBROUTINE logger_add_mid_section

    SUBROUTINE logger_add_start_entry(this,new_entry)
        CLASS(simulation_logger), INTENT(INOUT) :: this
        CLASS(simlog_entry), TARGET, INTENT(IN) :: new_entry
        call this%start_report%add_entry(new_entry)
    END SUBROUTINE logger_add_start_entry

    SUBROUTINE logger_add_start_section(this,new_section)
        CLASS(simulation_logger), INTENT(INOUT) :: this
        CLASS(simlog_subsection), TARGET, INTENT(IN) :: new_section
        call this%start_report%add_subsection(new_section)
    END SUBROUTINE logger_add_start_section

    SUBROUTINE logger_lap_frame(this)
        CLASS(simulation_logger), INTENT(INOUT) :: this

        ! set a lap for the global simulation timer which will track the frames
        call this%global_timer%lap(.FALSE.)
        this%total_time = this%global_timer%get_total()

        ! only log the frame data if the user has requested
        IF (this%log_frames) call this%print_frame_report()

        IF (this%global_timer%ilaps == this%total_frames) RETURN

        ! check if the mid report should be triggered this frame
        IF (MODULO(this%global_timer%ilaps, &
                this%mid_report_trigger) == 0) THEN
            call this%mid_report%print()
        ENDIF

    END SUBROUTINE logger_lap_frame

    SUBROUTINE logger_finalize(this)
        CLASS(simulation_logger), INTENT(INOUT) :: this

        CHARACTER(len=80) :: log_line, log_endline

        ! print the end message
        call line_section("End Simulation",14,log_line,log_endline)
        WRITE(*,*) log_line

        ! update the timer data
        this%avg_time = this%global_timer%get_avg()
        this%total_time = this%global_timer%get_total()
    END SUBROUTINE logger_finalize

    SUBROUTINE logger_force_sync(this)
        CLASS(simulation_logger), INTENT(IN) :: this

        num = this%irank
    END SUBROUTINE logger_force_sync

    SUBROUTINE logger_init(this,irank,inum_procs,total_frames)
        CLASS(simulation_logger), TARGET, INTENT(INOUT) :: this
        INTEGER, INTENT(IN) :: irank, inum_procs, total_frames

        TYPE(simlog_entry_count) :: frame_counter
        TYPE(simlog_entry_time) :: sim_timer
        TYPE(simlog_entry_time) :: total_timer

        this%irank = irank
        this%inum_procs = inum_procs
        this%total_frames = total_frames

        IF (irank /= 0) RETURN

        ! initialize the simulation start report
        this%start_report%title = "Simulation Start Report"

        ! initialize the intermediate report
        this%mid_report%title = "Intermediate Report"

        ! set the values to track for the frame counter and running timer
        frame_counter%label = "Frames:"
        frame_counter%value => this%global_timer%ilaps
        frame_counter%total = this%total_frames
        sim_timer%label = "Current Run-Time:"
        sim_timer%value => this%total_time
        sim_timer%use_newline = .TRUE.

        ! add the entries to the intermediate report
        call this%mid_report%add_entry(frame_counter)
        call this%mid_report%add_entry(sim_timer)

        this%end_report%title = "Completion Report"

        ! set the values to the total simulation time entry
        total_timer%label = "Simulation Time:"
        total_timer%value => this%total_time
        total_timer%use_newline = .TRUE.
        
        ! add the entries to the end report
        call this%end_report%add_entry(total_timer)

    END SUBROUTINE logger_init

    SUBROUTINE logger_print_frame_report(this)
        CLASS(simulation_logger), INTENT(INOUT) :: this

        CHARACTER(len=80) :: log_line
        CHARACTER(len=60) :: log_text
        CHARACTER(len=15) :: time_text

        IF (this%irank /= 0) RETURN

        ! format the header with basic frame data
        WRITE(log_text,'("- Frame ",I6.4," - Time:")') &
        this%global_timer%ilaps

        call format_time(this%global_timer%last_interval,time_text)
        WRITE(log_text,'(A,T28,A," / ")') log_text, time_text
        call format_time(this%global_timer%total_time,time_text)
        WRITE(log_text,'(A,T46,A)') log_text, time_text

        call line_log(log_text,60,log_line)

        WRITE(*,*) log_line
        
        ! print the data within the report
        call this%frame_report%print()

    END SUBROUTINE logger_print_frame_report

    SUBROUTINE logger_run_end_report(this)
        CLASS(simulation_logger), INTENT(INOUT) :: this
        IF(this%irank.EQ.0) call this%end_report%print()
    END SUBROUTINE logger_run_end_report

    SUBROUTINE logger_start(this)
        CLASS(simulation_logger), INTENT(INOUT) :: this

        CHARACTER(len=80) :: log_line, log_endline

        IF(this%irank /= 0) THEN
            call this%global_timer%init()
            RETURN
        ENDIF

        ! print the start report just before beginning the simulation
        call this%start_report%print()

        ! initialize the global simulation timer
        call this%global_timer%init()

        call line_section("Start Simulation",16,log_line,log_endline)
        WRITE(*,*) log_line
        IF(this%log_frames) call this%print_frame_report()
    END SUBROUTINE logger_start

! node-tracker procedures =====================================================

    SUBROUTINE node_tracker_add_entry(this,new_entry)
        CLASS(simlog_node_tracker), INTENT(INOUT) :: this
        CLASS(simlog_entry), TARGET, INTENT(IN) :: new_entry

        TYPE(simlog_entry_mpi) :: mpi_entry

        IF (.NOT.this%binitialized) RETURN
        IF (this%num_elements >= SIMLOG_MAX_ELEMENTS) RETURN
        this%num_elements = this%num_elements + 1

        IF (this%irank.EQ.0) THEN
            DO i = 1, this%inum_procs
                call this%node_entries(i,this%num_elements)%init(this%irank, &
                        this%inum_procs+1,i)
                call this%node_entries(i,this%num_elements)%link(new_entry)
            END DO
        ELSE
            call this%node_entries(1,this%num_elements)%init(this%irank, &
                    this%inum_procs+1,this%irank)
            call this%node_entries(1,this%num_elements)%link(new_entry)
            call this%node_entries(1,this%num_elements)%update()
        ENDIF
    END SUBROUTINE

    SUBROUTINE node_tracker_init(this,irank,inum_procs)
        CLASS(simlog_node_tracker), INTENT(INOUT) :: this
        INTEGER, INTENT(IN) :: irank, inum_procs

        CHARACTER(len=30) :: label = " "

        this%irank = irank
        this%inum_procs = inum_procs - 1

        IF (irank == 0) THEN
            ALLOCATE(this%node_labels(this%inum_procs))
            DO i = 1, this%inum_procs
                WRITE(label,*) i
                label = ADJUSTL(label)
                this%node_labels(i)%label = " "
                WRITE(this%node_labels(i)%label,'("Node ",A,":")') TRIM(label)
            END DO

            ALLOCATE(this%node_entries(this%inum_procs,SIMLOG_MAX_ELEMENTS))
        ELSE
            ALLOCATE(this%node_entries(1,SIMLOG_MAX_ELEMENTS))
        ENDIF

        this%binitialized = .TRUE.
    END SUBROUTINE

    SUBROUTINE node_tracker_print(this)
        CLASS(simlog_node_tracker), INTENT(INOUT) :: this

        IF (this%irank.NE.0) RETURN
        IF (this%num_elements <= 0) RETURN

        DO i = 1, this%inum_procs
            call this%node_labels(i)%print()
            DO j = 1, this%num_elements
                call this%node_entries(i,j)%print()
            END DO
        END DO
    END SUBROUTINE node_tracker_print

    SUBROUTINE node_tracker_sync(this)
        CLASS(simlog_node_tracker), INTENT(INOUT) :: this

        IF (.NOT.this%binitialized) RETURN
        IF (this%num_elements <= 0) RETURN

        IF (irank.EQ.0) THEN
            DO i = 1, this%num_elements
                DO j = 1, this%inum_procs
                    call this%node_entries(j,i)%sync()
                END DO
            END DO
        ELSE
            DO i = 1, this%num_elements
                call this%node_entries(1,i)%sync()
            END DO
        ENDIF
    END SUBROUTINE node_tracker_sync

    SUBROUTINE node_tracker_update(this)
        CLASS(simlog_node_tracker), INTENT(INOUT) :: this

        IF (this%irank.EQ.0) THEN
            DO i = 1, this%num_elements
                DO j = 1, this%inum_procs
                    call this%node_entries(j,i)%update()
                END DO
            END DO
        ELSE
            DO i = 1, this%num_elements
                call this%node_entries(1,i)%update()
            END DO
        ENDIF
    END SUBROUTINE node_tracker_update

! report procedures ===========================================================
    
    SUBROUTINE report_add_entry(this,log_entry)
        CLASS(simlog_report), INTENT(INOUT) :: this
        CLASS(simlog_entry), TARGET, INTENT(IN) :: log_entry

        IF (this%num_elements >= SIMLOG_MAX_ELEMENTS) RETURN
        call this%add_element(log_entry)
    END SUBROUTINE report_add_entry

    SUBROUTINE report_add_subsection(this,log_subsection)
        CLASS(simlog_report), INTENT(INOUT) :: this
        CLASS(simlog_subsection), TARGET, INTENT(IN) :: log_subsection
        
        IF (this%num_elements >= SIMLOG_MAX_ELEMENTS) RETURN
        
        call this%add_element(log_subsection)
    END SUBROUTINE report_add_subsection

    SUBROUTINE report_print(this)
        CLASS(simlog_report), INTENT(INOUT) :: this

        CHARACTER(len=80) :: line, endline

        WRITE(*,*)
        call line_section(this%title,58,line,endline)
        WRITE(*,*) line

        call this%simlog_collection%print()

        WRITE(*,*) endline
        WRITE(*,*)

    END SUBROUTINE report_print

! subsection procedures =======================================================

    SUBROUTINE subsection_add_entry(this,log_entry)
        CLASS(simlog_subsection), INTENT(INOUT) :: this
        CLASS(simlog_entry), TARGET, INTENT(IN) :: log_entry

        IF (this%num_elements >= SIMLOG_MAX_ELEMENTS) RETURN
        call this%simlog_collection%add_element(log_entry)
    END SUBROUTINE subsection_add_entry

    SUBROUTINE subsection_print(this)
        CLASS(simlog_subsection), INTENT(INOUT) :: this

        CHARACTER(len=80) :: line

        line = " "
        call line_subsection(this%title,58,line)
        WRITE(*,*) line

        call this%simlog_collection%print()

        WRITE(*,*)

    END SUBROUTINE subsection_print

! timer procedures ============================================================

    SUBROUTINE timer_init(this)
        CLASS(simlog_timer), INTENT(INOUT) :: this

        call CPU_TIME(this%start_time)
        this%last_time = this%start_time

        this%binitialized = .TRUE.

    END SUBROUTINE timer_init

    SUBROUTINE timer_lap(this,bpause)
        CLASS(simlog_timer), INTENT(INOUT) :: this
        LOGICAL, INTENT(IN) :: bpause

        IF (this%bpaused) RETURN

        ! increment the frame counter for the timer
        this%ilaps = this%ilaps + 1

        ! update the timer
        this%last_interval = this%last_time
        call CPU_TIME(this%last_time)
        this%last_interval = this%last_time - this%last_interval

        ! update the total and average times
        this%total_time = this%total_time + this%last_interval
        this%avg_time = this%total_time / this%ilaps

        IF (bpause) this%bpaused = .TRUE.

    END SUBROUTINE timer_lap

    FUNCTION timer_get_avg(this) RESULT(avg)
        CLASS(simlog_timer), INTENT(IN) :: this
        REAL :: avg
        avg = this%avg_time
    END FUNCTION timer_get_avg

    FUNCTION timer_get_current(this) RESULT(current)
        CLASS(simlog_timer), INTENT(IN) :: this
        REAL :: current

        call CPU_TIME(current)
        current = current - this%last_time 
    END FUNCTION

    FUNCTION timer_get_last(this) RESULT(last)
        CLASS(simlog_timer), INTENT(IN) :: this
        REAL :: last
        last = this%last_interval
    END FUNCTION timer_get_last

    FUNCTION timer_get_total(this) RESULT(total)
        CLASS(simlog_timer), INTENT(IN) :: this
        REAL :: total

        total = this%last_time - this%start_time
    END FUNCTION timer_get_total

    SUBROUTINE timer_pause(this)
        CLASS(simlog_timer), INTENT(INOUT) :: this

        IF (this%bpaused) RETURN

        this%last_interval = this%last_time
        call CPU_TIME(this%last_time)
        this%last_interval = this%last_time - this%last_interval

        this%total_time = this%total_time + this%last_interval
        this%avg_time = this%total_time / this%ilaps

        this%bpaused = .TRUE.
        
    END SUBROUTINE timer_pause

    SUBROUTINE timer_start(this)
        CLASS(simlog_timer), INTENT(INOUT) :: this

        IF (.NOT.this%binitialized) call this%init()
        IF (.NOT.this%bpaused) RETURN

        call CPU_TIME(this%last_time)

        this%bpaused = .TRUE.
    END SUBROUTINE

! utility procedures ==========================================================

    ! format the given time (in seconds) to a formatted character array
    !
    !   arguments (in) --------------------------------------------------------
    !     time : REAL
    !       the time (in seconds)
    !
    !   arguments (out)
    !     out : CHARACTER(len=15)
    !       the formatted time in the format HH:MM:SS.SSS (hour can be 5 digits)
    SUBROUTINE format_time(time,out)
        REAL, INTENT(IN) :: time
        CHARACTER(len=15), INTENT(OUT) :: out

        INTEGER :: hour, minute, second
        REAL :: millis

        ! calculate the time intervals
        millis = MODULO(time,1.0)
        second = MODULO(INT(time), 60)
        minute = MODULO(INT(time/60.0),60)
        hour = INT(time/3600)

        WRITE(out,'(I5,":",I2.2,":",I2.2,F4.3)') hour, minute, second, millis

    END SUBROUTINE format_time

    ! creates a log line with the time stamp and given text
    !
    !   arguments (in) --------------------------------------------------------
    !     text : CHARACTER(len=ilen)
    !       the title text to insert into the section divider
    !     ilen : INTEGER
    !       the length of the text input. anything over 58 characters is cut
    !
    !   arguments (out) -----------------------------------------------------
    !     line : CHARACTER(len=80)
    !       the 80-character print line with the timestamp and text filled
    SUBROUTINE line_log(text,ilen,line)
        CHARACTER(len=ilen), INTENT(IN) :: text
        CHARACTER(len=80), INTENT(OUT) :: line

        ! add the timestamp to the log line
        line = " "
        call line_timestamp(line)

        ! add the text to the line
        IF (ilen > 60) THEN
            WRITE(line,'(A,T20,A60)') line, text
        ELSE
            WRITE(line,'(A,T20,A)') line, text
        ENDIF

    END SUBROUTINE line_log

    ! creates the section divider print line with timestamp and a given title
    !
    !   arguments (in) --------------------------------------------------------
    !     text : CHARACTER(len=ilen)
    !       the title text to insert into the section divider
    !     ilen : INTEGER
    !       the length of the text input. anything over 58 characters is cut
    !
    !   arguments (out) -----------------------------------------------------
    !     line : CHARACTER(len=80)
    !       the 80-character print line that will be filled with the section
    !       divider including a timestamp, a dividing line ('===') and the title
    !     endline : CHARACTER(len=80)
    !       the 80-character print line that will be filled with the section
    !       divider for the end of the section including just the dividing line
    !       and " End " plus the title
    SUBROUTINE line_section(text,ilen,line,endline)
        CHARACTER(len=ilen), INTENT(IN) :: text
        CHARACTER(len=80), INTENT(OUT) :: line, endline

        ! start by filling the lines with the dividers
        WRITE(line,'(79A)') ("=", i=1,79)
        WRITE(endline,'(T20,60A)') ("=", i=1,60)
        WRITE(endline,'(A,T21,A)') endline, " End "

        ! add the timestamp to the section header line
        call line_timestamp(line)

        ! insert the title based on its length
        IF (ilen > 58) THEN
            WRITE(line,'(A,T21," ",A58," ")') line, TRIM(text)
            WRITE(endline,'(A,T26,A,A54," ")') endline, TRIM(text)
        ELSE IF (ilen > 54) THEN
            WRITE(line,'(A,T21," ",A," ")') line, TRIM(text)
            WRITE(endline,'(A,T26,A,A54," ")') endline, TRIM(text)
        ELSE
            WRITE(line,'(A,T21," ",A," ")') line, TRIM(text)
            WRITE(endline,'(A,T26,A," ")') endline, TRIM(text)
        ENDIF

    END SUBROUTINE line_section

    ! creates the subsection divider print line with timestamp and a given title
    !
    !   arguments (in) --------------------------------------------------------
    !     text : CHARACTER(len=ilen)
    !       the title text to insert into the subsection divider
    !     ilen : INTEGER
    !       the length of the text input. anything over 58 characters is cut
    !
    !   arguments (out) -----------------------------------------------------
    !     line : CHARACTER(len=80)
    !       the 80-character print line that will be filled with the subsection
    !       divider using a positioned dividing line ('---') and the title
    SUBROUTINE line_subsection(text,ilen,line)
        CHARACTER(len=ilen), INTENT(IN) :: text
        CHARACTER(len=80), INTENT(OUT) :: line

        WRITE(line,'(T20,60A)') ("-", i=1,60)
        IF (ilen > 58) THEN
            WRITE(line,'(A,T21," ",A58," ")') line, TRIM(text)
        ELSE
            WRITE(line,'(A,T21," ",A," ")') line, TRIM(text)
        ENDIF

    END SUBROUTINE line_subsection

    ! inserts the timestamp at the beginning of the character array given
    !
    !   arguments (inout) -------------------------------------------------------
    !     line : CHARACTER(len=80)
    !       the 80-character line that receives the timestamp at character 1
    !       in the format [YYYY-MM-DD HH:MM]
    SUBROUTINE line_timestamp(line)
        CHARACTER(len=80), INTENT(INOUT) :: line

        INTEGER, DIMENSION(8) :: date_values

        ! retrieve the date and time and format it for output
        call date_and_time(values=date_values)
        WRITE(line,'(A,T1"[",I4.4,"-",I2.2,"-",I2.2," ",I2.2,":",I2.2,"] ")') &
            line, date_values(1), date_values(2), date_values(3), &
            date_values(5), date_values(6)

    END SUBROUTINE line_timestamp

END MODULE sim_log