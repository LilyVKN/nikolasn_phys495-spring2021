PROGRAM main
    use sph_eos
    use sph_kernel
    
    INTEGER, PARAMETER :: PCOUNT = 1000, TIMESTEPS = 240
    REAL, PARAMETER :: gamma = 5.0 / 3.0

    REAL, DIMENSION(PCOUNT,3) :: pos, vel, acc
    REAL, DIMENSION(PCOUNT,PCOUNT,3) :: submat
    REAL, DIMENSION(PCOUNT) :: m, rho, u, P, du_dt

    CLASS(cubic_spline_kernel), ALLOCATABLE :: kern

    ! file management variables
    INTEGER, DIMENSION(8) :: date_values
    CHARACTER(len=15) :: timestamp
    CHARACTER(len=27) :: dir_name
    CHARACTER(len=80) :: filename

    REAL :: e_i, e_f, energy
    INTEGER :: istat

    ALLOCATE(kern)

    call date_and_time(values=date_values)
    WRITE(timestamp,'(I2.2,"-",I2.2,"-",I4.4,"_",I4.4)') date_values(3), &
        date_values(2), date_values(1), &
        ((date_values(5) * 100) + date_values(6))
    dir_name = 'SPH_fortran_' // timestamp
    call execute_command_line('mkdir -p ' // dir_name)

    ! create a directory for temporary files
    call execute_command_line('mkdir -p ./tmp/fortran')

    m = 1.0
    call part_init_default(pos,vel,PCOUNT,1.0)
    vel = vel * 0.01
    u = SUM(vel**2,2)
    e_i = SUM(u)

    call sph_rho(m,pos,PCOUNT,kern,rho)
    call ideal_eos(gamma,rho,u,P)

    DO i = 1, TIMESTEPS
        call sph_du_dt(m,rho,P,pos,vel,icount,kern,du_dt)
        u = u + (du_dt / 24.0)

        call sph_dv_dt(m,rho,P,pos,PCOUNT,kern,acc)

        call force_grav_submat(pos,1,1,PCOUNT,submat,energy,istat)
        acc = acc + SUM(submat,1)

        ! kick-drift-kick
        pos = pos + (vel / 48.0)

        vel = vel + (acc / 24.0)

        pos = pos + (vel / 48.0)

        call sph_rho(m,pos,PCOUNT,kern,rho)

        call ideal_eos(gamma,rho,u,P)

        WRITE(filename,'(A,"/SPH_fortran_",I4.4,".fcache")') dir_name, i
        OPEN(file=filename,unit=16)
        DO j = 1, PCOUNT
            WRITE(16,*) pos(j,:), vel(j,:), m(j), rho(j)
        END DO
        CLOSE(16)
    END DO
    e_f = SUM(u)

END PROGRAM main