! calculates the basic energy equation for updating the internal energy of each
! particle.
!
!   arguments (in) ------------------------------------------------------------
!     m : REAL, DIMENSION(icount)
!       the array of masses for each particle
!     rho : REAL, DIMENSION(icount)
!       the array of densities at each particle's position
!     P : REAL, DIMENSION(icount)
!       the array of pressure values at each particle's position
!     pos : REAL, DIMENSION(icount,3)
!       the array of particle positions
!     vel : REAL, DIMENSION(icount,3)
!       the array of particle velocities
!     icount : INTEGER
!       the number of particles
!     kern : CLASS(kernel)
!       the SPH kernel for calculating weights (see sph_kernel.f90 for details)
!
!   arguments (out) -----------------------------------------------------------
!     du_dt : REAL, DIMENSION(icount)
!       the change in internal energy for each particle
!
!   notes ---------------------------------------------------------------------
!     since the kernel is of a polymorphic type, the argument passed to it
!     should be a pointer or allocatable; otherwise, there will be a segfault
!     when attempting to access its members
SUBROUTINE sph_du_dt(m,rho,P,pos,vel,icount,kern,du_dt)
    use sph_kernel

    REAL, DIMENSION(icount), INTENT(IN) :: m, rho, P
    REAL, DIMENSION(icount,3), INTENT(IN) :: pos, vel
    CLASS(kernel), INTENT(IN) :: kern
    REAL, DIMENSION(icount), INTENT(OUT) :: du_dt

    REAL :: temp

    du_dt = 0.0
    DO i = 1, icount
        DO j = i+1, icount
            temp = dot_product(vel(i,:)-vel(j,:),kern%dW(pos(i,:)-pos(j,:)))

            ! mirror the additional term along pairwise symmetry
            du_dt(i) = du_dt(i) + (m(j) * temp * P(i) / (rho(i)**2))
            du_dt(j) = du_dt(j) + (m(i) * temp * P(j) / (rho(j)**2))
        END DO
    END DO
END SUBROUTINE sph_du_dt