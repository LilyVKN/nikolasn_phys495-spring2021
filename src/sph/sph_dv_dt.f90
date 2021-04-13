! calculates the momentum contribution of the internal pressure among particles
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
!     icount : INTEGER
!       the number of particles
!     kern : CLASS(kernel)
!       the SPH kernel for calculating weights (see sph_kernel.f90 for details)
!
!   arguments (out) -----------------------------------------------------------
!     dv_dt : REAL, DIMENSION(icount)
!       the acceleration due to internal pressure for each particle
!
!   notes ---------------------------------------------------------------------
!     since the kernel is of a polymorphic type, the argument passed to it
!     should be a pointer or allocatable; otherwise, there will be a segfault
!     when attempting to access its members
SUBROUTINE sph_dv_dt(m,rho,P,pos,icount,kern,dv_dt)
    use sph_kernel

    REAL, DIMENSION(icount), INTENT(IN) :: m, rho, P
    REAL, DIMENSION(icount,3), INTENT(IN) :: pos
    CLASS(kernel), INTENT(IN) :: kern
    REAL, DIMENSION(icount,3), INTENT(OUT) :: dv_dt

    REAL, DIMENSION(3) :: temp

    dv_dt = 0.0
    DO i = 1, icount
        DO j = i+1, icount
            temp = kern%dW(pos(i,:) - pos(j,:))
            temp = (P(i) / (rho(i)**2)) + (P(j) / (rho(j)**2)) * temp

            ! mirror the additional term along pairwise symmetry
            dv_dt(i,:) = dv_dt(i,:) - (m(j) * temp)
            dv_dt(j,:) = dv_dt(j,:) + (m(i) * temp)
        END DO
    END DO

END SUBROUTINE sph_dv_dt