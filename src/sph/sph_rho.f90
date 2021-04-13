
! calculates the density at each particle position as determined by the entire
! system
!
!   arguments (in) ------------------------------------------------------------
!     m : REAL, DIMENSION(icount)
!       the array of masses for the particles
!     pos : REAL, DIMENSION(icount)
!       the array of positions for the particles
!     icount : INTEGER
!       the number of particles to check
!     kern : CLASS(kernel)
!       the SPH kernel for calculating weights (see sph_kernel.f90 for details)
!
!   arguments (out) -----------------------------------------------------------
!     rho : REAL
!       the density at the designated position
!
!   notes ---------------------------------------------------------------------
!     since the kernel is of a polymorphic type, the argument passed to it
!     should be a pointer or allocatable; otherwise, there will be a segfault
!     when attempting to access its members
SUBROUTINE sph_rho(m,pos,icount,kern,rho)
    use sph_kernel

    REAL, DIMENSION(icount), INTENT(IN) :: m
    REAL, DIMENSION(icount,3), INTENT(IN) :: pos
    CLASS(kernel), INTENT(IN) :: kern
    REAL, DIMENSION(icount), INTENT(OUT) :: rho

    REAL :: temp

    rho = 0.0
    DO i = 1, icount
        DO j = i, icount
            temp = kern%W(pos(j,:)-pos(i,:))

            ! mirror the mass contributions along pairwise symmetry
            rho(i) = rho(i) + m(i) * temp
            rho(j) = rho(j) + m(j) * temp
        END DO
    END DO

END SUBROUTINE sph_rho