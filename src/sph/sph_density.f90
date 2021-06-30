! calculates the density at the given coordinates, r, as determined by the
! particles in the vicinity
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
!     r : REAL, DIMENSION(3)
!       the coordinates at which to evaluate the density
!
!   arguments (out) -----------------------------------------------------------
!     rho : REAL
!       the density at the designated position
!
!   notes ---------------------------------------------------------------------
!     since the kernel is of a polymorphic type, the argument passed to it
!     should be a pointer or allocatable; otherwise, there will be a segfault
!     when attempting to access its members
SUBROUTINE sph_density(m,pos,icount,kern,r,rho)
    use sph_kernel

    REAL, DIMENSION(icount), INTENT(IN) :: m
    REAL, DIMENSION(icount,3), INTENT(IN) :: pos
    CLASS(kernel), INTENT(IN) :: kern
    REAL, DIMENSION(3), INTENT(IN) :: r
    REAL, INTENT(OUT) :: rho

    rho = 0.0
    DO i = 1, icount
        rho = rho + m(i) * kern%W(r-pos(i,:))
    END DO

END SUBROUTINE sph_density