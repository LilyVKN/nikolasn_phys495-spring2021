! calculates the mono-, di-, and quadrupole moments of the set of particles
! given using Cartesian values of the multipole expansion. this allows for a
! faster, cleaner, and more accurate calculation of the expansions when limited
! multipole precision is needed.
!
!   arguments (in) ------------------------------------------------------------
!     pos : REAL, DIMENSION(icount,3)
!       the Cartesian position vector of the particles relative to the point
!       around which the moments are calculated
!     icount : INTEGER
!       the number of particles in the array to contribute to the moments
!     vals : REAL, DIMENSION(icount)
!       the list of weights for the particles. these can represent masses,
!       charges, etc.
!
!   arguments (out) -----------------------------------------------------------
!     mono : REAL
!       the monopole moment of the particles, a sum over all their values
!     di : REAL, DIMENSION(3)
!       the dipole moments which is represented by a vector, a sum over the
!       displacement vectors weighted by each particle's value
!     quad : REAL, DIMENSION(3,3)
!       the quadrupole moments which is represented by a rank-2 tensor (here as
!       a 3x3 matrix)
!     moments : COMPLEX, DIMENSION(9)
!       the spherical-form multipole moments, q_ml
!
!   references ----------------------------------------------------------------
!   .. [1] Pfalzner, S., & Gibbon, P. (2005). Appendix 2 Spherical Harmonics. In
!       Many-body tree methods in physics (pp. 152-154). Cambridge: Cambridge
!       University Press.
SUBROUTINE quadrupole(pos,icount,vals,mono,di,quad,moments)
    REAL, DIMENSION(icount,3), INTENT(IN) :: pos    ! particle positions
    REAL, DIMENSION(icount), INTENT(IN) :: vals     ! particle values
    REAL, INTENT(OUT) :: mono                       ! monopole value
    REAL, DIMENSION(3), INTENT(OUT) :: di           ! dipole values
    REAL, DIMENSION(3,3), INTENT(OUT) :: quad       ! quadrupole values
    COMPLEX, DIMENSION(9), INTENT(OUT) :: moments   ! multipole moments

    REAL, PARAMETER :: PI_FACTOR = SQRT(1.0 / (16.0 * ATAN(1.0)))

    ! the mono- and dipoles are trivial sum calculations
    mono = SUM(vals)
    di = SUM(SPREAD(vals,2,3) * pos,1)

    ! iterate over the matrix elements
    DO i = 1, 3
        DO j = 1, 3
            quad(i,j) = SUM(3.0 * vals * pos(:,i) * pos(:,j),1)
            IF (i.EQ.j) quad(i,j) = quad(i,j) - SUM(vals * SUM(pos**2,2),1)
        END DO
    END DO

    moments(1) = PI_FACTOR * mono

    moments(2) = SQRT(1.5) * PI_FACTOR * COMPLEX(di(1),di(2))
    moments(3) = SQRT(3.0) * PI_FACTOR * di(3)
    moments(4) = -SQRT(1.5) * PI_FACTOR * COMPLEX(di(1),-di(2))

    moments(5) = COMPLEX(quad(1,1) - quad(2,2),2.0 * quad(1,2))
    moments(5) = SQRT(5.0/24.0) * PI_FACTOR * moments(5)
    moments(6) = SQRT(5.0/6.0) * PI_FACTOR * COMPLEX(quad(1,3),quad(2,3))
    moments(7) = SQRT(1.25) * PI_FACTOR * quad(3,3)
    moments(8) = -SQRT(5.0/6.0) * PI_FACTOR * COMPLEX(quad(1,3),-quad(2,3))
    moments(9) = COMPLEX(quad(1,1) - quad(2,2),-2.0 * quad(1,2))
    moments(9) = SQRT(5.0/24.0) * PI_FACTOR * moments(9)

END SUBROUTINE quadrupole