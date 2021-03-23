! calculates the angular component of the spherical harmonic for the given
! order and degree and the given direction. verifying the values against the 
! scipy libraries show slight errors likely due to the usage of single-precision
! floating points.
!
!   arguments (in) ------------------------------------------------------------
!     m : INTEGER
!       the order of the angular component, must be in range [-l,l]
!     l : INTEGER
!       the degree of the angular component, must be in the range [0,inf)
!     theta : REAL
!       the polar angle (colatitude), should be in range [0, pi]
!     phi : REAL
!       the azimuthal angle, should be in range [0, 2*pi]
!
!   arguments (out) -----------------------------------------------------------
!     result : COMPLEX
!       the value of the angular component, both real and complex parts
!
!   notes ---------------------------------------------------------------------
!     i'm so sorry about the angular variable naming. the papers i reference are
!     all written by physicists not mathematicians, and it's just easier to keep
!     their notation...
SUBROUTINE yml(m,l,theta,phi,result)
    INTEGER, INTENT(IN) :: m, l
    REAL, INTENT(IN) :: theta, phi
    COMPLEX, INTENT(OUT) :: result

    ! define the parameters, j = sqrt(-1) (i like saving 'i' for indices)
    COMPLEX, PARAMETER :: j = (0,1.0)
    REAL, PARAMETER :: PI_FACTOR = 1.0 / (16.0 * ATAN(1.0))

    ! local variables for intermediate calculations
    REAL :: pml
    REAL :: factorial    ! the factorial can easily exceed the REAL(4) max

    ! for l of 0, the result is simply 1.0
    IF (l == 0) THEN
        result = SQRT(REAL(PI_FACTOR))
        RETURN
    ENDIF

    ! calculate the phi component
    result = EXP(j * m * phi)

    ! calculate the theta component (associated Legendre of cos(theta))
    call pmlcos(m,l,theta,pml)
    result = result * pml

    ! calculate the normalization prefactor (l - m)! / (l + m)!
    factorial = 1.0
    IF (m == 0) THEN
        factorial = 1.0
    ELSE IF (m > 0) THEN
        DO i = l - m + 1, l + m         ! for positive m, the divisor dominates
            factorial = factorial * SQRT(REAL(i))
        END DO
        factorial = 1.0 / factorial
    ELSE
        DO i = l + m + 1, l - m         ! for negative m, the dividend dominates
            factorial = factorial * SQRT(REAL(i))
        END DO
    ENDIF

    ! apply the total normalization [1/(4PI) here is equivalent to ATAN(1.0)]
    result = SQRT((2 * l + 1) * PI_FACTOR) * factorial * result
    
END SUBROUTINE yml