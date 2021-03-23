! calculates the angular component of the spherical harmonic for the given
! order and degree and the given direction. verifying the values against the 
! scipy libraries show slight errors likely due to the usage of single-precision
! floating points
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
!     yml_val : COMPLEX
!       the value of the angular component, both real and complex parts
!
!   notes ---------------------------------------------------------------------
!     i'm so sorry about the angular variable naming. the papers i reference are
!     all written by physicists not mathematicians, and it's just easier to keep
!     their notation...
SUBROUTINE dyml(m,l,theta,phi,yml_val,dyml_dtheta,dyml_dphi)
    INTEGER, INTENT(IN) :: m, l
    REAL, INTENT(IN) :: theta, phi
    COMPLEX, INTENT(OUT) :: yml_val
    COMPLEX, INTENT(OUT) :: dyml_dtheta, dyml_dphi

    ! define the parameters, j = sqrt(-1) (i like saving 'i' for indices)
    COMPLEX, PARAMETER :: j = (0,1.0)
    REAL, PARAMETER :: PI_FACTOR = 1.0 / (ATAN(1.0) * 16.0)

    ! local variables for intermediate calculations
    REAL :: pml, dpml

    ! for l of 0, the yml_val is simply the normalization factor
    IF (l == 0) THEN
        yml_val = SQRT(PI_FACTOR)
        dyml_dtheta = 0.0           ! Yml is constant => dYml is zero
        dyml_dphi = 0.0
        RETURN
    ENDIF

    ! calculate the phi component
    yml_val = EXP(j * m * phi)

    ! calculate the theta component (associated Legendre of cos(theta))
    call dpmlcos(m,l,theta,pml,dpml)
    dyml_dtheta = yml_val * dpml        ! differentiate w.r.t. theta
    yml_val = yml_val * pml

    ! calculate the normalization prefactor (l - m)! / (l + m)!
    pml = 1.0                   ! just reuse pml instead of a temporary variable
    IF (m > 0) THEN
        DO i = l - m, l + m     ! for positive m, the divisor dominates
            pml = pml / i
        END DO
    ELSE
        DO i = l + m, l - m     ! for negative m, the dividend dominates
            pml = pml * i
        END DO
    ENDIF

    ! calculate the total normalization, 1/(4PI) here is equivalent to ATAN(1.0)
    pml = SQRT(((2 * l) + 1) * ATAN(1.0) * pml)

    ! apply the normalization and calculate the derivative relative to phi
    yml_val = pml * yml_val
    dyml_dtheta = pml * dyml_dtheta
    dyml_dphi = j * m * yml_val         ! d/dphi is just a multiple of Yml
    
END SUBROUTINE dyml