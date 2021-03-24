! calculates the entire list of angular components of the spherical harmonics
! for all l in [0, l_max] and valid m in [-l, l]. therefore, the values are
! stored in an array of complex numbers of length (l+1)^2. verifying the values 
! against the scipy libraries and against other implementations in this source
! show rather large errors (~5%) for negative m values of large l, but only for
! some angle ranges. this may be due to the use of single-precision values which
! cannot handle the high recursion values of factorials. for this reason values
! of l_max higher than 10 are not recommended.
!
!   arguments (in) ------------------------------------------------------------
!     l_max : INTEGER
!       the highest degree of the polynomial to calculate to. it must be
!       non-negative, [0,inf), for providing values to the infinite series of
!       spherical harmonics
!     theta : REAL
!       the polar angle (colatitude), should be in range [0, pi]
!     phi : REAL
!       the azimuthal angle, should be in range [0, 2*pi]
!
!   arguments (out) -----------------------------------------------------------
!     result : COMPLEX, DIMENSION(l_max**2)
!       the value of the angular component, both real and complex parts
!
!   notes ---------------------------------------------------------------------
!     i'm so sorry about the angular variable naming. the papers i reference are
!     all written by physicists not mathematicians, and it's just easier to keep
!     their notation...
SUBROUTINE yml_full(l_max,theta,phi,result)
    INTEGER, INTENT(IN) :: l_max
    REAL, INTENT(IN) :: theta, phi
    COMPLEX, DIMENSION((l_max+1)**2), INTENT(OUT) :: result

    ! define the parameters, j = sqrt(-1) (i like saving 'i' for indices)
    COMPLEX, PARAMETER :: j = (0,1.0)
    REAL, PARAMETER :: PI_FACTOR = 1.0 / (16.0 * ATAN(1.0))

    ! local variables for intermediate calculations
    REAL, DIMENSION((l_max+1)**2) :: pml
    REAL :: temp, prefactor

    ! for l of 0, the result is simply 1.0
    IF (l == 0) THEN
        RETURN
    ENDIF

    ! apply the normalization factor for the l = 0, m = 0 value
    result(1) = SQRT(PI_FACTOR)

    ! calculate the phi components
    i = 1
    DO l = 1, l_max
        ! increment the index to point to the next m = l
        i = i + (2 * l)

        ! set the common l-dependent prefactor
        prefactor = SQRT((2 * l + 1) * PI_FACTOR)
        result(i) = prefactor
        
        ! reset the temporary value for the m-dependent prefactor, (l-m)!/(l+m)!
        temp = 1.0
        DO m = 1, l
            temp = temp * SQRT((l + m) * (l - m + 1.0))
            
            ! apply all the normalization constants to the results
            result(i - m) = prefactor * temp
            result(i + m) = prefactor * (1.0 / temp)

            ! apply the phi-dependent component
            result(i - m) = result(i - m) * EXP(-j * m * phi)
            result(i + m) = result(i + m) * EXP( j * m * phi)
        END DO
    END DO

    ! calculate the theta component (associated Legendre of cos(theta))
    call pmlcos_full(l_max,theta,pml)
    result = result * pml
    
END SUBROUTINE yml_full