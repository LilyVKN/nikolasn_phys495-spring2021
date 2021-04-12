! calculates the array of values used for the multipole- and local-expansion
! translation matrics a_lm as described in [1]. the values are stored in an 
! (l+1)^2 array ordered by increasing m, then by increasing l (i.e. / a_0,0, 
! a_-1,1, a_0,1, a_1,1, a_-2,2, ... , a_l,l /). the values can then be stored
! for reuse in any of the translation matrix calculations. it was chosen to
! include the 1 / (2l + 1) factor here to avoid performing it multiple times
! later and to help avoid transcription errors.
!
!   arguments (in) ------------------------------------------------------------
!     l_max : INTEGER
!       the highest degree of the polynomial to calculate to. it must be
!       non-negative, [0,inf), for providing values to the infinite series of
!       spherical harmonics
!
!   arguments (out) -----------------------------------------------------------
!     result : REAL, DIMENSION((l_max+1)**2)
!       the array of values of a_lm to the maximum degree
!
!   references ----------------------------------------------------------------
!   .. [1] Schmidt, K.E. & Lee, M.A. (1991) Implementing the Fast Multipole
!          Method in Three Dimensions. Journal of Statistical Physics, 63(5/6).
!
!   notes ---------------------------------------------------------------------
!     both the reference [1] and the Pfalzner & Gibbon (2005) book which cites
!     it have a few errors in the translation matrices including the factor of
!     1 / (4*pi) and the pairing of aml with a factor of 1 / (2*l+1)
SUBROUTINE aml(l_max,result)
    REAL, DIMENSION((l_max+1)**2) :: result

    ! precalculate the pi-based prefactor as a parameter
    REAL, PARAMETER :: PI_FACTOR = 1.0 / SQRT(16.0 * ATAN(1.0))

    ! initialize the a_00 value
    result(1) = PI_FACTOR

    ! start the iteration over l
    i = 1
    temp_fact = 1.0
    DO l = 1, l_max
        ! increment the index to point to the next m = 0
        i = i + (2 * l)

        ! contribute to the initial factorial of l!
        temp_fact = temp_fact / l

        ! combine the factorial part with the common prefactor for m = 1
        result(i) =  PI_FACTOR * temp_fact / SQRT((2.0 * l + 1.0))
        IF (MODULO(l,2).EQ.1) result(i) = -result(i)

        ! iterate over the remaining |m| values
        DO m = 1, l
            ! adjust the recursive value to [(l + m)!(l - m)!]^(-1/2)
            result(i + m) = -result(i + m - 1) * SQRT((l - m + 1.0) / (l + m))

            ! mirror the value
            result(i - m) = result(i + m)
        END DO
    END DO
END SUBROUTINE aml