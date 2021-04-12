! translates a local expansion by a vector (r, theta, phi) using the given
! pre-computed spherical harmonic values. local expansion translation has a
! small region of convergence for small r. anything above that is highly error
! prone due to the exponential factor of r.
!
!   arguments (in) ------------------------------------------------------------
!     l_in : COMPLEX, DIMENSION((l_max+1)**2)
!       the local expansion to be translated. filled from l in [0,l_max] and
!       m in [-l,l]
!     l_max : INTEGER
!       the maximum l-value to use for the local expansions
!     r : REAL
!       the distance to translate the local expansion, the magnitude of the
!       vector in spherical coordinates
!     yml : COMPLEX, DIMENSION((l_max+1)**2)
!       the pre-computed spherical harmonic values for the direction of the
!       translation vector, (theta, phi)
!     aml : REAL, DIMENSION((l_max+1)**2)
!       the pre-computed aml factor values for use in multipole and local
!       expansion translations combined with the 1/(2l+1) factor to reduce
!       repeated calculations
!
!   arguments (out) -----------------------------------------------------------
!     l_out : COMPLEX, DIMENSION((l_max+1)**2)
!       the translated local expansion up to the l-value of l_max
!
!   notes ---------------------------------------------------------------------
!     since the suite of translation subroutines this belongs to is intended for
!     use on a grid, the values of Yml, and aml are expected to be calculated
!     before the call to this subroutine. that way, the values can be stored and
!     reused. for example, in subdividing grids, there will be the same eight
!     Yml values for moving from parent to child and vice versa. the only value
!     necessary to change then would be the radius.
SUBROUTINE tll(l_in,l_max,r,yml,aml,l_out)
    IMPLICIT NONE

    COMPLEX, DIMENSION((l_max+1)**2), INTENT(IN) :: l_in
    INTEGER, INTENT(IN) :: l_max
    REAL, INTENT(IN) :: r
    COMPLEX, DIMENSION((l_max+1)**2), INTENT(IN) :: yml
    REAL, DIMENSION((l_max+1)**2), INTENT(IN) :: aml
    COMPLEX, DIMENSION((l_max+1)**2), INTENT(OUT) :: l_out

    REAL, PARAMETER :: PI_FACTOR = 16.0 * ATAN(1.0)     ! parameter 4*pi
    
    ! indices for tracking ls, ms, and array positions
    INTEGER :: ind_p, l_p, m_p, ind, l, m, ind_del

    ! factors for storing subcalculation
    COMPLEX :: sum_factor, factor, r_factor

    ! start the outer loop (l_p,m_p) over the values of m_out
    ind_p = 1
    DO l_p = 0, l_max
        DO m_p = -l_p, l_p

            ! store the factors dependent on only l_p and m_p first
            l_out(ind_p) = PI_FACTOR * aml(ind_p)

            ! clear the sum for this entry
            sum_factor = 0.0

            ! clear the r^(l-l_p) factor for the new inner loop
            r_factor = 1.0

            ! start the inner loop for the sum over the input values
            DO l = l_p, l_max
                ! calculate the indices pointing to the m=0 entry for the
                ! degrees l and (l - l_p)
                ind = l**2 + l + 1
                ind_del = (l - l_p)**2 + l - l_p + 1

                DO m = -l, l
                    ! perform a check to exclude the zero entries for any |m|>l
                    IF (ABS(m - m_p).GT.(l - l_p)) CYCLE

                    ! calculate the factor for this entry of the input
                    factor = aml(ind + m)
                    factor = r_factor / factor
                    factor = factor * yml(ind_del + m - m_p)
                    factor = factor * aml(ind_del + m - m_p)
                    factor = factor * l_in(ind + m)

                    ! add the value to the sum
                    sum_factor = sum_factor + factor
                END DO

                ! update the factor of r^(l-l_p)
                r_factor = r_factor * r
            END DO

            ! merge the sum with the pre-factor
            l_out(ind_p) = l_out(ind_p) * sum_factor

            ! increment the index to track the current output entry
            ind_p = ind_p + 1

        END DO
    END DO

END SUBROUTINE