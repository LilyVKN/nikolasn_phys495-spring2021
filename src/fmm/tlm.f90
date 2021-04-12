! translates a multipole expansion into a local expansion displaced by a vector
! (r, theta, phi) using the given pre-computed spherical harmonic values. the
! local expansion is valid for evaluations R << r (the particle positions).
! therefore the magnitude of the translation should be sufficiently large to
! avoid divergence.
!
!   arguments (in) ------------------------------------------------------------
!     m_in : COMPLEX, DIMENSION((l_max+1)**2)
!       the multipole expansion to be transformed. filled from l in [0,l_max]
!       and m in [-l,l].
!     l_max : INTEGER
!       the maximum l-value to use for the multipole expansions
!     r : REAL
!       the distance to translate the local expansion, the magnitude of the
!       vector in spherical coordinates
!     yml : COMPLEX, DIMENSION(((2*l_max)+1)**2)
!       the pre-computed spherical harmonic values for the direction of the
!       translation vector, (theta, phi)
!     aml : REAL, DIMENSION(((2*l_max)+1)**2)
!       the pre-computed aml factor values for use in multipole and local
!       expansion translations combined with the 1/(2l+1) factor to reduce
!       repeated calculations. must be to a maximum of 2*l_max to accomodate the
!       factor of (l_p + l)
!
!   arguments (out) -----------------------------------------------------------
!     l_out : COMPLEX, DIMENSION((l_max+1)**2)
!       the local expansion up to the l-value of l_max at the new origin
!
!   notes ---------------------------------------------------------------------
!     since the suite of translation subroutines this belongs to is intended for
!     use on a grid, the values of Yml, and aml are expected to be calculated
!     before the call to this subroutine. that way, the values can be stored and
!     reused. for example, in subdividing grids, there will be the same eight
!     Yml values for moving from parent to child and vice versa. the only value
!     necessary to change then would be the radius.
!
!     the transformation requires values of up to 2*l_max in order to attain
!     tolerable accuracy.
SUBROUTINE tlm(m_in,l_max,r,yml,aml,l_out)
    IMPLICIT NONE

    COMPLEX, DIMENSION((l_max+1)**2), INTENT(IN) :: m_in
    INTEGER, INTENT(IN) :: l_max
    REAL, INTENT(IN) :: r
    COMPLEX, DIMENSION(((2*l_max)+1)**2), INTENT(IN) :: yml
    REAL, DIMENSION(((2*l_max)+1)**2), INTENT(IN) :: aml
    COMPLEX, DIMENSION((l_max+1)**2), INTENT(OUT) :: l_out

    REAL, PARAMETER :: PI_FACTOR = 16.0 * ATAN(1.0)     ! parameter 4*pi

    ! indices for tracking ls, ms, and array positions
    INTEGER :: ind_p, l_p, m_p, ind, l, m, ind_del

    ! factors for storing subcalculations
    REAL :: r_factor, r_lfactor
    COMPLEX :: sum_factor, factor

    r_factor = 1.0

    ! start the outer loop over all the output entries
    ind_p = 1
    DO l_p = 0, l_max

        ! update the factor of r^-l_p-1
        r_factor = r_factor / r

        DO m_p = -l_p, l_p

            ! store the factors dependent on only l_p and m_p with the sign
            l_out(ind_p) = PI_FACTOR * r_factor * aml(ind_p)
            IF(MODULO(l_p+m_p,2).EQ.1) l_out(ind_p) = -l_out(ind_p)

            ! clear the sum for this entry
            sum_factor = 0.0

            ! clear the factor of r^-l
            r_lfactor = 1.0

            ! start the inner loop for the sum over the input values
            DO l = 0, l_max
                ! calculate the indices pointing to the m=0 entry for the
                ! degrees l and (l_p + 1)
                ind = l**2 + l + 1
                ind_del = (l_p + l)**2 + l_p + l + 1

                DO m = -l, l
                    ! perform a check to exclude the zero entries for any m > l
                    IF (ABS(m_p - m).GT.(l_p + l)) CYCLE

                    ! calculate the factor for this entry of the input
                    factor = aml(ind_del + m_p - m) * (2 * (l_p + l) + 1.0)
                    factor = r_lfactor / factor
                    factor = factor * CONJG(yml(ind_del + m_p - m))
                    factor = factor * aml(ind + m)
                    factor = factor * m_in(ind + m)

                    ! apply the sign factor
                    IF (MODULO(m_p-m,2).EQ.1) factor = -factor

                    ! add the value tot the sum
                    sum_factor = sum_factor + factor
                END DO

                ! update the factor of r^-l
                r_lfactor = r_lfactor / r

            END DO

            ! merge the sum with the pre-factor
            l_out(ind_p) = l_out(ind_p) * sum_factor

            ! increment the index to track the current output entry
            ind_p = ind_p + 1
        END DO
    END DO

END SUBROUTINE