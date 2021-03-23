! calculates the entire list of P_ml(cos(theta)) values and their derivatives
! for l in [0, l_max] and all valid m in [-l, l] using a combination of
! recursions provided by [1] and [2]. the values are stored in (l+1)^2 size
! arrays ordered by increasing m, then by increasing l (i.e. / P_0,0, P_-1,1, 
! P_0,1, P_1,1, P_-2,2, ... , P_l,l /). this algorithm was chosen as the most
! efficient recursion as described in [1]. it is meant only for calculating
! values for spherical harmonics and thus has stricter rules than a general
! solution. the derivative recursion was implemented based on [2]. verifying the
! values against the scipy libraries show slight errors likely due to the usage
! of single-precision floating points
!
!   arguments (in) ------------------------------------------------------------
!     l_max : INTEGER
!       the highest degree of the polynomial to calculate to. it must be
!       non-negative, [0,inf), for providing values to the infinite series of
!       spherical harmonics
!     theta : REAL
!       the value of theta to be used for the function
!
!   arguments (out) -----------------------------------------------------------
!     result : REAL, DIMENSION((l_max+1)**2)
!       the value output by the all polynomial up to order l_max
!     deriv : REAL, DIMENSION((l_max+1)**2)
!       the value of the derivative at all polynomials up to order l_max
!
!   references ----------------------------------------------------------------
!   .. [1] Wiggins, R., & Saito, M. (1971). Evaluation of Computational
!       Algorithms for the Associated Legendre Polynomiasl By Interval Analysis.
!       Bulletin of the Seismological Society of America, 61(2), pp.375-381.
!   .. [2] Bosch, W. (2000) On the Computation of Derivatives of Legendre
!       Functions. Physics and Chemistry of the Earth (A), 25(9-11), pp.655-659.
SUBROUTINE dpmlcos_full(l_max,theta,result,deriv)
    ! argument variables
    INTEGER, INTENT(IN) :: l_max
    REAL, INTENT(IN) :: theta
    REAL, DIMENSION((l_max+1)**2), INTENT(OUT) :: result, deriv

    ! variables to store reused trigonometric calculations
    REAL :: sin_val, cos_val, cot_val, csc_val

    ! initialize result to zero which becomes useful for recursion
    result = 0.0

    ! add the initial values of P and exit if these satisfy the arguments
    result(1) = 1.0             ! P_0,0(cos(theta))
    IF (l == 0) RETURN

    ! calculate the basic trig values needed for l >= 1
    sin_val = SIN(THETA)
    cos_val = COS(THETA)

    result(2) = 0.5 * sin_val   ! P_(-1),1(cos(theta))
    result(3) = cos_val         ! P_0,1(cos(theta))
    result(4) = -sin_val        ! P_1,1(cos(theta))
    IF (l_max == 1) RETURN

    ! calculate some convenient, reused trig values to reduce calculations
    csc_val = 1.0 / sin_val
    cot_val = cos_val / sin_val

    ! set the start index for the current value of l
    istart_ind = 5
    DO l = 2, l_max
        ! set an index for the value of m = l
        ind = istart_ind + (2 * l)

        ! calculate P_ll(cos(theta)), the initial value for the current l
        result(ind) = -result(istart_ind - 1) * sin_val * (2 * l - 1)

        ! define the first derivative value before recursion
        deriv(ind) = l * result(ind) * cot_val

        ! set the next start index to be after this l value
        istart_ind = ind + 1

        ! start the recursion of P_m-1,n based on P_m,n and P_m+1,n
        DO m = l, -l + 1, -1
            result(ind - 1) = (2 * m * cot_val * result(ind)) + result(ind + 1)
            result(ind - 1) = -result(ind - 1) / ((l + m) * (l - m + 1))

            ! use the values of P_m,n and P_m+1,n to get the derivative 
            deriv(ind - 1) = (m - 1) * cot_val * result(ind - 1)
            deriv(ind - 1) = deriv(ind - 1) - result(ind)

            ! move to the next value
            ind = ind - 1
        END DO

    END DO

END SUBROUTINE dpmlcos_full