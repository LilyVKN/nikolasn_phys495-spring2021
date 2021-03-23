! calculates the associated Legendre polynomial P_ml(cos(theta)) using a
! recursive formula with the derivative also output. this algorithm was chosen
! as the most efficient recursion as described in [1]. it is meant only for
! calculating values for spherical harmonics and thus has stricter rules than a
! general solution. the derivative recursion was implemented based on [2]. 
! verifying the values against the scipy libraries show slight errors likely due
! to the usage of single-precision floating points
!
!   arguments (in) ------------------------------------------------------------
!     m : INTEGER
!       the order of the polynomial. it must be within the range [-l,l]
!     l : INTEGER
!       the degree of the polynomial. it must be non-negative, [0,inf), for
!       providing values to the infinite series of spherical harmonics
!     theta : REAL
!       the value of theta to be used for the function
!
!   arguments (out) -----------------------------------------------------------
!     result : REAL
!       the value output by the requested polynomial order and degree
!     deriv : REAL
!       the value of the derivative at the requested polynomial order and degree
!
!   references ----------------------------------------------------------------
!   .. [1] Wiggins, R., & Saito, M. (1971). Evaluation of Computational
!       Algorithms for the Associated Legendre Polynomiasl By Interval Analysis.
!       Bulletin of the Seismological Society of America, 61(2), pp.375-381.
!   .. [2] Bosch, W. (2000) On the Computation of Derivatives of Legendre
!       Functions. Physics and Chemistry of the Earth (A), 25(9-11), pp.655-659.
SUBROUTINE dpmlcos(m,l,theta,result,deriv)
    ! argument variables
    INTEGER, INTENT(IN) :: m, l
    REAL, DIMENSION(:), INTENT(IN) :: theta
    REAL, DIMENSION(SIZE(theta)), INTENT(OUT) :: result, deriv

    ! variables to store intermediate calculations
    REAL, DIMENSION(SIZE(theta)) :: sin_val, csc_val, cos_val, cot_val, &
                                    P_next, P_curr, P_last

    ! ensure the order and degree are valid (l in [0,inf); m in [-l,l])
    IF (l < 0) THEN
        result = 0.0
        RETURN
    ELSE IF (ABS(m) > l) THEN
        result = 0.0
        RETURN
    ENDIF

    ! cover the trivial cases of l = 0, 1
    IF (l == 0) THEN
        result = 1.0
        deriv = 0.0
    ELSE IF (l == 1) THEN
        IF (m == 0) THEN
            result = COS(theta)
            deriv = -SIN(theta)
        ELSE IF (m == 1) THEN
            result = -SIN(theta)
            deriv = -COS(theta)
        ELSE
            result = 0.5 * SIN(theta)
            deriv = 0.5 * COS(theta)
        ENDIF        
    ELSE
        ! determine the trig values of the given input
        sin_val = SIN(THETA)
        cos_val = COS(THETA)
        csc_val = 1.0 / sin_val         ! store csc and cot to reduce divisions
        cot_val = cos_val / sin_val

        ! calculate P_ll(cos(theta)), the initial value
        P_curr = ((0.5 * sin_val)**l)
        DO i = l + 1, (2 * l)           ! multiply by the factor (2l)! / l!
            P_curr = P_curr * i
        END DO
        
        ! start the recursion of P_m-1,n based on P_m,n and P_m+1,n
        P_next = 0.0
        P_last = 0.0
        DO m_curr = l, m + 1, -1
            P_next = (2 * m_curr * cot_val * P_curr) - P_last
            P_next = P_next / ((l + m_curr) * (l - m_curr + 1))

            ! move the values forward
            P_last = P_curr
            P_curr = P_next

            ! update the current derivative using the current Pnm values
            deriv = (m_curr + 1) * cos_val * P_curr
            deriv = deriv - (m_curr + l + 1) * P_last
            deriv = deriv * csc_val
        END DO

        ! set the result to the final iteration
        result = P_curr

    ENDIF
END SUBROUTINE dpmlcos