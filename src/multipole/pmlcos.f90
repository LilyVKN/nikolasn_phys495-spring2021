! calculates the associated Legendre polynomial P_ml(cos(theta)) using a
! recursive formula. this algorithm was chosen as the most efficient as
! described in [1]. it is meant only for calculating values for spherical
! harmonics and thus has stricter rules than a general solution. verifying the
! values against the scipy libraries show slight errors likely due to the usage
! of single-precision floating points
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
!
!   references ----------------------------------------------------------------
!   .. [1] Wiggins, R., & Saito, M. (1971). Evaluation of Computational
!       Algorithms for the Associated Legendre Polynomiasl By Interval Analysis.
!       Bulletin of the Seismological Society of America, 61(2), pp.375-381.
SUBROUTINE pmlcos(m,l,theta,result)
    ! argument variables
    INTEGER, INTENT(IN) :: m, l
    REAL, INTENT(IN) :: theta
    REAL, INTENT(OUT) :: result

    ! variables to store intermediate calculations
    REAL :: sin_val, cos_val, cot_val, P_next, P_curr, P_last

    ! ensure the order and degree are valid (l in [0,inf); m in [-l,l])
    IF (l < 0) THEN
        result = 0.0
        RETURN
    ELSE IF (ABS(m) > l) THEN
        result = 0.0
        RETURN
    ENDIF

    ! pre-check easy values of P
    IF (l == 0) THEN                ! for l = 0, P is 1.0
        result = 1.0
    ELSE IF (l == 1) THEN           ! for l = 1, the values are basic trig
        IF (m == 0) THEN
            result = COS(theta)
        ELSE IF (m == 1) THEN
            result = -SIN(theta)
        ELSE
            result = 0.5 * SIN(theta)
        ENDIF
    ELSE                            ! for higher l, user recursion
        ! calculate some convenient reused trig values
        sin_val = SIN(theta)
        cos_val = COS(theta)
        cot_val = cos_val / sin_val

        ! calculate P_ll(cos(theta)), the initial value
        P_curr = ((0.5 * sin_val)**l)
        DO i = l + 1, (2 * l)
            P_curr = P_curr * i
        END DO

        ! apply the standard sign convention for the P_ll value
        IF (MODULO(l,2).EQ.1) P_curr = -P_curr
        
        ! start the recursion of P_m-1,n based on P_m,n and P_m+1,n
        P_next = 0.0
        P_last = 0.0
        DO m_curr = l, m + 1, -1
            P_next = (2 * m_curr * cot_val * P_curr) + P_last
            P_next = -P_next / ((l + m_curr) * (l - m_curr + 1))

            ! move the recursion values forward
            P_last = P_curr
            P_curr = P_next
        END DO

        ! set the result to the final recursion value
        result = P_curr
    ENDIF
END SUBROUTINE pmlcos