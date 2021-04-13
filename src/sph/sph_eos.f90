! this is a place to define a few common equations-of-state i've encountered
MODULE sph_eos
CONTAINS

    ! applies the equation of state for an ideal gas
    !
    !   arguments (in) --------------------------------------------------------
    !     gamma : REAL
    !       the specific heat ratio for this gas
    !     rho : REAL
    !       the gas density at this position
    !     u : REAL
    !       the internal energy of the gas at this position
    !
    !   arguments (out) ----------------------------------------------------
    !     P : REAL
    !       the updated pressure state for this value
    ELEMENTAL SUBROUTINE ideal_eos(gamma,rho,u,P)
        REAL, INTENT(IN) :: gamma
        REAL, INTENT(IN) :: rho, u
        REAL, INTENT(OUT) :: P

        P = (gamma - 1.0) * rho * u

    END SUBROUTINE ideal_eos

    ! applies the equation of state for a polytropic process
    !
    !   arguments (in) --------------------------------------------------------
    !     K : REAL
    !       the entropy-dependent factor of this EOS
    !     GAMMA : REAL
    !       the polytropic exponent for this process
    !     rho : REAL
    !       the gas density at this position
    !
    !   arguments (out) ----------------------------------------------------
    !     P : REAL
    !       the updated pressure state
    ELEMENTAL SUBROUTINE polytropic_eos(K,GAMMA,rho,P)
        REAL, INTENT(IN) :: K, GAMMA
        REAL, INTENT(IN) :: rho
        REAL, INTENT(OUT) :: P

        P = K * (rho**GAMMA)

    END SUBROUTINE polytropic_eos

END MODULE