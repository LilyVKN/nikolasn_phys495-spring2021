! initialize the particles using the default method: even distribution within a
! sphere with tangential velocities
!
!   arguments (inout) =========================================================
!     pos : REAL, DIMENSION(icount,3)
!       the particle position array to be filled
!     vel : REAL, DIMENSION(icount,3)
!       the particle velocity array to be filled
!
!   arguments (in) ============================================================
!     icount : INTEGER
!       the number of particles to be initialized in the array
!     r : REAL
!       the maximum distance the particles can be from the center 
!
! TODO: Create a better velocity initialization function
SUBROUTINE part_init_default(pos,vel,icount,r)
    ! load the shared simulation parameters
    REAL, DIMENSION(icount,3), INTENT(INOUT) :: pos
    REAL, DIMENSION(icount,3), INTENT(INOUT) :: vel

    REAL, PARAMETER :: PI = 4.0 * ATAN(1.0)

    ! initialize the particles
    DO i = 1, icount
        rTheta = 2.0 * PI * RAND()          ! random azimuthal angle [0, 2*pi]
        rPhi = PI * (RAND() - 0.5)          ! random polar angle [-pi/2, +pi/2]
        rRadius = RAND() * r                ! random distance from center

        ! set position in a spherically even distribution
        pos(i,1) = COS(rTheta) * SIN(rPhi)
        pos(i,2) = SIN(rTheta) * SIN(rPhi)
        pos(i,3) = COS(rPhi)
        pos(i,:) = rRadius * pos(i,:)

        ! set velocity tangential to z-axis
        vel(i,1) = -SIN(rTheta)
        vel(i,2) = COS(rTheta)
        IF (rRadius /= 0.0) THEN
            vel(i,1:2) = vel(i,1:2) * (0.0125 / rRadius)
        END IF
        vel(i,3) = 0.0
    END DO
END SUBROUTINE part_init_default