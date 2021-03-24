! calculates the multipole moments for a list of points with a common weighting
! values (val can represent mass, charge, etc.) up to a given maximum degree. 
! the values are stored in (l+1)^2 size arrays ordered by increasing m, then by
! increasing l (i.e. / M_0,0, M_-1,1, M_0,1, M_1,1, M_-2,2, ... , M_l,l /).
! because the high l, m < 0 values are increasingly larger, the values can
! overflow fairly easily. the mono-, di-, and quadrupoles were verified using
! the derived equations in Cartesian coordinates.
!
!   arguments (in) ------------------------------------------------------------
!     pos : REAL, DIMENSION(icount,3)
!       the 3-dimensional Cartesian positions of each particle relative to the
!       origin around which the multipole is calculated
!     icount : INTEGER
!       the number of particles contributing to the moments
!     val : REAL
!       the weighting value (e.g. mass or charge), only supports scalar for now
!     l_max : INTEGER
!       the total number of degrees to calculate the moments to
!
!   arguments (out) -----------------------------------------------------------
!     moments : REAL, DIMENSION((l_max+1)**2)
!       the list of multipole moments from l = [0,l_max] and m = [-l,l] which
!       amounts to (l_max+1)^2 total elements. they are sorted by rising m then
!       rising l.
!
!   references ----------------------------------------------------------------
!   .. [1] Pfalzner, S., & Gibbon, P. (2005). The Fast Multipole Method. In
!       Many-body tree methods in physics (pp. 126-148). Cambridge: Cambridge
!       University Press.
SUBROUTINE multipole(pos,icount,val,l_max,moments)
    REAL, DIMENSION(icount,3), INTENT(IN) :: pos
    REAL, INTENT(IN) :: val
    COMPLEX, DIMENSION((l_max+1)**2), INTENT(OUT) :: moments

    ! a temporary array for holding the current moments to add to the total
    COMPLEX, DIMENSION((l_max+1)**2) :: curr_moment

    ! clear all the values
    moments = (0.0,0.0)

    ! start iterating over each particle
    DO i = 1, icount
        ! translate the values from Cartesian to spherical coordinates
        r = SQRT(SUM(pos(i,:)**2))      ! r = SQRT(x**2 + y**2 + z**2)
        theta = ACOS(pos(i,3)/r)        ! theta = ACOS(z/r)
        phi = ATAN(pos(i,2)/pos(i,1))   ! phi = ATAN(y/x)

        ! calculate the angular component
        call yml_full(l_max,theta,phi,curr_moment)

        ! apply the radial component of r^l by cumulative multiplication over l
        istart_ind = 2
        DO l = 1, l_max
            curr_moment(istart_ind:) = curr_moment(istart_ind:) * r
            istart_ind = istart_ind + (2 * l + 1)
        END DO

        ! add the multipole moments from this particle to the total
        moments = moments + curr_moment
    END DO

    ! multiply the weighting value for all the particles
    moments = val * moments

    ! apply the complex conjugate of the Yml (only complex contribution to Mml)
    moments = CONJG(moments)
END SUBROUTINE multipole