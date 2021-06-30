! calculates the local_expansion locals for a list of points with a common weighting
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
!       origin around which the local_expansion is calculated
!     icount : INTEGER
!       the number of particles contributing to the locals
!     val : REAL
!       the weighting value (e.g. mass or charge), only supports scalar for now
!     l_max : INTEGER
!       the total number of degrees to calculate the locals to
!
!   arguments (out) -----------------------------------------------------------
!     locals : REAL, DIMENSION((l_max+1)**2)
!       the list of local_expansion locals from l = [0,l_max] and m = [-l,l] which
!       amounts to (l_max+1)^2 total elements. they are sorted by rising m then
!       rising l.
!
!   references ----------------------------------------------------------------
!   .. [1] Pfalzner, S., & Gibbon, P. (2005). The Fast local_expansion Method. In
!       Many-body tree methods in physics (pp. 126-148). Cambridge: Cambridge
!       University Press.
SUBROUTINE local_expansion(pos,icount,val,l_max,locals)
    REAL, DIMENSION(icount,3), INTENT(IN) :: pos
    REAL, INTENT(IN) :: val
    COMPLEX, DIMENSION((l_max+1)**2), INTENT(OUT) :: locals

    ! a temporary array for holding the current locals to add to the total
    COMPLEX, DIMENSION((l_max+1)**2) :: curr_moment

    ! clear all the values
    locals = (0.0,0.0)

    ! start iterating over each particle
    DO i = 1, icount
        ! translate the values from Cartesian to spherical coordinates
        r = SQRT(SUM(pos(i,:)**2))      ! r = SQRT(x**2 + y**2 + z**2)
        theta = ACOS(pos(i,3)/r)        ! theta = ACOS(z/r)
        phi = ATAN2(pos(i,2),pos(i,1))   ! phi = ATAN(y/x)

        ! calculate the complex conjugate of the angular component
        call yml_full(l_max,theta,phi,curr_moment)
        curr_moment = CONJG(curr_moment)

        ! apply the radial component of r^-l-1 by cumulative multiplication
        istart = 1
        iend = 1
        r_factor = 1.0 / r
        DO l = 0, l_max
            iend = istart + (2 * l)

            curr_moment(istart:iend) = curr_moment(istart:iend) * &
                (r_factor / (2.0 * l + 1.0))

            r_factor = r_factor / r

            istart = iend + 1
        END DO

        ! add the local_expansion locals from this particle to the total
        locals = locals + curr_moment
    END DO

    ! multiply the weighting value for all the particles
    locals = val * locals
END SUBROUTINE local_expansion