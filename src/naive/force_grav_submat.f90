! calculate an MxM submatrix representing the gravitational forces. the MxM
! submatrix (i,j) stores the gravitational force between particles
! ((i-1)*M, i*M] and ((j-1)*M, j*M]
!
!   arguments (in) ============================================================
!     pos : REAL, DIMENSION(PCOUNT)
!       the 3-D positions of the particles interacting with each other
!     i : INTEGER
!       the 1st index of the submatrix, must be in the range of 1 to PCOUNT / M
!     j : INTEGER
!       the 2nd index of the submatrix, must be in the range of 1 to PCOUNT / M
!     M : INTEGER
!       the size of the submatrix (MxM)
!
!   arguments (out) ===========================================================
!     submat : REAL, DIMENSION(M,M)
!       the submatrix output of size MxM which store the force between the
!       particles indicated by the submatrix indices
!
!   notes =====================================================================
!     the force interaction matrix is symmetric, which means any submatrix where
!     i == j is also symmetric and the transpose of submatrix (i,j) is equal to
!     submatrix (j,i)
SUBROUTINE force_grav_submat(pos,i,j,M,submat)
    ! load the shared simulation parameters
    use SIM_PARAMS
    ! arguments
    REAL, DIMENSION(PCOUNT,3), INTENT(IN) :: pos   ! particle positions
    INTEGER, INTENT(IN) :: i                        ! submatrix 1st index
    INTEGER, INTENT(IN) :: j                        ! submatrix 2nd index
    INTEGER, INTENT(IN) :: M                        ! submatrix size (MxM)
    REAL, DIMENSION(M,M,3), INTENT(OUT) :: submat   ! submatrix to return
    
    REAL, DIMENSION(M,M,3) :: temp_mat  ! a temporary matrix for reducing memory

    ! fill the submatrix with the distance vector r_i - r_j
    temp_mat = SPREAD(pos(((i-1)*M)+1:i*M,:),1,M)
    submat = SPREAD(pos(((j-1)*M)+1:j*M,:),2,M)
    submat = temp_mat - submat
    
    ! perform the inverse square avoiding zeros and taking advantage of the
    ! force matrix skew-symmetry (i.e. for submatrices on the diagonal, i == j,
    ! the submatrix is also skew-symmetric)
    IF (i == j) THEN
        ! iterate over the upper triangle
        DO k = 1,M
            DO l = k + 1,M
                ! calculate the norm and ignore zero values
                r = NORM2(submat(k,l,:))
                IF (r /= 0) THEN
                    ! apply the softening, inverse square, and normalization
                    r = ((r ** 2) + SOFTENING) * r
                    r = 1.0 / r

                    ! calculate the gravitational force and apply skew-symmetry
                    submat(k,l,:) = -G * MASS * r * submat(k,l,:)
                    submat(l,k,:) = -submat(k,l,:)
                END IF
            END DO
        END DO
    ELSE
        ! iterate over all elements
        DO k = 1,M
            DO l = 1,M
                ! calculate the norm and ignore zero values
                r = NORM2(submat(k,l,:))
                IF (r /= 0) THEN
                    ! apply the softening, inverse square, and normalization
                    r = ((r**2) + SOFTENING) * r
                    r = 1.0 / r

                    ! calculate the gravitational force for this element
                    submat(k,l,:) = -G * MASS * r * submat(k,l,:)
                END IF
            END DO
        END DO
    ENDIF
END SUBROUTINE