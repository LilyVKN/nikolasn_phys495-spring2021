! calculate an MxM submatrix representing the gravitational forces. the MxM
! submatrix (i,j) stores the gravitational force between particles
! ((i-1)*M, i*M] and ((j-1)*M, j*M]. this function has caching capabilities to
! make use of the force matrix' skew-symmetry. it may be possible that multiple
! nodes try to create a lock file at the same time, but only if the simulation
! is incredibly fast.
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
!     energy : REAL
!       the output of the potential contributed by these particles. Only the
!       first instance of the matrix will add to the total energy since the
!       cached data is stored as force vectors. this also helps prevent double
!       counting particle pairs.
!     istat : INTEGER
!       cache status flag: the output is 0 for a cache miss, 1 indicates a
!       cache-hit
!
!   notes =====================================================================
!     the force interaction matrix is symmetric, which means any submatrix where
!     i == j is also symmetric and the transpose of submatrix (i,j) is equal to
!     submatrix (j,i)
SUBROUTINE force_grav_submat(pos,i,j,M,submat,energy,istat)
    ! load the shared simulation parameters
    use sim_params
    ! arguments
    REAL, DIMENSION(PCOUNT,3), INTENT(IN) :: pos   ! particle positions
    INTEGER, INTENT(IN) :: i                        ! submatrix 1st index
    INTEGER, INTENT(IN) :: j                        ! submatrix 2nd index
    INTEGER, INTENT(IN) :: M                        ! submatrix size (MxM)
    REAL, DIMENSION(M,M,3), INTENT(OUT) :: submat   ! submatrix to return
    REAL, INTENT(OUT) :: energy                     ! potential from the forces
    INTEGER, INTENT(OUT) :: istat
    
    REAL, DIMENSION(M,M,3) :: temp_mat  ! a temporary matrix for reducing memory

    ! cache file management
    LOGICAL :: lexist
    CHARACTER(len=80) :: filename, lockname

    ! initialize the energy to zero
    energy = 0.0
    istat = 0

    ! perform the inverse square avoiding zeros and taking advantage of the
    ! force matrix skew-symmetry (i.e. for submatrices on the diagonal, 
    ! i == j, the submatrix is also skew-symmetric) - no caching is needed
    IF (i == j) THEN
        ! fill the submatrix with the distance vector r_i - r_j
        temp_mat = SPREAD(pos(((i-1)*M)+1:i*M,:),1,M)
        submat = SPREAD(pos(((j-1)*M)+1:j*M,:),2,M)
        submat = temp_mat - submat

        ! iterate over the upper triangle
        DO k = 1,M
            DO l = k + 1,M
                ! calculate the norm and ignore zero values
                r = NORM2(submat(k,l,:))
                IF (r /= 0) THEN
                    ! add this particle's contribution to the potential
                    energy = energy + (1.0 / r)

                    ! apply the softening, inverse square, and normalization
                    r = ((r ** 2) + SOFTENING) * r
                    r = 1.0 / r

                    ! calculate the gravitational force and skew-symmetry
                    submat(k,l,:) = -G * MASS * r * submat(k,l,:)
                    submat(l,k,:) = -submat(k,l,:)
                END IF
            END DO
        END DO

        ! multiply the energy by the appropriate common factor
        energy = - G * (MASS**2) * energy

        RETURN
    ENDIF

    ! create the filename grav_i-j_M.submat for upper triangle submatrices
    WRITE(filename,'("./tmp/fortran/grav_",I2.2,"-",I2.2,"_",I4.4,".submat")') &
        MIN(i,j), MAX(i,j), M

    ! create the cache file lock name grav_i-j_M.lock
    WRITE(lockname,'("./tmp/fortran/grav_",I2.2,"-",I2.2,"_",I4.4,".lock")') &
        MIN(i,j), MAX(i,j), M
    
    ! check if the file exists, so we can make use of the cached data if so
    INQUIRE(file=filename,exist=lexist)
    IF (.not.lexist) THEN

        ! create the files so that the next node knows another node is
        ! performing the calculations. technically, this will stall the other
        ! node until this node is done (a fairly bad bottleneck), but if the
        ! calculations are so fast that the other node catches up to the
        ! skew-symmetric submatrix, the simulation is very small
        OPEN(file=filename,unit=17)     ! create the cache file
        OPEN(file=lockname,unit=18)     ! create the lock

        ! fill the submatrix with the displacement vector r_i - r_j
        temp_mat = SPREAD(pos(((i-1)*M)+1:i*M,:),1,M)
        submat = SPREAD(pos(((j-1)*M)+1:j*M,:),2,M)
        submat = temp_mat - submat
        
        ! iterate over all elements
        DO k = 1,M
            DO l = 1,M
                ! calculate the norm and ignore zero values
                r = NORM2(submat(k,l,:))
                IF (r /= 0) THEN
                    ! add this particle's contribution to the potential
                    energy = energy + (1.0 / r)

                    ! apply the softening, inverse square, and normalization
                    r = ((r**2) + SOFTENING) * r
                    r = 1.0 / r

                    ! calculate the gravitational force for this element
                    submat(k,l,:) = -G * MASS * r * submat(k,l,:)
                END IF
            END DO
        END DO

        ! write the submatrix to the file for submatrix caching
        IF (i > j) THEN                        ! for upper triangle submatrices
            DO k = 1, M                         ! just save the values as normal
                DO l = 1, M
                    WRITE(17,*) submat(k,l,:)
                END DO
            END DO
        ELSE IF (i < j) THEN                    ! for lower triangle submatrices
            DO k = 1, M                         ! save the skew-symmetric values
                DO l = 1, M                     ! to make it the upper triangle
                    WRITE(17,*) -submat(l,k,:)  ! submatrix
                END DO
            END DO
        ENDIF                                   ! ignore diagonal submatrices
        
        CLOSE(17)
        CLOSE(18,status="DELETE",iostat=ios)    ! delete the lock - we're done

        IF (ios /= 0) THEN 
            WRITE(*,'("Error on lock deletion: ",I4)') ios
            WRITE(*,*) char(9), "File: ", lockname
        ENDIF

    ELSE    ! if the cache already exists, read the data
        
        istat = 1   ! set the status flag to indicate this is a cache-hit

        DO  ! keep testing to see if there's a lock on the file before accessing
            INQUIRE(file=lockname,exist=lexist)
            IF (.not.lexist) EXIT
        END DO

        ! read the submatrix data from the file for caching
        OPEN(file=filename,unit=17)

        IF (i > j) THEN                         ! for upper triangle submatrices
            DO k = 1, M                         ! read the values as normal
                DO l = 1, M
                    READ(17,*) submat(k,l,:)
                END DO
            END DO
        ELSE                                    ! for lower triangle submatrices
            DO k = 1, M                         ! read as transpose values for
                DO l = 1, M                     ! implementing skew-symmetry
                    READ(17,*) submat(l,k,:)
                END DO
            END DO
            submat = -submat        ! for skew-symmetric data, make it negative
        ENDIF

        ! close the cache file
        CLOSE(17)
        
    ENDIF

    ! multiply the energy by the appropriate common factor
    energy = - G * (MASS**2) * energy

END SUBROUTINE