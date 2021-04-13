! this module stores the SPH kernel framework for determining the weights and
! weight gradients and applying them across the entire SPH system cleanly
!
!   notes ---------------------------------------------------------------------
!     the kernel type is passed as a polymorphic class to the relevant SPH
!     procedures. therefore, any instantiations of kernel classes including new,
!     user-derived classes should be declared as allocatable or pointers to
!     allow the functions to access their member variables and functions
MODULE sph_kernel
    
    ! this acts as a wrapper for the weight and weight gradient functions using
    ! the given smoothing factor, h.
    TYPE, PUBLIC, ABSTRACT :: kernel
        REAL :: h = 1.0
    CONTAINS
        PROCEDURE(func_kernel), DEFERRED :: W
        PROCEDURE(func_del_kernel), DEFERRED :: dW
    END TYPE kernel

    ! the abstract interface for defining weight and weight gradient functions
    ABSTRACT INTERFACE
        FUNCTION func_kernel(this,pos) RESULT(W)        ! weight function
            import kernel
            CLASS(kernel), INTENT(IN) :: this
            REAL, DIMENSION(3), INTENT(IN) :: pos
            REAL :: W
        END FUNCTION func_kernel

        FUNCTION func_del_kernel(this,pos) RESULT(dW)   ! gradient func
            import kernel
            CLASS(kernel), INTENT(IN) :: this
            REAL, DIMENSION(3), INTENT(IN) :: pos
            REAL, DIMENSION(3) :: dW
        END FUNCTION func_del_kernel
    END INTERFACE

    ! a default, cubic spline SPH kernel. to use, you must declare any
    ! instantiantion of these classes as pointers or allocatables
    TYPE, PUBLIC, EXTENDS(kernel) :: cubic_spline_kernel
    CONTAINS
        PROCEDURE :: W => func_cubic_spline_kernel
        PROCEDURE :: dW => func_del_cubic_spline_kernel
    END TYPE

CONTAINS

    ! the common cubic spline kernel
    FUNCTION func_cubic_spline_kernel(this,pos) RESULT(W)
        CLASS(cubic_spline_kernel), INTENT(IN) :: this
        REAL, DIMENSION(3), INTENT(IN) :: pos
        
        REAL :: W

        REAL, PARAMETER :: PI_FACTOR = 1.0 / (4.0 * ATAN(1.0))

        q = SQRT(SUM(pos**2)) / this%h
        IF (q.LE.1.0) THEN
            w = 1.0 - 1.5 * (q**2) + 0.75 * (q**3)
            w = PI_FACTOR * w / (this%h**3)
        ELSE IF (q.LE.2.0) THEN
            w = 0.25 * ((2.0 - q)**3)
            w = PI_FACTOR * w / (this%h**3)
        ELSE
            w = 0.0
        ENDIF
    END FUNCTION func_cubic_spline_kernel

    ! the gradient of the common cubic spline kernel
    FUNCTION func_del_cubic_spline_kernel(this,pos) RESULT(dW)
        CLASS(cubic_spline_kernel), INTENT(IN) :: this
        REAL, DIMENSION(3), INTENT(IN) :: pos

        REAL, DIMENSION(3) :: dW

        REAL, PARAMETER :: PI_FACTOR  = 1.0 / (4.0 * ATAN(1.0))

        q = SQRT(SUM(pos**2)) / this%h
        IF (q.LE.1.0) THEN
            dW = -(3.0 * pos) + (2.25 * pos * q)
            dW = PI_FACTOR * dW / (this%h**5)
        ELSE IF (q.LE.2.0) THEN
            dW = -0.75 * ((2.0 - q)**2) * (pos / q)
            dW = PI_FACTOR * dW / (this%h**5)
        ELSE
            dW = 0.0
        ENDIF
    END FUNCTION func_del_cubic_spline_kernel

END MODULE sph_kernel