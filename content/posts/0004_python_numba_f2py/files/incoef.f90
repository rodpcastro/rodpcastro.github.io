module incoef
! Calculation of boundary element influence coefficients G and Q.
!
! Subroutines
! -----------
! analytical_fortran(field_local, element_length, G, Q)
!   Get influence coefficients at a field point using analytical integration.
! gauss_fortran(field_local, element_length, G, Q)
!   Get influence coefficients at a field point using Gauss-Legendre quadrature.
!
! Author
! ------
! Rodrigo Castro (GitHub: rodpcastro)
!
! History
! -------
! 14-06-2025 - Rodrigo Castro - Original code
!
! References
! ----------
! 1. Rodrigo Castro. 2025. 2D constant boundary element influence coefficients.
!    https://rodpcastro.github.io/posts/3_2d_constant_boundary_element/

  implicit none
  private
  public :: analytical_fortran, gauss_fortran

contains

  logical function ismall(x, ref)
    ! Evaluates the smallness of a variable compared to a reference value.
    !
    ! .true. if x is small compared to ref according to eps, and .false. otherwise.
    !
    ! Parameters
    ! ----------
    ! x : real(dp)
    !   Variable to be checked for smallness.
    !    
    ! Returns
    ! -------
    ! ref : real(dp) 
    !   Reference value.

    implicit none
    integer, parameter :: dp = kind(1.0d0)
    real(dp), intent(in) :: x
    real(dp), intent(in), optional :: ref
    real(dp), parameter :: eps = epsilon(1.0_dp)
    real(dp) :: ref_

    ref_ = 1.0_dp
    if (present(ref)) ref_ = ref

    ismall = abs(x) < eps * abs(ref_)
  end function ismall

  subroutine analytical_fortran(field_local, element_length, G, Q)
    ! Get influence coefficients at a field point using analytical integration.
    !
    ! Parameters
    ! ----------
    ! field_local : real(dp), dimension(2)
    !   Array with field point's coordinates in the local system.
    ! element_length : real(dp)
    !   Boundary element's length.
    !    
    ! Returns
    ! -------
    ! G : real(dp)
    !   Integral of the Green's function over the element.
    ! Q : real(dp)
    !   Integral of the Green's function normal derivative over the element.

    implicit none
    integer, parameter :: dp = kind(1.0d0)
    real(dp), intent(in) :: field_local(2)
    real(dp), intent(in) :: element_length
    real(dp), intent(out) :: G, Q
    real(dp), parameter :: pi = 3.141592653589793_dp
    real(dp) :: a, x, y
    real(dp) :: xpa, xma
    real(dp) :: r1, r2, t1, t2

    a = 0.5_dp * element_length
    x = field_local(1)
    y = field_local(2)

    if (ismall(y, element_length) .and. ismall(abs(x) - a, element_length)) then
      G = a / pi * (log(2.0_dp * a) - 1.0_dp)
      Q = 0.0_dp
    else
      xpa = x + a
      xma = x - a

      r1 = sqrt(xma**2 + y**2)
      r2 = sqrt(xpa**2 + y**2)
      t1 = atan2(y, xma)
      t2 = atan2(y, xpa)

      G = 0.5_dp / pi * (y * (t1 - t2) - xma * log(r1) + xpa * log(r2) - 2.0_dp * a)

      if (ismall(y, element_length)) then
        ! Q is discontinuous in |x| < a and y = 0.
        Q = 0.0_dp
      else
        Q = -0.5_dp / pi * (t1 - t2)
      end if
    end if
  end subroutine analytical_fortran

  subroutine gauss_fortran(field_local, element_length, G, Q)
    ! Get influence coefficients at a field point using Gauss-Legendre quadrature.
    !
    ! Parameters
    ! ----------
    ! field_local : real(dp), dimension(2)
    !   Array with field point's coordinates in the local system.
    ! element_length : real(dp)
    !   Boundary element's length.
    !    
    ! Returns
    ! -------
    ! G : real(dp)
    !   Integral of the Green's function over the element.
    ! Q : real(dp)
    !   Integral of the Green's function normal derivative over the element.

    implicit none
    integer, parameter :: dp = kind(1.0d0)
    real(dp), intent(in) :: field_local(2)
    real(dp), intent(in) :: element_length
    real(dp), intent(out) :: G, Q
    real(dp), parameter :: pi = 3.141592653589793_dp
    real(dp), parameter :: roots(2) = [0.3399810435848563_dp, 0.8611363115940526_dp]
    real(dp), parameter :: weights(2) = [0.6521451548625461_dp, 0.3478548451374538_dp]
    real(dp) :: a, x, y
    integer :: i

    a = 0.5_dp * element_length
    x = field_local(1)
    y = field_local(2)

    G = 0.0_dp
    Q = 0.0_dp
    do i = 1, size(roots)
      G = G + weights(i) * (intG(a * roots(i)) + intG(-a * roots(i)))
      Q = Q + weights(i) * (intQ(a * roots(i)) + intQ(-a * roots(i)))
    end do

    G = G * 0.25_dp * a / pi
    Q = -Q * 0.5_dp * a / pi

  contains
    
    real(dp) function intG(t)
      real(dp), intent(in) :: t
      intG = log((x - t)**2 + y**2)
    end function intG

    real(dp) function intQ(t)
      real(dp), intent(in) :: t
      intQ = y / ((x - t)**2 + y**2)
    end function intQ

  end subroutine gauss_fortran

end module incoef
