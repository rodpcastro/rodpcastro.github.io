! Extract from incoef.f90

module incoef

  implicit none
  private
  public :: gauss_fortran

contains

  subroutine gauss_fortran(field_local, element_length, G, Q)
    ! Get influence coefficients at a field point using Gauss-Legendre quadrature.

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
