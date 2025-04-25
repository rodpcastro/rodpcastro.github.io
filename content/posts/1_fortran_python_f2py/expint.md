```fortran
subroutine expi(x, ei)
  ! Exponential integral Ei(x).
  !
  ! Parameters
  ! ----------
  ! x : real(real64)
  !   Real number ≥ 0.
  !    
  ! Returns
  ! -------
  ! ei : real(real64) 
  !   Exponential integral of x.
  !
  ! Author
  ! ------
  ! Rodrigo Castro (GitHub: rodpcastro)
  !
  ! History
  ! -------
  ! 24-04-2025 - Rodrigo Castro - Original code
  !
  ! Reference
  ! ----------
  ! Shanjie Zhang, Jianming Jin (1996). Computation of Special Functions.

  use, intrinsic :: iso_fortran_env, only: int16, real64

  implicit none
  real(real64), parameter :: pi = 3.141592653589793d0
  real(real64), parameter :: gm = 0.5772156649015329d0 ! Euler's constant
  real(real64), intent(in) :: x
  real(real64), intent(out) :: ei
  real(real64) :: r
  integer(int16) :: n

  if (x == 0.0d0) then
    ei = -1.0d+300
  else if (x <= 40.0d0) then
    ei = x
    r = x
    do n = 2, 100
      r = r * x * (n-1) / n**2
      ei = ei + r
      if (abs(r) <= 1.0d-15*abs(ei)) exit
    end do
    ei = ei + gm + log(x)
  else
    ei = 1.0d0
    r = 1.0d0
    do n = 1, 20
      r = r * n / x
      ei = ei + r
    end do
    ei = exp(x) * ei / x 
  end if

end subroutine expi
```
