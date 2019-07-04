module const_mod

  implicit none

  real(8), parameter :: pi = atan(1.0d0) * 4.0d0
  real(8), parameter :: rad = pi / 180.0d0
  real(8), parameter :: deg = 180.0d0 / pi
  real(8), parameter :: missing_value = 1.0d20

end module const_mod
