module LinAlgParameters
  use, intrinsic :: iso_fortran_env
  implicit none
  integer, parameter :: dp = real64
  integer, parameter :: sp = real32
  integer, parameter :: kp = int32
  integer :: iseed(4) = (/3239, 4241, 1903, 1093/)
end module LinAlgParameters
