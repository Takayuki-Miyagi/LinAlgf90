module SingleDouble
  use VectorSingle, only: SVec
  use MatrixSingle, only: SMat
  use VectorDouble, only: DVec
  use MatrixDouble, only: DMat
  implicit none
  public :: DVec2SVec, SVec2DVec, SMat2DMat, DMat2SMat
contains
  subroutine DVec2SVec(a,b)
    type(SVec), intent(out) :: a
    type(DVec), intent(in) :: b
    integer :: n
    n = size(b%v)
    allocate(a%v(n))
    a%v = real(b%v)
  end subroutine DVec2SVec

  subroutine SVec2DVec(a,b)
    type(DVec), intent(out) :: a
    type(SVec), intent(in) :: b
    integer :: n
    n = size(b%v)
    allocate(a%v(n))
    a%v = dble(b%v)
  end subroutine SVec2DVec

  subroutine DMat2SMat(a,b)
    type(SMat), intent(out) :: a
    type(DMat), intent(in) :: b
    integer :: n, m
    n = size(b%m,1)
    m = size(b%m,2)
    allocate(a%m(n,m))
    a%m = real(b%m)
  end subroutine DMat2SMat

  subroutine SMat2DMat(a,b)
    type(DMat), intent(out) :: a
    type(SMat), intent(in) :: b
    integer :: n, m
    n = size(b%m,1)
    m = size(b%m,2)
    allocate(a%m(n,m))
    a%m = dble(b%m)
  end subroutine SMat2DMat
end module SingleDouble
