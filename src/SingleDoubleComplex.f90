module SingleDoubleComplex
  use VectorSingle, only: SVec
  use MatrixSingle, only: SMat
  use VectorDouble, only: DVec
  use MatrixDouble, only: DMat
  use VectorComplex, only: CVec
  use MatrixComplex, only: CMat
  implicit none
  public :: DVec2SVec
  public :: SVec2DVec
  public :: SVec2CVec
  public :: DVec2CVec
  public :: CVec2DVecReal
  public :: CVec2SVecReal

  public :: SMat2DMat
  public :: DMat2SMat
  public :: SMat2CMat
  public :: DMat2CMat
  public :: CMat2DMatReal
  public :: CMat2SMatReal
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

  subroutine SVec2CVec(a,b)
    type(CVec), intent(out) :: a
    type(SVec), intent(in) :: b
    integer :: n
    n = size(b%v)
    call a%zeros(n)
    a%v = dble(b%v)
  end subroutine SVec2CVec

  subroutine DVec2CVec(a,b)
    type(CVec), intent(out) :: a
    type(DVec), intent(in) :: b
    integer :: n
    n = size(b%v)
    call a%zeros(n)
    a%v = b%v
  end subroutine DVec2CVec

  subroutine CVec2DVecReal(a,b)
    type(DVec), intent(out) :: a
    type(CVec), intent(in) :: b
    integer :: n
    n = size(b%v)
    call a%ini(n)
    a%v = dble(b%v)
  end subroutine CVec2DVecReal

  subroutine CVec2SVecReal(a,b)
    type(SVec), intent(out) :: a
    type(CVec), intent(in) :: b
    integer :: n
    n = size(b%v)
    call a%ini(n)
    a%v = real(dble(b%v))
  end subroutine CVec2SVecReal

  subroutine DMat2SMat(a,b)
    type(SMat), intent(out) :: a
    type(DMat), intent(in) :: b
    integer :: n, m
    n = size(b%m,1)
    m = size(b%m,2)
    call a%ini(n,m)
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

  subroutine SMat2CMat(a,b)
    type(CMat), intent(out) :: a
    type(SMat), intent(in) :: b
    integer :: n, m
    n = size(b%m,1)
    m = size(b%m,2)
    call a%zeros(n,m)
    a%m = dble(b%m)
  end subroutine SMat2CMat

  subroutine DMat2CMat(a,b)
    type(CMat), intent(out) :: a
    type(DMat), intent(in) :: b
    integer :: n, m
    n = size(b%m,1)
    m = size(b%m,2)
    call a%zeros(n,m)
    a%m = dble(b%m)
  end subroutine DMat2CMat

  subroutine CMat2SMatReal(a,b)
    type(SMat), intent(out) :: a
    type(CMat), intent(in) :: b
    integer :: n, m
    n = size(b%m,1)
    m = size(b%m,2)
    call a%ini(n,m)
    a%m = real(dble(b%m))
  end subroutine CMat2SMatReal

  subroutine CMat2DMatReal(a,b)
    type(DMat), intent(out) :: a
    type(CMat), intent(in) :: b
    integer :: n, m
    n = size(b%m,1)
    m = size(b%m,2)
    call a%ini(n,m)
    a%m = dble(b%m)
  end subroutine CMat2DMatReal
end module SingleDoubleComplex
