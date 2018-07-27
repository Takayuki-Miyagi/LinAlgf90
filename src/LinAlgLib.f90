module LinAlgLib
  use VectorDouble, only: DVec, VectorCopyD, VectorSumD, VectorSubtractD, &
    & VectorScaleRD, VectorScaleLD, InnerProductD, VectorDivideD
  use VectorComplex, only: CVec, VectorCopyC, VectorSumC, VectorSubtractC, &
    & VectorScaleRC, VectorScaleLC, InnerProductC, VectorDivideC
  use MatrixDouble, only: DMat, MatrixCopyD, MatrixSumD, MatrixSubtractD, &
    & MatrixScaleRD, MatrixScaleLD, MatrixProductD, MatrixScaleDivideD
  use MatrixComplex, only: CMat, MatrixCopyC, MatrixSumC, MatrixSubtractC, &
    & MatrixScaleRC, MatrixScaleLC, MatrixProductC, MatrixScaleDivideC
  use MatVecDouble, only: OuterProductD, VMProductD, MVProductD
  use MatVecComplex, only: OuterProductC, VMProductC, MVProductC
  implicit none
  private :: VectorCopyD, VectorCopyC, MatrixCopyD, MatrixCopyC, &
    & VectorSumD, VectorSumC, MatrixSumD, MatrixSumC, &
    & VectorSubtractD, VectorSubtractC, MatrixSubtractD, MatrixSubtractC, &
    & VectorScaleRD, VectorScaleRC, VectorScaleLD, VectorScaleLC, InnerProductD, &
    & InnerProductC, MatrixScaleLC, MatrixScaleLD, MatrixScaleRC, MatrixScaleRD, &
    & MatrixProductC, MatrixProductD, MVProductD,  VMProductD, MVProductC, VMProductC, &
    & VectorDivideD, VectorDivideC, MatrixScaleDivideD, MatrixScaleDivideC, &
    & OuterProductD, OuterProductC, InitEigenSolSymD, FinEigenSolSymD, DiagSym, Eigenval
  public :: assignment(=), operator(+), operator(-), operator(*), &
    & operator(/), operator(.x.), EigenSolSymD

  interface assignment(=)
    module procedure :: VectorCopyD, &
      & VectorCopyC, &
      & MatrixCopyD, &
      & MatrixCopyC
  end interface assignment(=)

  interface operator(+)
    module procedure :: VectorSumD, &
      & VectorSumC, &
      & MatrixSumD, &
      & MatrixSumC
  end interface operator(+)

  interface operator(-)
    module procedure :: VectorSubtractD, &
      & VectorSubtractC, &
      & MatrixSubtractD, &
      & MatrixSubtractC
  end interface operator(-)

  interface operator(*)
    module procedure :: VectorScaleRD, &
      & VectorScaleRC, &
      & VectorScaleLD, &
      & VectorScaleLC, &
      & InnerProductD, &
      & InnerProductC, &
      & MatrixScaleLC, &
      & MatrixScaleLD, &
      & MatrixScaleRC, &
      & MatrixScaleRD, &
      & MatrixProductC,&
      & MatrixProductD,&
      & MVProductD   , &
      & VMProductD   , &
      & MVProductC   , &
      & VMProductC
  end interface operator(*)

  interface operator(/)
    module procedure :: VectorDivideD, &
      & VectorDivideC, &
      & MatrixScaleDivideD, &
      & MatrixScaleDivideC
  end interface operator(/)

  interface operator(.x.)
    module procedure :: OuterProductD, &
      & OuterProductC
  end interface operator(.x.)

  interface exp
    module procedure :: ExpD, ExpC
  end interface exp

  type :: EigenSolSymD
    type(DVec) :: eig
    type(DMat) :: vec
  contains
    procedure :: init => InitEigenSolSymD
    procedure :: fin => FinEigenSolSymD
    procedure :: DiagSym ! eigen values and eigen vectors
    procedure :: Eigenval! only eigen values
  end type EigenSolSymD
contains

  subroutine InitEigenSolSymD(this, A)
    class(EigenSolSymD) :: this
    type(DMat), intent(in) :: A
    integer :: n
    n = size(A%m, 1)
    call this%eig%ini(n)
    call this%vec%ini(n,n)
  end subroutine InitEigenSolSymD

  subroutine FinEigenSolSymD(this)
    class(EigenSolSymD) :: this
    call this%eig%fin()
    call this%vec%fin()
  end subroutine FinEigenSolSymD

  subroutine DiagSym(this, A, qmin, qmax, m, error)
    use LinAlgParameters, only: eps
    class(EigenSolSymD) :: this
    type(DMat), intent(in) :: A
    real(8), intent(in), optional :: qmin, qmax
    integer, intent(in), optional :: m
    integer, intent(in), optional :: error
    integer :: num
    real(8), allocatable :: work(:), rcondz(:), zerrbd(:)
    integer, allocatable :: iwork(:), ifailv(:)
    integer :: info, lwork, n, i
    real(8) :: lw, dlamch, e, eerbd
    n = size(A%M, 1)
    this%vec = A
    if(.not. present(m) .and. .not. present(qmin) .and. &
        & .not. present(qmax)) then
      call dsyev('v', 'u', n, A%m, n, this%eig%v, lw, -1, info)
      lwork = int(lw)
      allocate(work(lwork))
      call dsyev('v', 'u', n, this%vec%m, n, this%eig%v, work, lwork, info)
      if(info /= 0) then
        write(*,'(a, i6)') 'Error in DiagSym: info = ', info
        stop
      end if
      deallocate(work)
      if(present(error)) then
        allocate(rcondz(n), zerrbd(n))
        e = epsilon(1.d0)
        eerbd = e * max(abs(this%eig%v(1)), abs(this%eig%v(n)))
        call ddisna('Eigenvectors', n, n, this%eig%v, rcondz, info)
        do i = 1, n
          zerrbd(i) = eerbd / rcondz(i)
        end do
        write(*,'(a)') 'Error estimate for eigen values'
        write(*, '(1es12.4)') eerbd
        write(*,'(a)') 'Error estimate for eigen vectors'
        write(*, '(10es12.4)') zerrbd
        deallocate(rcondz, zerrbd)
      end if

    elseif(present(m)) then
      allocate(iwork(5*n), ifailv(n))
      call dsyevx('v', 'i', 'u', n, A%m, n, -1.d100, 1.d100, 1, n, dlamch('S'), &
          &  num, this%eig%v, this%vec%m, n, lw, -1, iwork, ifailv, info)
      deallocate( iwork, ifailv)

    else
      allocate(iwork(5*n), ifailv(n))
      call dsyevx('v', 'i', 'u', n, A%m, n, -1.d100, 1.d100, 1, n, dlamch('S'), &
          &  num, this%eig%v, this%vec%m, n, lw, -1, iwork, ifailv, info)
      lwork = int(lw)
      allocate(work(1:lwork))
      call dsyevx('v', 'v', 'u', n, A%m, n, qmin, qmax, 1, n, eps, &
          &  num, this%eig%v, this%vec%m, n, work, lwork, iwork, ifailv, info)
      deallocate( iwork, ifailv)
    end if
  end subroutine DiagSym

  subroutine Eigenval(this, A, m)
    class(EigenSolSymD) :: this
    type(DMat), intent(in) :: A
    integer, intent(in) :: m
    integer, allocatable :: iwork(:), iblock(:), isplit(:)
    real(8), allocatable :: work(:), d(:), e(:), tau(:), w(:)
    real(8) :: dlamch, lw
    integer :: n, info, lwork, nsplit, i

    n = size(A%M, 1)
    allocate(d(n), e(max(1, n-1)), tau(max(1, n-1)), w(n))
    allocate(iblock(n), isplit(n))
    call dsytrd('u',n,A%m,n,d,e,tau,lw,-1,info)
    lwork = int(lw)
    allocate(work(lwork))
    call dsytrd('u',n,A%m,n,d,e,tau,work,lwork,info)
    call dstebz('i','e',n,0.d0,0.d0,1,min(n,m),dlamch('S'), &
        & d, e, m, nsplit, w, iblock, isplit, work, iwork, info)
    do i = 1, min(n,m)
      this%eig%v(i) = w(i)
    end do
    deallocate(work)
    deallocate(d, e, tau,iblock, isplit)
  end subroutine Eigenval

  type(DMat) function ExpD(a, ord) result(r)
    type(DMat), intent(in) :: a
    type(DMat) :: b
    integer, intent(in), optional :: ord
    integer :: i
    integer :: iord = 12
    if(present(ord)) iord = ord
    call r%eye(size(a%m, 1))
    b = r
    do i = 1, iord
      b = b * a / dble(i)
      r = r + b
    end do
  end function ExpD

  type(CMat) function ExpC(a, ord) result(r)
    type(CMat), intent(in) :: a
    type(CMat) :: b
    integer, intent(in), optional :: ord
    integer :: i
    integer :: iord = 12
    if(present(ord)) iord = ord
    call r%eye(size(a%m, 1))
    b = r
    do i = 1, iord
      b = b * a / dble(i)
      r = r + b
    end do
  end function ExpC
end module LinAlgLib
