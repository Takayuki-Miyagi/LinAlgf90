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
  private :: VectorCopyD, VectorCopyC, MatrixCopyD, MatrixCopyC
  private :: VectorSumD, VectorSumC, MatrixSumD, MatrixSumC
  private :: VectorSubtractD, VectorSubtractC, MatrixSubtractD, MatrixSubtractC
  private :: VectorScaleRD, VectorScaleRC, VectorScaleLD, VectorScaleLC, InnerProductD
  private :: InnerProductC, MatrixScaleLC, MatrixScaleLD, MatrixScaleRC, MatrixScaleRD
  private :: MatrixProductC, MatrixProductD, MVProductD,  VMProductD, MVProductC, VMProductC
  private :: VectorDivideD, VectorDivideC, MatrixScaleDivideD, MatrixScaleDivideC
  private :: OuterProductD, OuterProductC, InitEigenSolSymD, FinEigenSolSymD, DiagSymD, EigenvalD
  private :: InitEigenSolHermite, FinEigenSolHermite, DiagHermite!, EigenvalHermite

  public :: assignment(=), operator(+), operator(-), operator(*)
  public :: operator(/), operator(.x.), EigenSolSymD, EigenSolHermite, exp

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
    procedure :: DiagSym => DiagSymD   ! eigen values and eigen vectors
    procedure :: Eigenval => EigenvalD ! only eigen values
  end type EigenSolSymD

  type :: EigenSolHermite
    type(DVec) :: eig
    type(CMat) :: vec
  contains
    procedure :: init => InitEigenSolHermite
    procedure :: fin => FinEigenSolHermite
    procedure :: DiagSym => DiagHermite      ! eigen values and eigen vectors
    !procedure :: Eigenval => EigenvalHermite ! only eigen values
  end type EigenSolHermite
contains

  subroutine InitEigenSolSymD(this, A)
    class(EigenSolSymD) :: this
    type(DMat), intent(in) :: A
    integer(4) :: n
    n = size(A%m, 1)
    call this%eig%ini(n)
    call this%vec%ini(n,n)
  end subroutine InitEigenSolSymD

  subroutine FinEigenSolSymD(this)
    class(EigenSolSymD) :: this
    call this%eig%fin()
    call this%vec%fin()
  end subroutine FinEigenSolSymD

  subroutine DiagSymD(this, A, qmin, qmax, m, error)
    class(EigenSolSymD) :: this
    type(DMat), intent(in) :: A
    real(8), intent(in), optional :: qmin, qmax
    integer(4), intent(in), optional :: m
    integer(4), intent(in), optional :: error
    real(8), allocatable :: work(:), rcondz(:), zerrbd(:), mat(:,:)
    integer(4), allocatable :: iwork(:), ifailv(:)
    integer(4) :: info, lwork, n, i, num
    real(8) :: lw, dlamch, e, eerbd
    n = size(A%M, 1)
    this%vec = A

    if(.not. present(m) .and. .not. present(qmin) .and. .not. present(qmax)) then
      !
      ! solve all eigen values and eigen vectors
      !

      call dsyev('v', 'u', n, this%vec%m, n, this%eig%v, lw, -1, info)
      lwork = int(lw)
      allocate(work(lwork))
      call dsyev('v', 'u', n, this%vec%m, n, this%eig%v, work, lwork, info)
      if(info /= 0) then
        write(*,'(a, i6)') 'Error in DiagSym: info = ', info
        stop
      end if
      deallocate(work)
      do i = 1, n
        if(this%vec%m(1,i) < 0.d0) this%vec%m(:,i) = this%vec%m(:,i) * (-1.d0)
      end do
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
      !
      ! solve lowest m eigen values and eigen vectors
      !
      this%eig%v(:) = 0.d0
      allocate(mat(n,n))
      allocate(iwork(5*n), ifailv(n))
      mat = A%m
      call dsyevx('v', 'i', 'u', n, mat, n, -1.d100, 1.d100, 1, m, dlamch('S'), &
          &  num, this%eig%v, this%vec%m, n, lw, -1, iwork, ifailv, info)
      lwork = int(lw)
      allocate(work(1:lwork))
      call dsyevx('v', 'i', 'u', n, mat, n, -1.d100, 1.d100, 1, m, dlamch('S'), &
          &  num, this%eig%v, this%vec%m, n, work, lwork, iwork, ifailv, info)
      this%vec%m(:,num+1:n) = 0.d0
      deallocate( iwork, ifailv, work, mat)
      do i = 1, num
        if(this%vec%m(1,i) < 0.d0) this%vec%m(:,i) = this%vec%m(:,i) * (-1.d0)
      end do

    elseif(present(qmin) .and. present(qmax)) then
      !
      ! solve eigen values in (qmin, qmax) and
      ! corrsponding eigen vectors
      !

      allocate(mat(n,n))
      allocate(iwork(5*n), ifailv(n))
      mat = A%m
      this%eig%v(:) = 0.d0
      call dsyevx('v', 'v', 'u', n, mat, n, qmin, qmax, 1, n, dlamch('S'), &
          &  num, this%eig%v, this%vec%m, n, lw, -1, iwork, ifailv, info)
      lwork = int(lw)
      allocate(work(1:lwork))
      call dsyevx('v', 'v', 'u', n, mat, n, qmin, qmax, 1, n, dlamch('S'), &
          &  num, this%eig%v, this%vec%m, n, work, lwork, iwork, ifailv, info)
      this%vec%m(:,num+1:n) = 0.d0
      deallocate( iwork, ifailv, work, mat )
      do i = 1, num
        if(this%vec%m(1,i) < 0.d0) this%vec%m(:,i) = this%vec%m(:,i) * (-1.d0)
      end do
    end if
  end subroutine DiagSymD

  subroutine EigenvalD(this, A, m)
    class(EigenSolSymD) :: this
    type(DMat), intent(in) :: A
    integer(4), intent(in) :: m
    integer(4), allocatable :: iwork(:), iblock(:), isplit(:)
    real(8), allocatable :: work(:), d(:), e(:), tau(:), w(:)
    real(8) :: dlamch, lw
    integer(4) :: n, info, lwork, nsplit, i

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
  end subroutine EigenvalD

  subroutine InitEigenSolHermite(this, A)
    class(EigenSolHermite) :: this
    type(CMat), intent(in) :: A
    integer(4) :: n
    n = size(A%m, 1)
    call this%eig%ini(n)
    call this%vec%ini(n,n)
  end subroutine InitEigenSolHermite

  subroutine FinEigenSolHermite(this)
    class(EigenSolHermite) :: this
    call this%eig%fin()
    call this%vec%fin()
  end subroutine FinEigenSolHermite

  subroutine DiagHermite(this, A, qmin, qmax, m, error)
    class(EigenSolHermite) :: this
    type(CMat), intent(in) :: A
    real(8), intent(in), optional :: qmin, qmax
    integer(4), intent(in), optional :: m
    integer(4), intent(in), optional :: error
    complex(8), allocatable :: work(:), mat(:,:)
    real(8), allocatable :: rcondz(:), zerrbd(:), rwork(:)
    integer(4), allocatable :: iwork(:), ifailv(:)
    integer(4) :: info, lwork, n, i, num
    real(8) :: lw, dlamch, e, eerbd
    n = size(A%M, 1)
    this%vec = A

    if(.not. present(m) .and. .not. present(qmin) .and. .not. present(qmax)) then
      !
      ! solve all eigen values and eigen vectors
      !

      allocate(rwork(3*n - 2))
      call zheev('v', 'u', n, this%vec%m, n, this%eig%v, lw, -1, rwork, info)
      lwork = int(lw)
      allocate(work(lwork))
      call zheev('v', 'u', n, this%vec%m, n, this%eig%v, work, lwork, rwork, info)
      if(info /= 0) then
        write(*,'(a, i6)') 'Error in DiagSym: info = ', info
        stop
      end if
      deallocate(work, rwork)
      !if(present(error)) then
      !  allocate(rcondz(n), zerrbd(n))
      !  e = epsilon(1.d0)
      !  eerbd = e * max(abs(this%eig%v(1)), abs(this%eig%v(n)))
      !  call ddisna('Eigenvectors', n, n, this%eig%v, rcondz, info)
      !  do i = 1, n
      !    zerrbd(i) = eerbd / rcondz(i)
      !  end do
      !  write(*,'(a)') 'Error estimate for eigen values'
      !  write(*, '(1es12.4)') eerbd
      !  write(*,'(a)') 'Error estimate for eigen vectors'
      !  write(*, '(10es12.4)') zerrbd
      !  deallocate(rcondz, zerrbd)
      !end if

    !elseif(present(m)) then
    !  !
    !  ! solve lowest m eigen values and eigen vectors
    !  !
    !  this%eig%v(:) = 0.d0
    !  allocate(mat(n,n))
    !  allocate(iwork(5*n), ifailv(n))
    !  mat = A%m
    !  call dsyevx('v', 'i', 'u', n, mat, n, -1.d100, 1.d100, 1, m, dlamch('S'), &
    !      &  num, this%eig%v, this%vec%m, n, lw, -1, iwork, ifailv, info)
    !  lwork = int(lw)
    !  allocate(work(1:lwork))
    !  call dsyevx('v', 'i', 'u', n, mat, n, -1.d100, 1.d100, 1, m, dlamch('S'), &
    !      &  num, this%eig%v, this%vec%m, n, work, lwork, iwork, ifailv, info)
    !  this%vec%m(:,num+1:n) = 0.d0
    !  deallocate( iwork, ifailv, work, mat)

    !elseif(present(qmin) .and. present(qmax)) then
    !  !
    !  ! solve eigen values in (qmin, qmax) and
    !  ! corrsponding eigen vectors
    !  !

    !  allocate(mat(n,n))
    !  allocate(iwork(5*n), ifailv(n))
    !  mat = A%m
    !  this%eig%v(:) = 0.d0
    !  call dsyevx('v', 'v', 'u', n, mat, n, qmin, qmax, 1, n, dlamch('S'), &
    !      &  num, this%eig%v, this%vec%m, n, lw, -1, iwork, ifailv, info)
    !  lwork = int(lw)
    !  allocate(work(1:lwork))
    !  call dsyevx('v', 'v', 'u', n, mat, n, qmin, qmax, 1, n, dlamch('S'), &
    !      &  num, this%eig%v, this%vec%m, n, work, lwork, iwork, ifailv, info)
    !  this%vec%m(:,num+1:n) = 0.d0
    !  deallocate( iwork, ifailv, work, mat )
    end if
  end subroutine DiagHermite

  !subroutine EigenvalHermite(this, A, m)
  !  class(EigenSolHermite) :: this
  !  type(DMat), intent(in) :: A
  !  integer(4), intent(in) :: m
  !  integer(4), allocatable :: iwork(:), iblock(:), isplit(:)
  !  real(8), allocatable :: work(:), d(:), e(:), tau(:), w(:)
  !  real(8) :: dlamch, lw
  !  integer(4) :: n, info, lwork, nsplit, i

  !  n = size(A%M, 1)
  !  allocate(d(n), e(max(1, n-1)), tau(max(1, n-1)), w(n))
  !  allocate(iblock(n), isplit(n))
  !  call dsytrd('u',n,A%m,n,d,e,tau,lw,-1,info)
  !  lwork = int(lw)
  !  allocate(work(lwork))
  !  call dsytrd('u',n,A%m,n,d,e,tau,work,lwork,info)
  !  call dstebz('i','e',n,0.d0,0.d0,1,min(n,m),dlamch('S'), &
  !      & d, e, m, nsplit, w, iblock, isplit, work, iwork, info)
  !  do i = 1, min(n,m)
  !    this%eig%v(i) = w(i)
  !  end do
  !  deallocate(work)
  !  deallocate(d, e, tau,iblock, isplit)
  !end subroutine EigenvalHermite

  type(DMat) function ExpD(a, ord) result(r)
    type(DMat), intent(in) :: a
    type(DMat) :: b
    integer(4), intent(in), optional :: ord
    integer(4) :: i
    integer(4) :: iord = 12
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
    integer(4), intent(in), optional :: ord
    integer(4) :: i
    integer(4) :: iord = 12
    if(present(ord)) iord = ord
    call r%eye(size(a%m, 1))
    b = r
    do i = 1, iord
      b = b * a / dble(i)
      r = r + b
    end do
  end function ExpC
end module LinAlgLib
