module LinAlgLib
  use LinAlgParameters
  use SingleDoubleComplex
  use VectorSingle
  use VectorDouble
  use VectorComplex
  use MatrixSingle
  use MatrixDouble
  use MatrixComplex
  use MatVecSingle
  use MatVecDouble
  use MatVecComplex

  implicit none

  public :: assignment(=)
  public :: operator(+)
  public :: operator(-)
  public :: operator(*)
  public :: operator(/)
  public :: operator(.x.)
  public :: EigenSolSymD
  public :: EigenSolHermite
  public :: exp

  ! SVec methods
  private :: VectorCopyS
  private :: VectorSumS
  private :: VectorSubtractS
  private :: VectorScaleRS
  private :: VectorScaleLS
  private :: InnerProductS
  private :: VectorDivideS
  private :: OuterProductS

  ! DVec methods
  private :: VectorCopyD
  private :: VectorSumD
  private :: VectorSubtractD
  private :: VectorScaleRD
  private :: VectorScaleLD
  private :: InnerProductD
  private :: VectorDivideD
  private :: OuterProductD

  ! CVec methods
  private :: VectorCopyC
  private :: VectorSumC
  private :: VectorSubtractC
  private :: VectorScaleRC
  private :: VectorScaleLC
  private :: InnerProductC
  private :: VectorDivideC
  private :: OuterProductC

  ! SMat methods
  private :: MatrixCopyS
  private :: MatrixSumS
  private :: MatrixSubtractS
  private :: MatrixScaleLS
  private :: MatrixScaleRS
  private :: MatrixProductS
  private :: MatrixScaleDivideS

  ! DMat methods
  private :: MatrixCopyD
  private :: MatrixSumD
  private :: MatrixSubtractD
  private :: MatrixScaleLD
  private :: MatrixScaleRD
  private :: MatrixProductD
  private :: MatrixScaleDivideD

  ! CMat methods
  private :: MatrixCopyC
  private :: MatrixSumC
  private :: MatrixSubtractC
  private :: MatrixScaleLC
  private :: MatrixScaleRC
  private :: MatrixProductC
  private :: MatrixScaleDivideC

  private :: DVec2SVec
  private :: SVec2DVec
  private :: SMat2DMat
  private :: DMat2SMat

  ! MatVec single
  private :: MVProductS
  private :: VMProductS
  ! MatVec double
  private :: MVProductD
  private :: VMProductD
  ! MatVec complex
  private :: MVProductC
  private :: VMProductC

  ! diagonalization methods
  private :: InitEigenSolSymD
  private :: FinEigenSolSymD
  private :: DiagSymD
  private :: EigenvalD
  private :: InitEigenSolHermite
  private :: FinEigenSolHermite
  private :: DiagHermite
  !private :: EigenvalHermite

  interface assignment(=)
    module procedure :: VectorCopyD
    module procedure :: VectorCopyS
    module procedure :: VectorCopyC
    module procedure :: MatrixCopyS
    module procedure :: MatrixCopyD
    module procedure :: MatrixCopyC
    ! single <-> double <-> complex vector
    module procedure :: DVec2SVec
    module procedure :: SVec2DVec
    module procedure :: SVec2CVec
    module procedure :: DVec2CVec
    module procedure :: CVec2DVecReal
    module procedure :: CVec2SVecReal
    ! single <-> double <-> complex matrix
    module procedure :: SMat2DMat
    module procedure :: DMat2SMat
    module procedure :: SMat2CMat
    module procedure :: DMat2CMat
    module procedure :: CMat2DMatReal
    module procedure :: CMat2SMatReal
  end interface assignment(=)

  interface operator(+)
    module procedure :: VectorSumD
    module procedure :: VectorSumS
    module procedure :: VectorSumC
    module procedure :: MatrixSumS
    module procedure :: MatrixSumD
    module procedure :: MatrixSumC
  end interface operator(+)

  interface operator(-)
    module procedure :: VectorSubtractD
    module procedure :: VectorSubtractS
    module procedure :: VectorSubtractC
    module procedure :: MatrixSubtractS
    module procedure :: MatrixSubtractD
    module procedure :: MatrixSubtractC
  end interface operator(-)

  interface operator(*)
    module procedure :: VectorScaleRD
    module procedure :: VectorScaleRS
    module procedure :: VectorScaleRC
    module procedure :: VectorScaleLS
    module procedure :: VectorScaleLD
    module procedure :: VectorScaleLC
    module procedure :: InnerProductS
    module procedure :: InnerProductD
    module procedure :: InnerProductC
    module procedure :: MatrixScaleLS
    module procedure :: MatrixScaleLC
    module procedure :: MatrixScaleLD
    module procedure :: MatrixScaleRS
    module procedure :: MatrixScaleRC
    module procedure :: MatrixScaleRD
    module procedure :: MatrixProductS
    module procedure :: MatrixProductC
    module procedure :: MatrixProductD
    module procedure :: MVProductS
    module procedure :: VMProductS
    module procedure :: MVProductD
    module procedure :: VMProductD
    module procedure :: MVProductC
    module procedure :: VMProductC
  end interface operator(*)

  interface operator(/)
    module procedure :: VectorDivideD
    module procedure :: VectorDivideS
    module procedure :: VectorDivideC
    module procedure :: MatrixScaleDivideS
    module procedure :: MatrixScaleDivideD
    module procedure :: MatrixScaleDivideC
  end interface operator(/)

  interface operator(.x.)
    module procedure :: OuterProductD
    module procedure :: OuterProductS
    module procedure :: OuterProductC
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
    integer(kp) :: n
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
    real(dp), intent(in), optional :: qmin, qmax
    integer(kp), intent(in), optional :: m
    integer(kp), intent(in), optional :: error
    real(dp), allocatable :: work(:), rcondz(:), zerrbd(:), mat(:,:)
    integer(kp), allocatable :: iwork(:), ifailv(:)
    integer(kp) :: info, lwork, n, i, num
    real(dp) :: lw, dlamch, e, eerbd
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
    integer(kp), intent(in) :: m
    integer(kp), allocatable :: iwork(:), iblock(:), isplit(:)
    real(dp), allocatable :: work(:), d(:), e(:), tau(:), w(:)
    real(dp) :: dlamch, lw
    integer(kp) :: n, info, lwork, nsplit, i

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
    integer(kp) :: n
    n = size(A%m, 1)
    call this%eig%ini(n)
    call this%vec%ini(n,n)
  end subroutine InitEigenSolHermite

  subroutine FinEigenSolHermite(this)
    class(EigenSolHermite) :: this
    call this%eig%fin()
    call this%vec%fin()
  end subroutine FinEigenSolHermite

  !subroutine DiagHermite(this, A, qmin, qmax, m, error)
  subroutine DiagHermite(this, A, qmin, qmax, m)
    class(EigenSolHermite) :: this
    type(CMat), intent(in) :: A
    real(dp), intent(in), optional :: qmin, qmax
    integer(kp), intent(in), optional :: m
    !integer(kp), intent(in), optional :: error
    complex(dp), allocatable :: work(:)
    real(dp), allocatable :: rwork(:)
    integer(kp) :: info, lwork, n
    real(dp) :: lw
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
  !  integer(kp), intent(in) :: m
  !  integer(kp), allocatable :: iwork(:), iblock(:), isplit(:)
  !  real(dp), allocatable :: work(:), d(:), e(:), tau(:), w(:)
  !  real(dp) :: dlamch, lw
  !  integer(kp) :: n, info, lwork, nsplit, i

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
    integer(kp), intent(in), optional :: ord
    integer(kp) :: i
    integer(kp) :: iord = 12
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
    integer(kp), intent(in), optional :: ord
    integer(kp) :: i
    integer(kp) :: iord = 12
    if(present(ord)) iord = ord
    call r%eye(size(a%m, 1))
    b = r
    do i = 1, iord
      b = b * a / dble(i)
      r = r + b
    end do
  end function ExpC
end module LinAlgLib
