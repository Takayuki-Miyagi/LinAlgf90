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
  public :: DSingularValueDecomposition
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

  type :: GenEigenSolSymD
    type(DVec) :: eig
    type(DMat) :: vec
    integer :: itype = 1
  contains
    procedure :: init => InitGenEigenSolSymD
    procedure :: fin => FinGenEigenSolSymD
    procedure :: DiagSym => DiagGenSymD   ! eigen values and eigen vectors
  end type GenEigenSolSymD

  type :: EigenSolHermite
    type(DVec) :: eig
    type(CMat) :: vec
  contains
    procedure :: init => InitEigenSolHermite
    procedure :: fin => FinEigenSolHermite
    procedure :: DiagSym => DiagHermite      ! eigen values and eigen vectors
    !procedure :: Eigenval => EigenvalHermite ! only eigen values
  end type EigenSolHermite

  type :: DSingularValueDecomposition
    type(DMat) :: U, V
    type(DVec) :: Sigma
  contains
    procedure :: init => InitDSVD
    procedure :: fin => FinDSVD
    procedure :: SVD => DSVD
  end type DSingularValueDecomposition
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
    real(dp) :: dlamch, e, eerbd
    n = size(A%M, 1)
    this%vec = A

    if(.not. present(m) .and. .not. present(qmin) .and. .not. present(qmax)) then
      !
      ! solve all eigen values and eigen vectors
      !

      !call dsyev('v', 'u', n, this%vec%m, n, this%eig%v, lw, -1, info)
      !lwork = int(lw)
      allocate(work(1))
      call dsyev('v', 'u', n, this%vec%m, n, this%eig%v, work, -1, info)
      lwork = int(work(1))
      deallocate(work)
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
      !call dsyevx('v', 'i', 'u', n, mat, n, -1.d100, 1.d100, 1, m, dlamch('S'), &
      !    &  num, this%eig%v, this%vec%m, n, lw, -1, iwork, ifailv, info)
      !lwork = int(lw)
      allocate(work(1))
      call dsyevx('v', 'i', 'u', n, mat, n, -1.d100, 1.d100, 1, m, dlamch('S'), &
          &  num, this%eig%v, this%vec%m, n, work, -1, iwork, ifailv, info)
      lwork = int(work(1))
      deallocate(work)
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
      !call dsyevx('v', 'v', 'u', n, mat, n, qmin, qmax, 1, n, dlamch('S'), &
      !    &  num, this%eig%v, this%vec%m, n, lw, -1, iwork, ifailv, info)
      !lwork = int(lw)
      allocate(work(1))
      call dsyevx('v', 'v', 'u', n, mat, n, qmin, qmax, 1, n, dlamch('S'), &
          &  num, this%eig%v, this%vec%m, n, work, -1, iwork, ifailv, info)
      lwork = int(work(1))
      deallocate(work)
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
    real(dp) :: dlamch
    integer(kp) :: n, info, lwork, nsplit, i

    n = size(A%M, 1)
    allocate(d(n), e(max(1, n-1)), tau(max(1, n-1)), w(n))
    allocate(iblock(n), isplit(n))
    !call dsytrd('u',n,A%m,n,d,e,tau,lw,-1,info)
    !lwork = int(lw)
    allocate(work(1))
    call dsytrd('u',n,A%m,n,d,e,tau,work,-1,info)
    lwork = int(work(1))
    deallocate(work)
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

  subroutine InitGenEigenSolSymD(this, A, B, itype)
    class(GenEigenSolSymD) :: this
    type(DMat), intent(in) :: A, B
    integer(kp), intent(in), optional :: itype
    integer(kp) :: n
    if(present(itype)) this%itype = itype
    n = size(A%m, 1)
    if(this%itype == 2) n = size(A%m, 1)
    if(this%itype == 3) n = size(B%m, 1)
    call this%eig%ini(n)
    call this%vec%ini(n,n)
  end subroutine InitGenEigenSolSymD

  subroutine FinGenEigenSolSymD(this)
    class(GenEigenSolSymD) :: this
    call this%eig%fin()
    call this%vec%fin()
  end subroutine FinGenEigenSolSymD

  subroutine DiagGenSymD(this, A, B)
    class(GenEigenSolSymD) :: this
    type(DMat), intent(in) :: A, B
    integer(kp) :: n, lda, ldb, lwork, liwork, info
    integer(kp), allocatable :: iwork(:)
    real(dp), allocatable :: work(:)
    this%vec = A
    n = size(A%M, 1)
    lda = size(A%m,1)
    ldb = size(B%m,1)
    !call dsygvd(this%itype, 'V', 'U', n, A%m, lda, B%m, ldb, this%eig%v, dummy, -1, idummy, -1, info)
    !lwork = int(dummy)
    allocate(work(1), iwork(1))
    call dsygvd(this%itype, 'V', 'U', n, A%m, lda, B%m, ldb, this%eig%v, work, -1, iwork, -1, info)
    lwork = int(work(1))
    liwork = iwork(1)
    deallocate(work,iwork)
    allocate(work(lwork), iwork(liwork))
    call dsygvd(this%itype, 'V', 'U', n, A%m, lda, B%m, ldb, this%eig%v, work, lwork, iwork, liwork, info)
    deallocate(work, iwork)
    this%vec = A
  end subroutine DiagGenSymD

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
    n = size(A%M, 1)
    this%vec = A

    if(.not. present(m) .and. .not. present(qmin) .and. .not. present(qmax)) then

      !
      ! solve all eigen values and eigen vectors
      !

      allocate(rwork(3*n - 2))
      !call zheev('v', 'u', n, this%vec%m, n, this%eig%v, lw, -1, rwork, info)
      !lwork = int(lw)
      allocate(work(1))
      call zheev('v', 'u', n, this%vec%m, n, this%eig%v, work, -1, rwork, info)
      lwork = int(work(1))
      deallocate(work)
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

  type(DMat) function ExpD(a, ord, tol_in, show_process) result(r)
    type(DMat), intent(in) :: a
    type(DMat) :: b
    integer(kp), intent(in), optional :: ord
    real(dp), intent(in), optional :: tol_in
    logical, intent(in), optional :: show_process
    real(dp) :: tol = 1.d-8
    logical :: show = .false.
    integer(kp) :: i
    integer(kp) :: iord = 12
    if(present(ord)) iord = ord
    if(present(tol_in)) tol = tol_in
    if(present(show_process)) show = show_process
    call r%eye(size(a%m, 1))
    b = r
    do i = 1, iord
      b = (b * a) / dble(i)
      r = r + b
      if(show) then
        write(*,"(a,i4,a,es18.8)") "max value of matrix: order=",i,&
            & " maxval(A)",maxval(b%m)
      end if
      if((maxval(b%m)) < tol) return
    end do
    write(*,"(a,i4)") "warning: increase ord, current value is ", iord
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

  subroutine InitDSVD(this, A)
    class(DSingularValueDecomposition) :: this
    type(DMat), intent(in) :: A
    integer(kp) :: n
    call this%U%ini(A%n_row, min(A%n_row, A%n_col))
    call this%Sigma%ini(min(A%n_row, A%n_col))
    call this%V%ini(min(A%n_row, A%n_col), A%n_col)
  end subroutine InitDSVD

  subroutine FinDSVD(this)
    class(DSingularValueDecomposition) :: this
    call this%U%fin()
    call this%Sigma%fin()
    call this%V%fin()
  end subroutine FinDSVD

  subroutine DSVD(this, A)
    class(DSingularValueDecomposition) :: this
    type(DMat), intent(in) :: A
    real(dp), allocatable :: work(:)
    integer(kp) :: m, n, lwork, info
    type(DMat) :: Atmp
    Atmp = A
    m = Atmp%n_row
    n = Atmp%n_col

    ! Full decomposition
    allocate(work(1))
    lwork = -1
    call dgesvd("S","S",m,n,Atmp%m,m,this%Sigma%v,this%U%m,m,this%V%m,min(m,n),work,lwork,info)
    lwork = int(work(1))
    deallocate(work)
    allocate(work(lwork))
    call dgesvd("S","S",m,n,Atmp%m,m,this%Sigma%v,this%U%m,m,this%V%m,min(m,n),work,lwork,info)
    deallocate(work)
  end subroutine DSVD
end module LinAlgLib
