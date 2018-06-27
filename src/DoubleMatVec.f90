module DMatVec
  implicit none
  integer :: iseed(4) = (/3239, 4241, 1903, 1093/)
  type :: Vec
    real(8), allocatable :: V(:)
  contains
    procedure :: Ini => iniV
    procedure :: Fin => FinV
    procedure :: prt => VectorPrint
    procedure :: Random => GetRandomVector
    procedure :: Nrm
    procedure :: Nrm2
!    procedure :: DVec => DiagMatVec
  end type Vec

!  type :: MO
!    real(8), allocatable :: M(:,:)
!  contains
!    procedure :: Ini => IniM
!    procedure :: Fin => FinM
!    procedure :: eye
!    procedure :: T => Transepose
!    procedure :: Random => GetRandomMatrix
!    procedure :: prt => MatrixPrint
!    procedure :: DMat => VecDiagMat
!  end type MO

  interface assignment(=)
    procedure :: VectorCopy
!    procedure :: MatrixCopy
  end interface assignment(=)

  interface operator(+)
    procedure :: VectorSum
!    procedure :: MatrixSum
  end interface operator(+)

  interface operator(-)
    procedure :: VectorSubtract
!    procedure :: MatrixSubtract
  end interface operator(-)

  interface operator(*)
    procedure :: VectorScaleR
    procedure :: VectorScaleL
!    procedure :: MatrixScaleR
!    procedure :: MatrixScaleL
    procedure :: InnerProduct
!    procedure :: MatrixProduct
!    procedure :: MVProduct
!    procedure :: VMProduct
  end interface operator(*)

  interface operator(/)
    procedure :: VectorDivide
  end interface operator(/)

!  interface operator(.x.)
!    procedure :: OuterProduct
!  end interface operator(.x.)
contains
  subroutine IniV(a, n)
  class(Vec), intent(inout) :: a
    integer, intent(in) :: n
    if(allocated(a%V)) deallocate(a%V)
    allocate(a%V(n))
    a%V(:) = 0.d0
  end subroutine IniV

  subroutine FinV(a)
  class(Vec), intent(inout) :: a
    if(allocated(a%V)) deallocate(a%V)
  end subroutine FinV

  subroutine VectorCopy(b, a)
    type(Vec), intent(inout) :: b
    type(Vec), intent(in) :: a
    integer :: n
    n = size(a%V)
    call b%Ini(n)
    call dcopy(n, a%v, 1, b%v, 1)
  end subroutine VectorCopy

  type(Vec) function VectorSum(a, b) result(c)
  class(Vec), intent(in) :: a, b
    integer :: n
    if(size(a%v) /= size(b%v)) then
      write(*,'(a)') 'Error in VectorSum'
      stop
    end if
    n = size(a%v)
    c = a
    call daxpy(n, 1.d0, b%v, 1, c%v, 1)
  end function VectorSum

  type(Vec) function VectorSubtract(a, b) result(c)
    type(Vec), intent(in) :: a, b
    integer :: n
    if(size(a%v) /= size(b%v)) then
      write(*,'(a)') 'Error in VectorSubtract'
      stop
    end if
    n = size(a%v)
    c = a
    call daxpy(n, -1.d0, b%v, 1, c%v, 1)
  end function VectorSubtract

  type(Vec) function VectorScaleR(a, b) result(c)
    type(Vec), intent(in) :: a
    real(8), intent(in) :: b
    integer :: n
    n = size(a%v)
    c = a
    call dscal(n, b, c%v, 1)
  end function VectorScaleR


  type(Vec) function VectorScaleL(b, a) result(c)
    type(Vec), intent(in) :: a
    real(8), intent(in) :: b
    integer :: n
    n = size(a%v)
    c = a
    call dscal(n, b, c%v, 1)
  end function VectorScaleL

  type(Vec) function VectorDivide(a, b) result(c)
    type(Vec), intent(in) :: a
    real(8), intent(in) :: b
    integer :: n
    n = size(a%v)
    c = a
    call dscal(n, 1.d0 / b, c%v, 1)
  end function VectorDivide

  real(8) function InnerProduct(a, b) result(c)
    type(Vec), intent(in) :: a, b
    integer :: n
    real(8) :: ddot
    if(size(a%v) /= size(b%v)) then
      write(*,'(a)') 'Error in InnerProduct'
      stop
    end if
    n = size(a%v)
    c = ddot(n, a%v, 1, b%v, 1)
  end function InnerProduct

  real(8) function Nrm(a) result(b)
  class(Vec), intent(in) :: a
    integer :: n
    real(8) :: dnrm2
    n = size(a%v)
    b = dnrm2(n, a%v, 1)
  end function Nrm

  real(8) function Nrm2(a) result(b)
  class(Vec), intent(in) :: a
    integer :: n
    real(8) :: ddot
    n = size(a%v)
    b = ddot(n, a%v, 1, a%v, 1)
  end function Nrm2

  subroutine VectorPrint(this, string)
  class(Vec),intent(in)::this
    integer :: n
    character(*), intent(in), optional :: string
    n = size(this%v, 1)
    write(*,*)
    if(present(string)) write(*,*) string
    write(*,'(10f10.4)') this%v(:)
  end subroutine VectorPrint

  subroutine GetRandomVector(v, n, dist)
    ! idist = 1: uniform (0, 1)
    ! idist = 2: uniform (-1, 1)
    ! idist = 3: normal (-1, 1)
  class(Vec), intent(inout) :: v
    integer, intent(in) :: n
    integer, intent(in), optional :: dist
    integer :: idist = 3
    if(present(dist)) idist = dist
    call v%ini(n)
    call dlarnv(idist, iseed, n, v%v)
  end subroutine GetRandomVector

!  type(MO) function OuterProduct(a, b) result(c)
!    type(Vec), intent(in) :: a, b
!    integer :: n, m
!    real(8) :: ddot
!    n = size(a%v)
!    m = size(b%v)
!    call c%ini(n,m)
!    call dger(n, m, 1.d0, a%v, 1, b%v, 1, c%m, n)
!  end function OuterProduct

!  subroutine iniM(a, m, n)
!  class(MO), intent(inout) :: a
!    integer, intent(in) :: m,n
!    if(allocated(a%m)) deallocate(a%m)
!    allocate(a%m(m,n))
!    a%m = 0.d0
!  end subroutine iniM
!
!  subroutine eye(a, n)
!  class(MO), intent(inout) :: a
!    integer, intent(in) :: n
!    integer :: i
!    if(allocated(a%m)) deallocate(a%m)
!    allocate(a%m(n,n))
!    a%m = 0.d0
!    do i = 1, n
!      a%m(i,i) = 1.d0
!    end do
!  end subroutine eye
!
!  subroutine FinM(a)
!  class(MO), intent(inout) :: a
!    if(allocated(a%m)) deallocate(a%m)
!  end subroutine FinM
!
!  subroutine MatrixCopy(b, a)
!    type(MO), intent(inout) :: b
!    type(MO), intent(in) :: a
!    integer :: m, n
!    m = size(a%m, 1)
!    n = size(a%m, 2)
!    call b%Ini(m,n)
!    b%m(:,:) = a%m(:,:)
!  end subroutine MatrixCopy
!
!  type(MO) function MatrixProduct(a, b) result(c)
!    type(MO), intent(in) :: a, b
!    integer :: m, k, n
!    m = size(a%m, 1)
!    k = size(a%m, 2)
!    if(size(a%m, 2) /= size(b%m, 1)) then
!      write(*, '(a)') 'Error in MatrixProduct'
!      stop
!    end if
!    n = size(b%m, 2)
!    call c%Ini(m,n)
!    call dgemm('n','n',m,n,k,1.d0,a%m,m,b%m,k,0.d0,c%m,m)
!  end function MatrixProduct
!
!  type(Vec) function MVProduct(a, b) result(c)
!    type(MO), intent(in) :: a
!    type(Vec), intent(in) :: b
!    integer :: m, k, n
!    m = size(a%m, 1)
!    k = size(a%m, 2)
!    if(size(a%m, 2) /= size(b%v)) then
!      write(*, '(a)') 'Error in MVProduct'
!      stop
!    end if
!    call c%Ini(m)
!    call dgemv('n',m,k,1.d0,a%m,m,b%v,1,0.d0,c%v,1)
!  end function MVProduct
!
!  type(Vec) function VMProduct(a, b) result(c)
!    type(MO), intent(in) :: b
!    type(Vec), intent(in) :: a
!    integer :: m, k, n
!    m = size(b%m, 1)
!    k = size(b%m, 2)
!    n = size(a%v, 1)
!    if(m /= n) then
!      write(*, '(a)') 'Error in VMProduct'
!      stop
!    end if
!    call c%Ini(k)
!    call dgemv('t',m,k,1.d0,b%m,m,a%v,1,0.d0,c%v,1)
!  end function VMProduct
!
!  subroutine VecDiagMat(b, a)
!  class(MO), intent(inout) :: b
!    type(Vec), intent(in) :: a
!    integer :: n, i
!    n = size(a%V)
!    call b%Ini(n,n)
!    do i = 1, n
!      b%M(i,i) = a%V(i)
!    end do
!  end subroutine VecDiagMat
!
!  subroutine DiagMatVec(b, a)
!  class(Vec), intent(inout) :: b
!    type(MO), intent(in) :: a
!    integer :: n, i
!    n = size(a%M,1)
!    call b%Ini(n)
!    do i = 1, n
!      b%V(i) = a%M(i,i)
!    end do
!  end subroutine DiagMatVec
!
!  type(MO) function MatrixSum(a, b) result(c)
!  class(MO), intent(in) :: a, b
!    integer :: m, n
!    if(size(a%m, 1) /= size(b%m, 1) .or. size(a%m, 2) /= size(b%m, 2)) then
!      write(*, '(a)') 'Error in MatrixSum'
!      stop
!    end if
!    m = size(a%m, 1)
!    n = size(a%m, 2)
!    call c%Ini(m,n)
!    c%m = a%m + b%m
!  end function MatrixSum
!
!  type(MO) function MatrixSubtract(a, b) result(c)
!  class(MO), intent(in) :: a, b
!    integer :: m, n
!    if(size(a%m, 1) /= size(b%m, 1) .or. size(a%m, 2) /= size(b%m, 2)) then
!      write(*, '(a)') 'Error in MatrixSum'
!      stop
!    end if
!    m = size(a%m, 1)
!    n = size(a%m, 2)
!    call c%Ini(m,n)
!    c%m = a%m - b%m
!  end function MatrixSubtract
!
!  type(MO) function MatrixScaleL(b, a) result(c)
!  class(MO), intent(in) :: b
!    real(8), intent(in) :: a
!    integer :: m, n
!    m = size(b%m, 1)
!    n = size(b%m, 2)
!    call c%Ini(m,n)
!    c%m(:,:) = a * b%m(:,:)
!  end function MatrixScaleL
!
!  type(MO) function MatrixScaleR(a, b) result(c)
!  class(MO), intent(in) :: b
!    real(8), intent(in) :: a
!    integer :: m, n
!    m = size(b%m, 1)
!    n = size(b%m, 2)
!    call c%Ini(m,n)
!    c%m(:,:) = a * b%m(:,:)
!  end function MatrixScaleR
!
!  type(MO) function Transepose(a) result(b)
!  class(MO), intent(in) :: a
!    integer :: n, m, i, j
!    m = size(a%m, 1)
!    n = size(a%m, 2)
!    call b%Ini(n,m)
!    b%M = transpose(a%M)
!  end function Transepose
!
!  subroutine DiagonalizationSymmetric(r, eval, evec, qmin, qmax, m, tol)
!  class(MO), intent(in) :: r
!    type(Vec), intent(out) :: eval
!    type(MO), intent(out) :: evec
!    real(8), optional, intent(in) :: qmin, qmax, tol
!    integer, optional, intent(in) :: m
!    integer :: num
!    real(8), allocatable :: work(:), mat(:,:), eig(:)
!    integer, allocatable :: iwork(:), ifailv(:)
!    integer :: info, lwork, n, i
!    real(8) :: lw, dlamch
!
!    n = size(r%M, 1)
!    call eval%Ini(n)
!    call evec%Ini(n,n)
!    allocate(eig(n)); eig(:) = 0.d0
!
!    if(.not. present(m) .and. .not. present(qmin) .and. &
!      & .not. present(qmax) .and. .not. present(tol)) then
!      allocate(mat(n,n))
!      mat = r%m
!      call dsyev('v', 'u', n, mat, n, eig, lw, -1, info)
!      lwork = int(lw)
!      allocate(work(lwork))
!      call dsyev('v', 'u', n, mat, n, eig, work, lwork, info)
!      evec%m = mat
!      deallocate(mat)
!    elseif(present(m)) then
!      allocate(iwork(5*n), ifailv(n))
!      call dsyevx('v','i','u',n,r%m,n,-1.d100,1.d100,1,n,dlamch('S'), &
!        &  num,eig,evec%m,n,lw,-1,iwork,ifailv,info)
!      deallocate( iwork, ifailv)
!    else
!      allocate(iwork(5*n), ifailv(n))
!      call dsyevx('v','i','u',n,r%m,n,-1.d100,1.d100,1,n,dlamch('S'), &
!        &  num,eig,evec%m,n,lw,-1,iwork,ifailv,info)
!      lwork = int(lw)
!      allocate(work(1:lwork))
!      call dsyevx('v','v','u',n,r%m,n,qmin,qmax,1,n,tol, &
!        &  num,eig,evec%m,n,work,lwork,iwork,ifailv,info)
!      deallocate( iwork, ifailv)
!    end if
!    eval%v = eig
!    deallocate(eig)
!  end subroutine DiagonalizationSymmetric
!
!  subroutine EigenValuesSymmetric(r, eval, m)
!  class(MO), intent(inout) :: r
!    type(Vec), intent(out) :: eval
!    integer, intent(in) :: m
!    integer, allocatable :: iwork(:), iblock(:), isplit(:)
!    real(8), allocatable :: work(:), mat(:,:), d(:), e(:), tau(:), w(:)
!    real(8) :: dlamch, lw
!    integer :: n, info, lwork, nsplit, i
!    n = size(r%m)
!    call eval%Ini(n)
!    allocate(mat(n,n))
!    allocate(d(n), e(max(1, n-1)), tau(max(1, n-1)), w(n))
!    allocate(iblock(n), isplit(n))
!    mat = r%m
!    call dsytrd('u',n,mat,n,d,e,tau,lw,-1,info)
!    lwork = int(lw)
!    allocate(work(lwork))
!    call dsytrd('u',n,mat,n,d,e,tau,work,lwork,info)
!    call dstebz('i','e',n,0.d0,0.d0,1,min(n,m),dlamch('S'), &
!      & d, e, m, nsplit, w, iblock, isplit, work, iwork, info)
!    do i = 1, min(n,m)
!      eval%v(i) = w(i)
!    end do
!    deallocate(work)
!    deallocate(mat, d, e, tau,iblock, isplit)
!  end subroutine EigenValuesSymmetric
!
!  type(MO) function inverse(r) result(s)
!  class(MO), intent(in) :: r
!    real(8), allocatable :: a(:,:)
!    real(8), allocatable :: work(:)
!    integer, allocatable :: ipvt(:)
!    integer :: info, n
!    n = size(r%m, 1)
!    call s%Ini(n,n)
!    allocate(work(n*n),ipvt(n))
!    allocate(a(n,n))
!    a = r%m
!    call dgetrf(n,n,a,n,ipvt,info)
!    call dgetri(n,a,n,ipvt,work,n**2,info)
!    s%m = a
!    deallocate(a,work,ipvt)
!  end function inverse
!
!  real(8) function Det(r) result(d)
!  class(MO), intent(in) :: r
!    integer :: n, i, info
!    real(8), allocatable :: a(:,:)
!    integer, allocatable :: ipiv(:)
!    n = size(r%m, 1)
!    allocate(ipiv(n), a(n,n))
!    a = r%m
!    call dgetrf(n, n, a, n, ipiv, info)
!    if(info /= 0) then
!      write(*,'(a, i3)') "error in det: info = ", info
!      stop
!    end if
!    d = 1.d0
!    do i = 1, n
!      if(ipiv(i) .ne. i) then
!        d = -d * a(i, i)
!      else
!        d = d * a(i, i)
!      end if
!    end do
!    deallocate(ipiv, a)
!  end function Det
!
!  subroutine MatrixPrint(this,char)
!  class(mo),intent(in)::this
!    character(12)::cfmt
!    integer::i,n,m
!    character(*),intent(in),optional::char
!    cfmt='( xf10.4)'
!    n=size(this%m,1)
!    m=size(this%m,2)
!    write(cfmt(2:3),'(I2)') m
!    write(*,*)
!    if(present(char))write(*,*)char
!    do i=1,n
!      write(*,cfmt) this%m(i,:)
!    end do
!  end subroutine MatrixPrint
!
!  subroutine GetRandomMatrix(mat, m, n)
!  class(MO), intent(inout) :: mat
!    integer, intent(in) :: m, n
!    integer :: i
!    type(Vec) :: v
!    call mat%ini(m,n)
!    do i = 1, n
!      call v%GetRandomVector(m)
!      mat%m(:,i) = v%v(:)
!      call v%fin()
!    end do
!  end subroutine GetRandomMatrix

end module DMatVec
