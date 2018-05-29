module LinearAlgebra
  implicit none
  integer :: iseed(4) = (/3239, 4241, 1903, 1093/)
  type :: VO
    real(8), allocatable :: V(:)
  contains
    procedure :: Ini => iniV, Fin => FinV
    procedure :: InnerProduct   ! .x.   c = a .x. b  (a, b: Vector, c: real)
    procedure :: VectorSum      !  +    c = a  +  b  (a, b: Matrix)
    procedure :: VectorScale    !  *    c = a  *  b  (a: Matrix, b: scalar)
    procedure :: VectorSubtract !  -    c = a  -  b  (a, b: Matrix)
    procedure :: VectorCopy     !  =    b = a        (a: Matrix)
    procedure :: prt => VectorPrint
    procedure :: GetRandomVector
  end type VO

  type :: MO
    real(8), allocatable :: M(:,:)
  contains
    procedure :: Ini => IniM, Fin => FinM
    procedure :: eye
    procedure :: T => Transepose
    procedure :: DiagS => DiagonalizationSymmetric
    procedure :: EvalS => EigenValuesSymmetric
    procedure :: Inverse
    procedure :: Det
    procedure :: GetRandomMatrix
    procedure :: prt => MatrixPrint

    procedure :: MatrixProduct ! .x.   c = a .x. b  (a, b: Matrix)
    procedure :: MatrixSum     !  +    c = a  +  b  (a, b: Matrix)
    procedure :: MatrixScale   !  *    c = a  *  b  (a: Matrix, b: scalar)
    procedure :: MatrixSubtract!  -    c = a  -  b  (a, b: Matrix)
    procedure :: MatrixCopy    !  =    b = a        (a: Matrix)
  end type MO

  interface assignment(=)
    procedure :: VectorCopy
    procedure :: MatrixCopy
    procedure :: VecDiagMat
    procedure :: DiagMatVec
  end interface assignment(=)

  interface operator(+)
    procedure :: VectorSum
    procedure :: MatrixSum
  end interface operator(+)

  interface operator(-)
    procedure :: VectorSubtract
    procedure :: MatrixSubtract
  end interface operator(-)

  interface operator(*)
    procedure :: VectorScale
    procedure :: MatrixScale
    procedure :: InnerProduct
  end interface operator(*)

  interface operator(.x.)
    procedure :: OuterProduct
    procedure :: MatrixProduct
    procedure :: MVProduct
    procedure :: VMProduct
  end interface operator(.x.)
contains
  subroutine IniV(a, n)
  class(VO), intent(inout) :: a
    integer, intent(in) :: n
    if(allocated(a%V)) deallocate(a%V)
    allocate(a%V(n))
    a%V(:) = 0.d0
  end subroutine IniV

  subroutine FinV(a)
  class(VO), intent(inout) :: a
    if(allocated(a%V)) deallocate(a%V)
  end subroutine FinV

  subroutine VectorCopy(b, a)
  class(VO), intent(inout) :: b
    type(VO), intent(in) :: a
    integer :: n
    n = size(a%V)
    call b%Ini(n)
    b%V = a%V
  end subroutine VectorCopy

  type(VO) function VectorSum(a, b) result(c)
  class(VO), intent(in) :: a, b
    integer :: n
    if(size(a%v) /= size(b%v)) then
      write(*,'(a)') 'Error in VectorSum'
      stop
    end if
    n = size(a%v)
    call c%ini(n)
    c%v = a%v + b%v
  end function VectorSum

  type(VO) function VectorSubtract(a, b) result(c)
  class(VO), intent(in) :: a, b
    integer :: n
    if(size(a%v) /= size(b%v)) then
      write(*,'(a)') 'Error in VectorSubtract'
      stop
    end if
    n = size(a%v)
    call c%ini(n)
    c%v = a%v - b%v
  end function VectorSubtract

  type(VO) function VectorScale(a, b) result(c)
  class(VO), intent(in) :: a
    real(8), intent(in) :: b
    integer :: n
    n = size(a%v)
    call c%ini(n)
    c%v = b * a%v
  end function VectorScale

  real(8) function InnerProduct(a, b) result(c)
  class(VO), intent(in) :: a, b
    integer :: n
    real(8) :: ddot
    if(size(a%v) /= size(b%v)) then
      write(*,'(a)') 'Error in InnerProduct'
      stop
    end if
    n = size(a%v)
    c = ddot(n, a%v, 1, b%v, 1)
  end function InnerProduct

  type(MO) function OuterProduct(a, b) result(c)
  class(VO), intent(in) :: a, b
    integer :: n, m
    real(8) :: ddot
    n = size(a%v)
    m = size(b%v)
    call c%ini(n,m)
    call dger(n, m, 1.d0, a%v, 1, b%v, 1, c%m, n)
  end function OuterProduct

  subroutine iniM(a, m, n)
  class(MO), intent(inout) :: a
    integer, intent(in) :: m,n
    if(allocated(a%m)) deallocate(a%m)
    allocate(a%m(m,n))
    a%m = 0.d0
  end subroutine iniM

  subroutine eye(a, n)
  class(MO), intent(inout) :: a
    integer, intent(in) :: n
    integer :: i
    if(allocated(a%m)) deallocate(a%m)
    allocate(a%m(n,n))
    a%m = 0.d0
    do i = 1, n
      a%m(i,i) = 1.d0
    end do
  end subroutine eye

  subroutine FinM(a)
  class(MO), intent(inout) :: a
    if(allocated(a%m)) deallocate(a%m)
  end subroutine FinM

  subroutine MatrixCopy(b, a)
  class(MO), intent(inout) :: b
    type(MO), intent(in) :: a
    integer :: m, n
    m = size(a%m, 1)
    n = size(a%m, 2)
    call b%Ini(m,n)
    b%m(:,:) = a%m(:,:)
  end subroutine MatrixCopy

  type(MO) function MatrixProduct(a, b) result(c)
  class(MO), intent(in) :: a, b
    integer :: m, k, n
    m = size(a%m, 1)
    k = size(a%m, 2)
    if(size(a%m, 2) /= size(b%m, 1)) then
      write(*, '(a)') 'Error in MatrixProduct'
      stop
    end if
    n = size(b%m, 2)
    call c%Ini(m,n)
    call dgemm('n','n',m,n,k,1.d0,a%m,m,b%m,k,0.d0,c%m,m)
  end function MatrixProduct

  type(VO) function MVProduct(a, b) result(c)
    type(MO), intent(in) :: a
    type(VO), intent(in) :: b
    integer :: m, k, n
    m = size(a%m, 1)
    k = size(a%m, 2)
    if(size(a%m, 2) /= size(b%v)) then
      write(*, '(a)') 'Error in MVProduct'
      stop
    end if
    call c%Ini(m)
    call dgemv('n',m,k,1.d0,a%m,m,b%v,1,0.d0,c%v,1)
  end function MVProduct

  type(VO) function VMProduct(a, b) result(c)
    type(MO), intent(in) :: b
    type(VO), intent(in) :: a
    integer :: m, k, n
    m = size(b%m, 1)
    k = size(b%m, 2)
    n = size(a%v, 1)
    if(m /= n) then
      write(*, '(a)') 'Error in VMProduct'
      stop
    end if
    call c%Ini(k)
    call dgemv('t',m,k,1.d0,b%m,m,a%v,1,0.d0,c%v,1)
  end function VMProduct

  subroutine VecDiagMat(b, a)
  class(MO), intent(inout) :: b
    type(VO), intent(in) :: a
    integer :: n, i
    n = size(a%V)
    call b%Ini(n,n)
    do i = 1, n
      b%M(i,i) = a%V(i)
    end do
  end subroutine VecDiagMat

  subroutine DiagMatVec(b, a)
  class(VO), intent(inout) :: b
    type(MO), intent(in) :: a
    integer :: n, i
    n = size(a%M,1)
    call b%Ini(n)
    do i = 1, n
      b%V(i) = a%M(i,i)
    end do
  end subroutine DiagMatVec

  type(MO) function MatrixSum(a, b) result(c)
  class(MO), intent(in) :: a, b
    integer :: m, n
    if(size(a%m, 1) /= size(b%m, 1) .or. size(a%m, 2) /= size(b%m, 2)) then
      write(*, '(a)') 'Error in MatrixSum'
      stop
    end if
    m = size(a%m, 1)
    n = size(a%m, 2)
    call c%Ini(m,n)
    c%m = a%m + b%m
  end function MatrixSum

  type(MO) function MatrixSubtract(a, b) result(c)
  class(MO), intent(in) :: a, b
    integer :: m, n
    if(size(a%m, 1) /= size(b%m, 1) .or. size(a%m, 2) /= size(b%m, 2)) then
      write(*, '(a)') 'Error in MatrixSum'
      stop
    end if
    m = size(a%m, 1)
    n = size(a%m, 2)
    call c%Ini(m,n)
    c%m = a%m - b%m
  end function MatrixSubtract

  type(MO) function MatrixScale(b, a) result(c)
  class(MO), intent(in) :: b
    real(8), intent(in) :: a
    integer :: m, n
    m = size(b%m, 1)
    n = size(b%m, 2)
    call c%Ini(m,n)
    c%m(:,:) = a * b%m(:,:)
  end function MatrixScale

  type(MO) function Transepose(a) result(b)
  class(MO), intent(in) :: a
    integer :: n, m, i, j
    m = size(a%m, 1)
    n = size(a%m, 2)
    call b%Ini(n,m)
    do i = 1, n
      do j = 1, m
        b%M(i,j) = a%M(j,i)
      end do
    end do
  end function Transepose

  subroutine DiagonalizationSymmetric(r, eval, evec, qmin, qmax, m, tol)
  class(MO), intent(in) :: r
    type(VO), intent(out) :: eval
    type(MO), intent(out) :: evec
    real(8), optional, intent(in) :: qmin, qmax, tol
    integer, optional, intent(in) :: m
    integer :: num
    real(8), allocatable :: work(:), mat(:,:), eig(:)
    integer, allocatable :: iwork(:), ifailv(:)
    integer :: info, lwork, n, i
    real(8) :: lw, dlamch

    n = size(r%M, 1)
    call eval%Ini(n)
    call evec%Ini(n,n)
    allocate(eig(n)); eig(:) = 0.d0

    if(.not. present(m) .and. .not. present(qmin) .and. &
      & .not. present(qmax) .and. .not. present(tol)) then
      allocate(mat(n,n))
      mat = r%m
      call dsyev('v', 'u', n, mat, n, eig, lw, -1, info)
      lwork = int(lw)
      allocate(work(lwork))
      call dsyev('v', 'u', n, mat, n, eig, work, lwork, info)
      evec%m = mat
      deallocate(mat)
    elseif(present(m)) then
      allocate(iwork(5*n), ifailv(n))
      call dsyevx('v','i','u',n,r%m,n,-1.d100,1.d100,1,n,dlamch('S'), &
        &  num,eig,evec%m,n,lw,-1,iwork,ifailv,info)
      deallocate( iwork, ifailv)
    else
      allocate(iwork(5*n), ifailv(n))
      call dsyevx('v','i','u',n,r%m,n,-1.d100,1.d100,1,n,dlamch('S'), &
        &  num,eig,evec%m,n,lw,-1,iwork,ifailv,info)
      lwork = int(lw)
      allocate(work(1:lwork))
      call dsyevx('v','v','u',n,r%m,n,qmin,qmax,1,n,tol, &
        &  num,eig,evec%m,n,work,lwork,iwork,ifailv,info)
      deallocate( iwork, ifailv)
    end if
    eval%v = eig
    deallocate(eig)
  end subroutine DiagonalizationSymmetric

  subroutine EigenValuesSymmetric(r, eval, m)
  class(MO), intent(inout) :: r
    type(VO), intent(out) :: eval
    integer, intent(in) :: m
    integer, allocatable :: iwork(:), iblock(:), isplit(:)
    real(8), allocatable :: work(:), mat(:,:), d(:), e(:), tau(:), w(:)
    real(8) :: dlamch, lw
    integer :: n, info, lwork, nsplit, i
    n = size(r%m)
    call eval%Ini(n)
    allocate(mat(n,n))
    allocate(d(n), e(max(1, n-1)), tau(max(1, n-1)), w(n))
    allocate(iblock(n), isplit(n))
    mat = r%m
    call dsytrd('u',n,mat,n,d,e,tau,lw,-1,info)
    lwork = int(lw)
    allocate(work(lwork))
    call dsytrd('u',n,mat,n,d,e,tau,work,lwork,info)
    call dstebz('i','e',n,0.d0,0.d0,1,min(n,m),dlamch('S'), &
      & d, e, m, nsplit, w, iblock, isplit, work, iwork, info)
    do i = 1, min(n,m)
      eval%v(i) = w(i)
    end do
    deallocate(work)
    deallocate(mat, d, e, tau,iblock, isplit)
  end subroutine EigenValuesSymmetric

  type(MO) function inverse(r) result(s)
  class(MO), intent(in) :: r
    real(8), allocatable :: a(:,:)
    real(8), allocatable :: work(:)
    integer, allocatable :: ipvt(:)
    integer :: info, n
    n = size(r%m, 1)
    call s%Ini(n,n)
    allocate(work(n*n),ipvt(n))
    allocate(a(n,n))
    a = r%m
    call dgetrf(n,n,a,n,ipvt,info)
    call dgetri(n,a,n,ipvt,work,n**2,info)
    s%m = a
    deallocate(a,work,ipvt)
  end function inverse

  real(8) function Det(r) result(d)
  class(MO), intent(in) :: r
    integer :: n, i, info
    real(8), allocatable :: a(:,:)
    integer, allocatable :: ipiv(:)
    n = size(r%m, 1)
    allocate(ipiv(n), a(n,n))
    a = r%m
    call dgetrf(n, n, a, n, ipiv, info)
    if(info /= 0) then
      write(*,'(a, i3)') "error in det: info = ", info
      stop
    end if
    d = 1.d0
    do i = 1, n
      if(ipiv(i) .ne. i) then
        d = -d * a(i, i)
      else
        d = d * a(i, i)
      end if
    end do
    deallocate(ipiv, a)
  end function Det

  subroutine VectorPrint(this,char)
  class(vo),intent(in)::this
    integer::i,n,m
    character(*),intent(in),optional::char
    n=size(this%v,1)
    write(*,*)
    if(present(char))write(*,*)char
    write(*,'(10f10.4)') this%v(:)
  end subroutine VectorPrint

  subroutine MatrixPrint(this,char)
  class(mo),intent(in)::this
    character(12)::cfmt
    integer::i,n,m
    character(*),intent(in),optional::char
    cfmt='( xf10.4)'
    n=size(this%m,1)
    m=size(this%m,2)
    write(cfmt(2:3),'(I2)') m
    write(*,*)
    if(present(char))write(*,*)char
    do i=1,n
      write(*,cfmt) this%m(i,:)
    end do
  end subroutine MatrixPrint

  subroutine GetRandomVector(v, n)
  ! idist = 1: uniform (0, 1)
  ! idist = 2: uniform (-1, 1)
  ! idist = 3: normal (-1, 1)
  class(VO), intent(inout) :: v
    integer, intent(in) :: n
    integer, parameter :: idist = 3
    call v%ini(n)
    call dlarnv(idist, iseed, n, v%v)
  end subroutine GetRandomVector

  subroutine GetRandomMatrix(mat, m, n)
  class(MO), intent(inout) :: mat
    integer, intent(in) :: m, n
    integer :: i
    type(VO) :: v
    call mat%ini(m,n)
    do i = 1, n
      call v%GetRandomVector(m)
      mat%m(:,i) = v%v(:)
      call v%fin()
    end do
  end subroutine GetRandomMatrix

end module LinearAlgebra

!program test
!  use LinearAlgebra, only: MO, VO,  &
!      &   assignment(=), &
!      &   operator(+), &
!      &   operator(-), &
!      &   operator(*), &
!      &   operator(.x.)
!  type(MO) :: mat1, mat2, mat3
!  type(VO) :: vec1, vec2, vec3
!  integer :: m = 2, n = 3
!  call vec1%ini(n)
!  call vec2%ini(m)
!  vec1%v(:) = (/2.d0, 1.d0, 3.d0/)
!  vec2%v(:) = (/2.d0, 1.d0/)
!
!  mat1 = vec1 .x. vec2
!  mat2 = vec2 .x. vec1
!
!  call mat1%prt('mat1')
!  call mat2%prt('mat2')
!
!end program test
