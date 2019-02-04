module MatrixSingle
  use LinAlgParameters
  implicit none

  public :: SMat
  public :: MatrixCopyS
  public :: MatrixProductS
  public :: MatrixSumS
  public :: MatrixSubtractS
  public :: MatrixScaleLS
  public :: MatrixScaleRS
  public :: MatrixScaleDivideS

  private :: IniM
  private :: zeros
  private :: FinM
  private :: eye
  private :: Trans
  private :: Inverse
  private :: Det
  private :: GetRandomMatrix
  private :: MatrixPrint
  private :: DiagMat
  private :: block_SMat

  type :: SMat
    real(sp), allocatable :: M(:,:)
  contains
    procedure :: Ini => IniM
    procedure :: Fin => FinM
    procedure :: zeros
    procedure :: eye
    procedure :: T => Trans
    procedure :: Inv => Inverse
    procedure :: Det
    procedure :: Random => GetRandomMatrix
    procedure :: blk => block_SMat
    procedure :: prt => MatrixPrint
    procedure :: DiagMat
  end type SMat
contains
  subroutine iniM(a, m, n)
    class(SMat), intent(inout) :: a
    integer(kp), intent(in) :: m,n
    if(m < 1 .or. n < 1) return
    if(allocated(a%m)) deallocate(a%m)
    allocate(a%m(m,n))
  end subroutine iniM

  subroutine zeros(a, m, n)
    class(SMat), intent(inout) :: a
    integer(kp), intent(in) :: m,n
    if(m < 1 .or. n < 1) return
    call a%ini(m,n)
    a%m = 0.0
  end subroutine zeros

  subroutine eye(a, n)
    class(SMat), intent(inout) :: a
    integer(kp), intent(in) :: n
    integer(kp) :: i
    if(n < 1) return
    if(allocated(a%m)) deallocate(a%m)
    allocate(a%m(n,n))
    a%m = 0.0
    do i = 1, n
      a%m(i,i) = 1.0
    end do
  end subroutine eye

  subroutine FinM(a)
    class(SMat), intent(inout) :: a
    if(allocated(a%m)) deallocate(a%m)
  end subroutine FinM

  subroutine MatrixCopyS(b, a)
    type(SMat), intent(inout) :: b
    type(SMat), intent(in) :: a
    integer(kp) :: m, n, i
    m = size(a%m, 1)
    n = size(a%m, 2)
    if(m < 1 .or. n < 1) return
    call b%Ini(m,n)
    do i = 1, n
      call dcopy(m, a%m(:,i), 1, b%m(:,i), 1)
    end do
  end subroutine MatrixCopyS

  type(SMat) function MatrixProductS(a, b) result(c)
    type(SMat), intent(in) :: a, b
    integer(kp) :: m, k, n
    m = size(a%m, 1)
    k = size(a%m, 2)
    if(size(a%m, 2) /= size(b%m, 1)) then
      write(*, '(a)') 'Error in MatrixProduct'
      stop
    end if
    n = size(b%m, 2)
    if(m < 1 .or. n < 1) return
    call c%Ini(m,n)
    call sgemm('n','n',m,n,k,1.0,a%m,m,b%m,k,0.0,c%m,m)
  end function MatrixProductS

  type(SMat) function MatrixSumS(a, b) result(c)
    type(SMat), intent(in) :: a, b
    integer(kp) :: m, n, i
    if(size(a%m, 1) /= size(b%m, 1) .or. size(a%m, 2) /= size(b%m, 2)) then
      write(*, '(a)') 'Error in MatrixSum'
      stop
    end if
    m = size(a%m, 1)
    n = size(a%m, 2)
    if(m < 1 .or. n < 1) return
    call MatrixCopyS(c, a)
    do i = 1, n
      call saxpy(m, 1.0, b%m(:,i), 1, c%m(:,i), 1)
    end do
  end function MatrixSumS

  type(SMat) function MatrixSubtractS(a, b) result(c)
    type(SMat), intent(in) :: a, b
    integer(kp) :: m, n, i
    if(size(a%m, 1) /= size(b%m, 1) .or. size(a%m, 2) /= size(b%m, 2)) then
      write(*, '(a)') 'Error in MatrixSum'
      stop
    end if
    m = size(a%m, 1)
    n = size(a%m, 2)
    if(m < 1 .or. n < 1) return
    call MatrixCopyS(c, a)
    do i = 1, n
      call saxpy(m, -1.0, b%m(:,i), 1, c%m(:,i), 1)
    end do
  end function MatrixSubtractS

  type(SMat) function MatrixScaleLS(b, a) result(c)
    type(SMat), intent(in) :: b
    real(sp), intent(in) :: a
    integer(kp) :: m, n, i
    m = size(b%m, 1)
    n = size(b%m, 2)
    if(m < 1 .or. n < 1) return
    call MatrixCopyS(c, b)
    do i = 1, n
      call sscal(m, a, c%m(:,i), 1)
    end do
  end function MatrixScaleLS

  type(SMat) function MatrixScaleRS(a, b) result(c)
    type(SMat), intent(in) :: b
    real(sp), intent(in) :: a
    integer(kp) :: m, n, i
    m = size(b%m, 1)
    n = size(b%m, 2)
    if(m < 1 .or. n < 1) return
    call MatrixCopyS(c, b)
    do i = 1, n
      call sscal(m, a, c%m(:,i), 1)
    end do
  end function MatrixScaleRS

  type(SMat) function MatrixScaleDivideS(b, a) result(c)
    type(SMat), intent(in) :: b
    real(sp), intent(in) :: a
    integer(kp) :: m, n, i
    m = size(b%m, 1)
    n = size(b%m, 2)
    if(m < 1 .or. n < 1) return
    call MatrixCopyS(c, b)
    do i = 1, n
      call sscal(m, 1.0 / a, c%m(:,i), 1)
    end do
  end function MatrixScaleDivideS

  type(SMat) function Trans(a) result(b)
    class(SMat), intent(in) :: a
    integer(kp) :: n, m
    m = size(a%m, 1)
    n = size(a%m, 2)
    if(m < 1 .or. n < 1) return
    call b%Ini(n,m)
    b%M = transpose(a%M)
  end function Trans

  type(SMat) function inverse(r) result(s)
    class(SMat), intent(in) :: r
    real(sp), allocatable :: a(:,:)
    real(sp), allocatable :: work(:)
    integer(kp), allocatable :: ipvt(:)
    integer(kp) :: info, n
    n = size(r%m, 1)
    if(n < 1) return
    call s%Ini(n,n)
    allocate(work(n*n),ipvt(n))
    allocate(a(n,n))
    a = r%m
    call sgetrf(n,n,a,n,ipvt,info)
    call sgetri(n,a,n,ipvt,work,n**2,info)
    s%m = a
    deallocate(a,work,ipvt)
  end function inverse

  real(sp) function Det(r) result(d)
    class(SMat), intent(in) :: r
    integer(kp) :: n, i, info
    real(sp), allocatable :: a(:,:)
    integer(kp), allocatable :: ipiv(:)
    n = size(r%m, 1)
    d = 0.0
    if(n < 1) return
    allocate(ipiv(n), a(n,n))
    a = r%m
    call sgetrf(n, n, a, n, ipiv, info)
    if(info /= 0) then
      write(*,'(a, i3)') "error in det: info = ", info
      stop
    end if
    d = 1.0
    do i = 1, n
      if(ipiv(i) .ne. i) then
        d = -d * a(i, i)
      else
        d = d * a(i, i)
      end if
    end do
    deallocate(ipiv, a)
  end function Det

  subroutine MatrixPrint(this, msg, iunit, binary)
    class(SMat), intent(in) :: this
    integer, intent(in), optional :: iunit
    character(*), intent(in), optional :: msg
    logical, intent(in), optional :: binary
    logical :: bin
    character(12) :: cfmt
    integer(kp) :: i, j, n, m
    integer :: unt
    if(present(iunit)) then; unt = iunit
    else; unt = 6; end if

    if(present(binary)) then; bin = binary
    else; bin = .false.; end if

    if(bin) then
      write(unt) this%m
    else
      if(present(msg)) write(unt,*) msg
      n = size(this%m, 1)
      m = size(this%m, 2)
      if(unt == 6) then
        cfmt = '( xf10.4)'
        write(cfmt(2:3), '(I2)') m
        do i=1,n
          write(unt,cfmt) this%m(i,:)
        end do
      else
        do i=1,n
          do j=1,m
            write(unt,'(2i8,f14.6)') i,j,this%m(i,j)
          end do
        end do
      end if
    end if
  end subroutine MatrixPrint

  function block_SMat(this, m1, m2, n1, n2) result(r)
    class(SMat), intent(in) :: this
    type(SMat) :: r
    integer(kp), intent(in) :: m1, m2, n1, n2
    integer(kp) :: m, n
    m = m2 - m1 + 1
    n = n2 - n1 + 1
    if(m < 1 .or. n < 1) return
    call r%ini(m,n)
    r%m(:,:) = this%m(m1:m2,n1:n2)
  end function block_SMat

  subroutine GetRandomMatrix(mat, m, n)
    use VectorSingle, only: SVec
    class(SMat), intent(inout) :: mat
    integer(kp), intent(in) :: m, n
    integer(kp) :: i
    type(SVec) :: v
    call mat%ini(m,n)
    if(m < 1 .or. n < 1) return
    do i = 1, n
      call v%Random(m)
      mat%m(:,i) = v%v(:)
      call v%fin()
    end do
  end subroutine GetRandomMatrix

  subroutine DiagMat(b, a)
    use VectorSingle, only: SVec
    class(SMat), intent(inout) :: b
    type(SVec), intent(in) :: a
    integer(kp) :: n, i
    n = size(a%V)
    if(n < 1) return
    call b%Ini(n,n)
    do i = 1, n
      b%M(i,i) = a%V(i)
    end do
  end subroutine DiagMat
end module MatrixSingle
