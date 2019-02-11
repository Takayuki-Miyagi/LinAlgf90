module MatrixComplex
  use LinAlgParameters
  implicit none
  public :: CMat
  public :: MatrixCopyC
  public :: MatrixProductC
  public :: MatrixSumC
  public :: MatrixSubtractC
  public :: MatrixScaleLC
  public :: MatrixScaleRC
  public :: MatrixScaleDivideC

  private :: IniM
  private :: FinM
  private :: zeros
  private :: eye
  private :: Trans
  private :: ComplexConjugate
  private :: HermiteConjugate
  private :: Inverse
  private :: Det
  private :: GetRandomMatrix
  private :: MatrixPrint
  private :: DiagMat
  private :: block_cmat

  type :: CMat
    complex(dp), allocatable :: M(:,:)
  contains
    procedure :: Ini => IniM
    procedure :: Fin => FinM
    procedure :: zeros
    procedure :: eye
    procedure :: T => Trans
    procedure :: C => ComplexConjugate
    procedure :: H => HermiteConjugate
    procedure :: Inv => Inverse
    procedure :: blk => block_cmat
    procedure :: Det
    procedure :: Random => GetRandomMatrix
    procedure :: prt => MatrixPrint
    procedure :: DiagMat
  end type CMat
contains
  subroutine iniM(a, m, n)
    class(CMat), intent(inout) :: a
    integer(kp), intent(in) :: m,n
    if(m < 1 .or. n < 1) return
    if(allocated(a%m)) deallocate(a%m)
    allocate(a%m(m,n))
  end subroutine iniM

  subroutine zeros(a, m, n)
    class(CMat), intent(inout) :: a
    integer(kp), intent(in) :: m,n
    if(m < 1 .or. n < 1) return
    if(allocated(a%m)) deallocate(a%m)
    allocate(a%m(m,n))
    a%m = (0.d0,0.d0)
  end subroutine zeros

  subroutine eye(a, n)
    class(CMat), intent(inout) :: a
    integer(kp), intent(in) :: n
    integer(kp) :: i
    if(n < 1) return
    if(allocated(a%m)) deallocate(a%m)
    allocate(a%m(n,n))
    a%m = 0.d0
    do i = 1, n
      a%m(i,i) = (1.d0, 0.d0)
    end do
  end subroutine eye

  subroutine FinM(a)
    class(CMat), intent(inout) :: a
    if(allocated(a%m)) deallocate(a%m)
  end subroutine FinM

  subroutine MatrixCopyC(b, a)
    type(CMat), intent(inout) :: b
    type(CMat), intent(in) :: a
    integer(kp) :: m, n, i
    m = size(a%m, 1)
    n = size(a%m, 2)
    if(m < 1 .or. n < 1) return
    call b%Ini(m,n)
    do i = 1, n
      call zcopy(m, a%m(:,i), 1, b%m(:,i), 1)
    end do
  end subroutine MatrixCopyC

  type(CMat) function MatrixProductC(a, b) result(c)
    type(CMat), intent(in) :: a, b
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
    call zgemm('n','n',m,n,k,(1.d0,0.d0),a%m,m,b%m,k,(0.d0,0.d0),c%m,m)
  end function MatrixProductC

  type(CMat) function MatrixSumC(a, b) result(c)
    type(CMat), intent(in) :: a, b
    integer(kp) :: m, n, i
    if(size(a%m, 1) /= size(b%m, 1) .or. size(a%m, 2) /= size(b%m, 2)) then
      write(*, '(a)') 'Error in MatrixSum'
      stop
    end if
    m = size(a%m, 1)
    n = size(a%m, 2)
    if(m < 1 .or. n < 1) return
    call MatrixCopyC(c, a)
    do i = 1, n
      call zaxpy(m, (1.d0,0.d0), b%m(:,i), 1, c%m(:,i), 1)
    end do
  end function MatrixSumC

  type(CMat) function MatrixSubtractC(a, b) result(c)
    type(CMat), intent(in) :: a, b
    integer(kp) :: m, n, i
    if(size(a%m, 1) /= size(b%m, 1) .or. size(a%m, 2) /= size(b%m, 2)) then
      write(*, '(a)') 'Error in MatrixSum'
      stop
    end if
    m = size(a%m, 1)
    n = size(a%m, 2)
    if(m < 1 .or. n < 1) return
    call MatrixCopyC(c, a)
    do i = 1, n
      call zaxpy(m, (-1.d0,0.d0), b%m(:,i), 1, c%m(:,i), 1)
    end do
  end function MatrixSubtractC

  type(CMat) function MatrixScaleLC(b, a) result(c)
    type(CMat), intent(in) :: b
    complex(dp), intent(in) :: a
    integer(kp) :: m, n, i
    m = size(b%m, 1)
    n = size(b%m, 2)
    if(m < 1 .or. n < 1) return
    call MatrixCopyC(c, b)
    do i = 1, n
      call zscal(m, a, c%m(:,i), 1)
    end do
  end function MatrixScaleLC

  type(CMat) function MatrixScaleRC(a, b) result(c)
    type(CMat), intent(in) :: b
    complex(dp), intent(in) :: a
    integer(kp) :: m, n, i
    m = size(b%m, 1)
    n = size(b%m, 2)
    if(m < 1 .or. n < 1) return
    call MatrixCopyC(c, b)
    do i = 1, n
      call zscal(m, a, c%m(:,i), 1)
    end do
  end function MatrixScaleRC

  type(CMat) function MatrixScaleDivideC(b, a) result(c)
    type(CMat), intent(in) :: b
    real(dp), intent(in) :: a
    complex(dp) :: aa
    integer(kp) :: m, n, i
    m = size(b%m, 1)
    n = size(b%m, 2)
    if(m < 1 .or. n < 1) return
    call MatrixCopyC(c, b)
    aa = 1.d0 / a
    do i = 1, n
      call zscal(m, aa, c%m(:,i), 1)
    end do
  end function MatrixScaleDivideC

  type(CMat) function Trans(a) result(b)
    class(CMat), intent(in) :: a
    integer(kp) :: n, m
    m = size(a%m, 1)
    n = size(a%m, 2)
    if(m < 1 .or. n < 1) return
    call b%Ini(n,m)
    b%M = transpose(a%M)
  end function Trans

  type(CMat) function ComplexConjugate(a) result(b)
    class(CMat), intent(in) :: a
    integer(kp) :: n, m
    m = size(a%m, 1)
    n = size(a%m, 2)
    if(m < 1 .or. n < 1) return
    call b%Ini(n,m)
    b%M = conjg(a%M)
  end function ComplexConjugate

  type(CMat) function HermiteConjugate(a) result(b)
    class(CMat), intent(in) :: a
    b = a%C()
    b = b%T()
  end function HermiteConjugate

  type(CMat) function inverse(r) result(s)
    class(CMat), intent(in) :: r
    complex(dp), allocatable :: a(:,:)
    complex(dp), allocatable :: work(:)
    integer(kp), allocatable :: ipvt(:)
    integer(kp) :: info, n
    n = size(r%m, 1)
    if(n < 1) return
    call s%Ini(n,n)
    allocate(work(n*n),ipvt(n))
    allocate(a(n,n))
    a = r%m
    call zgetrf(n,n,a,n,ipvt,info)
    call zgetri(n,a,n,ipvt,work,n**2,info)
    s%m = a
    deallocate(a,work,ipvt)
  end function inverse

  complex(dp) function Det(r) result(d)
    class(CMat), intent(in) :: r
    integer(kp) :: n, i, info
    complex(dp), allocatable :: a(:,:)
    integer(kp), allocatable :: ipiv(:)
    n = size(r%m, 1)
    if(n < 1) return
    allocate(ipiv(n), a(n,n))
    a = r%m
    call zgetrf(n, n, a, n, ipiv, info)
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

  subroutine MatrixPrint(this, msg, iunit, binary)
    class(CMat), intent(in) :: this
    integer, intent(in), optional :: iunit
    character(*), intent(in), optional :: msg
    logical, intent(in), optional :: binary
    logical :: bin
    character(12) :: cfmt
    integer(kp) :: i, j, n, m, unt

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
        write(unt,'(a)') 'Real:'
        do i=1,n
          write(unt,cfmt) dble(this%m(i,:))
        end do
        write(unt,'(a)') 'Imag:'
        do i=1,n
          write(unt,cfmt) aimag(this%m(i,:))
        end do
      else
        do i = 1, n
          do j = 1, m
            write(unt,'(2i8,2f16.6)') i,j,this%m(i,j)
          end do
        end do
      end if
    end if
  end subroutine MatrixPrint

  function block_cmat(this, m1, m2, n1, n2) result(r)
    class(CMat), intent(in) :: this
    type(CMat) :: r
    integer(kp), intent(in) :: m1, m2, n1, n2
    integer(kp) :: m, n
    m = m2 - m1 + 1
    n = n2 - n1 + 1
    if(m < 1 .or. n < 1) return
    call r%ini(m,n)
    r%m(:,:) = this%m(m1:m2,n1:n2)
  end function block_cmat

  subroutine GetRandomMatrix(mat, m, n)
    use VectorComplex, only: CVec
    class(CMat), intent(inout) :: mat
    integer(kp), intent(in) :: m, n
    integer(kp) :: i
    type(CVec) :: v
    if(m < 1 .or. n < 1) return
    call mat%ini(m,n)
    do i = 1, n
      call v%Random(m)
      mat%m(:,i) = v%v(:)
      call v%fin()
    end do
  end subroutine GetRandomMatrix

  subroutine DiagMat(b, a)
    use VectorComplex, only: CVec
    class(CMat), intent(inout) :: b
    type(CVec), intent(in) :: a
    integer(kp) :: n, i
    n = size(a%V)
    if(n < 1) return
    call b%Ini(n,n)
    do i = 1, n
      b%M(i,i) = a%V(i)
    end do
  end subroutine DiagMat
end module MatrixComplex
