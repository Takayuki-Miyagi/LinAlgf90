module MatrixDouble
  use LinAlgParameters
  implicit none
  private :: IniM, FinM, eye, Trans, Inverse, Det
  private :: GetRandomMatrix, MatrixPrint
  private :: DiagMat, block_dmat

  public :: DMat, MatrixCopyD, MatrixProductD, MatrixSumD
  public :: MatrixSubtractD, MatrixScaleLD, MatrixScaleRD
  public :: MatrixScaleDivideD

  type :: DMat
    real(dp), allocatable :: M(:,:)
  contains
    procedure :: Ini => IniM
    procedure :: Fin => FinM
    procedure :: zeros
    procedure :: eye
    procedure :: T => Trans
    procedure :: Inv => Inverse
    procedure :: Det
    procedure :: Random => GetRandomMatrix
    procedure :: blk => block_dmat
    procedure :: prt => MatrixPrint
    procedure :: DiagMat
  end type DMat
contains
  subroutine iniM(a, m, n)
  class(DMat), intent(inout) :: a
    integer(kp), intent(in) :: m,n
    if(allocated(a%m)) deallocate(a%m)
    allocate(a%m(m,n))
  end subroutine iniM

  subroutine zeros(a, m, n)
  class(DMat), intent(inout) :: a
    integer(kp), intent(in) :: m,n
    if(allocated(a%m)) deallocate(a%m)
    allocate(a%m(m,n))
    a%m = 0.d0
  end subroutine zeros

  subroutine eye(a, n)
  class(DMat), intent(inout) :: a
    integer(kp), intent(in) :: n
    integer(kp) :: i
    if(allocated(a%m)) deallocate(a%m)
    allocate(a%m(n,n))
    a%m = 0.d0
    do i = 1, n
      a%m(i,i) = 1.d0
    end do
  end subroutine eye

  subroutine FinM(a)
  class(DMat), intent(inout) :: a
    if(allocated(a%m)) deallocate(a%m)
  end subroutine FinM

  subroutine MatrixCopyD(b, a)
    type(DMat), intent(inout) :: b
    type(DMat), intent(in) :: a
    integer(kp) :: m, n, i
    m = size(a%m, 1)
    n = size(a%m, 2)
    call b%Ini(m,n)
    do i = 1, n
      call dcopy(m, a%m(:,i), 1, b%m(:,i), 1)
    end do
  end subroutine MatrixCopyD

  type(DMat) function MatrixProductD(a, b) result(c)
    type(DMat), intent(in) :: a, b
    integer(kp) :: m, k, n
    m = size(a%m, 1)
    k = size(a%m, 2)
    if(size(a%m, 2) /= size(b%m, 1)) then
      write(*, '(a)') 'Error in MatrixProduct'
      stop
    end if
    n = size(b%m, 2)
    call c%Ini(m,n)
    call dgemm('n','n',m,n,k,1.d0,a%m,m,b%m,k,0.d0,c%m,m)
  end function MatrixProductD

  type(DMat) function MatrixSumD(a, b) result(c)
    type(DMat), intent(in) :: a, b
    integer(kp) :: m, n, i
    if(size(a%m, 1) /= size(b%m, 1) .or. size(a%m, 2) /= size(b%m, 2)) then
      write(*, '(a)') 'Error in MatrixSum'
      stop
    end if
    m = size(a%m, 1)
    n = size(a%m, 2)
    call MatrixCopyD(c, a)
    do i = 1, n
      call daxpy(m, 1.d0, b%m(:,i), 1, c%m(:,i), 1)
    end do
  end function MatrixSumD

  type(DMat) function MatrixSubtractD(a, b) result(c)
    type(DMat), intent(in) :: a, b
    integer(kp) :: m, n, i
    if(size(a%m, 1) /= size(b%m, 1) .or. size(a%m, 2) /= size(b%m, 2)) then
      write(*, '(a)') 'Error in MatrixSum'
      stop
    end if
    m = size(a%m, 1)
    n = size(a%m, 2)
    call MatrixCopyD(c, a)
    do i = 1, n
      call daxpy(m, -1.d0, b%m(:,i), 1, c%m(:,i), 1)
    end do
  end function MatrixSubtractD

  type(DMat) function MatrixScaleLD(b, a) result(c)
    type(DMat), intent(in) :: b
    real(dp), intent(in) :: a
    integer(kp) :: m, n, i
    m = size(b%m, 1)
    n = size(b%m, 2)
    call MatrixCopyD(c, b)
    do i = 1, n
      call dscal(m, a, c%m(:,i), 1)
    end do
  end function MatrixScaleLD

  type(DMat) function MatrixScaleRD(a, b) result(c)
    type(DMat), intent(in) :: b
    real(dp), intent(in) :: a
    integer(kp) :: m, n, i
    m = size(b%m, 1)
    n = size(b%m, 2)
    call MatrixCopyD(c, b)
    do i = 1, n
      call dscal(m, a, c%m(:,i), 1)
    end do
  end function MatrixScaleRD

  type(DMat) function MatrixScaleDivideD(b, a) result(c)
    type(DMat), intent(in) :: b
    real(dp), intent(in) :: a
    integer(kp) :: m, n, i
    m = size(b%m, 1)
    n = size(b%m, 2)
    call MatrixCopyD(c, b)
    do i = 1, n
      call dscal(m, 1.d0 / a, c%m(:,i), 1)
    end do
  end function MatrixScaleDivideD

  type(DMat) function Trans(a) result(b)
    class(DMat), intent(in) :: a
    integer(kp) :: n, m
    m = size(a%m, 1)
    n = size(a%m, 2)
    call b%Ini(n,m)
    b%M = transpose(a%M)
  end function Trans

  type(DMat) function inverse(r) result(s)
  class(DMat), intent(in) :: r
    real(dp), allocatable :: a(:,:)
    real(dp), allocatable :: work(:)
    integer(kp), allocatable :: ipvt(:)
    integer(kp) :: info, n
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

  real(dp) function Det(r) result(d)
  class(DMat), intent(in) :: r
    integer(kp) :: n, i, info
    real(dp), allocatable :: a(:,:)
    integer(kp), allocatable :: ipiv(:)
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

  subroutine MatrixPrint(this, msg, iunit, binary)
    class(DMat), intent(in) :: this
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

  function block_dmat(this, m1, m2, n1, n2) result(r)
  class(DMat), intent(in) :: this
    type(DMat) :: r
    integer(kp), intent(in) :: m1, m2, n1, n2
    integer(kp) :: m, n
    m = m2 - m1 + 1
    n = n2 - n1 + 1
    call r%ini(m,n)
    r%m(:,:) = this%m(m1:m2,n1:n2)
  end function block_dmat

  subroutine GetRandomMatrix(mat, m, n)
    use VectorDouble, only: DVec
  class(DMat), intent(inout) :: mat
    integer(kp), intent(in) :: m, n
    integer(kp) :: i
    type(DVec) :: v
    call mat%ini(m,n)
    do i = 1, n
      call v%Random(m)
      mat%m(:,i) = v%v(:)
      call v%fin()
    end do
  end subroutine GetRandomMatrix

  subroutine DiagMat(b, a)
    use VectorDouble, only: DVec
  class(DMat), intent(inout) :: b
    type(DVec), intent(in) :: a
    integer(kp) :: n, i
    n = size(a%V)
    call b%Ini(n,n)
    do i = 1, n
      b%M(i,i) = a%V(i)
    end do
  end subroutine DiagMat

end module MatrixDouble
