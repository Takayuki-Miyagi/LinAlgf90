module VectorDouble
  use LinAlgParameters
  implicit none

  public :: DVec
  public :: VectorCopyD
  public :: VectorSumD
  public :: VectorSubtractD
  public :: VectorScaleRD
  public :: VectorScaleLD
  public :: VectorDivideD
  public :: InnerProductD

  private :: IniV
  private :: zeros
  private :: FinV
  private :: VectorPrint
  private :: GetRandomVector
  private :: Nrm
  private :: Nrm2
  private :: block_dvec


  type :: DVec
    real(dp), allocatable :: V(:)
  contains
    procedure :: Ini => iniV
    procedure :: zeros
    procedure :: Fin => FinV
    procedure :: prt => VectorPrint
    procedure :: Random => GetRandomVector
    procedure :: blk => block_dvec
    procedure :: Nrm
    procedure :: Nrm2
  end type DVec
contains
  subroutine IniV(a, n)
    class(DVec), intent(inout) :: a
    integer(kp), intent(in) :: n
    if(allocated(a%V)) deallocate(a%V)
    allocate(a%V(n))
  end subroutine IniV

  subroutine zeros(a, n)
    class(DVec), intent(inout) :: a
    integer(kp), intent(in) :: n
    call a%ini(n)
    a%v(:) = 0.d0
  end subroutine zeros

  subroutine FinV(a)
    class(DVec), intent(inout) :: a
    if(allocated(a%V)) deallocate(a%V)
  end subroutine FinV

  subroutine VectorCopyD(b, a)
    type(DVec), intent(inout) :: b
    type(DVec), intent(in) :: a
    integer(kp) :: n
    n = size(a%V)
    if(n < 1) return
    call b%Ini(n)
    call dcopy(n, a%v, 1, b%v, 1)
  end subroutine VectorCopyD

  type(DVec) function VectorSumD(a, b) result(c)
    type(DVec), intent(in) :: a, b
    integer(kp) :: n
    if(size(a%v) /= size(b%v)) then
      write(*,'(a)') 'Error in DVectorSum'
      stop
    end if
    n = size(a%v)
    if(n < 1) return
    call VectorCopyD(c, a)
    call daxpy(n, 1.d0, b%v, 1, c%v, 1)
  end function VectorSumD

  type(DVec) function VectorSubtractD(a, b) result(c)
    type(DVec), intent(in) :: a, b
    integer(kp) :: n
    if(size(a%v) /= size(b%v)) then
      write(*,'(a)') 'Error in DVectorSubtract'
      stop
    end if
    n = size(a%v)
    if(n < 1) return
    call VectorCopyD(c, a)
    call daxpy(n, -1.d0, b%v, 1, c%v, 1)
  end function VectorSubtractD

  type(DVec) function VectorScaleRD(a, b) result(c)
    type(DVec), intent(in) :: a
    real(dp), intent(in) :: b
    integer(kp) :: n
    n = size(a%v)
    if(n < 1) return
    call VectorCopyD(c, a)
    call dscal(n, b, c%v, 1)
  end function VectorScaleRD

  type(DVec) function VectorScaleLD(b, a) result(c)
    type(DVec), intent(in) :: a
    real(dp), intent(in) :: b
    integer(kp) :: n
    n = size(a%v)
    if(n < 1) return
    call VectorCopyD(c, a)
    call dscal(n, b, c%v, 1)
  end function VectorScaleLD

  type(DVec) function VectorDivideD(a, b) result(c)
    type(DVec), intent(in) :: a
    real(dp), intent(in) :: b
    integer(kp) :: n
    n = size(a%v)
    if(n < 1) return
    call VectorCopyD(c, a)
    call dscal(n, 1.d0 / b, c%v, 1)
  end function VectorDivideD

  real(dp) function InnerProductD(a, b) result(c)
    type(DVec), intent(in) :: a, b
    integer(kp) :: n
    real(dp) :: ddot
    if(size(a%v) /= size(b%v)) then
      write(*,'(a)') 'Error in InnerProduct'
      stop
    end if
    n = size(a%v)
    c = 0.d0
    if(n < 1) return
    c = ddot(n, a%v, 1, b%v, 1)
  end function InnerProductD

  real(dp) function Nrm(a) result(b)
    class(DVec), intent(in) :: a
    integer(kp) :: n
    real(dp) :: dnrm2
    b = 0.d0
    n = size(a%v)
    if(n < 1) return
    b = dnrm2(n, a%v, 1)
  end function Nrm

  real(dp) function Nrm2(a) result(b)
    class(DVec), intent(in) :: a
    integer(kp) :: n
    real(dp) :: ddot
    b = 0.d0
    n = size(a%v)
    if(n < 1) return
    b = ddot(n, a%v, 1, a%v, 1)
  end function Nrm2

  subroutine VectorPrint(this, msg, iunit, binary)
    class(DVec),intent(in)::this
    integer(kp) :: i, n, unt
    integer(kp), intent(in), optional :: iunit
    character(*), intent(in), optional :: msg
    logical, intent(in), optional :: binary
    logical :: bin

    if(present(iunit)) then; unt = iunit
    else; unt = 6; end if

    if(present(binary)) then; bin = binary
    else; bin = .false.; end if

    if(bin) then
      write(unt) this%v
    else
      if(present(msg)) write(unt,*) msg
      if(unt == 6) then
        write(unt,'(10f10.4)') this%v(:)
      else
        n = size(this%v, 1)
        do i = 1, n
          write(unt,'(10f10.4)') this%v(i)
        end do
      end if
    end if
  end subroutine VectorPrint

  subroutine GetRandomVector(v, n, dist)
    ! idist = 1: uniform (0, 1)
    ! idist = 2: uniform (-1, 1)
    ! idist = 3: normal (-1, 1)
    class(DVec), intent(inout) :: v
    integer(kp), intent(in) :: n
    integer(kp), intent(in), optional :: dist
    integer(kp) :: idist = 3
    if(n < 1) return
    if(present(dist)) idist = dist
    call v%ini(n)
    call dlarnv(idist, iseed, n, v%v)
  end subroutine GetRandomVector

  function block_dvec(this, n1, n2) result(r)
    class(DVec), intent(in) :: this
    type(DVec) :: r
    integer(kp), intent(in) :: n1, n2
    integer(kp) :: n
    n = n2 - n1 + 1
    if(n < 1) return
    call r%ini(n)
    r%v(:) = this%v(n1:n2)
  end function block_dvec
end module VectorDouble
