module VectorDouble
  use LinAlgParameters
  implicit none

  private :: IniV, FinV, VectorPrintAscii, GetRandomVector, Nrm, Nrm2
  private :: block_dvec, VectorPrintBinary

  public :: DVec, VectorCopyD, VectorSumD, VectorSubtractD
  public :: VectorScaleRD, VectorScaleLD, VectorDivideD, InnerProductD

  type :: DVec
    real(8), allocatable :: V(:)
  contains
    procedure :: Ini => iniV
    procedure :: zeros
    procedure :: Fin => FinV
    procedure :: prt => VectorPrintAscii
    procedure :: prtbin => VectorPrintBinary
    procedure :: Random => GetRandomVector
    procedure :: blk => block_dvec
    procedure :: Nrm
    procedure :: Nrm2
  end type DVec
contains
  subroutine IniV(a, n)
  class(DVec), intent(inout) :: a
    integer(4), intent(in) :: n
    if(allocated(a%V)) deallocate(a%V)
    allocate(a%V(n))
  end subroutine IniV

  subroutine zeros(a, n)
  class(DVec), intent(inout) :: a
    integer(4), intent(in) :: n
    if(allocated(a%V)) deallocate(a%V)
    allocate(a%V(n))
    a%v(:) = 0.d0
  end subroutine zeros

  subroutine FinV(a)
  class(DVec), intent(inout) :: a
    if(allocated(a%V)) deallocate(a%V)
  end subroutine FinV

  subroutine VectorCopyD(b, a)
    type(DVec), intent(inout) :: b
    type(DVec), intent(in) :: a
    integer(4) :: n
    n = size(a%V)
    call b%Ini(n)
    call dcopy(n, a%v, 1, b%v, 1)
  end subroutine VectorCopyD

  type(DVec) function VectorSumD(a, b) result(c)
    type(DVec), intent(in) :: a, b
    integer(4) :: n
    if(size(a%v) /= size(b%v)) then
      write(*,'(a)') 'Error in DVectorSum'
      stop
    end if
    n = size(a%v)
    call VectorCopyD(c, a)
    call daxpy(n, 1.d0, b%v, 1, c%v, 1)
  end function VectorSumD

  type(DVec) function VectorSubtractD(a, b) result(c)
    type(DVec), intent(in) :: a, b
    integer(4) :: n
    if(size(a%v) /= size(b%v)) then
      write(*,'(a)') 'Error in DVectorSubtract'
      stop
    end if
    n = size(a%v)
    call VectorCopyD(c, a)
    call daxpy(n, -1.d0, b%v, 1, c%v, 1)
  end function VectorSubtractD

  type(DVec) function VectorScaleRD(a, b) result(c)
    type(DVec), intent(in) :: a
    real(8), intent(in) :: b
    integer(4) :: n
    n = size(a%v)
    call VectorCopyD(c, a)
    call dscal(n, b, c%v, 1)
  end function VectorScaleRD

  type(DVec) function VectorScaleLD(b, a) result(c)
    type(DVec), intent(in) :: a
    real(8), intent(in) :: b
    integer(4) :: n
    n = size(a%v)
    call VectorCopyD(c, a)
    call dscal(n, b, c%v, 1)
  end function VectorScaleLD

  type(DVec) function VectorDivideD(a, b) result(c)
    type(DVec), intent(in) :: a
    real(8), intent(in) :: b
    integer(4) :: n
    n = size(a%v)
    call VectorCopyD(c, a)
    call dscal(n, 1.d0 / b, c%v, 1)
  end function VectorDivideD

  real(8) function InnerProductD(a, b) result(c)
    type(DVec), intent(in) :: a, b
    integer(4) :: n
    real(8) :: ddot
    if(size(a%v) /= size(b%v)) then
      write(*,'(a)') 'Error in InnerProduct'
      stop
    end if
    n = size(a%v)
    c = ddot(n, a%v, 1, b%v, 1)
  end function InnerProductD

  real(8) function Nrm(a) result(b)
  class(DVec), intent(in) :: a
    integer(4) :: n
    real(8) :: dnrm2
    n = size(a%v)
    b = dnrm2(n, a%v, 1)
  end function Nrm

  real(8) function Nrm2(a) result(b)
  class(DVec), intent(in) :: a
    integer(4) :: n
    real(8) :: ddot
    n = size(a%v)
    b = ddot(n, a%v, 1, a%v, 1)
  end function Nrm2

  subroutine VectorPrintAscii(this, iunit, msg)
  class(DVec),intent(in)::this
    integer :: i, n, unt
    integer, intent(in), optional :: iunit
    character(*), intent(in), optional :: msg
    if(present(iunit)) then; unt = iunit
    else; unt = 6; end if
    n = size(this%v, 1)
    if(unt == 6) then
      if(present(msg)) write(*,*) msg
      write(unt,'(10f10.4)') this%v(:)
    else
      do i = 1, n
        write(unt,'(10f10.4)') this%v(i)
      end do
    end if
  end subroutine VectorPrintAscii

  subroutine VectorPrintBinary(this, iunit)
    class(DVec),intent(in)::this
    integer, intent(in) :: iunit
    write(iunit) this%v(:)
  end subroutine VectorPrintBinary

  subroutine GetRandomVector(v, n, dist)
    ! idist = 1: uniform (0, 1)
    ! idist = 2: uniform (-1, 1)
    ! idist = 3: normal (-1, 1)
  class(DVec), intent(inout) :: v
    integer(4), intent(in) :: n
    integer(4), intent(in), optional :: dist
    integer(4) :: idist = 3
    if(present(dist)) idist = dist
    call v%ini(n)
    call dlarnv(idist, iseed, n, v%v)
  end subroutine GetRandomVector

  function block_dvec(this, n1, n2) result(r)
  class(DVec), intent(in) :: this
    type(DVec) :: r
    integer(4), intent(in) :: n1, n2
    integer(4) :: n
    n = n2 - n1 + 1
    call r%ini(n)
    r%v(:) = this%v(n1:n2)
  end function block_dvec
end module VectorDouble
