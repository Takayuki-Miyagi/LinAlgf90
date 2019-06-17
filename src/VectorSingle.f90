module VectorSingle
  use LinAlgParameters
  implicit none

  public :: SVec
  public :: VectorCopyS
  public :: VectorSumS
  public :: VectorSubtractS
  public :: VectorScaleRS
  public :: VectorScaleLS
  public :: VectorDivideS
  public :: InnerProductS

  private :: IniV
  private :: zeros
  private :: FinV
  private :: VectorPrint
  private :: GetRandomVector
  private :: Nrm
  private :: Nrm2
  private :: block_svec

  type :: SVec
    real(sp), allocatable :: V(:)
    integer :: n_size = 0
  contains
    procedure :: Ini => iniV
    procedure :: zeros
    procedure :: Fin => FinV
    procedure :: prt => VectorPrint
    procedure :: Random => GetRandomVector
    procedure :: blk => block_svec
    procedure :: Nrm
    procedure :: Nrm2
  end type SVec
contains
  subroutine IniV(a, n)
    class(SVec), intent(inout) :: a
    integer(kp), intent(in) :: n
    if(allocated(a%V)) deallocate(a%V)
    a%n_size = n
    allocate(a%V(a%n_size))
  end subroutine IniV

  subroutine zeros(a, n)
    class(SVec), intent(inout) :: a
    integer(kp), intent(in) :: n
    call a%ini(n)
    a%v(:) = 0.0
  end subroutine zeros

  subroutine FinV(a)
    class(SVec), intent(inout) :: a
    if(allocated(a%V)) deallocate(a%V)
    a%n_size = 0
  end subroutine FinV

  subroutine VectorCopyS(b, a)
    type(SVec), intent(inout) :: b
    type(SVec), intent(in) :: a
    integer(kp) :: n
    n = a%n_size
    if(n < 1) return
    call b%Ini(n)
    call scopy(n, a%v, 1, b%v, 1)
  end subroutine VectorCopyS

  type(SVec) function VectorSumS(a, b) result(c)
    type(SVec), intent(in) :: a, b
    integer(kp) :: n
    if(a%n_size /= b%n_size) then
      write(*,'(a)') 'Error in DVectorSum'
      stop
    end if
    n = a%n_size
    if(n < 1) return
    call VectorCopyS(c, a)
    call saxpy(n, 1.0, b%v, 1, c%v, 1)
  end function VectorSumS

  type(SVec) function VectorSubtractS(a, b) result(c)
    type(SVec), intent(in) :: a, b
    integer(kp) :: n
    if(a%n_size /= b%n_size) then
      write(*,'(a)') 'Error in SVectorSubtract'
      stop
    end if
    n = a%n_size
    if(n < 1) return
    call VectorCopyS(c, a)
    call saxpy(n, -1.0, b%v, 1, c%v, 1)
  end function VectorSubtractS

  type(SVec) function VectorScaleRS(a, b) result(c)
    type(SVec), intent(in) :: a
    real(sp), intent(in) :: b
    integer(kp) :: n
    n = a%n_size
    if(n < 1) return
    call VectorCopyS(c, a)
    call sscal(n, b, c%v, 1)
  end function VectorScaleRS

  type(SVec) function VectorScaleLS(b, a) result(c)
    type(SVec), intent(in) :: a
    real(sp), intent(in) :: b
    integer(kp) :: n
    n = a%n_size
    if(n < 1) return
    call VectorCopyS(c, a)
    call sscal(n, b, c%v, 1)
  end function VectorScaleLS

  type(SVec) function VectorDivideS(a, b) result(c)
    type(SVec), intent(in) :: a
    real(sp), intent(in) :: b
    integer(kp) :: n
    n = a%n_size
    if(n < 1) return
    call VectorCopyS(c, a)
    call sscal(n, 1.0 / b, c%v, 1)
  end function VectorDivideS

  real(sp) function InnerProductS(a, b) result(c)
    type(SVec), intent(in) :: a, b
    integer(kp) :: n
    real(sp) :: sdot
    if(a%n_size /= b%n_size) then
      write(*,'(a)') 'Error in InnerProduct'
      stop
    end if
    c = 0.0
    n = a%n_size
    if(n < 1) return
    c = sdot(n, a%v, 1, b%v, 1)
  end function InnerProductS

  real(sp) function Nrm(a) result(b)
    class(SVec), intent(in) :: a
    integer(kp) :: n
    real(sp) :: snrm2
    b = 0.0
    n = a%n_size
    if(n < 1) return
    b = snrm2(n, a%v, 1)
  end function Nrm

  real(sp) function Nrm2(a) result(b)
    class(SVec), intent(in) :: a
    integer(kp) :: n
    real(sp) :: sdot
    b = 0.0
    n = a%n_size
    if(n < 1) return
    b = sdot(n, a%v, 1, a%v, 1)
  end function Nrm2

  subroutine VectorPrint(this, msg, iunit, binary)
    class(SVec),intent(in)::this
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
    class(SVec), intent(inout) :: v
    integer(kp), intent(in) :: n
    integer(kp), intent(in), optional :: dist
    integer(kp) :: idist = 3
    if(n < 1) return
    if(present(dist)) idist = dist
    call v%ini(n)
    call slarnv(idist, iseed, n, v%v)
  end subroutine GetRandomVector

  function block_SVec(this, n1, n2) result(r)
    class(SVec), intent(in) :: this
    type(SVec) :: r
    integer(kp), intent(in) :: n1, n2
    integer(kp) :: n
    n = n2 - n1 + 1
    if(n < 1) return
    call r%ini(n)
    r%v(:) = this%v(n1:n2)
  end function block_SVec

end module VectorSingle
