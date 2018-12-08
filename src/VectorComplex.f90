module VectorComplex
  use LinAlgParameters
  implicit none

  private :: iniV, FinV, ComplexConjugate, VectorPrintAscii, GetRandomVector
  private :: Nrm, Nrm2, block_cvec, VectorPrintBinary

  public :: CVec, VectorCopyC, VectorSumC, VectorSubtractC, VectorScaleRC
  public :: VectorScaleLC, VectorDivideC, InnerProductC

  type :: CVec
    complex(8), allocatable :: V(:)
  contains
    procedure :: Ini => iniV
    procedure :: zeros
    procedure :: Fin => FinV
    procedure :: CC => ComplexConjugate
    procedure :: prt => VectorPrintAscii
    procedure :: prtbin => VectorPrintBinary
    procedure :: Random => GetRandomVector
    procedure :: blk => block_cvec
    procedure :: Nrm
    procedure :: Nrm2
  end type CVec
contains
  subroutine IniV(a, n)
  class(CVec), intent(inout) :: a
    integer(4), intent(in) :: n
    if(allocated(a%V)) deallocate(a%V)
    allocate(a%V(n))
  end subroutine IniV

  subroutine zeros(a, n)
  class(CVec), intent(inout) :: a
    integer(4), intent(in) :: n
    if(allocated(a%V)) deallocate(a%V)
    allocate(a%V(n))
    a%V(:) = 0.d0
  end subroutine zeros

  subroutine FinV(a)
  class(CVec), intent(inout) :: a
    if(allocated(a%V)) deallocate(a%V)
  end subroutine FinV

  subroutine VectorCopyC(b, a)
    type(CVec), intent(inout) :: b
    type(CVec), intent(in) :: a
    integer(4) :: n
    n = size(a%V)
    call b%Ini(n)
    call zcopy(n, a%v, 1, b%v, 1)
  end subroutine VectorCopyC

  type(CVec) function VectorSumC(a, b) result(c)
    type(CVec), intent(in) :: a, b
    integer(4) :: n
    if(size(a%v) /= size(b%v)) then
      write(*,'(a)') 'Error in DVectorSum'
      stop
    end if
    n = size(a%v)
    call VectorCopyC(c,a)
    call zaxpy(n, 1.d0, b%v, 1, c%v, 1)
  end function VectorSumC

  type(CVec) function VectorSubtractC(a, b) result(c)
    type(CVec), intent(in) :: a, b
    integer(4) :: n
    if(size(a%v) /= size(b%v)) then
      write(*,'(a)') 'Error in DVectorSubtract'
      stop
    end if
    n = size(a%v)
    call VectorCopyC(c,a)
    call daxpy(n, -1.d0, b%v, 1, c%v, 1)
  end function VectorSubtractC

  type(CVec) function VectorScaleRC(a, b) result(c)
    type(CVec), intent(in) :: a
    complex(8), intent(in) :: b
    integer(4) :: n
    n = size(a%v)
    call VectorCopyC(c,a)
    call zscal(n, b, c%v, 1)
  end function VectorScaleRC


  type(CVec) function VectorScaleLC(b, a) result(c)
    type(CVec), intent(in) :: a
    complex(8), intent(in) :: b
    integer(4) :: n
    n = size(a%v)
    call VectorCopyC(c,a)
    call dscal(n, b, c%v, 1)
  end function VectorScaleLC

  type(CVec) function VectorDivideC(a, b) result(c)
    type(CVec), intent(in) :: a
    real(8), intent(in) :: b
    integer(4) :: n
    n = size(a%v)
    call VectorCopyC(c,a)
    call zdscal(n, 1.d0 / b, c%v, 1)
  end function VectorDivideC

  type(CVec) function ComplexConjugate(a) result(b)
  class(CVec), intent(in) :: a
    integer(4) :: n
    n = size(a%v)
    call b%ini(n)
    b%v = conjg(a%v)
  end function ComplexConjugate

  real(8) function InnerProductC(a, b) result(c)
    type(CVec), intent(in) :: a, b
    integer(4) :: n
    real(8) :: zdotu
    if(size(a%v) /= size(b%v)) then
      write(*,'(a)') 'Error in InnerProduct'
      stop
    end if
    n = size(a%v)
    c = zdotu(n, a%v, 1, b%v, 1)
  end function InnerProductC

  real(8) function Nrm(a) result(b)
  class(CVec), intent(in) :: a
    integer(4) :: n
    real(8) :: dznrm2
    n = size(a%v)
    b = dznrm2(n, a%v, 1)
  end function Nrm

  real(8) function Nrm2(a) result(b)
  class(CVec), intent(in) :: a
    integer(4) :: n
    real(8) :: dznrm2
    n = size(a%v)
    b = dznrm2(n, a%v, 1) ** 2
  end function Nrm2

  subroutine VectorPrintAscii(this, iunit, msg)
    class(CVec),intent(in)::this
    integer, intent(in), optional :: iunit
    character(*), intent(in), optional :: msg
    integer(4) :: n, i, unt
    if(present(iunit)) then; unt = iunit
    else; unt = 6; end if
    n = size(this%v, 1)
    if(unt == 6) then
      if(present(msg)) write(unt,*) msg
      write(unt,'(a)',advance='no') 'Real:'
      write(unt,'(10f10.4)') dble(this%v(:))
      write(unt,'(a)',advance='no') 'Imag:'
      write(unt,'(10f10.4)') dimag(this%v(:))
    else
      do i = 1, n
        write(unt,'(10f10.4)') this%v(i)
      end do
    end if
  end subroutine VectorPrintAscii

  subroutine VectorPrintBinary(this, iunit)
    class(CVec),intent(in)::this
    integer, intent(in) :: iunit
    integer(4) :: n
    write(iunit) this%v
  end subroutine VectorPrintBinary

  subroutine GetRandomVector(v, n, dist)
    ! idist = = 1:  real and imaginary parts each uniform (0,1)
    ! idist = = 2:  real and imaginary parts each uniform (-1,1)
    ! idist = = 3:  real and imaginary parts each normal (0,1)
    ! idist = = 4:  uniformly distributed on the disc abs(z) < 1
    ! idist = = 5:  uniformly distributed on the circle abs(z) = 1
  class(CVec), intent(inout) :: v
    integer(4), intent(in) :: n
    integer(4), intent(in), optional :: dist
    integer(4) :: idist = 4
    if(present(dist)) idist = dist
    call v%ini(n)
    call zlarnv(idist, iseed, n, v%v)
  end subroutine GetRandomVector

  function block_cvec(this, n1, n2) result(r)
  class(CVec), intent(in) :: this
    type(CVec) :: r
    integer(4), intent(in) :: n1, n2
    integer(4) :: n
    n = n2 - n1 + 1
    call r%ini(n)
    r%v(:) = this%v(n1:n2)
  end function block_cvec

end module VectorComplex
