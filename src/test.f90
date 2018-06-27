program test
  use VectorDouble, only: DVec
  use VectorComplex, only: CVec
  use MatrixDouble, only: DMat
  use MatrixComplex, only: CMat
  use LinAlgLib

  !call check_dvector()
  call check_dmatrix()
  !call check_zvector()
  call check_zmatrix()

contains
  subroutine check_dvector()
    type(DVec) :: a, b, c
    integer :: n = 2
    call a%Random(n)
    call b%Random(n)
    c = a + b
    call c%prt()
    c = a - b
    call c%prt()
    write(*,*) a * b
  end subroutine check_dvector

  subroutine check_zvector()
    type(CVec) :: a, b, c
    integer :: n = 2
    call a%Random(n)
    call b%Random(n)
    c = a + b
    call c%prt()
    c = a - b
    call c%prt()
    write(*,*) a * b
  end subroutine check_zvector

  subroutine check_dmatrix()
    type(DMat) :: a, b, c
    integer :: n = 2
    call a%Random(n,n)
    call b%Random(n,n)
    a%m(1,:) = (/1.d0, 0.d0/)
    a%m(2,:) = (/1.d0, 0.d0/)
    c = a%T()
    b = a * c
    call a%prt()
    call c%prt()
    call b%prt()
  end subroutine check_dmatrix

  subroutine check_zmatrix()
    type(CMat) :: a, b, c
    integer :: n = 2
    call a%Random(n,n)
    call b%Random(n,n)
    c = a%H()
    call a%prt()
    call c%prt()
  end subroutine check_zmatrix

end program test
