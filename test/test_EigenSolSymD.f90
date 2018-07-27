program test_EigenSolSymD
  !
  ! sample matrix is taken from
  ! http://www.nag-j.co.jp/lapack/dsyev.htm
  !
  !     (1 2 3 4)
  ! A = (2 2 3 4)
  !     (3 3 3 4)
  !     (4 4 4 4)


  use MatrixDouble, only: DMat
  use LinAlgLib

  type(DMat) :: a
  type(EigenSolSymD) :: sol
  call a%ini(4,4)
  a%m(:,1) = (/1.d0, 2.d0, 3.d0, 4.d0/)
  a%m(:,2) = (/2.d0, 2.d0, 3.d0, 4.d0/)
  a%m(:,3) = (/3.d0, 3.d0, 3.d0, 4.d0/)
  a%m(:,4) = (/4.d0, 4.d0, 4.d0, 4.d0/)

  call sol%init(a)
  call sol%DiagSym(a) ! w/o error estimation
  call sol%DiagSym(a, error=1) ! w/ error estimation

  call sol%eig%prt('Eigen Values')
  call sol%vec%prt('Eigen Vectors')


end program test_EigenSolSymD
