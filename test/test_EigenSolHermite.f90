program test_EigenSolSymD
  !
  ! sample matrix is taken from
  ! http://www.nag-j.co.jp/lapack/zheev.htm
  !
  !     (1   2-i  3-i  4-i )
  ! A = (2+i 2    3-2i 4-2i)
  !     (3+i 3+2i 3    4-3i)
  !     (4+i 4+2i 4+3i 4   )

  use MatrixComplex, only: CMat
  use LinAlgLib

  type(CMat) :: a
  type(EigenSolHermite) :: sol

  call a%ini(4,4)
  a%m(1,:) = (/(1.d0,0.d0), (2.d0,-1.d0), (3.d0,-1.d0), (4.d0,-1.d0)/)
  a%m(2,:) = (/(2.d0,1.d0), (2.d0, 0.d0), (3.d0,-2.d0), (4.d0,-2.d0)/)
  a%m(3,:) = (/(3.d0,1.d0), (3.d0, 2.d0), (3.d0, 0.d0), (4.d0,-3.d0)/)
  a%m(4,:) = (/(4.d0,1.d0), (4.d0, 2.d0), (4.d0, 3.d0), (4.d0, 0.d0)/)
  call a%prt('a')

  call sol%init(a)
  call sol%DiagSym(a) ! w/o error estimation
  !call sol%DiagSym(a, error=1) ! w/ error estimation
  call sol%eig%prt('All Eigen Values')
  call sol%vec%prt('All Eigen Vectors')
  call sol%fin()

end program test_EigenSolSymD
