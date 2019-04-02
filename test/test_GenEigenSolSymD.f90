program test_EigenSolSymD
  !
  ! sample matrix is taken from
  ! http://www.nag-j.co.jp/lapack/dsygvd.htm
  !
  !     ( 0.24  0.39  0.42 -0.16)
  ! A = ( 0.39 -0.11  0.79  0.63)
  !     ( 0.42  0.79 -0.25  0.48)
  !     (-0.16  0.63  0.48 -0.03)
  !
  !     ( 4.16 -3.12  0.56 -0.10)
  ! B = (-3.12  5.03 -0.83  1.09)
  !     ( 0.56 -0.83  0.76  0.34)
  !     (-0.10  1.09  0.34  1.18)

  use LinAlgLib

  type(DMat) :: a, b
  type(GenEigenSolSymD) :: sol

  call a%ini(4,4)
  a%m(:,1) = [ 0.24d0,  0.39d0,  0.42d0, -0.16d0]
  a%m(:,2) = [ 0.39d0, -0.11d0,  0.79d0,  0.63d0]
  a%m(:,3) = [ 0.42d0,  0.79d0, -0.25d0,  0.48d0]
  a%m(:,4) = [-0.16d0,  0.63d0,  0.48d0, -0.03d0]

  call b%ini(4,4)
  b%m(:,1) = [ 4.16d0, -3.12d0,  0.56d0, -0.10d0]
  b%m(:,2) = [-3.12d0,  5.03d0, -0.83d0,  1.09d0]
  b%m(:,3) = [ 0.56d0, -0.83d0,  0.76d0,  0.34d0]
  b%m(:,4) = [-0.10d0,  1.09d0,  0.34d0,  1.18d0]

  call sol%init(a,b,2)
  call sol%DiagSym(a,b)
  call sol%eig%prt('All Eigen Values')
  call sol%vec%prt('All Eigen Vectors')
  call sol%fin()


  call a%fin()
  call b%fin()
end program test_EigenSolSymD
