program test
  use MatrixDouble, only: DMat
  use LinAlgLib

  type(DMat) :: a, b, c
  type(DSingularValueDecomposition) :: sol
  call a%Random(4,6)
  a%m(:,1) = [ 2.27d0,  -1.54d0,   1.15d0,  -1.94d0]
  a%m(:,2) = [ 0.28d0,  -1.67d0,   0.94d0,  -0.78d0]
  a%m(:,3) = [-0.48d0,  -3.09d0,   0.99d0,  -0.21d0]
  a%m(:,4) = [ 1.07d0,   1.22d0,   0.79d0,   0.63d0]
  a%m(:,5) = [-2.35d0,   2.93d0,  -1.45d0,   2.30d0]
  a%m(:,6) = [ 0.62d0,  -7.39d0,   1.03d0,  -2.57d0]

  a = a%T()
  call sol%init(a)
  call sol%SVD(a)
  call c%DiagMat(sol%Sigma)
  call a%prt("Original")
  call sol%U%prt("U")
  call c%prt("Sigma")
  call sol%V%prt("V")
  b = sol%U * c * sol%V
  call b%prt("U * Sigma * V")
  call sol%fin()
end program test
