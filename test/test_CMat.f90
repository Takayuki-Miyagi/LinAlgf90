program test
  use MatrixComplex, only: CMat
  use LinAlgLib

  type(CMat) :: a, b, c, d
  integer :: n = 2

  call a%Random(n,n)
  call b%Random(n,n)
  c = a + a%T()
  call c%prt(iunit=6,msg='a + a^t')

  c = a + b
  call c%prt(iunit=6,msg='a + b')

  c = a - b
  call c%prt(iunit=6,msg='a - b')

  c = a * a%inv()
  call c%prt(iunit=6,msg='a * a^-1')
  write(*,*) c%m

  c = b * a%inv()
  call c%prt(iunit=6,msg='b * a^-1')

  c = a%blk(1,1,1,n)
  call c%prt(iunit=6,msg='1st column vector')

  c = a%blk(1,n,1,1)
  call c%prt(iunit=6,msg='1st row vector')

  d = a%blk(1,1,1,n) * a%blk(1,n,1,1)
  call d%prt(iunit=15,msg='')

end program test
