program test
  use VectorDouble, only: DVec
  use LinAlgLib
  type(DVec) :: a, b, c
  integer :: n = 2
  call a%Random(n)
  call b%Random(n)

  call a%prt('a')
  call b%prt('b')
  c = a + b
  call c%prt('a + b')
  c = a - b
  call c%prt('a - b')

  write(*,'(a, f12.6)') 'a * b = ', a * b
  write(*,'(a, f12.6)') 'Norm a = ', a%Nrm()
  write(*,'(a, f12.6)') 'Norm b = ', b%Nrm()
  write(*,'(a, f12.6)') 'Norm**2 a = ', a%Nrm2()
  write(*,'(a, f12.6)') 'Norm**2 b = ', b%Nrm2()
end program test

