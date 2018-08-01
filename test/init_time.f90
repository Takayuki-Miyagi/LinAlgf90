program test
  use VectorDouble, only: DVec
  use VectorComplex, only: CVec
  use MatrixDouble, only: DMat
  use MatrixComplex, only: CMat
  use LinAlgLib
  !$ use omp_lib

  type(DVec) :: a
  type(CVec) :: b
  type(DMat) :: c
  type(CMat) :: d
  integer :: n = 1000
  integer :: i, imax = 10000
  !$ real(8) :: t

  !$ t = omp_get_wtime()

  do i = 1, imax
    call a%zeros(n)
    call a%fin()
  end do

  !$ write(*, '(a, i6, a, f12.6, a)') 'time for ', imax, ' times zeros double vector: ', &
  !$    & omp_get_wtime() - t, ' sec'

  !$ t = omp_get_wtime()

  do i = 1, imax
    call a%ini(n)
    call a%fin()
  end do

  !$ write(*, '(a, i6, a, f12.6, a)') 'time for ', imax, ' times ini double vecotor: ', &
  !$    & omp_get_wtime() - t, ' sec'

  !$ t = omp_get_wtime()

  do i = 1, imax
    call b%zeros(n)
    call b%fin()
  end do

  !$ write(*, '(a, i6, a, f12.6, a)') 'time for ', imax, ' times zeros complex vector: ', &
  !$    & omp_get_wtime() - t, ' sec'

  !$ t = omp_get_wtime()

  do i = 1, imax
    call b%ini(n)
    call b%fin()
  end do

  !$ write(*, '(a, i6, a, f12.6, a)') 'time for ', imax, ' times ini complex vecotor: ', &
  !$    & omp_get_wtime() - t, ' sec'

  !$ t = omp_get_wtime()

  do i = 1, imax
    call c%zeros(n,n)
    call c%fin()
  end do

  !$ write(*, '(a, i6, a, f12.6, a)') 'time for ', imax, ' times zeros double matrix: ', &
  !$    & omp_get_wtime() - t, ' sec'

  !$ t = omp_get_wtime()

  do i = 1, imax
    call c%ini(n,n)
    call c%fin()
  end do

  !$ write(*, '(a, i6, a, f12.6, a)') 'time for ', imax, ' times ini double matrix: ', &
  !$    & omp_get_wtime() - t, ' sec'

  !$ t = omp_get_wtime()

  do i = 1, imax
    call c%eye(n)
    call c%fin()
  end do

  !$ write(*, '(a, i6, a, f12.6, a)') 'time for ', imax, ' times I double matrix: ', &
  !$    & omp_get_wtime() - t, ' sec'

  !$ t = omp_get_wtime()

  do i = 1, imax
    call d%zeros(n,n)
    call d%fin()
  end do

  !$ write(*, '(a, i6, a, f12.6, a)') 'time for ', imax, ' times zeros complex matrix: ', &
  !$    & omp_get_wtime() - t, ' sec'

  !$ t = omp_get_wtime()

  do i = 1, imax
    call d%ini(n,n)
    call d%fin()
  end do

  !$ write(*, '(a, i6, a, f12.6, a)') 'time for ', imax, ' times ini complex matrix: ', &
  !$    & omp_get_wtime() - t, ' sec'

  !$ t = omp_get_wtime()

  do i = 1, imax
    call d%eye(n)
    call d%fin()
  end do

  !$ write(*, '(a, i6, a, f12.6, a)') 'time for ', imax, ' times I complex matrix: ', &
  !$    & omp_get_wtime() - t, ' sec'
end program test
