program test
  use LinAlgLib

  type(CVec) :: cv
  type(DVec) :: dv
  type(SVec) :: sv
  type(CMat) :: cm
  type(DMat) :: dm
  type(SMat) :: sm
  integer :: n = 3

  call cm%Random(n,n)
  dm = cm
  call dm%prt("real CMat => DMat")
  dm = (0.d0,-1.d0) * cm
  call dm%prt("img CMat => DMat")

  sm = cm
  call sm%prt("real CMat => SMat")
  sm = (0.d0,-1.d0) * cm
  call sm%prt("img CMat => SMat")

  cm = dm
  call cm%prt('DMat => CMat')
  cm = sm
  call cm%prt('SMat => CMat')

  call cv%Random(n)
  dv = cv
  call dv%prt("real CVec => DVec")
  dv = cv * (0.d0, -1.d0)
  call dv%prt("real CVec => DVec")
  sv = cv
  call sv%prt("real CVec => SVec")
  sv = cv * (0.d0, -1.d0)
  call sv%prt("real CVec => SVec")

  cv = dv
  call cv%prt('DVec => CVec')
  cv = sv
  call cv%prt('SVec => CVec')
end program test
