module LinAlgLib
  use VectorDouble, only: DVec, VectorCopyD, VectorSumD, VectorSubtractD, &
    & VectorScaleRD, VectorScaleLD, InnerProductD, VectorDivideD
  use VectorComplex, only: CVec, VectorCopyC, VectorSumC, VectorSubtractC, &
    & VectorScaleRC, VectorScaleLC, InnerProductC, VectorDivideC
  use MatrixDouble, only: DMat, MatrixCopyD, MatrixSumD, MatrixSubtractD, &
    & MatrixScaleRD, MatrixScaleLD, MatrixProductD, MatrixScaleDivideD
  use MatrixComplex, only: CMat, MatrixCopyC, MatrixSumC, MatrixSubtractC, &
    & MatrixScaleRC, MatrixScaleLC, MatrixProductC, MatrixScaleDivideC
  use MatVecDouble, only: OuterProductD, VMProductD, MVProductD
  use MatVecComplex, only: OuterProductC, VMProductC, MVProductC

  interface assignment(=)
    procedure :: VectorCopyD
    procedure :: VectorCopyC
    procedure :: MatrixCopyD
    procedure :: MatrixCopyC
  end interface assignment(=)

  interface operator(+)
    procedure :: VectorSumD
    procedure :: VectorSumC
    procedure :: MatrixSumD
    procedure :: MatrixSumC
  end interface operator(+)

  interface operator(-)
    procedure :: VectorSubtractD
    procedure :: VectorSubtractC
    procedure :: MatrixSubtractD
    procedure :: MatrixSubtractC
  end interface operator(-)

  interface operator(*)
    procedure :: VectorScaleRD
    procedure :: VectorScaleRC
    procedure :: VectorScaleLD
    procedure :: VectorScaleLC
    procedure :: InnerProductD
    procedure :: InnerProductC
    procedure :: MatrixScaleLC
    procedure :: MatrixScaleLD
    procedure :: MatrixScaleRC
    procedure :: MatrixScaleRD
    procedure :: MatrixProductC
    procedure :: MatrixProductD
    procedure :: MVProductD
    procedure :: VMProductD
    procedure :: MVProductC
    procedure :: VMProductC
  end interface operator(*)

  interface operator(/)
    procedure :: VectorDivideD
    procedure :: VectorDivideC
    procedure :: MatrixScaleDivideD
    procedure :: MatrixScaleDivideC
  end interface operator(/)

  interface operator(.x.)
    procedure :: OuterProductD
    procedure :: OuterProductC
  end interface operator(.x.)

  type :: EigenSolSymD
    type(DVec) :: eig
    type(DMat) :: vec
  contains
    procedure :: init
    procedure :: fin
    procedure :: DiagSym
  end type EigenSolSymD
contains

  subroutine init(this, A)
  class(EigenSolSymD) :: this
    type(DMat), intent(in) :: A
    integer :: n
    n = size(A%m)
    call this%eig%ini(n)
    call this%vec%ini(n,n)
  end subroutine init

  subroutine fin(this)
  class(EigenSolSymD) :: this
    call this%eig%fin()
    call this%vec%fin()
  end subroutine fin

  subroutine DiagSym(this, A, qmin, qmax, m)
    use parameters, only: eps
  class(EigenSolSymD) :: this
    type(DMat), intent(in) :: A
    real(8), intent(in), optional :: qmin, qmax
    integer, intent(in), optional :: m
    integer :: num
    real(8), allocatable :: work(:)
    integer, allocatable :: iwork(:), ifailv(:)
    integer :: info, lwork, n
    real(8) :: lw, dlamch
    n = size(A%M, 1)
    call this%eig%ini(n)
    call this%vec%ini(n,n)
    this%vec = A
    if(.not. present(m) .and. .not. present(qmin) .and. &
      & .not. present(qmax)) then
      call dsyev('v', 'u', n, A%m, n, eig, lw, -1, info)
      lwork = int(lw)
      allocate(work(lwork))
      call dsyev('v', 'u', n, this%vec%m, n, this%eig%v, work, lwork, info)
      deallocate(work)

    elseif(present(m)) then
      allocate(iwork(5*n), ifailv(n))
      call dsyevx('v', 'i', 'u', n, A%m, n, -1.d100, 1.d100, 1, n, dlamch('S'), &
        &  num, this%eig%v, this%vec%m, n, lw, -1, iwork, ifailv, info)
      deallocate( iwork, ifailv)

    else
      allocate(iwork(5*n), ifailv(n))
      call dsyevx('v', 'i', 'u', n, A%m, n, -1.d100, 1.d100, 1, n, dlamch('S'), &
        &  num, this%eig%v, this%vec%m, n, lw, -1, iwork, ifailv, info)
      lwork = int(lw)
      allocate(work(1:lwork))
      call dsyevx('v', 'v', 'u', n, A%m, n, qmin, qmax, 1, n, eps, &
        &  num, this%eig%v, this%vec%m, n, work, lwork, iwork, ifailv, info)
      deallocate( iwork, ifailv)
    end if
  end subroutine DiagSym

!  subroutine DiagonalizationSymmetric(r, eval, evec, qmin, qmax, m, tol)
!  class(MO), intent(in) :: r
!    type(Vec), intent(out) :: eval
!    type(MO), intent(out) :: evec
!    real(8), optional, intent(in) :: qmin, qmax, tol
!    integer, optional, intent(in) :: m
!    integer :: num
!    real(8), allocatable :: work(:), mat(:,:), eig(:)
!    integer, allocatable :: iwork(:), ifailv(:)
!    integer :: info, lwork, n, i
!    real(8) :: lw, dlamch
!
!    n = size(r%M, 1)
!    call eval%Ini(n)
!    call evec%Ini(n,n)
!    allocate(eig(n)); eig(:) = 0.d0
!
!    if(.not. present(m) .and. .not. present(qmin) .and. &
!      & .not. present(qmax) .and. .not. present(tol)) then
!      allocate(mat(n,n))
!      mat = r%m
!      call dsyev('v', 'u', n, mat, n, eig, lw, -1, info)
!      lwork = int(lw)
!      allocate(work(lwork))
!      call dsyev('v', 'u', n, mat, n, eig, work, lwork, info)
!      evec%m = mat
!      deallocate(mat)
!    elseif(present(m)) then
!      allocate(iwork(5*n), ifailv(n))
!      call dsyevx('v','i','u',n,r%m,n,-1.d100,1.d100,1,n,dlamch('S'), &
!        &  num,eig,evec%m,n,lw,-1,iwork,ifailv,info)
!      deallocate( iwork, ifailv)
!    else
!      allocate(iwork(5*n), ifailv(n))
!      call dsyevx('v','i','u',n,r%m,n,-1.d100,1.d100,1,n,dlamch('S'), &
!        &  num,eig,evec%m,n,lw,-1,iwork,ifailv,info)
!      lwork = int(lw)
!      allocate(work(1:lwork))
!      call dsyevx('v','v','u',n,r%m,n,qmin,qmax,1,n,tol, &
!        &  num,eig,evec%m,n,work,lwork,iwork,ifailv,info)
!      deallocate( iwork, ifailv)
!    end if
!    eval%v = eig
!    deallocate(eig)
!  end subroutine DiagonalizationSymmetric
!
!  subroutine EigenValuesSymmetric(r, eval, m)
!  class(MO), intent(inout) :: r
!    type(Vec), intent(out) :: eval
!    integer, intent(in) :: m
!    integer, allocatable :: iwork(:), iblock(:), isplit(:)
!    real(8), allocatable :: work(:), mat(:,:), d(:), e(:), tau(:), w(:)
!    real(8) :: dlamch, lw
!    integer :: n, info, lwork, nsplit, i
!    n = size(r%m)
!    call eval%Ini(n)
!    allocate(mat(n,n))
!    allocate(d(n), e(max(1, n-1)), tau(max(1, n-1)), w(n))
!    allocate(iblock(n), isplit(n))
!    mat = r%m
!    call dsytrd('u',n,mat,n,d,e,tau,lw,-1,info)
!    lwork = int(lw)
!    allocate(work(lwork))
!    call dsytrd('u',n,mat,n,d,e,tau,work,lwork,info)
!    call dstebz('i','e',n,0.d0,0.d0,1,min(n,m),dlamch('S'), &
!      & d, e, m, nsplit, w, iblock, isplit, work, iwork, info)
!    do i = 1, min(n,m)
!      eval%v(i) = w(i)
!    end do
!    deallocate(work)
!    deallocate(mat, d, e, tau,iblock, isplit)
!  end subroutine EigenValuesSymmetric

end module LinAlgLib
