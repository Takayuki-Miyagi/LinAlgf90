module LinAlgLib
  use VectorDouble
  use VectorComplex
  use MatrixDouble
  use MatrixComplex
  use MatVecDouble
  use MatVecComplex

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

end module LinAlgLib
