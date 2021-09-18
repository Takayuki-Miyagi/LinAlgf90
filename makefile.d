obj/LinAlgLib.o : src/LinAlgLib.f90 obj/MatVecComplex.o obj/MatVecDouble.o obj/MatVecSingle.o obj/MatrixComplex.o obj/MatrixDouble.o obj/MatrixSingle.o obj/VectorComplex.o obj/VectorDouble.o obj/VectorSingle.o obj/SingleDoubleComplex.o obj/LinAlgParameters.o 
obj/LinAlgParameters.o : src/LinAlgParameters.f90 
obj/MatVecComplex.o : src/MatVecComplex.f90 obj/MatrixComplex.o obj/VectorComplex.o obj/LinAlgParameters.o 
obj/MatVecDouble.o : src/MatVecDouble.f90 obj/MatrixDouble.o obj/VectorDouble.o obj/LinAlgParameters.o 
obj/MatVecSingle.o : src/MatVecSingle.f90 obj/MatrixSingle.o obj/VectorSingle.o obj/LinAlgParameters.o 
obj/MatrixComplex.o : src/MatrixComplex.f90 obj/VectorComplex.o obj/LinAlgParameters.o 
obj/MatrixDouble.o : src/MatrixDouble.f90 obj/VectorDouble.o obj/LinAlgParameters.o 
obj/MatrixSingle.o : src/MatrixSingle.f90 obj/VectorSingle.o obj/LinAlgParameters.o 
obj/SingleDoubleComplex.o : src/SingleDoubleComplex.f90 obj/MatrixComplex.o obj/VectorComplex.o obj/MatrixDouble.o obj/VectorDouble.o obj/MatrixSingle.o obj/VectorSingle.o 
obj/VectorComplex.o : src/VectorComplex.f90 obj/LinAlgParameters.o 
obj/VectorDouble.o : src/VectorDouble.f90 obj/LinAlgParameters.o 
obj/VectorSingle.o : src/VectorSingle.f90 obj/LinAlgParameters.o 
