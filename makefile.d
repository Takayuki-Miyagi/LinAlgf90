obj/LinAlgLib.o : src/LinAlgLib.f90 obj/MatVecComplex.o obj/MatVecDouble.o obj/MatrixComplex.o obj/MatrixDouble.o obj/VectorComplex.o obj/VectorDouble.o 
obj/MatVecComplex.o : src/MatVecComplex.f90 obj/MatrixComplex.o obj/VectorComplex.o 
obj/MatVecDouble.o : src/MatVecDouble.f90 obj/MatrixDouble.o obj/VectorDouble.o 
obj/MatrixComplex.o : src/MatrixComplex.f90 obj/VectorComplex.o 
obj/MatrixDouble.o : src/MatrixDouble.f90 obj/VectorDouble.o 
obj/Parameters.o : src/Parameters.f90 
obj/VectorComplex.o : src/VectorComplex.f90 obj/Parameters.o 
obj/VectorDouble.o : src/VectorDouble.f90 obj/Parameters.o 
