./DoubleMatVec.o : ./DoubleMatVec.f90 
./LinAlgLib.o : ./LinAlgLib.f90 ./MatVecComplex.o ./MatVecDouble.o ./MatrixComplex.o ./MatrixDouble.o ./VectorComplex.o ./VectorDouble.o 
./MatVecComplex.o : ./MatVecComplex.f90 ./MatrixComplex.o ./VectorComplex.o 
./MatVecDouble.o : ./MatVecDouble.f90 ./MatrixDouble.o ./VectorDouble.o 
./MatrixComplex.o : ./MatrixComplex.f90 ./VectorComplex.o 
./MatrixDouble.o : ./MatrixDouble.f90 ./VectorDouble.o 
./VectorComplex.o : ./VectorComplex.f90 ./parameters.o 
./VectorDouble.o : ./VectorDouble.f90 ./parameters.o 
./parameters.o : ./parameters.f90 
