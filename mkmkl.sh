#/bin/sh

icc -L$MKLROOT/lib/intel64 -lmkl_rt -Iinclude -o MKLSolver ./src/Crushing.cpp ./src/GlobalVar.cpp ./src/Hash.cpp ./src/Helper.cpp ./src/HSDSolver.cpp ./src/Init.cpp ./src/Main.cpp ./src/MPSRead.cpp ./src/Presolve.cpp ./src/UsePardiso.cpp

