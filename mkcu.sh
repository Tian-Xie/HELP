#/bin/sh
icc -L/usr/local/cuda/lib64 -lcuda -lcudart -lcusolver -lcusparse -L$MKLROOT/lib/intel64 -lmkl_rt -I/usr/local/cuda/include -Iinclude -o CUSolver ./src/Crushing.cpp ./src/GlobalVar.cpp ./src/Hash.cpp ./src/Helper.cpp ./src/HSDSolver.cpp ./src/Init.cpp ./src/Main.cpp ./src/MPSRead.cpp ./src/Presolve.cpp ./src/UsecuSolver_chol_PREVIEW.cpp
