#/bin/sh
icc -L/usr/local/cuda/lib64 -Llib -lblas -llapack -lamd -lblas -lcamd -lccolamd -lcholmod -lcolamd -lmetis -lsuitesparseconfig -Iinclude -o GPUSolver ./src/Crushing.cpp ./src/GlobalVar.cpp ./src/Hash.cpp ./src/Helper.cpp ./src/HSDSolver.cpp ./src/Init.cpp ./src/Main.cpp ./src/MPSRead.cpp ./src/Presolve.cpp ./src/Setting_GPU.cpp

