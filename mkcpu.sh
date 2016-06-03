#/bin/sh
icc -L/usr/local/cuda/lib64 -Llib -L$MKLROOT/lib/intel64 -lmkl_rt -llapack -lamd -lcamd -lccolamd -lcholmod -lcolamd -lmetis -lsuitesparseconfig -Iinclude -o CPUSolver ./src/Crushing.cpp ./src/GlobalVar.cpp ./src/Hash.cpp ./src/Helper.cpp ./src/HSDSolver.cpp ./src/Init.cpp ./src/Main.cpp ./src/MPSRead.cpp ./src/Presolve.cpp ./src/Setting_CPU.cpp ./src/UseCholmod.cpp

