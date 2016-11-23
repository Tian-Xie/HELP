CC = icc
CC_ARG = -O2 -openmp
MKL_LIB = $MKLROOT/lib/intel64
MKL_ARG = -lmkl_rt
CUDA_LIB = /usr/local/cuda/lib64
CUDA_INCLUDE = /usr/local/cuda/include
CUDA_ARG = -lcuda -lcudart -lcusolver -lcusparse
SUITESPARSE_LIB = lib
SUITESPARSE_INCLUDE = lib
SUITESPARSE_ARG = -lamd -lcamd -lccolamd -lcholmod -lcolamd -lmetis -lsuitesparseconfig

SRC_DIR = 
LUSOL_DIR = $(SRC_DIR)lusol/
LUSOL_SRC = $(LUSOL_DIR)commonlib.cpp $(LUSOL_DIR)lusol.cpp $(LUSOL_DIR)myblas.cpp $(LUSOL_DIR)sparselib.cpp
GENERAL_SRC = $(SRC_DIR)Crushing.cpp $(SRC_DIR)GlobalVar.cpp $(SRC_DIR)Hash.cpp $(SRC_DIR)Helper.cpp $(SRC_DIR)HSDSolver.cpp $(SRC_DIR)Init.cpp $(SRC_DIR)Main.cpp $(SRC_DIR)MPSRead.cpp $(SRC_DIR)Presolve.cpp $(SRC_DIR)Presolve_Lindep.cpp $(SRC_DIR)Report.cpp

cholmod_cpu: 
	$(CC) $(CC_ARG) -L$(CUDA_LIB) -L$(SUITESPARSE_LIB) -I$(SUITESPARSE_INCLUDE) $(SUITESPARSE_ARG) -L$(MKL_LIB) $(MKL_ARG) -o CHOLMOD_CPU $(LUSOL_SRC) $(GENERAL_SRC) $(SRC_DIR)UseCholmod.cpp $(SRC_DIR)Setting_CPU.cpp

cholmod_gpu: 
	$(CC) $(CC_ARG) -L$(CUDA_LIB) -L$(SUITESPARSE_LIB) -I$(SUITESPARSE_INCLUDE) $(SUITESPARSE_ARG) -L$(MKL_LIB) $(MKL_ARG) -o CHOLMOD_GPU $(LUSOL_SRC) $(GENERAL_SRC) $(SRC_DIR)UseCholmod.cpp $(SRC_DIR)Setting_GPU.cpp

cusolver: 
	$(CC) $(CC_ARG) -L$(CUDA_LIB) -I$(CUDA_INCLUDE) -L$(CUDA_ARG) -L$(MKL_LIB) $(MKL_ARG) -o CUSOLVER $(LUSOL_SRC) $(GENERAL_SRC) $(SRC_DIR)UsecuSolver_chol_PREVIEW.cpp

pardiso:
	$(CC) $(CC_ARG) -L$(MKL_LIB) $(MKL_ARG) -o PARDISO $(LUSOL_SRC) $(GENERAL_SRC) $(SRC_DIR)UsePardiso.cpp
