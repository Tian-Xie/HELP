# HELP (Homogeneous and self-dual Efficient Linear Programming)

An experimental implementation of Homogeneous and Self-Dual Algorithm for Linear Programming. 

### Author

Tian Xie (Research Center for Management Science and Information Analytics, Shanghai University of Finance and Economics)

### Target

To explore the possibility of parallelization on multi-core CPU and GPU.

### Progress

Currently working on J. Gondzio's Matrix-free regime, an indirect solver (PCG) for large-scale testcases that direct solver cannot handle. 

### Disclaimer

This program is currently not finished testing. No warrenty. Use it at your own risk! 

# Credits

1. Fundamental implementation idea originated from COPL_LP (Xiong Zhang and Yinyu Ye). See http://web.stanford.edu/~yyye/Col.html .

2. Erling D. Andersen, Knud D. Andersen. (2000). The Mosek Interior Point Optimizer for Linear Programming: An Implementation of the Homogeneous Algorithm. Applied Optimization 33, 197-232.

3. Jacek Gondzio, Matrix-free interior point method, Computational Optimization and Applications, (2010).

# Direct Solver: Third Party Libraries

##A. Sparse Cholesky Decomposition

Sparse Cholesky Decomposition is currently supported by three alternative choices: 

### 1. CHOLMOD in SuiteSparse (Timothy A. Davis). 

**Website**: <http://www.suitesparse.com>

**Reference**: Algorithm 887: CHOLMOD, Supernodal Sparse Cholesky Factorization and Update/Downdate, Yanqing Chen, Timothy A. Davis, William W. Hager, Sivasankaran Rajamanickam, ACM Transactions on Mathematical Software, Vol 35, Issue 3, 2008, pp 22:1 - 22:14.

Both CPU and CPU-GPU version. You should install SuiteSparse yourself. 

CPU version: `mkcpu.sh`. You can adjust MKL_NUM_THREADS by the some command like `export MKL_NUM_THREADS=16`.

CPU-GPU version: `mkgpu.sh`. If NVIDIA GPU is not detected, it will turn to CPU version automatically. Currently only for large enough problems, CPU-GPU version can outperform CPU version. 

### 2. cuSOLVER in CUDAToolkit 7.5 (NVIDIA Corp.)

**Website**: <https://developer.nvidia.com/cusolver>

Only CPU-GPU version is provided. You should install CUDAToolkit (>= 7.5) yourself. 

CPU-GPU version: `mkcu.sh`. You need an NVIDIA GPU to run this program. 

### 3. Intel MKL PARDISO (Intel Corp.)

**Website**: <https://software.intel.com/en-us/node/470282>

Only CPU version is provided. You should install Intel MKL (>= 11.3) yourself.

CPU version: `mkmkl.sh`. You can adjust MKL_NUM_THREADS. 

Note that we provide a new version of ADAt calculation (preliminarily optimized), which is accelerated with OpenMP. 
You can adjust `OMP_NUM_THREADS` by the some command like `export OMP_NUM_THREADS=16`. 

**PLEASE GUARANTEE THAT `OMP_NUM_THREADS` IS NO MORE THAN `OMP_THREADS_MAX` DEFINED IN `LP.h`!**

##B. Sparse LU Decomposition

### LUSOL (Fortran 77 version by Michael Saunders, C Translation by Kjell Eikland)

**Website**: <http://web.stanford.edu/group/SOL/software/lusol/>

**Reference**: P. E. Gill, W. Murray, M. A. Saunders and M. H. Wright (1987). Maintaining LU factors of a general sparse matrix, Linear Algebra and its Applications 88/89, 239-270.

LUSOL supports Bartels-Golub-Reid updates for column replacement. Note that C Translation is used in another MILP solver: lp_solve (<http://lpsolve.sourceforge.net/5.5/>).

## C. Iterative Methods

These methods are used to solve the normal equation, (A D A^T) x = b. However, due to the ill-conditionness, the performance is far from satisfactory and is only used for fun. For testing, users can 

(1) Include "UsePCG_CPU.cpp" and "PCG.cpp" for a CPU implementation, this is a preconditioned conjugate gradient written by myself. 
For preconditioner, the following is referred:

J. Gondzio, Matrix-Free Interior Point Method, Computational Optimization and Applications 51 (2012) pp. 457-480.

(2) Another choice is use CG_DESCENT by including "UseCgDescent.cpp" and all files in subdirectory "CG_DESCENT", which is credit to

[1] W. W. Hager and H. Zhang, A new conjugate gradient method with guaranteed descent and an efficient line search, SIAM Journal on Optimization, 16 (2005), 170-192.
[2] W. W. Hager and H. Zhang, Algorithm 851: CG_DESCENT, A conjugate gradient method with guaranteed descent, ACM Transactions on Mathematical Software, 32 (2006), 113-137.
[3] W. W. Hager and H. Zhang, A survey of nonlinear conjugate gradient methods, Pacific Journal of Optimization, 2 (2006), pp. 35-58.
[4] W. W. Hager and H. Zhang, Limited memory conjugate gradients, www.math.ufl.edu/~hager/papers/CG/lcg.pdf
