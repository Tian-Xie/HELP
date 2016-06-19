# LPSolver
An experimental implementation of Homogeneous and Self-Dual Algorithm for Linear Programming. 

**Author**: Tian Xie (Research Center for Management Science and Information Analytics, Shanghai University of Finance and Economics)

**Progress**: Currently working on Presolve and Cross-over. HSD Algorithm without presolving is preliminarily completed. Some numerical issues about LLt factorization is known. (LDLt is more numerically robust.)

**Disclaimer**: This program is currently not finished testing. No warrenty. Use it at your own risk! 

#### Credits ####

Fundamental implementation idea originated from COPL_LP (Xiong Zhang and Yinyu Ye). See http://web.stanford.edu/~yyye/Col.html .

#### Third Party Libraries ####

Matrix Decomposition is currently supported by three alternative choices: 

###### 1. CHOLMOD in SuiteSparse (Timothy A. Davis). ######

**Website**: <http://www.suitesparse.com>

**Reference**: Algorithm 887: CHOLMOD, Supernodal Sparse Cholesky Factorization and Update/Downdate, Yanqing Chen, Timothy A. Davis, William W. Hager, Sivasankaran Rajamanickam, ACM Transactions on Mathematical Software, Vol 35, Issue 3, 2008, pp 22:1 - 22:14.

Both CPU and CPU-GPU version. You should install SuiteSparse yourself. 

CPU version: `mkcpu.sh`. You can adjust MKL_NUM_THREADS by the some command like `export MKL_NUM_THREADS=16`.

CPU-GPU version: `mkgpu.sh`. If NVIDIA GPU is not detected, it will turn to CPU version automatically. Currently only for large enough problems, CPU-GPU version can outperform CPU version. 

###### 2. cuSOLVER in CUDAToolkit 7.5 (NVIDIA Corp.) ######

**Website**: <https://developer.nvidia.com/cusolver>

Only CPU-GPU version is provided. You should install CUDAToolkit (>= 7.5) yourself. 

CPU-GPU version: `mkcu.sh`. You need an NVIDIA GPU to run this program. 

###### 3. Intel MKL PARDISO (Intel Corp.) ######

**Website**: <https://software.intel.com/en-us/node/470282>

Only CPU version is provided. You should install Intel MKL (>= 11.3) yourself.

CPU version: `mkmkl.sh`. You can adjust MKL_NUM_THREADS. 

Note that we provide a new version of ADAt calculation, which is accelerated with OpenMP. 
You can adjust `OMP_NUM_THREADS` by the some command like `export OMP_NUM_THREADS=16`. 
+PLEASE GUARANTEE THAT `OMP_NUM_THREADS` IS NO MORE THAN `OMP_THREADS_MAX` DEFINED IN `LP.h`!+