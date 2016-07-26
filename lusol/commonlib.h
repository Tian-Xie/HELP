#ifndef HEADER_commonlib
#define HEADER_commonlib

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#define BIGNUMBER    1.0e+30
#define TINYNUMBER   1.0e-4
#define MACHINEPREC  2.22e-16
#define MATHPREC     1.0e-16
#define ERRLIMIT     1.0e-6

#ifndef LINEARSEARCH
  #define LINEARSEARCH 5
#endif

#if 0
  #define INTEGERTIME
#endif

/* ************************************************************************ */
/* Define sizes of standard number types                                    */
/* ************************************************************************ */
#ifndef LLONG
  #if defined __BORLANDC__
    #define LLONG __int64
  #elif !defined _MSC_VER || _MSC_VER >= 1310
    #define LLONG long long
  #else
    #define LLONG __int64
  #endif
#endif

#ifndef MYBOOL
  #if 0
    #define MYBOOL unsigned int
  #else
    #define MYBOOL unsigned char
  #endif
#endif

#ifndef REAL
  #define REAL     double
#endif
#ifndef BLAS_prec
  #define BLAS_prec "d" /* The BLAS precision prefix must correspond to the REAL type */
#endif

#ifndef REALXP
  #if 1
    #define REALXP long double  /* Set local accumulation variable as long double */
  #else
    #define REALXP REAL          /* Set local accumulation as default precision */
  #endif
#endif

#ifndef my_boolstr
  #define my_boolstr(x)          (!(x) ? "FALSE" : "TRUE")
#endif

#ifndef NULL
  #define NULL 	       0
#endif

#ifndef FALSE
  #define FALSE        0
  #define TRUE         1
#endif

#ifndef DOFASTMATH
  #define DOFASTMATH
#endif


#ifndef CALLOC
#define CALLOC(ptr, nr)\
  if(!((void *) ptr = calloc((size_t)(nr), sizeof(*ptr))) && nr) {\
    printf("calloc of %d bytes failed on line %d of file %s\n",\
           (size_t) nr * sizeof(*ptr), __LINE__, __FILE__);\
  }
#endif

#ifndef MALLOC
#define MALLOC(ptr, nr)\
  if(!((void *) ptr = malloc((size_t)((size_t) (nr) * sizeof(*ptr)))) && nr) {\
    printf("malloc of %d bytes failed on line %d of file %s\n",\
           (size_t) nr * sizeof(*ptr), __LINE__, __FILE__);\
  }
#endif

#ifndef REALLOC
#define REALLOC(ptr, nr)\
  if(!((void *) ptr = realloc(ptr, (size_t)((size_t) (nr) * sizeof(*ptr)))) && nr) {\
    printf("realloc of %d bytes failed on line %d of file %s\n",\
           (size_t) nr * sizeof(*ptr), __LINE__, __FILE__);\
  }
#endif

#ifndef FREE
#define FREE(ptr)\
  if((void *) ptr != NULL) {\
    free(ptr);\
    ptr = NULL; \
  }
#endif

#ifndef MEMCOPY
#define MEMCOPY(nptr, optr, nr)\
  memcpy((nptr), (optr), (size_t)((size_t)(nr) * sizeof(*(optr))))
#endif

#ifndef MEMMOVE
#define MEMMOVE(nptr, optr, nr)\
  memmove((nptr), (optr), (size_t)((size_t)(nr) * sizeof(*(optr))))
#endif

#ifndef MEMALLOCCOPY
#define MEMALLOCCOPY(nptr, optr, nr)\
  {MALLOC(nptr, (size_t)(nr));\
   MEMCOPY(nptr, optr, (size_t)(nr));}
#endif

#ifndef STRALLOCCOPY
#define STRALLOCCOPY(nstr, ostr)\
  {nstr = (char *) malloc((size_t) (strlen(ostr) + 1));\
   strcpy(nstr, ostr);}
#endif

#ifndef MEMCLEAR
/*#define useMMX*/
#ifdef useMMX
  #define MEMCLEAR(ptr, nr)\
    mem_set((ptr), '\0', (size_t)((size_t)(nr) * sizeof(*(ptr))))
#else
  #define MEMCLEAR(ptr, nr)\
    memset((ptr), '\0', (size_t)((size_t)(nr) * sizeof(*(ptr))))
#endif
#endif


#define MIN(x, y)         ((x) < (y) ? (x) : (y))
#define MAX(x, y)         ((x) > (y) ? (x) : (y))
#define SETMIN(x, y)      if((x) > (y)) x = y
#define SETMAX(x, y)      if((x) < (y)) x = y
#define LIMIT(lo, x, hi)  ((x < (lo) ? lo : ((x) > hi ? hi : x)))
#define IF(t, x, y)       ((t) ? (x) : (y))
#define SIGN(x)           ((x) < 0 ? -1 : 1)

#define DELTA_SIZE(newSize, oldSize) ((int) ((newSize) * pow(1.5, fabs((double)newSize)/((oldSize+newSize)+1))))

#ifndef CMP_CALLMODEL
#ifdef WIN32
# define CMP_CALLMODEL _cdecl
#else
# define CMP_CALLMODEL
#endif
#endif

#define CMP_ATTRIBUTES(item) (((char *) attributes)+(item)*recsize)
typedef int (CMP_CALLMODEL findCompare_func)(const void *current, const void *candidate);

#if 1

typedef struct _QSORTrec
{
  void     *self;  /* This could also hold a 32-bit integer */
  void     *prev;  /* prev+next could also hold a 64-bit double */
  void     *next;
} QSORTrec;
#define QSitem_double(QS) ((double *) &(QS.prev))
#define QSitem_float1(QS) ((float *) &(QS.prev))
#define QSitem_float2(QS) ((float *) &(QS.next))

#else

typedef struct _QSORTrec1
  void     *self;  /* This could also hold a 32-bit integer */
  void     *prev;  /* prev+next could also hold a 64-bit double */
  void     *next;
} QSORTrec1;
typedef struct _QSORTrec2
  int      intval;
  REAL     realval;
} QSORTrec2;
typedef struct _QSORTrec3
  int      intval;
  int      intpar1;
  int      intpar2;
} QSORTrec3;
typedef union _QSORTrec
{
  QSORTrec1 pointers;
  QSORTrec2 intreal;
  QSORTrec3 intintint;
} QSORTrec;

#endif


#ifdef __cplusplus
  extern "C" {
#endif

int mod(int n, int d);
int gcd(LLONG a, LLONG b, int *c, int *d);

int findIndex(int target, int *attributes, int count, int offset);
int findIndexEx(void *target, void *attributes, int count, int offset, int recsize, findCompare_func findCompare, MYBOOL ascending);

int CMP_CALLMODEL compareCHAR(const void *current, const void *candidate);
int CMP_CALLMODEL compareINT(const void *current, const void *candidate);
int CMP_CALLMODEL compareREAL(const void *current, const void *candidate);
void hpsort(void *attributes, int count, int offset, int recsize, MYBOOL descending, findCompare_func findCompare);
void hpsortex(void *attributes, int count, int offset, int recsize, MYBOOL descending, findCompare_func findCompare, int *tags);

int QS_addfirst(QSORTrec a[], void *mydata);
int QS_append(QSORTrec a[], int ipos, void *mydata);
void QS_replace(QSORTrec a[], int ipos, void *mydata);
void QS_insert(QSORTrec a[], int ipos, void *mydata, int epos);
void QS_delete(QSORTrec a[], int ipos, int epos);
void QS_swap(QSORTrec a[], int i, int j);
MYBOOL QS_execute(QSORTrec a[], int count, findCompare_func findCompare, MYBOOL islinkedlist, int *nswaps);

int sortByREAL(int *item, REAL *weight, int size, int offset, MYBOOL unique);
int sortByINT(int *item, int *weight, int size, int offset, MYBOOL unique);
REAL sortREALByINT(REAL *item, int *weight, int size, int offset, MYBOOL unique);


double timeNow(void);

void blockWriteBOOL(FILE *output, char *label, MYBOOL *myvector, int first, int last, MYBOOL asRaw);
void blockWriteINT(FILE *output, char *label, int *myvector, int first, int last);
void blockWriteREAL(FILE *output, char *label, REAL *myvector, int first, int last);

void printvec( int n, REAL *x, int modulo );
void printmatSQ( int size, int n, REAL *X, int modulo );
void printmatUT( int size, int n, REAL *U, int modulo );

#if defined _MSC_VER
int fileCount( char *filemask );
MYBOOL fileSearchPath( char *envvar, char *searchfile, char *foundpath );
#endif

#ifdef __cplusplus
  }
#endif

#endif /* HEADER_commonlib */

