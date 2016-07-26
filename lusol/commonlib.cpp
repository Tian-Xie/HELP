
#include <sys/types.h>

#ifdef INTEGERTIME
# include <time.h>
#else
# include <sys/timeb.h>
#endif

#include <stdlib.h>
#include <stdio.h>
#ifdef WIN32
# include <io.h>       /* Used in file search functions */
#endif
#include <string.h>
#include <math.h>
#include "commonlib.h"

#ifdef FORTIFY
# include "lp_fortify.h"
#endif


/* Math operator equivalence function */
int mod(int n, int d)
{
  return(n % d);
}


/* Return the greatest common divisor of a and b, or -1 if it is 
   not defined. Return through the pointer arguments the integers
   such that gcd(a,b) = c*a + b*d. */
int gcd(LLONG a, LLONG b, int *c, int *d)
{
  LLONG q,r,t;
  int   cret,dret,C,D,rval, sgn_a = 1,sgn_b = 1, swap = 0;

  if((a == 0) || (b == 0)) 
    return( -1 );

  /* Use local multiplier instances, if necessary */
  if(c == NULL)
    c = &cret;
  if(d == NULL)
    d = &dret;

  /* Normalize so that 0 < a <= b */
  if(a < 0){
    a = -a;
    sgn_a = -1;
  }
  if(b < 0){
    b = -b;
    sgn_b = -1;
  }
  if(b < a){
    t = b;
    b = a;
    a = t;
    swap = 1;
  }

  /* Now a <= b and both >= 1. */
  q = b/a;
  r = b - a*q;
  if(r == 0) {
    if(swap){
      *d = 1;
      *c = 0;
    }
    else {
      *c = 1;
      *d = 0;
    }
    *c = sgn_a*(*c);
    *d = sgn_b*(*d);
    return( (int) a );
  }

  rval = gcd(a,r,&C,&D);
  if(swap){
    *d = (int) (C-D*q);
    *c = D;
  }
  else {
    *d = D;
    *c = (int) (C-D*q);
  }
  *c = sgn_a*(*c);
  *d = sgn_b*(*d);
  return( rval );
}

/* Array search functions */
int findIndex(int target, int *attributes, int count, int offset)
{
  int focusPos, beginPos, endPos;
  int focusAttrib, beginAttrib, endAttrib;

 /* Set starting and ending index offsets */
  beginPos = offset;
  endPos = beginPos + count - 1;
  if(endPos < beginPos)
    return(-1);

 /* Do binary search logic based on a sorted (decending) attribute vector */
  focusPos = (beginPos + endPos) / 2;
  beginAttrib = attributes[beginPos];
  focusAttrib = attributes[focusPos];
  endAttrib   = attributes[endPos];

  while(endPos - beginPos > LINEARSEARCH) {
    if(beginAttrib == target) {
      focusAttrib = beginAttrib;
      endPos = beginPos;
    }
    else if(endAttrib == target) {
      focusAttrib = endAttrib;
      beginPos = endPos;
    }
    else if(focusAttrib < target) {
      beginPos = focusPos + 1;
      beginAttrib = attributes[beginPos];
      focusPos = (beginPos + endPos) / 2;
      focusAttrib = attributes[focusPos];
    }
    else if(focusAttrib > target) {
      endPos = focusPos - 1;
      endAttrib = attributes[endPos];
      focusPos = (beginPos + endPos) / 2;
      focusAttrib = attributes[focusPos];
    }
    else {
      beginPos = focusPos;
      endPos = focusPos;
    }
  }

 /* Do linear (unsorted) search logic */
  if(endPos - beginPos <= LINEARSEARCH) {

    /* CPU intensive loop; provide alternative evaluation models */
#if defined DOFASTMATH
    /* Do fast pointer arithmetic */
    int *attptr = attributes + beginPos;
    while((beginPos < endPos) && ((*attptr) < target)) {
      beginPos++;
      attptr++;
    }
    focusAttrib = (*attptr);
#else
    /* Do traditional indexed access */
    focusAttrib = attributes[beginPos];
    while((beginPos < endPos) && (focusAttrib < target)) {
      beginPos++;
      focusAttrib = attributes[beginPos];
    }
#endif
  }

 /* Return the index if a match was found, or signal failure with a -1        */
  if(focusAttrib == target)             /* Found; return retrieval index      */
    return(beginPos);
  else if(focusAttrib > target)         /* Not found; last item               */
    return(-beginPos);
  else if(beginPos > offset+count-1)
    return(-(endPos+1));                /* Not found; end of list             */
  else
    return(-(beginPos+1));              /* Not found; intermediate point      */

}
int findIndexEx(void *target, void *attributes, int count, int offset, int recsize, findCompare_func findCompare, MYBOOL ascending)
{
  int  focusPos, beginPos, endPos, compare, order;
  void *focusAttrib, *beginAttrib, *endAttrib;

 /* Set starting and ending index offsets */
  beginPos = offset;
  endPos = beginPos + count - 1;
  if(endPos < beginPos)
    return(-1);
  order = (ascending ? -1 : 1);

 /* Do binary search logic based on a sorted attribute vector */
  focusPos = (beginPos + endPos) / 2;
  beginAttrib = CMP_ATTRIBUTES(beginPos);
  focusAttrib = CMP_ATTRIBUTES(focusPos);
  endAttrib   = CMP_ATTRIBUTES(endPos);

  compare = 0;
  while(endPos - beginPos > LINEARSEARCH) {
    if(findCompare(target, beginAttrib) == 0) {
      focusAttrib = beginAttrib;
      endPos = beginPos;
    }
    else if(findCompare(target, endAttrib) == 0) {
      focusAttrib = endAttrib;
      beginPos = endPos;
    }
    else {
      compare = findCompare(target, focusAttrib)*order;
      if(compare < 0) {
        beginPos = focusPos + 1;
        beginAttrib = CMP_ATTRIBUTES(beginPos);
        focusPos = (beginPos + endPos) / 2;
        focusAttrib = CMP_ATTRIBUTES(focusPos);
      }
      else if(compare > 0) {
        endPos = focusPos - 1;
        endAttrib = CMP_ATTRIBUTES(endPos);
        focusPos = (beginPos + endPos) / 2;
        focusAttrib = CMP_ATTRIBUTES(focusPos);
      }
      else {
        beginPos = focusPos;
        endPos = focusPos;
      }
    }
  }

 /* Do linear (unsorted) search logic */
  if(endPos - beginPos <= LINEARSEARCH) {

    /* Do traditional indexed access */
    focusAttrib = CMP_ATTRIBUTES(beginPos);
    if(beginPos == endPos)
      compare = findCompare(target, focusAttrib)*order;
    else
    while((beginPos < endPos) &&
          ((compare = findCompare(target, focusAttrib)*order) < 0)) {
      beginPos++;
      focusAttrib = CMP_ATTRIBUTES(beginPos);
    }
  }

 /* Return the index if a match was found, or signal failure with a -1        */
  if(compare == 0)                      /* Found; return retrieval index      */
    return(beginPos);
  else if(compare > 0)                  /* Not found; last item               */
    return(-beginPos);
  else if(beginPos > offset+count-1)
    return(-(endPos+1));                /* Not found; end of list             */
  else
    return(-(beginPos+1));              /* Not found; intermediate point      */

}

/* Simple sorting and searching comparison "operators" */
int CMP_CALLMODEL compareCHAR(const void *current, const void *candidate)
{
  if(*(char *) current < *(char *) candidate)
    return( -1 );
  else if(*(char *) current > *(char *) candidate)
    return( 1 );
  else
    return( 0 );
}
int CMP_CALLMODEL compareINT(const void *current, const void *candidate)
{
  if(*(int *) current < *(int *) candidate)
    return( -1 );
  else if(*(int *) current > *(int *) candidate)
    return( 1 );
  else
    return( 0 );
}
int CMP_CALLMODEL compareREAL(const void *current, const void *candidate)
{
  if(*(REAL *) current < *(REAL *) candidate)
    return( -1 );
  else if(*(REAL *) current > *(REAL *) candidate)
    return( 1 );
  else
    return( 0 );
}

/* Heap sort function (procedurally based on the Numerical Recipes version, 
   but expanded and generalized to hande any object with the use of 
   qsort-style comparison operator).  An expanded version is also implemented,
   where interchanges are reflected in a caller-initialized integer "tags" list. */
void hpsort(void *attributes, int count, int offset, int recsize, MYBOOL descending, findCompare_func findCompare)
{
  register int  i, j, k, ir, order;
  register char *hold, *base;
  char          *save;

  if(count < 2)
    return;
  offset -= 1;
  attributes = CMP_ATTRIBUTES(offset);
  base = CMP_ATTRIBUTES(1);
  save = (char *) malloc(recsize);
  if(descending)
    order = -1;
  else
    order = 1;
  
  k = (count >> 1) + 1;
  ir = count;

  for(;;) {
    if(k > 1) {
      MEMCOPY(save, CMP_ATTRIBUTES(--k), recsize);
    }
    else {
      hold = CMP_ATTRIBUTES(ir);
      MEMCOPY(save, hold, recsize);
      MEMCOPY(hold, base, recsize);
      if(--ir == 1) {
        MEMCOPY(base, save, recsize);
        break;
      }
    }

    i = k;
    j = k << 1;
    while(j <= ir) {
      hold = CMP_ATTRIBUTES(j);
      if( (j < ir) && (findCompare(hold, CMP_ATTRIBUTES(j+1))*order < 0) ) {
        hold += recsize;
        j++;
      }
      if(findCompare(save, hold)*order < 0) {
        MEMCOPY(CMP_ATTRIBUTES(i), hold, recsize);
        i = j;
        j <<= 1;
	    }
      else
        break;
    }
    MEMCOPY(CMP_ATTRIBUTES(i), save, recsize);
  }
  
  FREE(save);
}
void hpsortex(void *attributes, int count, int offset, int recsize, MYBOOL descending, findCompare_func findCompare, int *tags)
{
  if(count < 2)
    return;
  if(tags == NULL) {
    hpsort(attributes, count, offset, recsize, descending, findCompare);
    return;
  }
  else {
    register int  i, j, k, ir, order;
    register char *hold, *base;
    char          *save;
    int           savetag;
  
    offset -= 1;
    attributes = CMP_ATTRIBUTES(offset);
    tags += offset;
    base = CMP_ATTRIBUTES(1);
    save = (char *) malloc(recsize);
    if(descending)
      order = -1;
    else
      order = 1;
  
    k = (count >> 1) + 1;
    ir = count;

    for(;;) {
      if(k > 1) {
        MEMCOPY(save, CMP_ATTRIBUTES(--k), recsize);
        savetag = k;
      }
      else {
        hold = CMP_ATTRIBUTES(ir);
        MEMCOPY(save, hold, recsize);
        MEMCOPY(hold, base, recsize);
        savetag = tags[ir];
        tags[ir] = tags[1];
        if(--ir == 1) {
          MEMCOPY(base, save, recsize);
          tags[1] = savetag;
          break;
        }
      }

      i = k;
      j = k << 1;
      while(j <= ir) {
        hold = CMP_ATTRIBUTES(j);
        if( (j < ir) && (findCompare(hold, CMP_ATTRIBUTES(j+1))*order < 0) ) {
          hold += recsize;
          j++;
        }
        if(findCompare(save, hold)*order < 0) {
          MEMCOPY(CMP_ATTRIBUTES(i), hold, recsize);
          tags[i] = tags[j];
          i = j;
          j <<= 1;
  	    }
        else
          break;
      }
      MEMCOPY(CMP_ATTRIBUTES(i), save, recsize);
      tags[i] = savetag;
    }
  
    FREE(save);
  }
}


/* This is a "specialized generic" version of C.A.R Hoare's Quick Sort algorithm. 
   It will handle arrays that are already sorted, and arrays with duplicate keys. 
   The implementation here requires the user to pass a comparison operator and
   assumes that the array passed has the QSORTrec format, which i.a. includes
   the ability for to do linked list sorting. If the passed comparison operator 
   is NULL, the comparison is assumed to be for integers. */
#define QS_IS_switch 4    /* Threshold for switching to insertion sort */
int QS_addfirst(QSORTrec a[], void *mydata)
{
  a[0].self = mydata;
  return( 0 );
}
int QS_append(QSORTrec a[], int ipos, void *mydata)
{
  if(ipos <= 0)
    ipos = QS_addfirst(a, mydata);
  else
    a[ipos].self = mydata;
  return( ipos );
}
void QS_replace(QSORTrec a[], int ipos, void *mydata)
{
  a[ipos].self = mydata;
}
void QS_insert(QSORTrec a[], int ipos, void *mydata, int epos)
{
  for(; epos > ipos; epos--)
    a[epos] = a[epos-1];
  a[ipos].self = mydata;
}
void QS_delete(QSORTrec a[], int ipos, int epos)
{
  for(; epos > ipos; epos--)
    a[epos] = a[epos-1];
}
MYBOOL QS_validate(QSORTrec a[], int count)
{
  MYBOOL iserror = FALSE;

  count--;
  iserror = (a[0].prev != NULL) || (a[count].next != NULL);
  while((count > 0) && !iserror) {
    iserror = (a[count].prev != a[count-1].self) ||
              (a[count-1].next != a[count].self);
    count--;
  }
  if(iserror)
    printf("QS_validate: Error in linked list consistency at position %d\n", count);
  return( ! iserror );
}
void QS_updatelink(QSORTrec a[], int i, int j)
{
  register void *Pi1, *Pi2, *Pj1, *Pj2;

  /* Make sure j has the lowest index, i.e. was the previous bottom-most item */
  if(i < j) {
    int k = i;
    i = j;
    j = k;
  }
  Pi1 = a[i].prev;
  Pi2 = a[i].next;
  Pj1 = a[j].prev;
  Pj2 = a[j].next;

  /* Update new top-most item - j */
  if(j > 0)
    a[j-1].next = a[j].self;
  a[j].prev = Pi1;
  if(i - j > 1)
    a[j+1].prev = a[j].self;
  a[j].next = Pi2;

  /* Update new bottom-most item - i */
  if(Pj2 != NULL)
    a[i+1].prev = a[i].self;
  a[i].next = Pj2;
  if(i - j > 1)
    a[i-1].next = a[i].self;
  a[i].prev = Pj1;

}
void QS_swap(QSORTrec a[], int i, int j)
{
  QSORTrec T = a[i];
  a[i] = a[j];
  a[j] = T;
}
int QS_sort(QSORTrec a[], int l, int r, findCompare_func findCompare)
{
  register int i, j, nmove = 0;
  QSORTrec       v;

  /* Perform the a fast QuickSort */
  if((r-l) > QS_IS_switch) {
    i = (r+l)/2;

    /* Tri-Median Method */
    if(findCompare((char *) &a[l], (char *) &a[i]) > 0)
      { nmove++; QS_swap(a,l,i); } 
    if(findCompare((char *) &a[l], (char *) &a[r]) > 0) 
      { nmove++; QS_swap(a,l,r); }
    if(findCompare((char *) &a[i], (char *) &a[r]) > 0) 
      { nmove++; QS_swap(a,i,r); }

    j = r-1;
    QS_swap(a,i,j);
    i = l;
    v = a[j];
    for(;;) {
      while(findCompare((char *) &a[++i], (char *) &v) < 0);
      while(findCompare((char *) &a[--j], (char *) &v) > 0);
      if(j < i) break;
      nmove++; QS_swap (a,i,j);
    }
    nmove++; QS_swap(a,i,r-1);
    nmove += QS_sort(a,l,j,findCompare);
    nmove += QS_sort(a,i+1,r,findCompare);
  }
  return( nmove );
}
int QS_finish(QSORTrec a[], int lo0, int hi0, findCompare_func findCompare)
{
  int      i, j, nmove = 0;
  QSORTrec v;

  /* This is actually InsertionSort, which is faster for local sorts */
  for(i = lo0+1; i <= hi0; i++) {

    /* Save bottom-most item */
    v = a[i];

    /* Shift down! */
    j = i;
    while ((j > lo0) && (findCompare((char *) &a[j-1], (char *) &v) > 0)) {
      a[j] = a[j-1];
      j--;
      nmove++;
    }

    /* Store bottom-most item at the top */
    a[j] = v;
  }
  return( nmove );
}
MYBOOL QS_execute(QSORTrec a[], int count, findCompare_func findCompare, MYBOOL islinkedlist, int *nswaps)
{
  int iswaps = 0;

  /* Check and initialize */
  if(count <= 1)
    goto Finish;
  count--;

  /* Perform sort */
  iswaps = QS_sort(a, 0, count, findCompare);
#if QS_IS_switch > 0
  iswaps += QS_finish(a, 0, count, findCompare);
#endif

  /* Update linked list, if required */
  if(islinkedlist) {
    a[count].next = NULL;
    for(; count > 0; count--) {
      a[count].prev = a[count-1].self;
      a[count-1].next = a[count].self;
    }
    a[0].prev = NULL;
  }
#if 1
  else if(islinkedlist)
    QS_validate(a, count+1);
#endif

Finish:
  if(nswaps != NULL)
    *nswaps = iswaps;
  return( TRUE );
}

 

/* Simple specialized bubble/insertion sort functions */
int sortByREAL(int *item, REAL *weight, int size, int offset, MYBOOL unique)
{
  int i, ii, saveI;
  REAL saveW;

  for(i = 1; i < size; i++) {
    ii = i+offset-1;
    while ((ii >= offset) && (weight[ii] >= weight[ii+1])) {
      if(weight[ii] == weight[ii+1]) {
        if(unique)
          return(item[ii]);
      }
      else {
        saveI = item[ii];
        saveW = weight[ii];
        item[ii] = item[ii+1];
        weight[ii] = weight[ii+1];
        item[ii+1] = saveI;
        weight[ii+1] = saveW;
      }
      ii--;
    }
  }
  return(0);
}
int sortByINT(int *item, int *weight, int size, int offset, MYBOOL unique)
{
  int i, ii, saveI;
  int saveW;

  for(i = 1; i < size; i++) {
    ii = i+offset-1;
    while ((ii >= offset) && (weight[ii] >= weight[ii+1])) {
      if(weight[ii] == weight[ii+1]) {
        if(unique)
          return(item[ii]);
      }
      else {
        saveI = item[ii];
        saveW = weight[ii];
        item[ii] = item[ii+1];
        weight[ii] = weight[ii+1];
        item[ii+1] = saveI;
        weight[ii+1] = saveW;
      }
      ii--;
    }
  }
  return(0);
}
REAL sortREALByINT(REAL *item, int *weight, int size, int offset, MYBOOL unique)
{
  int  i, ii, saveW;
  REAL saveI;

  for(i = 1; i < size; i++) {
    ii = i+offset-1;
    while ((ii >= offset) && (weight[ii] >= weight[ii+1])) {
      if(weight[ii] == weight[ii+1]) {
        if(unique)
          return(item[ii]);
      }
      else {
        saveI = item[ii];
        saveW = weight[ii];
        item[ii] = item[ii+1];
        weight[ii] = weight[ii+1];
        item[ii+1] = saveI;
        weight[ii+1] = saveW;
      }
      ii--;
    }
  }
  return(0);
}


/* Time and message functions */
double timeNow(void)
{
#ifdef INTEGERTIME
  return((double)time(NULL));
#elif defined CLOCKTIME
  return((double)clock()/CLOCKS_PER_SEC /* CLK_TCK */);
#else
  struct timeb buf;

  ftime(&buf);
  return((double)buf.time+((double) buf.millitm)/1000.0);
#endif
}


/* Miscellaneous reporting functions */

/* List a vector of INT values for the given index range */
void blockWriteINT(FILE *output, char *label, int *myvector, int first, int last)
{
  int i, k = 0;

  fprintf(output, label);
  fprintf(output, "\n");
  for(i = first; i <= last; i++) {
    fprintf(output, " %5d", myvector[i]);
    k++;
    if(k % 12 == 0) {
      fprintf(output, "\n");
      k = 0;
    }
  }
  if(k % 12 != 0)
    fprintf(output, "\n");
}

/* List a vector of MYBOOL values for the given index range */
void blockWriteBOOL(FILE *output, char *label, MYBOOL *myvector, int first, int last, MYBOOL asRaw)
{
  int i, k = 0;

  fprintf(output, label);
  fprintf(output, "\n");
  for(i = first; i <= last; i++) {
    if(asRaw)
      fprintf(output, " %1d", myvector[i]);
    else
      fprintf(output, " %5s", my_boolstr(myvector[i]));
    k++;
    if(k % 36 == 0) {
      fprintf(output, "\n");
      k = 0;
    }
  }
  if(k % 36 != 0)
    fprintf(output, "\n");
}

/* List a vector of REAL values for the given index range */
void blockWriteREAL(FILE *output, char *label, REAL *myvector, int first, int last)
{
  int i, k = 0;

  fprintf(output, label);
  fprintf(output, "\n");
  for(i = first; i <= last; i++) {
    fprintf(output, " %18g", myvector[i]);
    k++;
    if(k % 4 == 0) {
      fprintf(output, "\n");
      k = 0;
    }
  }
  if(k % 4 != 0)
    fprintf(output, "\n");
}


/* CONSOLE vector and matrix printing routines */
void printvec( int n, REAL *x, int modulo )
{
  int i;

  if (modulo <= 0) modulo = 5;
  for (i = 1; i<=n; i++) {
    if(mod(i, modulo) == 1)
      printf("\n%2d:%12g", i, x[i]);
    else
      printf(" %2d:%12g", i, x[i]);
  }
  if(i % modulo != 0) printf("\n");
}


void printmatUT( int size, int n, REAL *U, int modulo)
{
   int i, ll;
   ll = 0;
   for(i = 1; i<=n; i++) {
     printvec(n-i+1, &U[ll], modulo);
     ll += size-i+1;
   }
}


void printmatSQ( int size, int n, REAL *X, int modulo)
{
   int i, ll;
   ll = 0;
   for(i = 1; i<=n; i++) {
     printvec(n, &X[ll], modulo);
     ll += size;
   }
}

/* Miscellaneous file functions */
#if defined _MSC_VER
/* Check MS versions before 7 */
#if _MSC_VER < 1300
# define intptr_t long
#endif

int fileCount( char *filemask )
{
  struct   _finddata_t c_file;
  intptr_t hFile;
  int      count = 0;

  /* Find first .c file in current directory */
  if( (hFile = _findfirst( filemask, &c_file )) == -1L )
    ;
  /* Iterate over all matching names */
  else {
     while( _findnext( hFile, &c_file ) == 0 )
       count++;
    _findclose( hFile );
  }
  return( count );
}
MYBOOL fileSearchPath( char *envvar, char *searchfile, char *foundpath )
{
   char pathbuffer[_MAX_PATH];

   _searchenv( searchfile, envvar, pathbuffer );
   if(pathbuffer[0] == '\0')
     return( FALSE );
   else {
     if(foundpath != NULL)
       strcpy(foundpath, pathbuffer);
     return( TRUE );
   }
}
#endif

