
/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   LUSOL routines from the Stanford Optimization Laboratory
   The parts included are:
    lusol1      Factor a given matrix A from scratch (lu1fac).
    lusol2      Heap-management routines for lu1fac.
    lusol6      Solve with the current LU factors.
    lusol7      Utilities for all update routines.
    lusol8      Replace a column (Bartels-Golub update).
   ------------------------------------------------------------------
   26 Apr 2002: TCP implemented using heap data structure.
   01 May 2002: lu1DCP implemented.
   07 May 2002: lu1mxc must put 0.0 at top of empty columns.
   09 May 2002: lu1mCP implements Markowitz with cols searched
                in heap order.
                Often faster (searching 20 or 40 cols) but more dense.
   11 Jun 2002: TRP implemented.
                lu1mRP implements Markowitz with Threshold Rook
                Pivoting.
                lu1mxc maintains max col elements  (was lu1max.)
                lu1mxr maintains max row elements.
   12 Jun 2002: lu1mCP seems too slow on big problems (e.g. memplus).
                Disabled it for the moment.  (Use lu1mar + TCP.)
   14 Dec 2002: TSP implemented.
                lu1mSP implements Markowitz with TSP.
   07 Mar 2003: character*1, character*2 changed to f90 form.
                Comments changed from * in column to ! in column 1.
                Comments kept within column 72 to avoid compiler
                warning.
   06 Mar 2004: Translation to C by Kjell Eikland with the addition
                of data wrappers, parametric constants, various
                helper routines, and dynamic memory reallocation.
   26 May 2004: Added LUSOL_IP_UPDATELIMIT parameter and provided
                for dynamic memory expansion based on possible
                forward requirements.
   08 Jul 2004: Revised logic in lu6chk based on new code from
                Michael Saunders.
   ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
/* #include <varargs.h>  For UNIX 5 compatibility */
#include <string.h>
#include <float.h>
#include <math.h>
#include "lusol.h"
#include "myblas.h"


/* LUSOL Object creation and destruction */

void *clean_realloc(void *oldptr, int width, int newsize, int oldsize)
{
  newsize *= width;
  oldsize *= width;
  oldptr = realloc(oldptr, newsize);
  if(newsize > oldsize)
/*    MEMCLEAR(oldptr+oldsize, newsize-oldsize); */
    memset((char *)oldptr+oldsize, '\0', newsize-oldsize);
  return(oldptr);
}

MYBOOL LUSOL_realloc_a(LUSOLrec *LUSOL, int newsize)
{
  int oldsize;

  if(newsize < 0)
    newsize = LUSOL->lena + MAX(abs(newsize), LUSOL_MINDELTA_a);

  oldsize = LUSOL->lena;
  LUSOL->lena = newsize;
  if(newsize > 0)
    newsize++;
  if(oldsize > 0)
    oldsize++;

  LUSOL->a    = (REAL *) clean_realloc(LUSOL->a,    sizeof(*(LUSOL->a)),
                                                    newsize, oldsize);
  LUSOL->indc = (int *)  clean_realloc(LUSOL->indc, sizeof(*(LUSOL->indc)),
                                                    newsize, oldsize);
  LUSOL->indr = (int *)  clean_realloc(LUSOL->indr, sizeof(*(LUSOL->indr)),
                                                    newsize, oldsize);
  if((newsize == 0) ||
     ((LUSOL->a != NULL) && (LUSOL->indc != NULL) && (LUSOL->indr != NULL)))
    return( TRUE );
  else
    return( FALSE );
}

MYBOOL LUSOL_expand_a(LUSOLrec *LUSOL, int *delta_lena, int *right_shift)
{
#ifdef StaticMemAlloc
  return( FALSE );
#else
  int LENA, NFREE, LFREE;

  /* Add expansion factor to avoid having to resize too often/too much;
     (exponential formula suggested by Michael A. Saunders) */
  LENA = LUSOL->lena;
  *delta_lena = DELTA_SIZE(*delta_lena, LENA);

  /* Expand it! */
  if((*delta_lena <= 0) || !LUSOL_realloc_a(LUSOL, LENA+(*delta_lena)))
    return( FALSE );

  /* Make sure we return the actual memory increase of a */
  *delta_lena = LUSOL->lena-LENA;

  /* Shift the used memory area to the right */
  LFREE = *right_shift;
  NFREE = LFREE+*delta_lena;
  LENA  -= LFREE-1;
  MEMMOVE(LUSOL->a+NFREE,    LUSOL->a+LFREE,    LENA);
  MEMMOVE(LUSOL->indr+NFREE, LUSOL->indr+LFREE, LENA);
  MEMMOVE(LUSOL->indc+NFREE, LUSOL->indc+LFREE, LENA);

  /* Also return the new starting position for the used memory area of a */
  *right_shift  = NFREE;

  LUSOL->expanded_a++;
  return( TRUE );
#endif
}

MYBOOL LUSOL_realloc_r(LUSOLrec *LUSOL, int newsize)
{
  int oldsize;

  if(newsize < 0)
    newsize = LUSOL->maxm + MAX(abs(newsize), LUSOL_MINDELTA_rc);

  oldsize = LUSOL->maxm;
  LUSOL->maxm = newsize;
  if(newsize > 0)
    newsize++;
  if(oldsize > 0)
    oldsize++;

  LUSOL->lenr  = (int *) clean_realloc(LUSOL->lenr,  sizeof(*(LUSOL->lenr)),
                                                     newsize, oldsize);
  LUSOL->ip    = (int *) clean_realloc(LUSOL->ip,    sizeof(*(LUSOL->ip)),
                                                     newsize, oldsize);
  LUSOL->iqloc = (int *) clean_realloc(LUSOL->iqloc, sizeof(*(LUSOL->iqloc)),
                                                     newsize, oldsize);
  LUSOL->ipinv = (int *) clean_realloc(LUSOL->ipinv, sizeof(*(LUSOL->ipinv)),
                                                     newsize, oldsize);
  LUSOL->locr  = (int *) clean_realloc(LUSOL->locr,  sizeof(*(LUSOL->locr)),
                                                     newsize, oldsize);

  if((newsize == 0) ||
     ((LUSOL->lenr != NULL) &&
      (LUSOL->ip != NULL) && (LUSOL->iqloc != NULL) &&
      (LUSOL->ipinv != NULL) && (LUSOL->locr != NULL))) {

#ifndef ClassicHamaxR
#ifdef AlwaysSeparateHamaxR
    if(LUSOL->luparm[LUSOL_IP_PIVOTTYPE] == LUSOL_PIVMOD_TRP)
#endif
    {
      LUSOL->amaxr = (REAL *) clean_realloc(LUSOL->amaxr, sizeof(*(LUSOL->amaxr)),
                                                          newsize, oldsize);
      if((newsize > 0) && (LUSOL->amaxr == NULL))
        return( FALSE );
    }
#endif
    return( TRUE );
  }
  else
    return( FALSE );
}

MYBOOL LUSOL_realloc_c(LUSOLrec *LUSOL, int newsize)
{
  int oldsize;

  if(newsize < 0)
    newsize = LUSOL->maxn + MAX(abs(newsize), LUSOL_MINDELTA_rc);

  oldsize = LUSOL->maxn;
  LUSOL->maxn = newsize;
  if(newsize > 0)
    newsize++;
  if(oldsize > 0)
    oldsize++;

  LUSOL->lenc  = (int *)  clean_realloc(LUSOL->lenc,  sizeof(*(LUSOL->lenc)),
                                                      newsize, oldsize);
  LUSOL->iq    = (int *)  clean_realloc(LUSOL->iq,    sizeof(*(LUSOL->iq)),
                                                      newsize, oldsize);
  LUSOL->iploc = (int *)  clean_realloc(LUSOL->iploc, sizeof(*(LUSOL->iploc)),
                                                      newsize, oldsize);
  LUSOL->iqinv = (int *)  clean_realloc(LUSOL->iqinv, sizeof(*(LUSOL->iqinv)),
                                                      newsize, oldsize);
  LUSOL->locc  = (int *)  clean_realloc(LUSOL->locc,  sizeof(*(LUSOL->locc)),
                                                      newsize, oldsize);
  LUSOL->w     = (REAL *) clean_realloc(LUSOL->w,     sizeof(*(LUSOL->w)),
                                                      newsize, oldsize);
#ifdef LUSOLSafeFastUpdate
  LUSOL->vLU6L = (REAL *) clean_realloc(LUSOL->vLU6L, sizeof(*(LUSOL->vLU6L)),
                                                      newsize, oldsize);
#else
  LUSOL->vLU6L = LUSOL->w;
#endif

  if((newsize == 0) ||
     ((LUSOL->w != NULL) && (LUSOL->lenc != NULL) &&
      (LUSOL->iq != NULL) && (LUSOL->iploc != NULL) &&
      (LUSOL->iqinv != NULL) && (LUSOL->locc != NULL))) {

#ifndef ClassicHamaxR
    if(LUSOL->luparm[LUSOL_IP_PIVOTTYPE] == LUSOL_PIVMOD_TCP) {
      LUSOL->Ha = (REAL *) clean_realloc(LUSOL->Ha,   sizeof(*(LUSOL->Ha)),
                                                      newsize, oldsize);
      LUSOL->Hj = (int *)  clean_realloc(LUSOL->Hj,   sizeof(*(LUSOL->Hj)),
                                                      newsize, oldsize);
      LUSOL->Hk = (int *)  clean_realloc(LUSOL->Hk,   sizeof(*(LUSOL->Hk)),
                                                      newsize, oldsize);
      if((newsize > 0) &&
         ((LUSOL->Ha == NULL) || (LUSOL->Hj == NULL) || (LUSOL->Hk == NULL)))
        return( FALSE );
    }
#endif
#ifndef ClassicdiagU
    if(LUSOL->luparm[LUSOL_IP_KEEPLU] == FALSE) {
      LUSOL->diagU = (REAL *) clean_realloc(LUSOL->diagU, sizeof(*(LUSOL->diagU)),
                                                          newsize, oldsize);
      if((newsize > 0) && (LUSOL->diagU == NULL))
        return( FALSE );
    }
#endif

    return( TRUE );
  }
  else
    return( FALSE );
}

LUSOLrec *LUSOL_create(FILE *outstream, int msgfil, int pivotmodel, int updatelimit)
{
  LUSOLrec *newLU;

  newLU = (LUSOLrec *) calloc(1, sizeof(*newLU));
  if(newLU == NULL)
    return( newLU );

  newLU->luparm[LUSOL_IP_SCALAR_NZA]       = LUSOL_MULT_nz_a;
  newLU->outstream = outstream;
  newLU->luparm[LUSOL_IP_PRINTUNIT]        = msgfil;
  newLU->luparm[LUSOL_IP_PRINTLEVEL]       = LUSOL_MSG_SINGULARITY;

  LUSOL_setpivotmodel(newLU, pivotmodel, LUSOL_PIVTOL_DEFAULT);

  newLU->parmlu[LUSOL_RP_GAMMA]            = LUSOL_DEFAULT_GAMMA;

  newLU->parmlu[LUSOL_RP_ZEROTOLERANCE]    = 3.0e-13;

  newLU->parmlu[LUSOL_RP_SMALLDIAG_U]      = /*3.7e-11;*/
  newLU->parmlu[LUSOL_RP_EPSDIAG_U]        = 3.7e-11;

  newLU->parmlu[LUSOL_RP_COMPSPACE_U]      = 3.0e+0;

  newLU->luparm[LUSOL_IP_MARKOWITZ_MAXCOL] = 5;
  newLU->parmlu[LUSOL_RP_MARKOWITZ_CONLY]  = 0.3e+0;
  newLU->parmlu[LUSOL_RP_MARKOWITZ_DENSE]  = 0.5e+0;

  newLU->parmlu[LUSOL_RP_SMARTRATIO]       = LUSOL_DEFAULT_SMARTRATIO;
#ifdef ForceRowBasedL0
  newLU->luparm[LUSOL_IP_USEROWL0]         = LUSOL_BASEORDER;
#endif
  newLU->luparm[LUSOL_IP_KEEPLU]           = TRUE;
  newLU->luparm[LUSOL_IP_UPDATELIMIT]      = updatelimit;

  init_BLAS();

  return( newLU );
}

MYBOOL LUSOL_sizeto(LUSOLrec *LUSOL, int init_r, int init_c, int init_a)
{
  if(LUSOL_realloc_a(LUSOL, init_a) &&
     LUSOL_realloc_r(LUSOL, init_r) &&
     LUSOL_realloc_c(LUSOL, init_c))
    return( TRUE );
  else
    return( FALSE );
}

char *LUSOL_pivotLabel(LUSOLrec *LUSOL)
{
  static /*const*/ char *pivotText[LUSOL_PIVMOD_MAX+1] =
  {"TPP", "TRP", "TCP", "TSP"};
  return(pivotText[LUSOL->luparm[LUSOL_IP_PIVOTTYPE]]);
}

void LUSOL_setpivotmodel(LUSOLrec *LUSOL, int pivotmodel, int initlevel)
{
  REAL newFM, newUM;

  /* Set pivotmodel if specified */
  if(pivotmodel > LUSOL_PIVMOD_NOCHANGE) {
    if((pivotmodel <= LUSOL_PIVMOD_DEFAULT) || (pivotmodel > LUSOL_PIVMOD_MAX))
      pivotmodel = LUSOL_PIVMOD_TPP;
    LUSOL->luparm[LUSOL_IP_PIVOTTYPE]        = pivotmodel;
  }

  /* Check if we need bother about changing tolerances */
  if((initlevel <= LUSOL_PIVTOL_NOCHANGE) || (initlevel > LUSOL_PIVTOL_MAX))
    return;

  /* Set default pivot tolerances
     (note that UPDATEMAX should always be <= FACTORMAX) */
  if(initlevel == LUSOL_PIVTOL_BAGGY) {        /* Extra-loose pivot thresholds */
    newFM = 500.0;
    newUM = newFM / 20;
  }
  else if(initlevel == LUSOL_PIVTOL_LOOSE) {  /* Moderately tight pivot tolerances */
    newFM = 100.0;
    newUM = newFM / 10;
  }
  else if(initlevel == LUSOL_PIVTOL_NORMAL) { /* Standard pivot tolerances */
    newFM = 28.0;
    newUM = newFM / 4;
  }
  else if(initlevel == LUSOL_PIVTOL_SLIM) {   /* Better accuracy pivot tolerances */
    newFM = 10.0;
    newUM = newFM / 2;
  }
  else if(initlevel == LUSOL_PIVTOL_TIGHT) {  /* Enhanced accuracy pivot tolerances */
    newFM = 5.0;
    newUM = newFM / 2;
  }
  else {                                      /* Very tight pivot tolerances for extra accuracy */
    newFM = 2.5;
    newUM = newFM / 2;
  }

  /* Set the tolerances */
  LUSOL->parmlu[LUSOL_RP_FACTORMAX_Lij] = newFM;
  LUSOL->parmlu[LUSOL_RP_UPDATEMAX_Lij] = newUM;
}

MYBOOL LUSOL_tightenpivot(LUSOLrec *LUSOL)
{
  REAL newvalue;

  /* Give up tightening if we are already less than limit and we cannot change strategy */
  if(MIN(LUSOL->parmlu[LUSOL_RP_UPDATEMAX_Lij],
         LUSOL->parmlu[LUSOL_RP_UPDATEMAX_Lij]) < 1.1) {
    if(LUSOL->luparm[LUSOL_IP_PIVOTTYPE] >= LUSOL_PIVMOD_TRP)
      return( FALSE );
    LUSOL_setpivotmodel(LUSOL, LUSOL->luparm[LUSOL_IP_PIVOTTYPE]+1, LUSOL_PIVTOL_DEFAULT+1);
    return( 2 );
  }

  /* Otherwise tighten according to defined schedule */
#if 0   /* This is Michael Saunder's proposed tightening procedure */
  newvalue = sqrt(LUSOL->parmlu[LUSOL_RP_FACTORMAX_Lij]);
  LUSOL->parmlu[LUSOL_RP_FACTORMAX_Lij] = newvalue;
  LUSOL->parmlu[LUSOL_RP_UPDATEMAX_Lij] = MIN(newvalue,
                                              LUSOL->parmlu[LUSOL_RP_UPDATEMAX_Lij]);
#elif 0 /* This is Kjell Eikland's schedule #1 */
  newvalue = sqrt(LUSOL->parmlu[LUSOL_RP_FACTORMAX_Lij]);
  LUSOL->parmlu[LUSOL_RP_FACTORMAX_Lij] = newvalue;
  LUSOL->parmlu[LUSOL_RP_UPDATEMAX_Lij] = 1.0 + (newvalue - 1.0) / 2;
#else   /* This was Kjell Eikland's schedule #2 */
  LUSOL->parmlu[LUSOL_RP_FACTORMAX_Lij] = 1.0 + LUSOL->parmlu[LUSOL_RP_FACTORMAX_Lij]/3.0;
  LUSOL->parmlu[LUSOL_RP_UPDATEMAX_Lij] = 1.0 + LUSOL->parmlu[LUSOL_RP_UPDATEMAX_Lij]/3.0;
#endif
  return( TRUE );
}


char *LUSOL_informstr(LUSOLrec *LUSOL, int inform)
{
  static char *informText[LUSOL_INFORM_MAX-LUSOL_INFORM_MIN+1] =
  {"LUSOL_RANKLOSS: Lost rank",
   "LUSOL_LUSUCCESS: Success",
   "LUSOL_LUSINGULAR: Singular A",
   "LUSOL_LUUNSTABLE: Unstable factorization",
   "LUSOL_ADIMERR: Row or column count exceeded",
   "LUSOL_ADUPLICATE: Duplicate A matrix entry found",
   "",
   "",
   "LUSOL_ANEEDMEM: Insufficient memory for factorization",
   "LUSOL_FATALERR: Fatal internal error",
   "LUSOL_NOPIVOT: Found no suitable pivot",
   "LUSOL_NOMEMLEFT: Could not obtain more memory"};
  if(inform < LUSOL_INFORM_MIN || inform > LUSOL_INFORM_MAX)
    inform = /* LUSOL->luparm[LUSOL_IP_INFORM] */ LUSOL_INFORM_FATALERR;
  return(informText[inform-LUSOL_INFORM_MIN]);
}

void LUSOL_clear(LUSOLrec *LUSOL, MYBOOL nzonly)
{
  int len;

  LUSOL->nelem = 0;
  if(!nzonly) {

   /* lena arrays */
    len = LUSOL->lena + LUSOL_ARRAYOFFSET;
    MEMCLEAR(LUSOL->a,    len);
    MEMCLEAR(LUSOL->indc, len);
    MEMCLEAR(LUSOL->indr, len);

   /* maxm arrays */
    len = LUSOL->maxm + LUSOL_ARRAYOFFSET;
    MEMCLEAR(LUSOL->lenr,  len);
    MEMCLEAR(LUSOL->ip,    len);
    MEMCLEAR(LUSOL->iqloc, len);
    MEMCLEAR(LUSOL->ipinv, len);
    MEMCLEAR(LUSOL->locr,  len);

#ifndef ClassicHamaxR
    if((LUSOL->amaxr != NULL)
#ifdef AlwaysSeparateHamaxR
       && (LUSOL->luparm[LUSOL_IP_PIVOTTYPE] == LUSOL_PIVMOD_TRP)
#endif
      )
      MEMCLEAR(LUSOL->amaxr, len);
#endif

   /* maxn arrays */
    len = LUSOL->maxn + LUSOL_ARRAYOFFSET;
    MEMCLEAR(LUSOL->lenc,  len);
    MEMCLEAR(LUSOL->iq,    len);
    MEMCLEAR(LUSOL->iploc, len);
    MEMCLEAR(LUSOL->iqinv, len);
    MEMCLEAR(LUSOL->locc,  len);
    MEMCLEAR(LUSOL->w,     len);

    if(LUSOL->luparm[LUSOL_IP_PIVOTTYPE] == LUSOL_PIVMOD_TCP) {
      MEMCLEAR(LUSOL->Ha,  len);
      MEMCLEAR(LUSOL->Hj,  len);
      MEMCLEAR(LUSOL->Hk,  len);
    }
#ifndef ClassicdiagU
    if(LUSOL->luparm[LUSOL_IP_KEEPLU] == FALSE) {
      MEMCLEAR(LUSOL->diagU, len);
    }
#endif

  }
}


MYBOOL LUSOL_assign(LUSOLrec *LUSOL, int iA[], int jA[], REAL Aij[], int nzcount, MYBOOL istriplet)
{
  int k, m, n, ij, kol;

  /* Adjust the size of the a structure */
  if(nzcount > (LUSOL->lena/LUSOL->luparm[LUSOL_IP_SCALAR_NZA]) &&
     !LUSOL_realloc_a(LUSOL, nzcount*LUSOL->luparm[LUSOL_IP_SCALAR_NZA]))
    return( FALSE );

  m = 0;
  n = 0;
  kol = 1;
  for(k = 1; k <= nzcount; k++) {
    /* First the row indicator */
    ij = iA[k];
    if(ij > m) {
      m = ij;
      if(m > LUSOL->maxm &&
         !LUSOL_realloc_r(LUSOL, -(m / LUSOL_MINDELTA_FACTOR + 1)))
        return( FALSE );
    }
    LUSOL->indc[k] = ij;

    /* Then the column indicator;
       Handle both triplet and column count formats */
    if(istriplet)
      ij = jA[k];
    else {
      if(k >= jA[kol])
        kol++;
      ij = kol;
    }
    if(ij > n) {
      n = ij;
      if(n > LUSOL->maxn &&
         !LUSOL_realloc_c(LUSOL, -(n / LUSOL_MINDELTA_FACTOR + 1)))
        return( FALSE );
    }
    LUSOL->indr[k] = ij;

    /* Lastly the matrix value itself */
    LUSOL->a[k] = Aij[k];
  }
  LUSOL->m = m;
  LUSOL->n = n;
  LUSOL->nelem = nzcount;
  return( TRUE );
}

int LUSOL_loadColumn(LUSOLrec *LUSOL, int iA[], int jA, REAL Aij[], int nzcount, int offset1)
{
  int i, ii, nz, k;

  nz = LUSOL->nelem;
  i = nz + nzcount;
  if(i > (LUSOL->lena/LUSOL->luparm[LUSOL_IP_SCALAR_NZA]) &&
     !LUSOL_realloc_a(LUSOL, i*LUSOL->luparm[LUSOL_IP_SCALAR_NZA]))
  return( -1 );

  k = 0;
  for(ii = 1; ii <= nzcount; ii++) {
    i = ii + offset1;
    if(Aij[i] == 0)
      continue;
    if(iA[i] <= 0 || iA[i] > LUSOL->m ||
       jA <= 0 || jA > LUSOL->n) {
      LUSOL_report(LUSOL, 0, "Variable index outside of set bounds (r:%d/%d, c:%d/%d)\n",
                             iA[i], LUSOL->m, jA, LUSOL->n);
      continue;
    }
    k++;
    nz++;
    LUSOL->a[nz]    = Aij[i];
    LUSOL->indc[nz] = iA[i];
    LUSOL->indr[nz] = jA;
  }
  LUSOL->nelem = nz;
  return( k );
}

void LUSOL_free(LUSOLrec *LUSOL)
{
  LUSOL_realloc_a(LUSOL, 0);
  LUSOL_realloc_r(LUSOL, 0);
  LUSOL_realloc_c(LUSOL, 0);
  if(LUSOL->L0 != NULL)
    LUSOL_matfree(&(LUSOL->L0));
  if(!is_nativeBLAS())
    unload_BLAS();
  free(LUSOL);
}

void LUSOL_report(LUSOLrec *LUSOL, int msglevel, char *format, ...)
{
  va_list ap;

  va_start(ap, format);
  if(LUSOL == NULL) {
    vfprintf(stderr, format, ap);
  }
  else if(msglevel >= 0  /*LUSOL->luparm[2]*/) {
    if(LUSOL->writelog != NULL) {
      char buff[255];

      vsprintf(buff, format, ap);
      LUSOL->writelog(LUSOL, LUSOL->loghandle, buff);
    }
    if(LUSOL->outstream != NULL) {
      vfprintf(LUSOL->outstream, format, ap);
      fflush(LUSOL->outstream);
    }
  }
  va_end(ap);
}

void LUSOL_timer(LUSOLrec *LUSOL, int timerid, char *text)
{
  LUSOL_report(LUSOL, -1, "TimerID %d at %s - %s\n",
                          timerid, "", text);
}

int LUSOL_factorize(LUSOLrec *LUSOL)
{
  int inform;

  LU1FAC( LUSOL, &inform );
  return( inform );
}

int LUSOL_ftran(LUSOLrec *LUSOL, REAL b[], int NZidx[], MYBOOL prepareupdate)
{
  int  inform;
  REAL *vector;

  if(prepareupdate)
    vector = LUSOL->vLU6L;
  else
    vector = LUSOL->w;

  /* Copy RHS vector, but make adjustment for offset since this
     can create a memory error when the calling program uses
     a 0-base vector offset back to comply with LUSOL. */
  MEMCOPY(vector+1, b+1, LUSOL->n);
  vector[0] = 0;

  LU6SOL(LUSOL, LUSOL_SOLVE_Aw_v, vector, b, NZidx, &inform);
  LUSOL->luparm[LUSOL_IP_FTRANCOUNT]++;

  return(inform);
}


int LUSOL_btran(LUSOLrec *LUSOL, REAL b[], int NZidx[])
{
  int inform;

  /* Copy RHS vector, but make adjustment for offset since this
     can create a memory error when the calling program uses
     a 0-base vector offset back to comply with LUSOL. */
  MEMCOPY(LUSOL->w+1, b+1, LUSOL->m);
  LUSOL->w[0] = 0;

  LU6SOL(LUSOL, LUSOL_SOLVE_Atv_w, b, LUSOL->w, NZidx, &inform);
  LUSOL->luparm[LUSOL_IP_BTRANCOUNT]++;

  return(inform);
}


int LUSOL_replaceColumn(LUSOLrec *LUSOL, int jcol, REAL v[])
{
  int  inform;
  REAL DIAG, VNORM;

  LU8RPC(LUSOL, LUSOL_UPDATE_OLDNONEMPTY, LUSOL_UPDATE_NEWNONEMPTY,
                jcol, v, NULL,
                &inform, &DIAG, &VNORM);

  LUSOL->replaced_c++;
  return( inform );
}

int LUSOL_findColumnPosition(LUSOLrec *LUSOL, int jcol)
/* The purpose of this routine is to find the slack row/column in
   user-index that was singular in the last unsuccessful column
   update; zero is returned if the search was unsuccessful.
   (Source: Michael A. Saunders; private communication to KE) */
{
  int j;

#if 1 /* Michael Saunders version */
  for(j = LUSOL->m; j > 0; j--)
    if(LUSOL->iq[j] == jcol)
      break;
#else /* Kjell Eikland version (note that iqinv could be invalid) */
  j = LUSOL->iqinv[jcol];
#endif
  return( j );
}

REAL LUSOL_vecdensity(LUSOLrec *LUSOL, REAL V[])
{
  int I, N = 0;

  for(I = 1; I <= LUSOL->m; I++)
    if(fabs(V[I]) > 0)
      N++;
  return( (REAL) N / (REAL) LUSOL->m );
}

char relationChar(REAL left, REAL right)
{
  if(left > right)
    return('>');
  else if(left == right)
    return('=');
  else
    return('<');
}

/* Retrieve the core modules ordered by order of dependency */

//#include "lusol2.c"      /* Heap management */

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   File  lusol2 LUSOL heap management routines
   Hbuild   Hchange  Hdelete  Hdown    Hinsert  Hup
   Heap-management routines for LUSOL's lu1fac.
   May be useful for other applications.
   ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   For LUSOL, the heap structure involves three arrays of length N.
   N        is the current number of entries in the heap.
   Ha(1:N)  contains the values that the heap is partially sorting.
            For LUSOL they are double precision values -- the largest
            element in each remaining column of the updated matrix.
            The biggest entry is in Ha(1), the top of the heap.
   Hj(1:N)  contains column numbers j.
            Ha(k) is the biggest entry in column j = Hj(k).
   Hk(1:N)  contains indices within the heap.  It is the
            inverse of Hj(1:N), so  k = Hk(j)  <=>  j = Hj(k).
            Column j is entry k in the heap.
   hops     is the number of heap operations,
            i.e., the number of times an entry is moved
            (the number of "hops" up or down the heap).
   Together, Hj and Hk let us find values inside the heap
   whenever we want to change one of the values in Ha.
   For other applications, Ha may need to be some other data type,
   like the keys that sort routines operate on.
   ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   11 Feb 2002: MATLAB  version derived from "Algorithms" by
                R. Sedgewick
   03 Mar 2002: F77     version derived from MATLAB version.
   07 May 2002: Safeguard input parameters k, N, Nk.
                We don't want them to be output!
   07 May 2002: Current version of lusol2.f.
   ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

/* ==================================================================
   Hdown  updates heap by moving down tree from node k.
   ------------------------------------------------------------------
   01 May 2002: Need Nk for length of Hk.
   05 May 2002: Change input paramter k to kk to stop k being output.
   05 May 2002: Current version of Hdown.
   ================================================================== */
void HDOWN(REAL HA[], int HJ[], int HK[], int N, int K, int *HOPS)
{
  int  J, JJ, JV, N2;
  REAL V;

  *HOPS = 0;
  V = HA[K];
  JV = HJ[K];
 N2 = N/2;
/*      while 1
        break */
x100:
  if(K>N2)
    goto x200;
  (*HOPS)++;
  J = K+K;
  if(J<N) {
    if(HA[J]<HA[J+1])
      J++;
  }
/*      break */
  if(V>=HA[J])
    goto x200;
  HA[K] = HA[J];
  JJ = HJ[J];
  HJ[K] = JJ;
  HK[JJ] = K;
  K = J;
  goto x100;
/*      end while */
x200:
  HA[K] = V;
  HJ[K] = JV;
  HK[JV] = K;
}

/* ==================================================================
   Hup updates heap by moving up tree from node k.
   ------------------------------------------------------------------
   01 May 2002: Need Nk for length of Hk.
   05 May 2002: Change input paramter k to kk to stop k being output.
   05 May 2002: Current version of Hup.
   ================================================================== */
void HUP(REAL HA[], int HJ[], int HK[], int K, int *HOPS)
{
  int  J, JV, K2;
  REAL V;

  *HOPS = 0;
  V = HA[K];
 JV = HJ[K];
/*      while 1
        break */
x100:
  if(K<2)
    goto x200;
  K2 = K/2;
/*      break */
  if(V<HA[K2])
    goto x200;
  (*HOPS)++;
  HA[K] = HA[K2];
  J = HJ[K2];
  HJ[K] = J;
  HK[J] = K;
  K = K2;
  goto x100;
/*      end while */
x200:
  HA[K] = V;
  HJ[K] = JV;
  HK[JV] = K;
}

/* ==================================================================
   Hinsert inserts (v,jv) into heap of length N-1
   to make heap of length N.
   ------------------------------------------------------------------
   03 Apr 2002: First version of Hinsert.
   01 May 2002: Require N to be final length, not old length.
                Need Nk for length of Hk.
   07 May 2002: Protect input parameters N, Nk.
   07 May 2002: Current version of Hinsert.
   ================================================================== */
void HINSERT(REAL HA[], int HJ[], int HK[], int N,
             REAL V, int JV, int *HOPS)
{
  HA[N] = V;
  HJ[N] = JV;
  HK[JV] = N;
  HUP(HA,HJ,HK,N,HOPS);
}

/* ==================================================================
   Hchange changes Ha(k) to v in heap of length N.
   ------------------------------------------------------------------
   01 May 2002: Need Nk for length of Hk.
   07 May 2002: Protect input parameters N, Nk, k.
   07 May 2002: Current version of Hchange.
   ================================================================== */
void HCHANGE(REAL HA[], int HJ[], int HK[], int N, int K,
             REAL V, int JV, int *HOPS)
{
  REAL V1;

  V1 = HA[K];
  HA[K] = V;
  HJ[K] = JV;
  HK[JV] = K;
  if(V1<V)
    HUP  (HA,HJ,HK,  K,HOPS);
  else
    HDOWN(HA,HJ,HK,N,K,HOPS);
}

/* ==================================================================
   Hdelete deletes Ha(k) from heap of length N.
   ------------------------------------------------------------------
   03 Apr 2002: Current version of Hdelete.
   01 May 2002: Need Nk for length of Hk.
   07 May 2002: Protect input parameters N, Nk, k.
   07 May 2002: Current version of Hdelete.
   ================================================================== */
void HDELETE(REAL HA[], int HJ[], int HK[], int *N, int K, int *HOPS)
{

  int  JV, NX;
  REAL V;

  NX = *N;
  V = HA[NX];
  JV = HJ[NX];
  (*N)--;
  *HOPS = 0;
  if(K<NX)
    HCHANGE(HA,HJ,HK,NX,K,V,JV,HOPS);
}

/* ==================================================================
   Hbuild initializes the heap by inserting each element of Ha.
   Input:  Ha, Hj.
   Output: Ha, Hj, Hk, hops.
   ------------------------------------------------------------------
   01 May 2002: Use k for new length of heap, not k-1 for old length.
   05 May 2002: Use kk in call to stop loop variable k being altered.
                (Actually Hinsert no longer alters that parameter.)
   07 May 2002: ftnchek wants us to protect Nk, Ha(k), Hj(k) too.
   07 May 2002: Current version of Hbuild.
   ================================================================== */
void HBUILD(REAL HA[], int HJ[], int HK[], int N, int *HOPS)
{
  int  H, JV, K, KK;
  REAL V;

  *HOPS = 0;
  for(K = 1; K <= N; K++) {
    KK = K;
    V = HA[K];
    JV = HJ[K];
    HINSERT(HA,HJ,HK,KK,V,JV,&H);
    (*HOPS) += H;
  }
}



//#include "lusol6a.c"     /* Singularity checking and equation solving */

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   File  lusol6a
      lu6sol   lu6L     lu6Lt     lu6U     Lu6Ut   lu6LD   lu6chk
   ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   26 Apr 2002: lu6 routines put into a separate file.
   15 Dec 2002: lu6sol modularized via lu6L, lu6Lt, lu6U, lu6Ut.
                lu6LD implemented to allow solves with LDL' or L|D|L'.
   15 Dec 2002: Current version of lusol6a.f.
   ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

/* ==================================================================
   lu6chk  looks at the LU factorization  A = L*U.
   If mode = 1, lu6chk is being called by lu1fac.
   (Other modes not yet implemented.)
   ------------------------------------------------------------------
   The important input parameters are

                  lprint = luparm(2)
                  keepLU = luparm(8)
                  Utol1  = parmlu(4)
                  Utol2  = parmlu(5)

   and the significant output parameters are

                  inform = luparm(10)
                  nsing  = luparm(11)
                  jsing  = luparm(12)
                  jumin  = luparm(19)
                  Lmax   = parmlu(11)
                  Umax   = parmlu(12)
                  DUmax  = parmlu(13)
                  DUmin  = parmlu(14)
                  and      w(*).

   Lmax  and Umax  return the largest elements in L and U.
   DUmax and DUmin return the largest and smallest diagonals of U
                   (excluding diagonals that are exactly zero).
   In general, w(j) is set to the maximum absolute element in
   the j-th column of U.  However, if the corresponding diagonal
   of U is small in absolute terms or relative to w(j)
   (as judged by the parameters Utol1, Utol2 respectively),
   then w(j) is changed to - w(j).
   Thus, if w(j) is not positive, the j-th column of A
   appears to be dependent on the other columns of A.
   The number of such columns, and the position of the last one,
   are returned as nsing and jsing.
   Note that nrank is assumed to be set already, and is not altered.
   Typically, nsing will satisfy      nrank + nsing = n,  but if
   Utol1 and Utol2 are rather large,  nsing > n - nrank   may occur.
   If keepLU = 0,
   Lmax  and Umax  are already set by lu1fac.
   The diagonals of U are in the top of A.
   Only Utol1 is used in the singularity test to set w(*).
   inform = 0  if  A  appears to have full column rank  (nsing = 0).
   inform = 1  otherwise  (nsing .gt. 0).
   ------------------------------------------------------------------
   00 Jul 1987: Early version.
   09 May 1988: f77 version.
   11 Mar 2001: Allow for keepLU = 0.
   17 Nov 2001: Briefer output for singular factors.
   05 May 2002: Comma needed in format 1100 (via Kenneth Holmstrom).
   06 May 2002: With keepLU = 0, diags of U are in natural order.
                They were not being extracted correctly.
   23 Apr 2004: TRP can judge singularity better by comparing
                all diagonals to DUmax.
   27 Jun 2004: (PEG) Allow write only if nout .gt. 0.
   ================================================================== */
#ifdef UseOld_LU6CHK_20040510
void LU6CHK(LUSOLrec *LUSOL, int MODE, int LENA2, int *INFORM)
{
  MYBOOL KEEPLU;
  int    I, J, JSING, JUMIN, K, L, L1, L2, LENL, LPRINT, NDEFIC, NRANK, NSING;
  REAL   AIJ, DIAG, DUMAX, DUMIN, LMAX, UMAX, UTOL1, UTOL2;

  LPRINT = LUSOL->luparm[LUSOL_IP_PRINTLEVEL];
  KEEPLU = (MYBOOL) (LUSOL->luparm[LUSOL_IP_KEEPLU]!=0);
  NRANK = LUSOL->luparm[LUSOL_IP_RANK_U];
  LENL  = LUSOL->luparm[LUSOL_IP_NONZEROS_L];
  UTOL1 = LUSOL->parmlu[LUSOL_RP_SMALLDIAG_U];
  UTOL2 = LUSOL->parmlu[LUSOL_RP_EPSDIAG_U];
  *INFORM = LUSOL_INFORM_LUSUCCESS;
  LMAX  = ZERO;
  UMAX  = ZERO;
  NSING = 0;
  JSING = 0;
  JUMIN = 0;
  DUMAX = ZERO;
  DUMIN = LUSOL_BIGNUM;

#ifdef LUSOLFastClear
  MEMCLEAR(LUSOL->w, LUSOL->n + 1);
#else
  for(I = 1; I <= LUSOL->n; I++)
    LUSOL->w[I] = ZERO;
#endif

  if(KEEPLU) {
/*     --------------------------------------------------------------
        Find  Lmax.
       -------------------------------------------------------------- */
    for(L = (LENA2+1)-LENL; L <= LENA2; L++) {
      LMAX = MAX(LMAX,fabs(LUSOL->a[L]));
     }
/*     --------------------------------------------------------------
        Find Umax and set w(j) = maximum element in j-th column of U.
       -------------------------------------------------------------- */
    for(K = 1; K <= NRANK; K++) {
      I = LUSOL->ip[K];
      L1 = LUSOL->locr[I];
      L2 = (L1+LUSOL->lenr[I])-1;
      for(L = L1; L <= L2; L++) {
        J = LUSOL->indr[L];
        AIJ = fabs(LUSOL->a[L]);
        LUSOL->w[J] = MAX(LUSOL->w[J],AIJ);
        UMAX = MAX(UMAX,AIJ);
      }
    }
/*     --------------------------------------------------------------
        Negate w(j) if the corresponding diagonal of U is
        too small in absolute terms or relative to the other elements
        in the same column of  U.
        Also find DUmax and DUmin, the extreme diagonals of U.
       -------------------------------------------------------------- */
    for(K = 1; K <= LUSOL->n; K++) {
      J = LUSOL->iq[K];
      if(K>NRANK)
        DIAG = ZERO;
      else {
        I = LUSOL->ip[K];
        L1 = LUSOL->locr[I];
        DIAG = fabs(LUSOL->a[L1]);
        DUMAX = MAX(DUMAX,DIAG);
        if(DUMIN>DIAG) {
          DUMIN = DIAG;
          JUMIN = J;
        }
      }
      if(DIAG<=UTOL1 || DIAG<=UTOL2*LUSOL->w[J]) {
        NSING++;
        JSING = J;
        LUSOL->w[J] = -LUSOL->w[J];
      }
    }
    LUSOL->parmlu[LUSOL_RP_MAXMULT_L] = LMAX;
    LUSOL->parmlu[LUSOL_RP_MAXELEM_U] = UMAX;
  }
   else {
/*     --------------------------------------------------------------
        keepLU = 0.
        Only diag(U) is stored.  Set w(*) accordingly.
       -------------------------------------------------------------- */
    for(K = 1; K <= LUSOL->n; K++) {
      J = LUSOL->iq[K];
      if(K>NRANK)
        DIAG = ZERO;
      else {
/* !             diag   = abs( diagU(k) ) ! 06 May 2002: Diags are in natural order */
        DIAG = fabs(LUSOL->diagU[J]);
        LUSOL->w[J] = DIAG;
        DUMAX = MAX(DUMAX,DIAG);
        if(DUMIN>DIAG) {
          DUMIN = DIAG;
          JUMIN = J;
        }
      }
      if(DIAG<=UTOL1) {
        NSING++;
        JSING = J;
        LUSOL->w[J] = -LUSOL->w[J];
      }
    }
  }
/*     -----------------------------------------------------------------
        Set output parameters.
       ----------------------------------------------------------------- */
  if(JUMIN==0)
    DUMIN = ZERO;
  LUSOL->luparm[LUSOL_IP_SINGULARITIES]  = NSING;
  LUSOL->luparm[LUSOL_IP_SINGULARINDEX]  = JSING;
  LUSOL->luparm[LUSOL_IP_COLINDEX_DUMIN] = JUMIN;
  LUSOL->parmlu[LUSOL_RP_MAXELEM_DIAGU]  = DUMAX;
  LUSOL->parmlu[LUSOL_RP_MINELEM_DIAGU]  = DUMIN;
/*      The matrix has been judged singular. */
  if(NSING>0) {
    *INFORM = LUSOL_INFORM_LUSINGULAR;
    NDEFIC = LUSOL->n-NRANK;
    if(LPRINT>=LUSOL_MSG_SINGULARITY) {
      LUSOL_report(LUSOL, 0, "Singular(m%cn)  rank:%9d  n-rank:%8d  nsing:%9d\n",
                             relationChar(LUSOL->m, LUSOL->n),NRANK,NDEFIC,NSING);
    }
  }
/*      Exit. */
  LUSOL->luparm[LUSOL_IP_INFORM] = *INFORM;
}
#else
void LU6CHK(LUSOLrec *LUSOL, int MODE, int LENA2, int *INFORM)
{
  MYBOOL KEEPLU, TRP;
  int    I, J, JSING, JUMIN, K, L, L1, L2, LENL, LDIAGU, LPRINT, NDEFIC, NRANK, NSING;
  REAL   AIJ, DIAG, DUMAX, DUMIN, LMAX, UMAX, UTOL1, UTOL2;

  LPRINT = LUSOL->luparm[LUSOL_IP_PRINTLEVEL];
  KEEPLU = (MYBOOL) (LUSOL->luparm[LUSOL_IP_KEEPLU] != 0);
  TRP    = (MYBOOL) (LUSOL->luparm[LUSOL_IP_PIVOTTYPE] == LUSOL_PIVMOD_TRP);
  NRANK  = LUSOL->luparm[LUSOL_IP_RANK_U];
  LENL   = LUSOL->luparm[LUSOL_IP_NONZEROS_L];
  UTOL1  = LUSOL->parmlu[LUSOL_RP_SMALLDIAG_U];
  UTOL2  = LUSOL->parmlu[LUSOL_RP_EPSDIAG_U];
  *INFORM = LUSOL_INFORM_LUSUCCESS;
  LMAX   = ZERO;
  UMAX   = ZERO;
  NSING  = 0;
  JSING  = 0;
  JUMIN  = 0;
  DUMAX  = ZERO;
  DUMIN  = LUSOL_BIGNUM;

#ifdef LUSOLFastClear
  MEMCLEAR(LUSOL->w, LUSOL->n + 1);
#else
  for(I = 1; I <= LUSOL->n; I++)
    LUSOL->w[I] = ZERO;
#endif

  if(KEEPLU) {
/*     --------------------------------------------------------------
        Find  Lmax.
       -------------------------------------------------------------- */
    for(L = (LENA2+1)-LENL; L <= LENA2; L++) {
      LMAX = MAX(LMAX,fabs(LUSOL->a[L]));
     }
/*     --------------------------------------------------------------
        Find Umax and set w(j) = maximum element in j-th column of U.
       -------------------------------------------------------------- */
    for(K = 1; K <= NRANK; K++) {
      I = LUSOL->ip[K];
      L1 = LUSOL->locr[I];
      L2 = (L1+LUSOL->lenr[I])-1;
      for(L = L1; L <= L2; L++) {
        J = LUSOL->indr[L];
        AIJ = fabs(LUSOL->a[L]);
        LUSOL->w[J] = MAX(LUSOL->w[J],AIJ);
        UMAX = MAX(UMAX,AIJ);
      }
    }
    LUSOL->parmlu[LUSOL_RP_MAXMULT_L] = LMAX;
    LUSOL->parmlu[LUSOL_RP_MAXELEM_U] = UMAX;
/*     --------------------------------------------------------------
       Find DUmax and DUmin, the extreme diagonals of U.
       -------------------------------------------------------------- */
    for(K = 1; K <= NRANK; K++) {
      J     = LUSOL->iq[K];
      I     = LUSOL->ip[K];
      L1    = LUSOL->locr[I];
      DIAG  = fabs(LUSOL->a[L1]);
      DUMAX = MAX( DUMAX, DIAG );
      if(DUMIN > DIAG) {
        DUMIN  = DIAG;
        JUMIN  = J;
      }
    }
  }
  else {
/*     --------------------------------------------------------------
       keepLU = 0.
       Only diag(U) is stored.  Set w(*) accordingly.
       Find DUmax and DUmin, the extreme diagonals of U.
       -------------------------------------------------------------- */
    LDIAGU = LENA2 - LUSOL->n;
    for(K = 1; K <= NRANK; K++) {
      J           = LUSOL->iq[K];
      DIAG        = fabs( LUSOL->a[LDIAGU + J] ); /* are in natural order */
      LUSOL->w[J] = DIAG;
      DUMAX  = MAX( DUMAX, DIAG );
      if(DUMIN > DIAG) {
        DUMIN = DIAG;
        JUMIN = J;
      }
    }
  }
/*     --------------------------------------------------------------
       Negate w(j) if the corresponding diagonal of U is
       too small in absolute terms or relative to the other elements
       in the same column of  U.
      
       23 Apr 2004: TRP ensures that diags are NOT small relative to
                    other elements in their own column.
                    Much better, we can compare all diags to DUmax.
      -------------------------------------------------------------- */
  if((MODE == 1) && TRP) 
    UTOL1 = MAX( UTOL1, UTOL2*DUMAX );

  if(KEEPLU) {
    for(K = 1; K <= LUSOL->n; K++) {
      J = LUSOL->iq[K];
      if(K>NRANK)
        DIAG = ZERO;
      else {
        I = LUSOL->ip[K];
        L1 = LUSOL->locr[I];
        DIAG = fabs(LUSOL->a[L1]);
      }
      if((DIAG<=UTOL1) || (DIAG<=UTOL2*LUSOL->w[J])) {
        NSING++;
        JSING = J;
        LUSOL->w[J] = -LUSOL->w[J];
      }
    }
  }
  else { /* keepLU = 0 */
    for(K = 1; K <= LUSOL->n; K++) {
      J = LUSOL->iq[K];
      DIAG = LUSOL->w[J];
      if(DIAG<=UTOL1) {
        NSING++;
        JSING = J;
        LUSOL->w[J] = -LUSOL->w[J];
      }
    }
  }
/*     -----------------------------------------------------------------
        Set output parameters.
       ----------------------------------------------------------------- */
  if(JUMIN==0)
    DUMIN = ZERO;
  LUSOL->luparm[LUSOL_IP_SINGULARITIES]  = NSING;
  LUSOL->luparm[LUSOL_IP_SINGULARINDEX]  = JSING;
  LUSOL->luparm[LUSOL_IP_COLINDEX_DUMIN] = JUMIN;
  LUSOL->parmlu[LUSOL_RP_MAXELEM_DIAGU]  = DUMAX;
  LUSOL->parmlu[LUSOL_RP_MINELEM_DIAGU]  = DUMIN;
/*      The matrix has been judged singular. */
  if(NSING>0) {
    *INFORM = LUSOL_INFORM_LUSINGULAR;
    NDEFIC = LUSOL->n-NRANK;
    if((LUSOL->outstream!=NULL) && (LPRINT>=LUSOL_MSG_SINGULARITY)) {
      LUSOL_report(LUSOL, 0, "Singular(m%cn)  rank:%9d  n-rank:%8d  nsing:%9d\n",
                             relationChar(LUSOL->m, LUSOL->n),NRANK,NDEFIC,NSING);
    }
  }
/*      Exit. */
  LUSOL->luparm[LUSOL_IP_INFORM] = *INFORM;
}
#endif


/* ------------------------------------------------------------------
   Include routines for row-based L0.
   20 Apr 2005 Current version - KE.
   ------------------------------------------------------------------ */
//#include "lusol6l0.c"

/* Create a row-based version of L0.
   This makes it possible to solve L0'x=h (btran) faster for sparse h,
   since we only run down the columns of L0' (rows of LO) for which
   the corresponding entry in h is non-zero. */
MYBOOL LU1L0(LUSOLrec *LUSOL, LUSOLmat **mat, int *inform)
{
  MYBOOL status = FALSE;
  int    K, L, LL, L1, L2, LENL0, NUML0, I;
  int    *lsumr;

  /* Assume success */
  *inform = LUSOL_INFORM_LUSUCCESS;

  /* Check if there is anything worth doing */
  if(mat == NULL)
    return( status );
  if(*mat != NULL)
    LUSOL_matfree(mat);
  NUML0 = LUSOL->luparm[LUSOL_IP_COLCOUNT_L0];
  LENL0 = LUSOL->luparm[LUSOL_IP_NONZEROS_L0];
  if((NUML0 == 0) || (LENL0 == 0) || (LUSOL->luparm[LUSOL_IP_USEROWL0] == LUSOL_BASEORDER))
    return( status );

  /* Allocate temporary array */
  lsumr = (int *) calloc((LUSOL->m+1), sizeof(*lsumr));
  if(lsumr == NULL) {
    *inform = LUSOL_INFORM_NOMEMLEFT;
    return( status );
  }

  /* Compute non-zero counts by permuted row index (order is unimportant) */
  K = 0;
  L2 = LUSOL->lena;
  L1 = L2-LENL0+1;
  for(L = L1; L <= L2; L++) {
    I = LUSOL->indc[L];
    lsumr[I]++;
    if(lsumr[I] == 1)
      K++;
  }
  LUSOL->luparm[LUSOL_IP_ROWCOUNT_L0] = K;

  /* Check if we should apply "smarts" before proceeding to the row matrix creation */
  if((LUSOL->luparm[LUSOL_IP_USEROWL0] == LUSOL_AUTOORDER) &&
     ((REAL) LUSOL->luparm[LUSOL_IP_ROWCOUNT_L0] / 
#if 0
             LUSOL->luparm[LUSOL_IP_COLCOUNT_L0]
#else
             LUSOL->m
#endif
      > LUSOL->parmlu[LUSOL_RP_SMARTRATIO]))
    goto Finish;

  /* We are Ok to create the new matrix object */
  *mat = LUSOL_matcreate(LUSOL->m, LENL0);
  if(*mat == NULL) {
    *inform = LUSOL_INFORM_NOMEMLEFT;
    goto Finish;
  }

  /* Cumulate row counts to get vector offsets; first row is leftmost
     (stick with Fortran array offset for consistency) */
  (*mat)->lenx[0] = 1;
  for(K = 1; K <= LUSOL->m; K++) {
    (*mat)->lenx[K] = (*mat)->lenx[K-1] + lsumr[K];
    lsumr[K] = (*mat)->lenx[K-1];
  }

  /* Map the matrix into row order by permuted index;
     Note: The first permuted row is located leftmost in the array.
           The column order is irrelevant, since the indeces will
           refer to constant / resolved values of V[] during solve. */
  L2 = LUSOL->lena;
  L1 = L2-LENL0+1;
  for(L = L1; L <= L2; L++) {
    I = LUSOL->indc[L];
    LL = lsumr[I]++;
    (*mat)->a[LL] = LUSOL->a[L];
    (*mat)->indr[LL] = LUSOL->indr[L];
    (*mat)->indc[LL] = I;
  }

  /* Pack row starting positions, and set mapper from original index to packed */
  I = 0;
  for(L = 1; L <= LUSOL->m; L++) {
    K = LUSOL->ip[L];
    if((*mat)->lenx[K] > (*mat)->lenx[K-1]) {
	  I++;
      (*mat)->indx[I] = K;
	}
  }

  /* Confirm that everything went well */
  status = TRUE;

  /* Clean up */
Finish:
  FREE(lsumr);
  return( status );
}

/* Solve L0' v = v based on row-based version of L0, constructed by LU1L0 */
void LU6L0T_v(LUSOLrec *LUSOL, LUSOLmat *mat, REAL V[], int NZidx[], int *INFORM)
{
#ifdef DoTraceL0
  REAL TEMP;
#endif
  int  LEN, K, KK, L, L1, NUML0;
  REAL SMALL;
  register REAL VPIV;
#if (defined LUSOLFastSolve) && !(defined DoTraceL0)
  REAL *aptr;
  int  *jptr;
#else
  int  J;
#endif

  NUML0 = LUSOL->luparm[LUSOL_IP_ROWCOUNT_L0];
  SMALL = LUSOL->parmlu[LUSOL_RP_ZEROTOLERANCE];

  /* Loop over the nz columns of L0' - from the end, going forward. */
  for(K = NUML0; K > 0; K--) {
    KK = mat->indx[K];
    L  = mat->lenx[KK];
    L1 = mat->lenx[KK-1];
    LEN = L - L1;
    if(LEN == 0)
      continue;
    /* Get value of the corresponding active entry of V[] */
    VPIV = V[KK];
    /* Only process the column of L0' if the value of V[] is non-zero */
    if(fabs(VPIV)>SMALL) {
/*     ***** This loop could be coded specially. */
#if (defined LUSOLFastSolve) && !(defined DoTraceL0)
      L--;
      for(aptr = mat->a+L, jptr = mat->indr+L;
          LEN > 0; LEN--, aptr--, jptr--)
        V[*jptr] += (*aptr) * VPIV;
#else
      for(; LEN > 0; LEN--) {
        L--;
        J = mat->indr[L];
#ifndef DoTraceL0
        V[J] += mat->a[L]*VPIV;
#else
        TEMP = V[J];
        V[J] += mat->a[L]*VPIV;
        printf("V[%3d] = V[%3d] + L[%d,%d]*V[%3d]\n", J, J, KK,J, KK);
        printf("%6g = %6g + %6g*%6g\n", V[J], TEMP, mat->a[L], VPIV);
#endif
      }
#endif
    }
#ifdef SetSmallToZero
    else
      V[KK] = 0;
#endif
  }

}



/* ------------------------------------------------------------------
   lu6L   solves   L v = v(input).
   ------------------------------------------------------------------
   15 Dec 2002: First version derived from lu6sol.
   15 Dec 2002: Current version.
   ------------------------------------------------------------------ */
void LU6L(LUSOLrec *LUSOL, int *INFORM, REAL V[], int NZidx[])
{
  int  JPIV, K, L, L1, LEN, LENL, LENL0, NUML, NUML0;
  REAL SMALL;
  register REAL VPIV;
#ifdef LUSOLFastSolve
  REAL *aptr;
  int  *iptr, *jptr;
#else
  int  I, J;
#endif

  NUML0 = LUSOL->luparm[LUSOL_IP_COLCOUNT_L0];
  LENL0 = LUSOL->luparm[LUSOL_IP_NONZEROS_L0];
  LENL  = LUSOL->luparm[LUSOL_IP_NONZEROS_L];
  SMALL = LUSOL->parmlu[LUSOL_RP_ZEROTOLERANCE];
  *INFORM = LUSOL_INFORM_LUSUCCESS;
  L1 = LUSOL->lena+1;
  for(K = 1; K <= NUML0; K++) {
    LEN = LUSOL->lenc[K];
    L = L1;
    L1 -= LEN;
    JPIV = LUSOL->indr[L1];
    VPIV = V[JPIV];
    if(fabs(VPIV)>SMALL) {
/*     ***** This loop could be coded specially. */
#ifdef LUSOLFastSolve
      L--;
      for(aptr = LUSOL->a+L, iptr = LUSOL->indc+L;
          LEN > 0; LEN--, aptr--, iptr--)
        V[*iptr] += (*aptr) * VPIV;
#else
      for(; LEN > 0; LEN--) {
        L--;
        I = LUSOL->indc[L];
        V[I] += LUSOL->a[L]*VPIV;
      }
#endif
    }
#ifdef SetSmallToZero
    else
      V[JPIV] = 0;
#endif
  }
  L = (LUSOL->lena-LENL0)+1;
  NUML = LENL-LENL0;
/*     ***** This loop could be coded specially. */
#ifdef LUSOLFastSolve
  L--;
  for(aptr = LUSOL->a+L, jptr = LUSOL->indr+L, iptr = LUSOL->indc+L;
      NUML > 0; NUML--, aptr--, jptr--, iptr--) {
    if(fabs(V[*jptr])>SMALL)
      V[*iptr] += (*aptr) * V[*jptr];
#ifdef SetSmallToZero
    else
      V[*jptr] = 0;
#endif
  }
#else
  for(; NUML > 0; NUML--) {
    L--;
    J = LUSOL->indr[L];
    if(fabs(V[J])>SMALL) {
      I = LUSOL->indc[L];
      V[I] += LUSOL->a[L]*V[J];
    }
#ifdef SetSmallToZero
    else
      V[J] = 0;
#endif
  }
#endif
/*      Exit. */
  LUSOL->luparm[LUSOL_IP_INFORM] = *INFORM;
}

/* ==================================================================
   lu6LD  assumes lu1fac has computed factors A = LU of a
   symmetric definite or quasi-definite matrix A,
   using Threshold Symmetric Pivoting (TSP),   luparm(6) = 3,
   or    Threshold Diagonal  Pivoting (TDP),   luparm(6) = 4.
   It also assumes that no updates have been performed.
   In such cases,  U = D L', where D = diag(U).
   lu6LDL returns v as follows:

   mode
    1    v  solves   L D v = v(input).
    2    v  solves   L|D|v = v(input).
   ------------------------------------------------------------------
   15 Dec 2002: First version of lu6LD.
   15 Dec 2002: Current version.
   ================================================================== */
void LU6LD(LUSOLrec *LUSOL, int *INFORM, int MODE, REAL V[], int NZidx[])
{
  int  IPIV, K, L, L1, LEN, NUML0;
  REAL DIAG, SMALL;
  register REAL VPIV;
#ifdef LUSOLFastSolve
  REAL *aptr;
  int  *jptr;
#else
  int  J;
#endif

/*      Solve L D v(new) = v  or  L|D|v(new) = v, depending on mode.
        The code for L is the same as in lu6L,
        but when a nonzero entry of v arises, we divide by
        the corresponding entry of D or |D|. */
  NUML0 = LUSOL->luparm[LUSOL_IP_COLCOUNT_L0];
  SMALL = LUSOL->parmlu[LUSOL_RP_ZEROTOLERANCE];
  *INFORM = LUSOL_INFORM_LUSUCCESS;
  L1 = LUSOL->lena+1;
  for(K = 1; K <= NUML0; K++) {
    LEN = LUSOL->lenc[K];
    L = L1;
    L1 -= LEN;
    IPIV = LUSOL->indr[L1];
    VPIV = V[IPIV];
    if(fabs(VPIV)>SMALL) {
/*     ***** This loop could be coded specially. */
#ifdef LUSOLFastSolve
      L--;
      for(aptr = LUSOL->a+L, jptr = LUSOL->indc+L;
          LEN > 0; LEN--, aptr--, jptr--)
        V[*jptr] += (*aptr)*VPIV;
#else
      for(; LEN > 0; LEN--) {
        L--;
        J = LUSOL->indc[L];
        V[J] += LUSOL->a[L]*VPIV;
      }
#endif
/*      Find diag = U(ipiv,ipiv) and divide by diag or |diag|. */
      L = LUSOL->locr[IPIV];
      DIAG = LUSOL->a[L];
      if(MODE==2)
        DIAG = fabs(DIAG);
      V[IPIV] = VPIV/DIAG;
    }
#ifdef SetSmallToZero
    else
      V[IPIV] = 0;
#endif
  }
}


/* ==================================================================
   lu6Lt  solves   L'v = v(input).
   ------------------------------------------------------------------
   15 Dec 2002: First version derived from lu6sol.
   15 Dec 2002: Current version.
   ================================================================== */
void LU6LT(LUSOLrec *LUSOL, int *INFORM, REAL V[], int NZidx[])
{
#ifdef DoTraceL0
  REAL    TEMP;
#endif
  int     K, L, L1, L2, LEN, LENL, LENL0, NUML0;
  REAL    SMALL;
  register REALXP SUM;
  register REAL HOLD;
#if (defined LUSOLFastSolve) && !(defined DoTraceL0)
  REAL    *aptr;
  int     *iptr, *jptr;
#else
  int     I, J;
#endif

  NUML0 = LUSOL->luparm[LUSOL_IP_COLCOUNT_L0];
  LENL0 = LUSOL->luparm[LUSOL_IP_NONZEROS_L0];
  LENL  = LUSOL->luparm[LUSOL_IP_NONZEROS_L];
  SMALL = LUSOL->parmlu[LUSOL_RP_ZEROTOLERANCE];
  *INFORM = LUSOL_INFORM_LUSUCCESS;
  L1 = (LUSOL->lena-LENL)+1;
  L2 = LUSOL->lena-LENL0;
  
/*     ***** This loop could be coded specially. */
#if (defined LUSOLFastSolve) && !(defined DoTraceL0)
  for(L = L1, aptr = LUSOL->a+L1, iptr = LUSOL->indr+L1, jptr = LUSOL->indc+L1;
      L <= L2; L++, aptr++, iptr++, jptr++) {
    HOLD = V[*jptr];
    if(fabs(HOLD)>SMALL)
      V[*iptr] += (*aptr)*HOLD;
#ifdef SetSmallToZero
    else
      V[*jptr] = 0;
#endif
  }
#else
  for(L = L1; L <= L2; L++) {
    J = LUSOL->indc[L];
    HOLD = V[J];
    if(fabs(HOLD)>SMALL) {
      I = LUSOL->indr[L];
      V[I] += LUSOL->a[L]*HOLD;
    }
#ifdef SetSmallToZero
    else
      V[J] = 0;
#endif
  }
#endif

  /* Do row-based L0 version, if available */
  if((LUSOL->L0 != NULL) || 
     ((LUSOL->luparm[LUSOL_IP_BTRANCOUNT] == 0) && LU1L0(LUSOL, &(LUSOL->L0), INFORM))) {
    LU6L0T_v(LUSOL, LUSOL->L0, V, NZidx, INFORM);
  }

  /* Alternatively, do the standard column-based L0 version */
  else  {
    /* Perform loop over columns */
    for(K = NUML0; K >= 1; K--) {
      SUM = ZERO;
      LEN = LUSOL->lenc[K];
      L1 = L2+1;
      L2 += LEN;
/*     ***** This loop could be coded specially. */
#if (defined LUSOLFastSolve) && !(defined DoTraceL0)
      for(L = L1, aptr = LUSOL->a+L1, jptr = LUSOL->indc+L1;
          L <= L2; L++, aptr++, jptr++)
        SUM += (*aptr) * V[*jptr];
#else
      for(L = L1; L <= L2; L++) {
        J = LUSOL->indc[L];
#ifndef DoTraceL0
        SUM += LUSOL->a[L]*V[J];
#else
        TEMP = V[LUSOL->indr[L1]] + SUM;
        SUM += LUSOL->a[L]*V[J];
        printf("V[%3d] = V[%3d] + L[%d,%d]*V[%3d]\n", LUSOL->indr[L1], LUSOL->indr[L1], J,LUSOL->indr[L1], J);
        printf("%6g = %6g + %6g*%6g\n", V[LUSOL->indr[L1]] + SUM, TEMP, LUSOL->a[L], V[J]);
#endif
      }
#endif
      V[LUSOL->indr[L1]] += SUM;
    }
  }

/*      Exit. */
  LUSOL->luparm[LUSOL_IP_INFORM] = *INFORM;
}

void print_L0(LUSOLrec *LUSOL)
{
  int  I, J, K, L, L1, L2, LEN, LENL0, NUML0;
  REAL *denseL0 = (REAL*) calloc(LUSOL->m+1, (LUSOL->n+1)*sizeof(*denseL0));

  NUML0 = LUSOL->luparm[LUSOL_IP_COLCOUNT_L0];
  LENL0 = LUSOL->luparm[LUSOL_IP_NONZEROS_L0];

  L2 = LUSOL->lena-LENL0;
  for(K = NUML0; K >= 1; K--) {
    LEN = LUSOL->lenc[K];
    L1 = L2+1;
    L2 += LEN;
    for(L = L1; L <= L2; L++) {
      I = LUSOL->indc[L];
      I = LUSOL->ipinv[I]; /* Undo row mapping */
      J = LUSOL->indr[L];
      denseL0[(LUSOL->n+1)*(J-1) + I] = LUSOL->a[L];
    }
  }

  for(I = 1; I <= LUSOL->n; I++) {
    for(J = 1; J <= LUSOL->m; J++)
      fprintf(stdout, "%10g", denseL0[(LUSOL->n+1)*(J-1) + I]);
    fprintf(stdout, "\n");
  }
  FREE(denseL0);
}

/* ==================================================================
   lu6U   solves   U w = v.          v  is not altered.
   ------------------------------------------------------------------
   15 Dec 2002: First version derived from lu6sol.
   15 Dec 2002: Current version.
   ================================================================== */
void LU6U(LUSOLrec *LUSOL, int *INFORM, REAL V[], REAL W[], int NZidx[])
{
  int  I, J, K, KLAST, L, L1, L2, L3, NRANK, NRANK1;
  REAL SMALL;
  register REALXP T;
#ifdef LUSOLFastSolve
  REAL *aptr;
  int  *jptr;
#endif

  NRANK = LUSOL->luparm[LUSOL_IP_RANK_U];
  SMALL = LUSOL->parmlu[LUSOL_RP_ZEROTOLERANCE];
  *INFORM = LUSOL_INFORM_LUSUCCESS;
  NRANK1 = NRANK+1;
/*      Find the first nonzero in v(1:nrank), counting backwards. */
  for(KLAST = NRANK; KLAST >= 1; KLAST--) {
    I = LUSOL->ip[KLAST];
    if(fabs(V[I])>SMALL)
      break;
  }
  L = LUSOL->n;
#ifdef LUSOLFastSolve
  for(K = KLAST+1, jptr = LUSOL->iq+K; K <= L; K++, jptr++)
    W[*jptr] = ZERO;
#else
  for(K = KLAST+1; K <= L; K++) {
    J = LUSOL->iq[K];
    W[J] = ZERO;
  }
#endif
/*      Do the back-substitution, using rows 1:klast of U. */
  for(K = KLAST; K >= 1; K--) {
    I = LUSOL->ip[K];
    T = V[I];
    L1 = LUSOL->locr[I];
    L2 = L1+1;
    L3 = (L1+LUSOL->lenr[I])-1;
/*     ***** This loop could be coded specially. */
#ifdef LUSOLFastSolve
    for(L = L2, aptr = LUSOL->a+L2, jptr = LUSOL->indr+L2;
        L <= L3; L++, aptr++, jptr++)
      T -= (*aptr) * W[*jptr];
#else
    for(L = L2; L <= L3; L++) {
      J = LUSOL->indr[L];
      T -= LUSOL->a[L]*W[J];
    }
#endif
    J = LUSOL->iq[K];
    if(fabs(T)<=SMALL)
      T = ZERO;
    else
      T /= LUSOL->a[L1];
    W[J] = T;
  }
/*      Compute residual for overdetermined systems. */
  T = ZERO;
  for(K = NRANK1; K <= LUSOL->m; K++) {
    I = LUSOL->ip[K];
    T += fabs(V[I]);
  }
/*      Exit. */
  if(T>ZERO)
    *INFORM = LUSOL_INFORM_LUSINGULAR;
  LUSOL->luparm[LUSOL_IP_INFORM]     = *INFORM;
  LUSOL->parmlu[LUSOL_RP_RESIDUAL_U] = T;
}

/* ==================================================================
   lu6Ut  solves   U'v = w.          w  is destroyed.
   ------------------------------------------------------------------
   15 Dec 2002: First version derived from lu6sol.
   15 Dec 2002: Current version.
   ================================================================== */
void LU6UT(LUSOLrec *LUSOL, int *INFORM, REAL V[], REAL W[], int NZidx[])
{
  int  I, J, K, L, L1, L2, NRANK, NRANK1;
  REAL SMALL;
  register REAL T;
#ifdef LUSOLFastSolve
  REAL *aptr;
  int  *jptr;
#endif

  NRANK = LUSOL->luparm[LUSOL_IP_RANK_U];
  SMALL = LUSOL->parmlu[LUSOL_RP_ZEROTOLERANCE];
  *INFORM = LUSOL_INFORM_LUSUCCESS;
  NRANK1 = NRANK+1;
  L = LUSOL->m;
#ifdef LUSOLFastSolve
  for(K = NRANK1, jptr = LUSOL->ip+K; K <= L; K++, jptr++)
    V[*jptr] = ZERO;
#else
  for(K = NRANK1; K <= L; K++) {
    I = LUSOL->ip[K];
    V[I] = ZERO;
  }
#endif
/*      Do the forward-substitution, skipping columns of U(transpose)
        when the associated element of w(*) is negligible. */
  for(K = 1; K <= NRANK; K++) {
    I = LUSOL->ip[K];
    J = LUSOL->iq[K];
    T = W[J];
    if(fabs(T)<=SMALL) {
      V[I] = ZERO;
      continue;
    }
    L1 = LUSOL->locr[I];
    T /= LUSOL->a[L1];
    V[I] = T;
    L2 = (L1+LUSOL->lenr[I])-1;
    L1++;
/*     ***** This loop could be coded specially. */
#ifdef LUSOLFastSolve
    for(L = L1, aptr = LUSOL->a+L1, jptr = LUSOL->indr+L1;
        L <= L2; L++, aptr++, jptr++)
      W[*jptr] -= T * (*aptr);
#else
    for(L = L1; L <= L2; L++) {
      J = LUSOL->indr[L];
      W[J] -= T*LUSOL->a[L];
    }
#endif
  }
/*      Compute residual for overdetermined systems. */
  T = ZERO;
  for(K = NRANK1; K <= LUSOL->n; K++) {
    J = LUSOL->iq[K];
    T += fabs(W[J]);
  }
/*      Exit. */
  if(T>ZERO)
    *INFORM = LUSOL_INFORM_LUSINGULAR;
  LUSOL->luparm[LUSOL_IP_INFORM]     = *INFORM;
  LUSOL->parmlu[LUSOL_RP_RESIDUAL_U] = T;
}

/* ==================================================================
   lu6sol  uses the factorization  A = L U  as follows:
   ------------------------------------------------------------------
   mode
    1    v  solves   L v = v(input).   w  is not touched.
    2    v  solves   L'v = v(input).   w  is not touched.
    3    w  solves   U w = v.          v  is not altered.
    4    v  solves   U'v = w.          w  is destroyed.
    5    w  solves   A w = v.          v  is altered as in 1.
    6    v  solves   A'v = w.          w  is destroyed.

   If mode = 3,4,5,6, v and w must not be the same arrays.
   If lu1fac has just been used to factorize a symmetric matrix A
   (which must be definite or quasi-definite), the factors A = L U
   may be regarded as A = LDL', where D = diag(U).  In such cases,

   mode
    7    v  solves   A v = L D L'v = v(input).   w  is not touched.
    8    v  solves       L |D| L'v = v(input).   w  is not touched.

   ip(*), iq(*)      hold row and column numbers in pivotal order.
   lenc(k)           is the length of the k-th column of initial L.
   lenr(i)           is the length of the i-th row of U.
   locc(*)           is not used.
   locr(i)           is the start  of the i-th row of U.

   U is assumed to be in upper-trapezoidal form (nrank by n).
   The first entry for each row is the diagonal element
   (according to the permutations  ip, iq).  It is stored at
   location locr(i) in a(*), indr(*).

   On exit, inform = 0 except as follows.
     if(mode = 3,4,5,6 and if U (and hence A) is singular,)
     inform = 1 if there is a nonzero residual in solving the system
     involving U.  parmlu(20) returns the norm of the residual.
   ------------------------------------------------------------------
     July 1987: Early version.
   09 May 1988: f77 version.
   27 Apr 2000: Abolished the dreaded "computed go to".
                But hard to change other "go to"s to "if then else".
   15 Dec 2002: lu6L, lu6Lt, lu6U, lu6Ut added to modularize lu6sol.
   ================================================================== */
void LU6SOL(LUSOLrec *LUSOL, int MODE, REAL V[], REAL W[], int NZidx[], int *INFORM)
{
  if(MODE==LUSOL_SOLVE_Lv_v) {          /*      Solve  L v(new) = v. */
    LU6L(LUSOL, INFORM,V, NZidx);
  }
  else if(MODE==LUSOL_SOLVE_Ltv_v) {    /*      Solve  L'v(new) = v. */
    LU6LT(LUSOL, INFORM,V, NZidx);
  }
  else if(MODE==LUSOL_SOLVE_Uw_v) {     /*      Solve  U w = v. */
    LU6U(LUSOL, INFORM,V,W, NZidx);
  }
  else if(MODE==LUSOL_SOLVE_Utv_w) {    /*      Solve  U'v = w. */
    LU6UT(LUSOL, INFORM,V,W, NZidx);
  }
  else if(MODE==LUSOL_SOLVE_Aw_v) {     /*      Solve  A w      = v */
    LU6L(LUSOL, INFORM,V, NZidx);       /*      via     L v(new) = v */
    LU6U(LUSOL, INFORM,V,W, NULL);      /*      ... and U w = v(new). */
  }
  else if(MODE==LUSOL_SOLVE_Atv_w) {    /*      Solve  A'v = w */
    LU6UT(LUSOL, INFORM,V,W, NZidx);    /*      via      U'v = w */
    LU6LT(LUSOL, INFORM,V, NULL);       /*      ... and  L'v(new) = v. */
  }
  else if(MODE==LUSOL_SOLVE_Av_v) {     /*      Solve  LDv(bar) = v */
    LU6LD(LUSOL, INFORM,1,V, NZidx);    /*      and    L'v(new) = v(bar). */
    LU6LT(LUSOL, INFORM,V, NULL);
  }
  else if(MODE==LUSOL_SOLVE_LDLtv_v) {  /*      Solve  L|D|v(bar) = v */
    LU6LD(LUSOL, INFORM,2,V, NZidx);    /*      and    L'v(new) = v(bar). */
    LU6LT(LUSOL, INFORM,V, NULL);
  }
}



//#include "lusol1.c"      /* Factorization and core components */

/* ==================================================================
   lu1DCP factors a dense m x n matrix A by Gaussian elimination,
   using Complete Pivoting (row and column interchanges) for stability.
   This version also uses column interchanges if all elements in a
   pivot column are smaller than (or equal to) "small".  Such columns
   are changed to zero and permuted to the right-hand end.
   As in LINPACK's dgefa, ipvt(!) keeps track of pivot rows.
   Rows of U are interchanged, but we don't have to physically
   permute rows of L.  In contrast, column interchanges are applied
   directly to the columns of both L and U, and to the column
   permutation vector iq(*).
   ------------------------------------------------------------------
   On entry:
      a       Array holding the matrix A to be factored.
      lda     The leading dimension of the array  a.
      m       The number of rows    in  A.
      n       The number of columns in  A.
      small   A drop tolerance.  Must be zero or positive.

   On exit:
      a       An upper triangular matrix and the multipliers
              which were used to obtain it.
              The factorization can be written  A = L*U  where
              L  is a product of permutation and unit lower
              triangular matrices and  U  is upper triangular.
      nsing   Number of singularities detected.
      ipvt    Records the pivot rows.
      iq      A vector to which column interchanges are applied.
   ------------------------------------------------------------------
   01 May 2002: First dense Complete Pivoting, derived from lu1DPP.
   07 May 2002: Another break needed at end of first loop.
   07 May 2002: Current version of lu1DCP.
   ================================================================== */
void LU1DCP(LUSOLrec *LUSOL, REAL DA[], int LDA, int M, int N, REAL SMALL,
            int *NSING, int IPVT[], int IX[])
{

  int       I, J, K, KP1, L, LAST, LENCOL, IMAX, JMAX, JLAST, JNEW;
  REAL      AIJMAX, AJMAX;
  register REAL T;
  register int IDA1, IDA2;

  *NSING = 0;
  LENCOL = M+1;
  LAST = N;
/*     -----------------------------------------------------------------
        Start of elimination loop.
       ----------------------------------------------------------------- */
  for(K = 1; K <= N; K++) {
    KP1 = K+1;
    LENCOL--;
/*      Find the biggest aij in row imax and column jmax. */
    AIJMAX = ZERO;
    IMAX = K;
    JMAX = K;
    JLAST = LAST;
    for(J = K; J <= JLAST; J++) {
x10:
      L = idamax(LENCOL,DA+DAPOS(K,J)-LUSOL_ARRAYOFFSET,1)+K-1;
      AJMAX = fabs(DA[DAPOS(L,J)]);
      if(AJMAX<=SMALL) {
/*     ========================================================
        Do column interchange, changing old column to zero.
        Reduce  "last"  and try again with same j.
       ======================================================== */
        (*NSING)++;
        JNEW = IX[LAST];
        IX[LAST] = IX[J];
        IX[J] = JNEW;
        for(I = 1; I <= K-1; I++) {
          IDA1 = DAPOS(I,LAST);
          IDA2 = DAPOS(I,J);
          T = DA[IDA1];
          DA[IDA1] = DA[IDA2];
          DA[IDA2] = T;
        }
        for(I = K; I <= M; I++) {
          IDA1 = DAPOS(I,LAST);
          T = DA[IDA1];
          DA[IDA1] = ZERO;
          DA[DAPOS(I,J)] = T;
        }
        LAST--;
        if(J<=LAST)
          goto x10;
        break;
      }
/*      Check if this column has biggest aij so far. */
      if(AIJMAX<AJMAX) {
        AIJMAX = AJMAX;
        IMAX = L;
        JMAX = J;
      }
      if(J>=LAST)
        break;
    }
    IPVT[K] = IMAX;
    if(JMAX!=K) {
/*     ==========================================================
        Do column interchange (k and jmax).
       ========================================================== */
      JNEW = IX[JMAX];
      IX[JMAX] = IX[K];
      IX[K] = JNEW;
      for(I = 1; I <= M; I++) {
        IDA1 = DAPOS(I,JMAX);
        IDA2 = DAPOS(I,K);
        T = DA[IDA1];
        DA[IDA1] = DA[IDA2];
        DA[IDA2] = T;
      }
    }
    if(M>K) {
/*     ===========================================================
        Do row interchange if necessary.
       =========================================================== */
      if(IMAX!=K) {
        IDA1 = DAPOS(IMAX,K);
        IDA2 = DAPOS(K,K);
        T = DA[IDA1];
        DA[IDA1] = DA[IDA2];
        DA[IDA2] = T;
      }
/*     ===========================================================
        Compute multipliers.
        Do row elimination with column indexing.
       =========================================================== */
      T = -ONE/DA[DAPOS(K,K)];
      dscal(M-K,T,DA+DAPOS(KP1,K)-LUSOL_ARRAYOFFSET,1);
      for(J = KP1; J <= LAST; J++) {
        IDA1 = DAPOS(IMAX,J);
        T = DA[IDA1];
        if(IMAX!=K) {
          IDA2 = DAPOS(K,J);
          DA[IDA1] = DA[IDA2];
          DA[IDA2] = T;
        }
        daxpy(M-K,T,DA+DAPOS(KP1,K)-LUSOL_ARRAYOFFSET,1,
                    DA+DAPOS(KP1,J)-LUSOL_ARRAYOFFSET,1);
      }
    }
    else
      break;
    if(K>=LAST)
      break;
  }
/*      Set ipvt(*) for singular rows. */
  for(K = LAST+1; K <= M; K++)
    IPVT[K] = K;

}

/* ==================================================================
   lu1DPP factors a dense m x n matrix A by Gaussian elimination,
   using row interchanges for stability, as in dgefa from LINPACK.
   This version also uses column interchanges if all elements in a
   pivot column are smaller than (or equal to) "small".  Such columns
   are changed to zero and permuted to the right-hand end.
   As in LINPACK, ipvt(*) keeps track of pivot rows.
   Rows of U are interchanged, but we don't have to physically
   permute rows of L.  In contrast, column interchanges are applied
   directly to the columns of both L and U, and to the column
   permutation vector iq(*).
   ------------------------------------------------------------------
   On entry:
        a       Array holding the matrix A to be factored.
        lda     The leading dimension of the array  a.
        m       The number of rows    in  A.
        n       The number of columns in  A.
        small   A drop tolerance.  Must be zero or positive.

   On exit:
        a       An upper triangular matrix and the multipliers
                which were used to obtain it.
                The factorization can be written  A = L*U  where
                L  is a product of permutation and unit lower
                triangular matrices and  U  is upper triangular.
        nsing   Number of singularities detected.
        ipvt    Records the pivot rows.
        iq      A vector to which column interchanges are applied.
   ------------------------------------------------------------------
   02 May 1989: First version derived from dgefa
                in LINPACK (version dated 08/14/78).
   05 Feb 1994: Generalized to treat rectangular matrices
                and use column interchanges when necessary.
                ipvt is retained, but column permutations are applied
                directly to iq(*).
   21 Dec 1994: Bug found via example from Steve Dirkse.
                Loop 100 added to set ipvt(*) for singular rows.
   ================================================================== */
void LU1DPP(LUSOLrec *LUSOL, REAL DA[], int LDA, int M, int N, REAL SMALL,
            int *NSING, int IPVT[], int IX[])
{
  int            I, J, K, KP1, L, LAST, LENCOL;
  register REAL T;
  register int  IDA1, IDA2;

  *NSING = 0;
  K = 1;
  LAST = N;
/*      ------------------------------------------------------------------
        Start of elimination loop.
        ------------------------------------------------------------------ */
x10:
  KP1 = K+1;
  LENCOL = (M-K)+1;
/*      Find l, the pivot row. */
  L = (idamax(LENCOL,DA+DAPOS(K,K)-LUSOL_ARRAYOFFSET,1)+K)-1;
  IPVT[K] = L;
  if(fabs(DA[DAPOS(L,K)])<=SMALL) {
/*         ===============================================================
           Do column interchange, changing old pivot column to zero.
           Reduce  "last"  and try again with same k.
           =============================================================== */
    (*NSING)++;
    J = IX[LAST];
    IX[LAST] = IX[K];
    IX[K] = J;
    for(I = 1; I <= K-1; I++) {
      IDA1 = DAPOS(I,LAST);
      IDA2 = DAPOS(I,K);
      T = DA[IDA1];
      DA[IDA1] = DA[IDA2];
      DA[IDA2] = T;
    }
    for(I = K; I <= M; I++) {
      IDA1 = DAPOS(I,LAST);
      T = DA[IDA1];
      DA[IDA1] = ZERO;
      DA[DAPOS(I,K)] = T;
    }
    LAST = LAST-1;
    if(K<=LAST)
      goto x10;
  }
  else if(M>K) {
/*         ===============================================================
           Do row interchange if necessary.
           =============================================================== */
    if(L!=K) {
      IDA1 = DAPOS(L,K);
      IDA2 = DAPOS(K,K);
      T = DA[IDA1];
      DA[IDA1] = DA[IDA2];
      DA[IDA2] = T;
    }
/*         ===============================================================
           Compute multipliers.
           Do row elimination with column indexing.
           =============================================================== */
    T = -ONE/DA[DAPOS(K,K)];
    dscal(M-K,T,DA+DAPOS(KP1,K)-LUSOL_ARRAYOFFSET,1);
    for(J = KP1; J <= LAST; J++) {
      IDA1 = DAPOS(L,J);
      T = DA[IDA1];
      if(L!=K) {
        IDA2 = DAPOS(K,J);
        DA[IDA1] = DA[IDA2];
        DA[IDA2] = T;
      }
      daxpy(M-K,T,DA+DAPOS(KP1,K)-LUSOL_ARRAYOFFSET,1,
                  DA+DAPOS(KP1,J)-LUSOL_ARRAYOFFSET,1);
    }
    K++;
    if(K<=LAST)
      goto x10;
  }
/*      Set ipvt(*) for singular rows. */
  for(K = LAST+1; K <= M; K++)
    IPVT[K] = K;

}


/* ==================================================================
   lu1pq1  constructs a permutation  iperm  from the array  len.
   ------------------------------------------------------------------
   On entry:
   len(i)  holds the number of nonzeros in the i-th row (say)
           of an m by n matrix.
   num(*)  can be anything (workspace).

   On exit:
   iperm   contains a list of row numbers in the order
           rows of length 0,  rows of length 1,..., rows of length n.
   loc(nz) points to the first row containing  nz  nonzeros,
           nz = 1, n.
   inv(i)  points to the position of row i within iperm(*).
   ================================================================== */
void LU1PQ1(LUSOLrec *LUSOL, int M, int N, int LEN[],
            int IPERM[], int LOC[], int INV[], int NUM[])
{
  int NZERO, NZ, I, L;

/*      Count the number of rows of each length. */
  NZERO = 0;
  for(NZ = 1; NZ <= N; NZ++) {
    NUM[NZ] = 0;
    LOC[NZ] = 0;
  }
  for(I = 1; I <= M; I++) {
    NZ = LEN[I];
    if(NZ==0)
      NZERO++;
    else
      NUM[NZ]++;
  }
/*      Set starting locations for each length. */
  L = NZERO+1;
  for(NZ = 1; NZ <= N; NZ++) {
    LOC[NZ] = L;
    L += NUM[NZ];
    NUM[NZ] = 0;
  }
/*      Form the list. */
  NZERO = 0;
  for(I = 1; I <= M; I++) {
    NZ = LEN[I];
    if(NZ==0) {
      NZERO++;
      IPERM[NZERO] = I;
    }
    else {
      L = LOC[NZ]+NUM[NZ];
      IPERM[L] = I;
      NUM[NZ]++;
    }
  }
/*      Define the inverse of iperm. */
  for(L = 1; L <= M; L++) {
    I = IPERM[L];
    INV[I] = L;
  }
}

/* ==================================================================
   lu1pq2 frees the space occupied by the pivot row,
   and updates the column permutation iq.
   Also used to free the pivot column and update the row perm ip.
   ------------------------------------------------------------------
   nzpiv   (input)    is the length of the pivot row (or column).
   nzchng  (output)   is the net change in total nonzeros.
   ------------------------------------------------------------------
   14 Apr 1989  First version.
   ================================================================== */
void LU1PQ2(LUSOLrec *LUSOL, int NZPIV, int *NZCHNG,
            int IND[], int LENOLD[], int LENNEW[], int IXLOC[], int IX[], int IXINV[])
{
  int LR, J, NZ, NZNEW, L, NEXT, LNEW, JNEW;

  *NZCHNG = 0;
  for(LR = 1; LR <= NZPIV; LR++) {
    J = IND[LR];
    IND[LR] = 0;
    NZ = LENOLD[LR];
    NZNEW = LENNEW[J];
    if(NZ!=NZNEW) {
      L = IXINV[J];
      *NZCHNG = (*NZCHNG+NZNEW)-NZ;
/*            l above is the position of column j in iq  (so j = iq(l)). */
      if(NZ<NZNEW) {
/*               Column  j  has to move towards the end of  iq. */
x110:
        NEXT = NZ+1;
        LNEW = IXLOC[NEXT]-1;
        if(LNEW!=L) {
          JNEW = IX[LNEW];
          IX[L] = JNEW;
          IXINV[JNEW] = L;
        }
        L = LNEW;
        IXLOC[NEXT] = LNEW;
        NZ = NEXT;
        if(NZ<NZNEW)
          goto x110;
      }
      else {
/*               Column  j  has to move towards the front of  iq. */
x120:
        LNEW = IXLOC[NZ];
        if(LNEW!=L) {
          JNEW = IX[LNEW];
          IX[L] = JNEW;
          IXINV[JNEW] = L;
        }
        L = LNEW;
        IXLOC[NZ] = LNEW+1;
        NZ = NZ-1;
        if(NZ>NZNEW)
          goto x120;
      }
      IX[LNEW] = J;
      IXINV[J] = LNEW;
    }
  }
}

/* ==================================================================
   lu1pq3  looks at the permutation  iperm(*)  and moves any entries
   to the end whose corresponding length  len(*)  is zero.
   ------------------------------------------------------------------
   09 Feb 1994: Added work array iw(*) to improve efficiency.
   ================================================================== */
void LU1PQ3(LUSOLrec *LUSOL, int MN, int LEN[], int IPERM[], int IW[], int *NRANK)
{
  int NZERO, K, I;

  *NRANK = 0;
  NZERO = 0;
  for(K = 1; K <= MN; K++) {
    I = IPERM[K];
    if(LEN[I]==0) {
      NZERO++;
      IW[NZERO] = I;
    }
    else {
      (*NRANK)++;
      IPERM[*NRANK] = I;
    }
  }
  for(K = 1; K <= NZERO; K++)
    IPERM[(*NRANK)+K] = IW[K];
}

/* ==================================================================
   lu1rec
   ------------------------------------------------------------------
   On exit:
   ltop         is the length of useful entries in ind(*), a(*).
   ind(ltop+1)  is "i" such that len(i), loc(i) belong to the last
                item in ind(*), a(*).
   ------------------------------------------------------------------
   00 Jun 1983: Original version of lu1rec followed John Reid's
                compression routine in LA05.  It recovered
                space in ind(*) and optionally a(*)
                by eliminating entries with ind(l) = 0.
                The elements of ind(*) could not be negative.
                If len(i) was positive, entry i contained
                that many elements, starting at  loc(i).
                Otherwise, entry i was eliminated.
   23 Mar 2001: Realised we could have len(i) = 0 in rare cases!
                (Mostly during TCP when the pivot row contains
                a column of length 1 that couldn't be a pivot.)
                Revised storage scheme to
                   keep        entries with       ind(l) >  0,
                   squeeze out entries with -n <= ind(l) <= 0,
                and to allow len(i) = 0.
                Empty items are moved to the end of the compressed
                ind(*) and/or a(*) arrays are given one empty space.
                Items with len(i) < 0 are still eliminated.
   27 Mar 2001: Decided to use only ind(l) > 0 and = 0 in lu1fad.
                Still have to keep entries with len(i) = 0.
   ================================================================== */
void LU1REC(LUSOLrec *LUSOL, int N, MYBOOL REALS, int *LTOP,
                             int IND[], int LEN[], int LOC[])
{
  int  NEMPTY, I, LENI, L, K, KLAST, ILAST, LPRINT;

  NEMPTY = 0;
  for(I = 1; I <= N; I++) {
    LENI = LEN[I];
    if(LENI>0) {
      L = (LOC[I]+LENI)-1;
      LEN[I] = IND[L];
      IND[L] = -(N+I);
    }
    else if(LENI==0)
      NEMPTY++;
  }
  K = 0;
/*      Previous k */
  KLAST = 0;
/*      Last entry moved. */
  ILAST = 0;
  for(L = 1; L <= *LTOP; L++) {
    I = IND[L];
    if(I>0) {
      K++;
      IND[K] = I;
      if(REALS)
        LUSOL->a[K] = LUSOL->a[L];
    }
    else if(I<-N) {
/*            This is the end of entry  i. */
      I = -(N+I);
      ILAST = I;
      K++;
      IND[K] = LEN[I];
      if(REALS)
        LUSOL->a[K] = LUSOL->a[L];
      LOC[I] = KLAST+1;
      LEN[I] = K-KLAST;
      KLAST = K;
    }
  }
/*      Move any empty items to the end, adding 1 free entry for each. */
  if(NEMPTY>0) {
    for(I = 1; I <= N; I++) {
      if(LEN[I]==0) {
        K++;
        LOC[I] = K;
        IND[K] = 0;
        ILAST = I;
      }
    }
  }
  LPRINT = LUSOL->luparm[LUSOL_IP_PRINTLEVEL];
  if(LPRINT>=LUSOL_MSG_PIVOT)
    LUSOL_report(LUSOL, 0, "lu1rec.  File compressed from %d to %d\n",
                        *LTOP,K,REALS,NEMPTY);
/*      ncp */
  LUSOL->luparm[LUSOL_IP_COMPRESSIONS_LU]++;
/*      Return ilast in ind(ltop + 1). */
  *LTOP = K;
  IND[(*LTOP)+1] = ILAST;
}

/* ==================================================================
   lu1slk  sets w(j) > 0 if column j is a unit vector.
   ------------------------------------------------------------------
   21 Nov 2000: First version.  lu1fad needs it for TCP.
                Note that w(*) is nominally an integer array,
                but the only spare space is the double array w(*).
   ================================================================== */
void LU1SLK(LUSOLrec *LUSOL)
{
  int J, LC1, LQ, LQ1, LQ2;

  for(J = 1; J <= LUSOL->n; J++) {
    LUSOL->w[J] = 0;
  }
  LQ1 = LUSOL->iqloc[1];
  LQ2 = LUSOL->n;
  if(LUSOL->m>1)
    LQ2 = LUSOL->iqloc[2]-1;
  for(LQ = LQ1; LQ <= LQ2; LQ++) {
    J = LUSOL->iq[LQ];
    LC1 = LUSOL->locc[J];
    if(fabs(LUSOL->a[LC1])==1) {
      LUSOL->w[J] = 1;
    }
  }
}

/* ==================================================================
   lu1gau does most of the work for each step of Gaussian elimination.
   A multiple of the pivot column is added to each other column j
   in the pivot row.  The column list is fully updated.
   The row list is updated if there is room, but some fill-ins may
   remain, as indicated by ifill and jfill.
   ------------------------------------------------------------------
   Input:
      ilast    is the row    at the end of the row    list.
      jlast    is the column at the end of the column list.
      lfirst   is the first column to be processed.
      lu + 1   is the corresponding element of U in au(*).
      nfill    keeps track of pending fill-in.
      a(*)     contains the nonzeros for each column j.
      indc(*)  contains the row indices for each column j.
      al(*)    contains the new column of L.  A multiple of it is
               used to modify each column.
      mark(*)  has been set to -1, -2, -3, ... in the rows
               corresponding to nonzero 1, 2, 3, ... of the col of L.
      au(*)    contains the new row of U.  Each nonzero gives the
               required multiple of the column of L.

   Workspace:
      markl(*) marks the nonzeros of L actually used.
               (A different mark, namely j, is used for each column.)

   Output:
      ilast     New last row    in the row    list.
      jlast     New last column in the column list.
      lfirst    = 0 if all columns were completed,
                > 0 otherwise.
      lu        returns the position of the last nonzero of U
                actually used, in case we come back in again.
      nfill     keeps track of the total extra space needed in the
                row file.
      ifill(ll) counts pending fill-in for rows involved in the new
                column of L.
      jfill(lu) marks the first pending fill-in stored in columns
                involved in the new row of U.
   ------------------------------------------------------------------
   16 Apr 1989: First version of lu1gau.
   23 Apr 1989: lfirst, lu, nfill are now input and output
                to allow re-entry if elimination is interrupted.
   23 Mar 2001: Introduced ilast, jlast.
   27 Mar 2001: Allow fill-in "in situ" if there is already room
                up to but NOT INCLUDING the end of the
                row or column file.
                Seems safe way to avoid overwriting empty rows/cols
                at the end.  (May not be needed though, now that we
                have ilast and jlast.)
   ================================================================== */
void LU1GAU(LUSOLrec *LUSOL, int MELIM, int NSPARE,
            REAL SMALL, int LPIVC1, int LPIVC2, int *LFIRST, int LPIVR2,
            int LFREE, int MINFRE, int ILAST, int *JLAST, int *LROW, int *LCOL,
            int *LU, int *NFILL,
            int MARK[],  REAL AL[], int MARKL[], REAL AU[], int IFILL[], int JFILL[])
{
  MYBOOL ATEND;
  int    LR, J, LENJ, NFREE, LC1, LC2, NDONE, NDROP, L, I, LL, K,
         LR1, LAST, LREP, L1, L2, LC, LENI;
  register REAL UJ;
  REAL   AIJ;

  for(LR = *LFIRST; LR <= LPIVR2; LR++) {
    J = LUSOL->indr[LR];
    LENJ = LUSOL->lenc[J];
    NFREE = LFREE - *LCOL;
    if(NFREE<MINFRE)
      goto x900;
/*         ---------------------------------------------------------------
           Inner loop to modify existing nonzeros in column  j.
           Loop 440 performs most of the arithmetic involved in the
           whole LU factorization.
           ndone counts how many multipliers were used.
           ndrop counts how many modified nonzeros are negligibly small.
           --------------------------------------------------------------- */
    (*LU)++;
    UJ = AU[*LU];
    LC1 = LUSOL->locc[J];
    LC2 = (LC1+LENJ)-1;
    ATEND = (MYBOOL) (J==*JLAST);
    NDONE = 0;
    if(LENJ==0)
      goto x500;
    NDROP = 0;
    for(L = LC1; L <= LC2; L++) {
      I = LUSOL->indc[L];
      LL = -MARK[I];
      if(LL>0) {
        NDONE++;
        MARKL[LL] = J;
        LUSOL->a[L] += AL[LL]*UJ;
        if(fabs(LUSOL->a[L])<=SMALL) {
          NDROP++;
        }
      }
    }
/*         ---------------------------------------------------------------
           Remove any negligible modified nonzeros from both
           the column file and the row file.
           --------------------------------------------------------------- */
    if(NDROP==0)
      goto x500;
    K = LC1;
    for(L = LC1; L <= LC2; L++) {
      I = LUSOL->indc[L];
      if(fabs(LUSOL->a[L])<=SMALL)
        goto x460;
      LUSOL->a[K] = LUSOL->a[L];
      LUSOL->indc[K] = I;
      K++;
      continue;
/*            Delete the nonzero from the row file. */
x460:
      LENJ--;
      LUSOL->lenr[I]--;
      LR1 = LUSOL->locr[I];
      LAST = LR1+LUSOL->lenr[I];
      for(LREP = LR1; LREP <= LAST; LREP++) {
        if(LUSOL->indr[LREP]==J)
          break;
      }
      LUSOL->indr[LREP] = LUSOL->indr[LAST];
      LUSOL->indr[LAST] = 0;
      if(I==ILAST)
        (*LROW)--;
    }
/*         Free the deleted elements from the column file. */
#ifdef LUSOLFastClear
    MEMCLEAR(LUSOL->indc+K, LC2-K+1);
#else
    for(L = K; L <= LC2; L++)
      LUSOL->indc[L] = ZERO;
#endif
    if(ATEND)
      *LCOL = K-1;
/*         ---------------------------------------------------------------
           Deal with the fill-in in column j.
           --------------------------------------------------------------- */
x500:
    if(NDONE==MELIM)
      goto x590;
/*         See if column j already has room for the fill-in. */
    if(ATEND)
      goto x540;
    LAST = (LC1+LENJ)-1;
    L1 = LAST+1;
    L2 = (LAST+MELIM)-NDONE;
/*      27 Mar 2001: Be sure it's not at or past end of the col file. */
    if(L2>=*LCOL)
      goto x520;
    for(L = L1; L <= L2; L++) {
      if(LUSOL->indc[L]!=0)
        goto x520;
    }
    goto x540;
/*         We must move column j to the end of the column file.
           First, leave some spare room at the end of the
           current last column. */
x520:
#if 1
    L1 = (*LCOL)+1;
    L2 = (*LCOL)+NSPARE;
    *LCOL = L2;
    for(L = L1; L <= L2; L++) {
#else
    for(L = (*LCOL)+1; L <= (*LCOL)+NSPARE; L++) {
      *LCOL = L;  /* ****** ERROR ???? */
#endif
/*      Spare space is free. */
      LUSOL->indc[L] = 0;
    }
    ATEND = TRUE;
    *JLAST = J;
    L1 = LC1;
    LC1 = (*LCOL)+1;
    LUSOL->locc[J] = LC1;
    for(L = L1; L <= LAST; L++) {
      (*LCOL)++;
      LUSOL->a[*LCOL] = LUSOL->a[L];
      LUSOL->indc[*LCOL] = LUSOL->indc[L];
/*      Free space. */
      LUSOL->indc[L] = 0;
    }
/*         ---------------------------------------------------------------
           Inner loop for the fill-in in column j.
           This is usually not very expensive.
           --------------------------------------------------------------- */
x540:
    LAST = (LC1+LENJ)-1;
    LL = 0;
    for(LC = LPIVC1; LC <= LPIVC2; LC++) {
      LL++;
      if(MARKL[LL]==J)
        continue;
      AIJ = AL[LL]*UJ;
      if(fabs(AIJ)<=SMALL)
        continue;
      LENJ++;
      LAST++;
      LUSOL->a[LAST] = AIJ;
      I = LUSOL->indc[LC];
      LUSOL->indc[LAST] = I;
      LENI = LUSOL->lenr[I];
/*            Add 1 fill-in to row i if there is already room.
              27 Mar 2001: Be sure it's not at or past the }
                           of the row file. */
      L = LUSOL->locr[I]+LENI;
      if(L>=*LROW)
        goto x550;
      if(LUSOL->indr[L]>0)
        goto x550;
      LUSOL->indr[L] = J;
      LUSOL->lenr[I] = LENI+1;
      continue;
/*            Row i does not have room for the fill-in.
              Increment ifill(ll) to count how often this has
              happened to row i.  Also, add m to the row index
              indc(last) in column j to mark it as a fill-in that is
              still pending.
              If this is the first pending fill-in for row i,
              nfill includes the current length of row i
              (since the whole row has to be moved later).
              If this is the first pending fill-in for column j,
              jfill(lu) records the current length of column j
              (to shorten the search for pending fill-ins later). */
x550:
      if(IFILL[LL]==0)
        (*NFILL) += LENI+NSPARE;
      if(JFILL[*LU]==0)
        JFILL[*LU] = LENJ;
      (*NFILL)++;
      IFILL[LL]++;
      LUSOL->indc[LAST] = LUSOL->m+I;
    }
    if(ATEND)
      *LCOL = LAST;
/*         End loop for column  j.  Store its final length. */
x590:
    LUSOL->lenc[J] = LENJ;
  }
/*      Successful completion. */
  *LFIRST = 0;
  return;
/*      Interruption.  We have to come back in after the
        column file is compressed.  Give lfirst a new value.
        lu and nfill will retain their current values. */
x900:
  *LFIRST = LR;
}

/* ==================================================================
   lu1mar  uses a Markowitz criterion to select a pivot element
   for the next stage of a sparse LU factorization,
   subject to a Threshold Partial Pivoting stability criterion (TPP)
   that bounds the elements of L.
   ------------------------------------------------------------------
   gamma  is "gamma" in the tie-breaking rule TB4 in the LUSOL paper.
   ------------------------------------------------------------------
   Search cols of length nz = 1, then rows of length nz = 1,
   then   cols of length nz = 2, then rows of length nz = 2, etc.
   ------------------------------------------------------------------
   00 Jan 1986  Version documented in LUSOL paper:
                Gill, Murray, Saunders and Wright (1987),
                Maintaining LU factors of a general sparse matrix,
                Linear algebra and its applications 88/89, 239-270.
   02 Feb 1989  Following Suhl and Aittoniemi (1987), the largest
                element in each column is now kept at the start of
                the column, i.e. in position locc(j) of a and indc.
                This should speed up the Markowitz searches.
   26 Apr 1989  Both columns and rows searched during spars1 phase.
                Only columns searched during spars2 phase.
                maxtie replaced by maxcol and maxrow.
   05 Nov 1993  Initializing  "mbest = m * n"  wasn't big enough when
                m = 10, n = 3, and last column had 7 nonzeros.
   09 Feb 1994  Realised that "mbest = maxmn * maxmn" might overflow.
                Changed to    "mbest = maxmn * 1000".
   27 Apr 2000  On large example from Todd Munson,
                that allowed  "if (mbest .le. nz1**2) go to 900"
                to exit before any pivot had been found.
                Introduced kbest = mbest / nz1.
                Most pivots can be rejected with no integer multiply.
                TRUE merit is evaluated only if it's as good as the
                best so far (or better).  There should be no danger
                of integer overflow unless A is incredibly
                large and dense.
   10 Sep 2000  TCP, aijtol added for Threshold Complete Pivoting.
   ================================================================== */
void LU1MAR(LUSOLrec *LUSOL, int MAXMN, MYBOOL TCP, REAL AIJTOL, REAL LTOL,
            int MAXCOL, int MAXROW, int *IBEST, int *JBEST, int *MBEST)
{
  int  KBEST, NCOL, NROW, NZ1, NZ, LQ1, LQ2, LQ, J, LC1, LC2, LC, I, LEN1, MERIT, LP1,
       LP2, LP, LR1, LR2, LR;
  REAL ABEST, LBEST, AMAX, AIJ, CMAX;

  ABEST = ZERO;
  LBEST = ZERO;
  *IBEST = 0;
  *MBEST = -1;
  KBEST = MAXMN+1;
  NCOL = 0;
  NROW = 0;
  NZ1 = 0;
  for(NZ = 1; NZ <= MAXMN; NZ++) {
/*         nz1    = nz - 1
           if (mbest .le. nz1**2) go to 900 */
    if(KBEST<=NZ1)
      goto x900;
    if(*IBEST>0) {
      if(NCOL>=MAXCOL)
        goto x200;
    }
    if(NZ>LUSOL->m)
      goto x200;
/*         ---------------------------------------------------------------
           Search the set of columns of length  nz.
           --------------------------------------------------------------- */
    LQ1 = LUSOL->iqloc[NZ];
    LQ2 = LUSOL->n;
    if(NZ<LUSOL->m)
      LQ2 = LUSOL->iqloc[NZ+1]-1;
    for(LQ = LQ1; LQ <= LQ2; LQ++) {
      NCOL = NCOL+1;
      J = LUSOL->iq[LQ];
      LC1 = LUSOL->locc[J];
      LC2 = LC1+NZ1;
      AMAX = fabs(LUSOL->a[LC1]);
/*            Test all aijs in this column.
              amax is the largest element (the first in the column).
              cmax is the largest multiplier if aij becomes pivot. */
      if(TCP) {
/*      Nothing in whole column */
        if(AMAX<AIJTOL)
          continue;
      }
      for(LC = LC1; LC <= LC2; LC++) {
        I = LUSOL->indc[LC];
        LEN1 = LUSOL->lenr[I]-1;
/*               merit  = nz1 * len1
                 if (merit > mbest) continue; */
        if(LEN1>KBEST)
          continue;
/*               aij  has a promising merit.
                 Apply the stability test.
                 We require  aij  to be sufficiently large compared to
                 all other nonzeros in column  j.  This is equivalent
                 to requiring cmax to be bounded by Ltol. */
        if(LC==LC1) {
/*                  This is the maximum element, amax.
                    Find the biggest element in the rest of the column
                    and hence get cmax.  We know cmax .le. 1, but
                    we still want it exactly in order to break ties.
                    27 Apr 2002: Settle for cmax = 1. */
          AIJ = AMAX;
          CMAX = ONE;
/*                  cmax   = zero
                    for (l = lc1 + 1; l <= lc2; l++)
                       cmax  = max( cmax, abs( a(l) ) );
                    cmax   = cmax / amax; */
        }
        else {
/*                  aij is not the biggest element, so cmax .ge. 1.
                    Bail out if cmax will be too big. */
          AIJ = fabs(LUSOL->a[LC]);
/*      Absolute test for Complete Pivoting */
          if(TCP) {
            if(AIJ<AIJTOL)
              continue;
/*      TPP */
          }
          else {
            if(AIJ*LTOL<AMAX)
              continue;
          }
          CMAX = AMAX/AIJ;
        }
/*               aij  is big enough.  Its maximum multiplier is cmax. */
        MERIT = NZ1*LEN1;
        if(MERIT==*MBEST) {
/*                  Break ties.
                    (Initializing mbest < 0 prevents getting here if
                    nothing has been found yet.)
                    In this version we minimize cmax
                    but if it is already small we maximize the pivot. */
          if(LBEST<=LUSOL->parmlu[LUSOL_RP_GAMMA] &&
             CMAX<=LUSOL->parmlu[LUSOL_RP_GAMMA]) {
            if(ABEST>=AIJ)
              continue;
          }
          else {
            if(LBEST<=CMAX)
              continue;
          }
        }
/*               aij  is the best pivot so far. */
        *IBEST = I;
        *JBEST = J;
        KBEST = LEN1;
        *MBEST = MERIT;
        ABEST = AIJ;
        LBEST = CMAX;
        if(NZ==1)
          goto x900;
      }
/*            Finished with that column. */
      if(*IBEST>0) {
        if(NCOL>=MAXCOL)
          goto x200;
      }
    }
/*         ---------------------------------------------------------------
           Search the set of rows of length  nz.
           --------------------------------------------------------------- */
x200:
/*    if (mbest .le. nz*nz1) go to 900 */
    if(KBEST<=NZ)
      goto x900;
    if(*IBEST>0) {
      if(NROW>=MAXROW)
        goto x290;
    }
    if(NZ>LUSOL->n)
      goto x290;
    LP1 = LUSOL->iploc[NZ];
    LP2 = LUSOL->m;
    if(NZ<LUSOL->n)
      LP2 = LUSOL->iploc[NZ+1]-1;
    for(LP = LP1; LP <= LP2; LP++) {
      NROW++;
      I = LUSOL->ip[LP];
      LR1 = LUSOL->locr[I];
      LR2 = LR1+NZ1;
      for(LR = LR1; LR <= LR2; LR++) {
        J = LUSOL->indr[LR];
        LEN1 = LUSOL->lenc[J]-1;
/*               merit  = nz1 * len1
                 if (merit .gt. mbest) continue */
        if(LEN1>KBEST)
          continue;
/*               aij  has a promising merit.
                 Find where  aij  is in column  j. */
        LC1 = LUSOL->locc[J];
        LC2 = LC1+LEN1;
        AMAX = fabs(LUSOL->a[LC1]);
        for(LC = LC1; LC <= LC2; LC++) {
          if(LUSOL->indc[LC]==I)
            break;
        }
/*               Apply the same stability test as above. */
        AIJ = fabs(LUSOL->a[LC]);
/*      Absolute test for Complete Pivoting */
        if(TCP) {
          if(AIJ<AIJTOL)
            continue;
        }
        if(LC==LC1) {
/*                  This is the maximum element, amax.
                    Find the biggest element in the rest of the column
                    and hence get cmax.  We know cmax .le. 1, but
                    we still want it exactly in order to break ties.
                    27 Apr 2002: Settle for cmax = 1. */
          CMAX = ONE;
/*                  cmax   = zero
                    for(l = lc1 + 1; l <= lc2; l++)
                       cmax  = max( cmax, fabs( a(l) ) )
                    cmax   = cmax / amax */
        }
        else {
/*                  aij is not the biggest element, so cmax .ge. 1.
                    Bail out if cmax will be too big. */
          if(TCP) {
/*      relax */
          }
          else {
            if(AIJ*LTOL<AMAX)
              continue;
          }
          CMAX = AMAX/AIJ;
        }
/*               aij  is big enough.  Its maximum multiplier is cmax. */
        MERIT = NZ1*LEN1;
        if(MERIT==*MBEST) {
/*                  Break ties as before.
                    (Initializing mbest < 0 prevents getting here if
                    nothing has been found yet.) */
          if(LBEST<=LUSOL->parmlu[LUSOL_RP_GAMMA] &&
             CMAX<=LUSOL->parmlu[LUSOL_RP_GAMMA]) {
            if(ABEST>=AIJ)
              continue;
          }
          else {
            if(LBEST<=CMAX)
              continue;
          }
        }
/*               aij  is the best pivot so far. */
        *IBEST = I;
        *JBEST = J;
        *MBEST = MERIT;
        KBEST = LEN1;
        ABEST = AIJ;
        LBEST = CMAX;
        if(NZ==1)
          goto x900;
      }
/*            Finished with that row. */
      if(*IBEST>0) {
        if(NROW>=MAXROW)
          goto x290;
      }
    }
/*         See if it's time to quit. */
x290:
    if(*IBEST>0) {
      if(NROW>=MAXROW && NCOL>=MAXCOL)
        goto x900;
    }
/*         Press on with next nz. */
    NZ1 = NZ;
    if(*IBEST>0)
      KBEST = *MBEST/NZ1;
  }
x900:
;
}

/* ==================================================================
   lu1mCP  uses a Markowitz criterion to select a pivot element
   for the next stage of a sparse LU factorization,
   subject to a Threshold Complete Pivoting stability criterion (TCP)
   that bounds the elements of L and U.
   ------------------------------------------------------------------
   gamma  is "gamma" in the tie-breaking rule TB4 in the LUSOL paper.
   ------------------------------------------------------------------
   09 May 2002: First version of lu1mCP.
                It searches columns only, using the heap that
                holds the largest element in each column.
   09 May 2002: Current version of lu1mCP.
   ================================================================== */
void LU1MCP(LUSOLrec *LUSOL, REAL AIJTOL, int *IBEST, int *JBEST, int *MBEST,
            int HLEN, REAL HA[], int HJ[])
{
  int  J, KHEAP, LC, LC1, LC2, LENJ, MAXCOL, NCOL, NZ1, I, LEN1, MERIT;
  REAL ABEST, AIJ, AMAX, CMAX, LBEST;

/*      ------------------------------------------------------------------
        Search up to maxcol columns stored at the top of the heap.
        The very top column helps initialize mbest.
        ------------------------------------------------------------------ */
  ABEST = ZERO;
  LBEST = ZERO;
  *IBEST = 0;
/*      Column at the top of the heap */
  *JBEST = HJ[1];
  LENJ = LUSOL->lenc[*JBEST];
/*      Bigger than any possible merit */
  *MBEST = LENJ*HLEN;
/*      ??? Big question */
  MAXCOL = 40;
/*      No. of columns searched */
  NCOL = 0;
  for(KHEAP = 1; KHEAP <= HLEN; KHEAP++) {
    AMAX = HA[KHEAP];
    if(AMAX<AIJTOL)
      continue;
    NCOL++;
    J = HJ[KHEAP];
/*         ---------------------------------------------------------------
           This column has at least one entry big enough (the top one).
           Search the column for other possibilities.
           --------------------------------------------------------------- */
    LENJ = LUSOL->lenc[J];
    NZ1 = LENJ-1;
    LC1 = LUSOL->locc[J];
    LC2 = LC1+NZ1;
/* --      amax   = abs( a(lc1) )
           Test all aijs in this column.
           amax is the largest element (the first in the column).
           cmax is the largest multiplier if aij becomes pivot. */
    for(LC = LC1; LC <= LC2; LC++) {
      I = LUSOL->indc[LC];
      LEN1 = LUSOL->lenr[I]-1;
      MERIT = NZ1*LEN1;
      if(MERIT>*MBEST)
        continue;
/*            aij  has a promising merit. */
      if(LC==LC1) {
/*               This is the maximum element, amax.
                 Find the biggest element in the rest of the column
                 and hence get cmax.  We know cmax .le. 1, but
                 we still want it exactly in order to break ties.
                 27 Apr 2002: Settle for cmax = 1. */
        AIJ = AMAX;
        CMAX = ONE;
/*               cmax   = ZERO;
                 for(l = lc1 + 1; l <= lc2; l++)
                    cmax  = max( cmax, abs( a(l) ) )
                 cmax   = cmax / amax; */
      }
      else {
/*               aij is not the biggest element, so cmax .ge. 1.
                 Bail out if cmax will be too big. */
        AIJ = fabs(LUSOL->a[LC]);
        if(AIJ<AIJTOL)
          continue;
        CMAX = AMAX/AIJ;
      }
/*            aij  is big enough.  Its maximum multiplier is cmax. */
      if(MERIT==*MBEST) {
/*               Break ties.
                 (Initializing mbest "too big" prevents getting here if
                 nothing has been found yet.)
                 In this version we minimize cmax
                 but if it is already small we maximize the pivot. */
        if(LBEST<=LUSOL->parmlu[LUSOL_RP_GAMMA] &&
           CMAX<=LUSOL->parmlu[LUSOL_RP_GAMMA]) {
          if(ABEST>=AIJ)
            continue;
        }
        else {
          if(LBEST<=CMAX)
            continue;
        }
      }
/*            aij  is the best pivot so far. */
      *IBEST = I;
      *JBEST = J;
      *MBEST = MERIT;
      ABEST = AIJ;
      LBEST = CMAX;
/*      Col or row of length 1 */
      if(MERIT==0)
        goto x900;
    }
    if(NCOL>=MAXCOL)
      goto x900;
  }
x900:
;
}

/* ==================================================================
   lu1mRP  uses a Markowitz criterion to select a pivot element
   for the next stage of a sparse LU factorization,
   subject to a Threshold Rook Pivoting stability criterion (TRP)
   that bounds the elements of L and U.
   ------------------------------------------------------------------
   11 Jun 2002: First version of lu1mRP derived from lu1mar.
   11 Jun 2002: Current version of lu1mRP.
   ================================================================== */
void LU1MRP(LUSOLrec *LUSOL, int MAXMN, REAL LTOL, int MAXCOL, int MAXROW,
  int *IBEST, int *JBEST, int *MBEST, REAL AMAXR[])
{
  int  I, J, KBEST, LC, LC1, LC2, LEN1, LP, LP1, LP2, LQ, LQ1,
       LQ2, LR, LR1, LR2, MERIT, NCOL, NROW, NZ, NZ1;
  REAL ABEST, AIJ, AMAX, ATOLI, ATOLJ;

/*      ------------------------------------------------------------------
        Search cols of length nz = 1, then rows of length nz = 1,
        then   cols of length nz = 2, then rows of length nz = 2, etc.
        ------------------------------------------------------------------ */
  ABEST = ZERO;
  *IBEST = 0;
  KBEST = MAXMN+1;
  *MBEST = -1;
  NCOL = 0;
  NROW = 0;
  NZ1 = 0;
  for(NZ = 1; NZ <= MAXMN; NZ++) {
/*         nz1    = nz - 1
           if (mbest .le. nz1**2) go to 900 */
    if(KBEST<=NZ1)
      goto x900;
    if(*IBEST>0) {
      if(NCOL>=MAXCOL)
        goto x200;
    }
    if(NZ>LUSOL->m)
      goto x200;
/*         ---------------------------------------------------------------
           Search the set of columns of length  nz.
           --------------------------------------------------------------- */
    LQ1 = LUSOL->iqloc[NZ];
    LQ2 = LUSOL->n;
    if(NZ<LUSOL->m)
      LQ2 = LUSOL->iqloc[NZ+1]-1;
    for(LQ = LQ1; LQ <= LQ2; LQ++) {
      NCOL = NCOL+1;
      J = LUSOL->iq[LQ];
      LC1 = LUSOL->locc[J];
      LC2 = LC1+NZ1;
      AMAX = fabs(LUSOL->a[LC1]);
/*      Min size of pivots in col j */
      ATOLJ = AMAX/LTOL;
/*            Test all aijs in this column. */
      for(LC = LC1; LC <= LC2; LC++) {
        I = LUSOL->indc[LC];
        LEN1 = LUSOL->lenr[I]-1;
/*               merit  = nz1 * len1
                 if (merit .gt. mbest) continue; */
        if(LEN1>KBEST)
          continue;
/*               aij  has a promising merit.
                 Apply the Threshold Rook Pivoting stability test.
                 First we require aij to be sufficiently large
                 compared to other nonzeros in column j.
                 Then  we require aij to be sufficiently large
                 compared to other nonzeros in row    i. */
        AIJ = fabs(LUSOL->a[LC]);
        if(AIJ<ATOLJ)
          continue;
        if(AIJ*LTOL<AMAXR[I])
          continue;
/*               aij  is big enough. */
        MERIT = NZ1*LEN1;
        if(MERIT==*MBEST) {
/*                  Break ties.
                    (Initializing mbest < 0 prevents getting here if
                    nothing has been found yet.) */
          if(ABEST>=AIJ)
            continue;
        }
/*               aij  is the best pivot so far. */
        *IBEST = I;
        *JBEST = J;
        KBEST = LEN1;
        *MBEST = MERIT;
        ABEST = AIJ;
        if(NZ==1)
          goto x900;
      }
/*            Finished with that column. */
      if(*IBEST>0) {
        if(NCOL>=MAXCOL)
          goto x200;
      }
    }
/*         ---------------------------------------------------------------
           Search the set of rows of length  nz.
           --------------------------------------------------------------- */
x200:
/*    if (mbest .le. nz*nz1) go to 900 */
    if(KBEST<=NZ)
      goto x900;
    if(*IBEST>0) {
      if(NROW>=MAXROW)
        goto x290;
    }
    if(NZ>LUSOL->n)
      goto x290;
    LP1 = LUSOL->iploc[NZ];
    LP2 = LUSOL->m;
    if(NZ<LUSOL->n)
      LP2 = LUSOL->iploc[NZ+1]-1;
    for(LP = LP1; LP <= LP2; LP++) {
      NROW = NROW+1;
      I = LUSOL->ip[LP];
      LR1 = LUSOL->locr[I];
      LR2 = LR1+NZ1;
/*      Min size of pivots in row i */
      ATOLI = AMAXR[I]/LTOL;
      for(LR = LR1; LR <= LR2; LR++) {
        J = LUSOL->indr[LR];
        LEN1 = LUSOL->lenc[J]-1;
/*               merit  = nz1 * len1
                 if (merit .gt. mbest) continue; */
        if(LEN1>KBEST)
          continue;
/*               aij  has a promising merit.
                 Find where  aij  is in column j. */
        LC1 = LUSOL->locc[J];
        LC2 = LC1+LEN1;
        AMAX = fabs(LUSOL->a[LC1]);
        for(LC = LC1; LC <= LC2; LC++) {
          if(LUSOL->indc[LC]==I)
            break;
        }
/*               Apply the Threshold Rook Pivoting stability test.
                 First we require aij to be sufficiently large
                 compared to other nonzeros in row    i.
                 Then  we require aij to be sufficiently large
                 compared to other nonzeros in column j. */
        AIJ = fabs(LUSOL->a[LC]);
        if(AIJ<ATOLI)
          continue;
        if(AIJ*LTOL<AMAX)
          continue;
/*               aij  is big enough. */
        MERIT = NZ1*LEN1;
        if(MERIT==*MBEST) {
/*                  Break ties as before.
                    (Initializing mbest < 0 prevents getting here if
                    nothing has been found yet.) */
          if(ABEST>=AIJ)
            continue;
        }
/*               aij  is the best pivot so far. */
        *IBEST = I;
        *JBEST = J;
        KBEST = LEN1;
        *MBEST = MERIT;
        ABEST = AIJ;
        if(NZ==1)
          goto x900;
      }
/*            Finished with that row. */
      if(*IBEST>0) {
        if(NROW>=MAXROW)
          goto x290;
      }
    }
/*         See if it's time to quit. */
x290:
    if(*IBEST>0) {
      if(NROW>=MAXROW && NCOL>=MAXCOL)
        goto x900;
    }
/*         Press on with next nz. */
    NZ1 = NZ;
    if(*IBEST>0)
      KBEST = *MBEST/NZ1;
  }
x900:
;
}

/* ==================================================================
   lu1mSP  is intended for symmetric matrices that are either
   definite or quasi-definite.
   lu1mSP  uses a Markowitz criterion to select a pivot element for
   the next stage of a sparse LU factorization of a symmetric matrix,
   subject to a Threshold Symmetric Pivoting stability criterion
   (TSP) restricted to diagonal elements to preserve symmetry.
   This bounds the elements of L and U and should have rank-revealing
   properties analogous to Threshold Rook Pivoting for unsymmetric
   matrices.
   ------------------------------------------------------------------
   14 Dec 2002: First version of lu1mSP derived from lu1mRP.
                There is no safeguard to ensure that A is symmetric.
   14 Dec 2002: Current version of lu1mSP.
   ================================================================== */
void LU1MSP(LUSOLrec *LUSOL, int MAXMN, REAL LTOL, int MAXCOL,
            int *IBEST, int *JBEST, int *MBEST)
{
  int  I, J, KBEST, LC, LC1, LC2, LQ, LQ1, LQ2, MERIT, NCOL, NZ, NZ1;
  REAL ABEST, AIJ, AMAX, ATOLJ;

/*      ------------------------------------------------------------------
        Search cols of length nz = 1, then cols of length nz = 2, etc.
        ------------------------------------------------------------------ */
  ABEST = ZERO;
  *IBEST = 0;
  *MBEST = -1;
  KBEST = MAXMN+1;
  NCOL = 0;
  NZ1 = 0;
  for(NZ = 1; NZ <= MAXMN; NZ++) {
/*         nz1    = nz - 1
           if (mbest .le. nz1**2) go to 900 */
    if(KBEST<=NZ1)
      goto x900;
    if(*IBEST>0) {
      if(NCOL>=MAXCOL)
        goto x200;
    }
    if(NZ>LUSOL->m)
      goto x200;
/*         ---------------------------------------------------------------
           Search the set of columns of length  nz.
           --------------------------------------------------------------- */
    LQ1 = LUSOL->iqloc[NZ];
    LQ2 = LUSOL->n;
    if(NZ<LUSOL->m)
      LQ2 = LUSOL->iqloc[NZ+1]-1;
    for(LQ = LQ1; LQ <= LQ2; LQ++) {
      NCOL++;
      J = LUSOL->iq[LQ];
      LC1 = LUSOL->locc[J];
      LC2 = LC1+NZ1;
      AMAX = fabs(LUSOL->a[LC1]);
/*      Min size of pivots in col j */
      ATOLJ = AMAX/LTOL;
/*            Test all aijs in this column.
              Ignore everything except the diagonal. */
      for(LC = LC1; LC <= LC2; LC++) {
        I = LUSOL->indc[LC];
/*      Skip off-diagonals. */
        if(I!=J)
          continue;
/*               merit  = nz1 * nz1
                 if (merit .gt. mbest) continue; */
        if(NZ1>KBEST)
          continue;
/*               aij  has a promising merit.
                 Apply the Threshold Partial Pivoting stability test
                 (which is equivalent to Threshold Rook Pivoting for
                 symmetric matrices).
                 We require aij to be sufficiently large
                 compared to other nonzeros in column j. */
        AIJ = fabs(LUSOL->a[LC]);
        if(AIJ<ATOLJ)
          continue;
/*               aij  is big enough. */
        MERIT = NZ1*NZ1;
        if(MERIT==*MBEST) {
/*                  Break ties.
                    (Initializing mbest < 0 prevents getting here if
                    nothing has been found yet.) */
          if(ABEST>=AIJ)
            continue;
        }
/*               aij  is the best pivot so far. */
        *IBEST = I;
        *JBEST = J;
        KBEST = NZ1;
        *MBEST = MERIT;
        ABEST = AIJ;
        if(NZ==1)
          goto x900;
      }
/*            Finished with that column. */
      if(*IBEST>0) {
        if(NCOL>=MAXCOL)
          goto x200;
      }
    }
/*         See if it's time to quit. */
x200:
    if(*IBEST>0) {
      if(NCOL>=MAXCOL)
        goto x900;
    }
/*         Press on with next nz. */
    NZ1 = NZ;
    if(*IBEST>0)
      KBEST = *MBEST/NZ1;
  }
x900:
;
}

/* ==================================================================
   lu1mxc  moves the largest element in each of columns iq(k1:k2)
   to the top of its column.
   If k1 > k2, nothing happens.
   ------------------------------------------------------------------
   06 May 2002: (and earlier)
                All columns k1:k2 must have one or more elements.
   07 May 2002: Allow for empty columns.  The heap routines need to
                find 0.0 as the "largest element".
   ================================================================== */
void LU1MXC(LUSOLrec *LUSOL, int K1, int K2, int IX[])
{
  int  I, J, K, L, LC, LENJ;
  REAL AMAX;

  for(K = K1; K <= K2; K++) {
    J = IX[K];
    LC = LUSOL->locc[J];
    LENJ = LUSOL->lenc[J];
    if(LENJ==0)
      LUSOL->a[LC] = ZERO;
    else {
      L = idamax(LUSOL->lenc[J], LUSOL->a + LC - LUSOL_ARRAYOFFSET,1) + LC - 1;
      if(L>LC) {
        AMAX = LUSOL->a[L];
        LUSOL->a[L] = LUSOL->a[LC];
        LUSOL->a[LC] = AMAX;
        I = LUSOL->indc[L];
        LUSOL->indc[L] = LUSOL->indc[LC];
        LUSOL->indc[LC] = I;
      }
    }
  }
}

/* ==================================================================
   lu1mxr  finds the largest element in each of row ip(k1:k2)
   and stores it in Amaxr(*).  The nonzeros are stored column-wise
   in (a,indc,lenc,locc) and their structure is row-wise
   in (  indr,lenr,locr).
   If k1 > k2, nothing happens.
   ------------------------------------------------------------------
   11 Jun 2002: First version of lu1mxr.
                Allow for empty columns.
   ================================================================== */
void LU1MXR(LUSOLrec *LUSOL, int K1, int K2, int IX[], REAL AMAXR[])
{
#define FastMXR
#ifdef FastMXR
  static int  I, *J, *IC, K, LC, LC1, LC2, LR, LR1, LR2;
  static REAL AMAX;
#else
  int  I, J, K, LC, LC1, LC2, LR, LR1, LR2;
  REAL AMAX;
#endif

  for(K = K1; K <= K2; K++) {
    AMAX = ZERO;
    I = IX[K];
/*      Find largest element in row i. */
    LR1 = LUSOL->locr[I];
    LR2 = (LR1+LUSOL->lenr[I])-1;
#ifdef FastMXR
    for(LR = LR1, J = LUSOL->indr + LR1;
        LR <= LR2; LR++, J++) {
/*      Find where  aij  is in column  j. */
      LC1 = LUSOL->locc[*J];
      LC2 = LC1+LUSOL->lenc[*J];
      for(LC = LC1, IC = LUSOL->indc + LC1;
          LC < LC2; LC++, IC++) {
        if(*IC==I)
          break;
      }
      AMAX = MAX(AMAX,fabs(LUSOL->a[LC]));
    }
#else
    for(LR = LR1; LR <= LR2; LR++) {
      J = LUSOL->indr[LR];
/*      Find where  aij  is in column  j. */
      LC1 = LUSOL->locc[J];
      LC2 = (LC1+LUSOL->lenc[J])-1;
      for(LC = LC1; LC <= LC2; LC++) {
        if(LUSOL->indc[LC]==I)
          break;
      }
      AMAX = MAX(AMAX,fabs(LUSOL->a[LC]));
    }
#endif
    AMAXR[I] = AMAX;
  }
}


/* ==================================================================
   lu1ful computes a dense (full) LU factorization of the
   mleft by nleft matrix that remains to be factored at the
   beginning of the nrowu-th pass through the main loop of lu1fad.
   ------------------------------------------------------------------
   02 May 1989: First version.
   05 Feb 1994: Column interchanges added to lu1DPP.
   08 Feb 1994: ipinv reconstructed, since lu1pq3 may alter ip.
   ================================================================== */
void LU1FUL(LUSOLrec *LUSOL, int LEND, int LU1, MYBOOL TPP,
            int MLEFT, int NLEFT, int NRANK, int NROWU,
            int *LENL, int *LENU, int *NSING,
            MYBOOL KEEPLU, REAL SMALL, REAL D[], int IPVT[])
{
  int  L, I, J, IPBASE, LDBASE, LQ, LC1, LC2, LC, LD, LKK, LKN, LU, K, L1,
       L2, IBEST, JBEST, LA, LL, NROWD, NCOLD;
  REAL AI, AJ;

/*      ------------------------------------------------------------------
        If lu1pq3 moved any empty rows, reset ipinv = inverse of ip.
        ------------------------------------------------------------------ */
  if(NRANK<LUSOL->m) {
    for(L = 1; L <= LUSOL->m; L++) {
      I = LUSOL->ip[L];
      LUSOL->ipinv[I] = L;
    }
  }
/*      ------------------------------------------------------------------
        Copy the remaining matrix into the dense matrix D.
         ------------------------------------------------------------------ */
#ifdef LUSOLFastClear
  MEMCLEAR((D+1), LEND);
#else
/*   dload(LEND, ZERO, D, 1); */
  for(J = 1; J <= LEND; J++) 
    D[J] = ZERO;
#endif

  IPBASE = NROWU-1;
  LDBASE = 1-NROWU;
  for(LQ = NROWU; LQ <= LUSOL->n; LQ++) {
    J = LUSOL->iq[LQ];
    LC1 = LUSOL->locc[J];
    LC2 = (LC1+LUSOL->lenc[J])-1;
    for(LC = LC1; LC <= LC2; LC++) {
      I = LUSOL->indc[LC];
      LD = LDBASE+LUSOL->ipinv[I];
      D[LD] = LUSOL->a[LC];
    }
    LDBASE += MLEFT;
  }
/*      ------------------------------------------------------------------
        Call our favorite dense LU factorizer.
        ------------------------------------------------------------------ */
  if(TPP)
    LU1DPP(LUSOL, D,MLEFT,MLEFT,NLEFT,SMALL,NSING,IPVT,LUSOL->iq+NROWU-LUSOL_ARRAYOFFSET);
  else
    LU1DCP(LUSOL, D,MLEFT,MLEFT,NLEFT,SMALL,NSING,IPVT,LUSOL->iq+NROWU-LUSOL_ARRAYOFFSET);

/*      ------------------------------------------------------------------
        Move D to the beginning of A,
        and pack L and U at the top of a, indc, indr.
        In the process, apply the row permutation to ip.
        lkk points to the diagonal of U.
        ------------------------------------------------------------------ */
#ifdef LUSOLFastCopy
  MEMCOPY(LUSOL->a+1,D+1,LEND);
#else
  dcopy(LEND,D,1,LUSOL->a,1);
#endif
#ifdef ClassicdiagU
  LUSOL->diagU = LUSOL->a + (LUSOL->lena-LUSOL->n);
#endif
  LKK = 1;
  LKN = (LEND-MLEFT)+1;
  LU = LU1;
  for(K = 1; K <= MIN(MLEFT,NLEFT); K++) {
    L1 = IPBASE+K;
    L2 = IPBASE+IPVT[K];
    if(L1!=L2) {
      I = LUSOL->ip[L1];
      LUSOL->ip[L1] = LUSOL->ip[L2];
      LUSOL->ip[L2] = I;
    }
    IBEST = LUSOL->ip[L1];
    JBEST = LUSOL->iq[L1];
    if(KEEPLU) {
/*            ===========================================================
              Pack the next column of L.
              =========================================================== */
      LA = LKK;
      LL = LU;
      NROWD = 1;
      for(I = K+1; I <= MLEFT; I++) {
        LA++;
        AI = LUSOL->a[LA];
        if(fabs(AI)>SMALL) {
          NROWD = NROWD+1;
          LL--;
          LUSOL->a[LL] = AI;
          LUSOL->indc[LL] = LUSOL->ip[IPBASE+I];
          LUSOL->indr[LL] = IBEST;
        }
      }
/*            ===========================================================
              Pack the next row of U.
              We go backwards through the row of D
              so the diagonal ends up at the front of the row of  U.
              Beware -- the diagonal may be zero.
              =========================================================== */
      LA = LKN+MLEFT;
      LU = LL;
      NCOLD = 0;
      for(J = NLEFT; J >= K; J--) {
        LA = LA-MLEFT;
        AJ = LUSOL->a[LA];
        if(fabs(AJ)>SMALL || J==K) {
          NCOLD++;
          LU--;
          LUSOL->a[LU] = AJ;
          LUSOL->indr[LU] = LUSOL->iq[IPBASE+J];
        }
      }
      LUSOL->lenr[IBEST] = -NCOLD;
      LUSOL->lenc[JBEST] = -NROWD;
      *LENL = ((*LENL)+NROWD)-1;
      *LENU = (*LENU)+NCOLD;
      LKN++;
    }
    else {
/*            ===========================================================
              Store just the diagonal of U, in natural order.
              =========================================================== */
      LUSOL->diagU[JBEST] = LUSOL->a[LKK];
    }
    LKK += MLEFT+1;
  }
}


/* ==================================================================
   lu1or1  organizes the elements of an  m by n  matrix  A  as
   follows.  On entry, the parallel arrays   a, indc, indr,
   contain  nelem  entries of the form     aij,    i,    j,
   in any order.  nelem  must be positive.
   Entries not larger than the input parameter  small  are treated as
   zero and removed from   a, indc, indr.  The remaining entries are
   defined to be nonzero.  numnz  returns the number of such nonzeros
   and  Amax  returns the magnitude of the largest nonzero.
   The arrays  lenc, lenr  return the number of nonzeros in each
   column and row of  A.
   inform = 0  on exit, except  inform = 1  if any of the indices in
   indc, indr  imply that the element  aij  lies outside the  m by n
   dimensions of  A.
   ------------------------------------------------------------------
   xx Feb 1985: Original version.
   17 Oct 2000: a, indc, indr now have size lena to allow nelem = 0.
   ================================================================== */
void LU1OR1(LUSOLrec *LUSOL, REAL SMALL,
            REAL *AMAX, int *NUMNZ, int *LERR, int *INFORM)
{
  int I, J, L, LDUMMY;

#ifdef LUSOLFastClear
  MEMCLEAR((LUSOL->lenr+1), LUSOL->m);
  MEMCLEAR((LUSOL->lenc+1), LUSOL->n);
#else
  for(I = 1; I <= LUSOL->m; I++)
    LUSOL->lenr[I] = ZERO;
  for(I = 1; I <= LUSOL->n; I++)
    LUSOL->lenc[I] = ZERO;
#endif

  *AMAX = 0;
  *NUMNZ = LUSOL->nelem;
  L = LUSOL->nelem+1;
  for(LDUMMY = 1; LDUMMY <= LUSOL->nelem; LDUMMY++) {
    L--;
    if(fabs(LUSOL->a[L])>SMALL) {
      I = LUSOL->indc[L];
      J = LUSOL->indr[L];
      *AMAX = MAX(*AMAX,fabs(LUSOL->a[L]));
      if(I<1 || I>LUSOL->m)
        goto x910;
      if(J<1 || J>LUSOL->n)
        goto x910;
      LUSOL->lenr[I]++;
      LUSOL->lenc[J]++;
    }
    else {
/*            Replace a negligible element by last element.  Since
              we are going backwards, we know the last element is ok. */
      LUSOL->a[L] = LUSOL->a[*NUMNZ];
      LUSOL->indc[L] = LUSOL->indc[*NUMNZ];
      LUSOL->indr[L] = LUSOL->indr[*NUMNZ];
      (*NUMNZ)--;
    }
  }
  *LERR = 0;
  *INFORM = LUSOL_INFORM_LUSUCCESS;
  return;

x910:
  *LERR = L;
  *INFORM = LUSOL_INFORM_LUSINGULAR;
}

/* ==================================================================
   lu1or2  sorts a list of matrix elements  a(i,j)  into column
   order, given  numa  entries  a(i,j),  i,  j  in the parallel
   arrays  a, inum, jnum  respectively.  The matrix is assumed
   to have  n  columns and an arbitrary number of rows.
   On entry,  len(*)  must contain the length of each column.
   On exit,  a(*) and inum(*)  are sorted,  jnum(*) = 0,  and
   loc(j)  points to the start of column j.
   lu1or2  is derived from  mc20ad,  a routine in the Harwell
   Subroutine Library, author J. K. Reid.
   ------------------------------------------------------------------
   xx Feb 1985: Original version.
   17 Oct 2000: a, inum, jnum now have size lena to allow nelem = 0.
   ================================================================== */
void LU1OR2(LUSOLrec *LUSOL)
{
  REAL ACE, ACEP;
  int  L, J, I, JCE, ICE, ICEP, JCEP, JA, JB;

/*      Set  loc(j)  to point to the beginning of column  j. */
  L = 1;
  for(J = 1; J <= LUSOL->n; J++) {
    LUSOL->locc[J] = L;
    L += LUSOL->lenc[J];
  }
/*      Sort the elements into column order.
        The algorithm is an in-place sort and is of order  numa. */
  for(I = 1; I <= LUSOL->nelem; I++) {
/*         Establish the current entry. */
    JCE = LUSOL->indr[I];
    if(JCE==0)
      continue;
    ACE = LUSOL->a[I];
    ICE = LUSOL->indc[I];
    LUSOL->indr[I] = 0;
/*         Chain from current entry. */
    for(J = 1; J <= LUSOL->nelem; J++) {
/*            The current entry is not in the correct position.
              Determine where to store it. */
      L = LUSOL->locc[JCE];
      LUSOL->locc[JCE]++;
/*            Save the contents of that location. */
      ACEP = LUSOL->a[L];
      ICEP = LUSOL->indc[L];
      JCEP = LUSOL->indr[L];
/*            Store current entry. */
      LUSOL->a[L] = ACE;
      LUSOL->indc[L] = ICE;
      LUSOL->indr[L] = 0;
/*            If next current entry needs to be processed,
              copy it into current entry. */
      if(JCEP==0)
        break;
      ACE = ACEP;
      ICE = ICEP;
      JCE = JCEP;
    }
  }
/*      Reset loc(j) to point to the start of column j. */
  JA = 1;
  for(J = 1; J <= LUSOL->n; J++) {
    JB = LUSOL->locc[J];
    LUSOL->locc[J] = JA;
    JA = JB;
  }
}

/* ==================================================================
   lu1or3  looks for duplicate elements in an  m by n  matrix  A
   defined by the column list  indc, lenc, locc.
   iw  is used as a work vector of length  m.
   ------------------------------------------------------------------
   xx Feb 1985: Original version.
   17 Oct 2000: indc, indr now have size lena to allow nelem = 0.
   ================================================================== */
void LU1OR3(LUSOLrec *LUSOL, int *LERR, int *INFORM)
{
  int I, J, L1, L2, L;

#ifdef LUSOLFastClear
  MEMCLEAR((LUSOL->ip+1), LUSOL->m);
#else
  for(I = 1; I <= LUSOL->m; I++)
    LUSOL->ip[I] = ZERO;
#endif

  for(J = 1; J <= LUSOL->n; J++) {
    if(LUSOL->lenc[J]>0) {
      L1 = LUSOL->locc[J];
      L2 = (L1+LUSOL->lenc[J])-1;
      for(L = L1; L <= L2; L++) {
        I = LUSOL->indc[L];
        if(LUSOL->ip[I]==J)
          goto x910;
        LUSOL->ip[I] = J;
      }
    }
  }
  *INFORM = LUSOL_INFORM_LUSUCCESS;
  return;
x910:
  *LERR = L;
  *INFORM = LUSOL_INFORM_LUSINGULAR;
}

/* ==================================================================
   lu1or4 constructs a row list  indr, locr
   from a corresponding column list  indc, locc,
   given the lengths of both columns and rows in  lenc, lenr.
   ------------------------------------------------------------------
   xx Feb 1985: Original version.
   17 Oct 2000: indc, indr now have size lena to allow nelem = 0.
   ================================================================== */
void LU1OR4(LUSOLrec *LUSOL)
{
  int L, I, L2, J, JDUMMY, L1, LR;

/*      Initialize  locr(i)  to point just beyond where the
        last component of row  i  will be stored. */
  L = 1;
  for(I = 1; I <= LUSOL->m; I++) {
    L += LUSOL->lenr[I];
    LUSOL->locr[I] = L;
  }
/*      By processing the columns backwards and decreasing  locr(i)
        each time it is accessed, it will end up pointing to the
        beginning of row  i  as required. */
  L2 = LUSOL->nelem;
  J = LUSOL->n+1;
  for(JDUMMY = 1; JDUMMY <= LUSOL->n; JDUMMY++) {
    J = J-1;
    if(LUSOL->lenc[J]>0) {
      L1 = LUSOL->locc[J];
      for(L = L1; L <= L2; L++) {
        I = LUSOL->indc[L];
        LR = LUSOL->locr[I]-1;
        LUSOL->locr[I] = LR;
        LUSOL->indr[LR] = J;
      }
      L2 = L1-1;
    }
  }
}

/* ==================================================================
   lu1pen deals with pending fill-in in the row file.
   ------------------------------------------------------------------
   ifill(ll) says if a row involved in the new column of L
             has to be updated.  If positive, it is the total
             length of the final updated row.
   jfill(lu) says if a column involved in the new row of U
             contains any pending fill-ins.  If positive, it points
             to the first fill-in in the column that has yet to be
             added to the row file.
   ------------------------------------------------------------------
   16 Apr 1989: First version of lu1pen.
   23 Mar 2001: ilast used and updated.
   ================================================================== */
void LU1PEN(LUSOLrec *LUSOL, int NSPARE, int *ILAST,
            int LPIVC1, int LPIVC2, int LPIVR1, int LPIVR2,
            int *LROW, int IFILL[], int JFILL[])
{
  int  LL, LC, L, I, LR1, LR2, LR, LU, J, LC1, LC2, LAST;

  LL = 0;
  for(LC = LPIVC1; LC <= LPIVC2; LC++) {
    LL++;
    if(IFILL[LL]==0)
      continue;
/*      Another row has pending fill.
        First, add some spare space at the }
        of the current last row. */
#if 1
    LC1 = (*LROW)+1;
    LC2 = (*LROW)+NSPARE;
    *LROW = LC2;
    for(L = LC1; L <= LC2; L++) {
#else
    for(L = (*LROW)+1; L <= (*LROW)+NSPARE; L++) {
      *LROW = L;  /* ******* ERROR ???? */
#endif
      LUSOL->indr[L] = 0;
    }
/*      Now move row i to the end of the row file. */
    I = LUSOL->indc[LC];
    *ILAST = I;
    LR1 = LUSOL->locr[I];
    LR2 = (LR1+LUSOL->lenr[I])-1;
    LUSOL->locr[I] = (*LROW)+1;
    for(LR = LR1; LR <= LR2; LR++) {
      (*LROW)++;
      LUSOL->indr[*LROW] = LUSOL->indr[LR];
      LUSOL->indr[LR] = 0;
    }
    (*LROW) += IFILL[LL];
  }
/*         Scan all columns of  D  and insert the pending fill-in
           into the row file. */
  LU = 1;
  for(LR = LPIVR1; LR <= LPIVR2; LR++) {
    LU++;
    if(JFILL[LU]==0)
      continue;
    J = LUSOL->indr[LR];
    LC1 = (LUSOL->locc[J]+JFILL[LU])-1;
    LC2 = (LUSOL->locc[J]+LUSOL->lenc[J])-1;
    for(LC = LC1; LC <= LC2; LC++) {
      I = LUSOL->indc[LC]-LUSOL->m;
      if(I>0) {
        LUSOL->indc[LC] = I;
        LAST = LUSOL->locr[I]+LUSOL->lenr[I];
        LUSOL->indr[LAST] = J;
        LUSOL->lenr[I]++;
      }
    }
  }
}


/* ==================================================================
   lu1fad  is a driver for the numerical phase of lu1fac.
   At each stage it computes a column of  L  and a row of  U,
   using a Markowitz criterion to select the pivot element,
   subject to a stability criterion that bounds the elements of  L.
   ------------------------------------------------------------------
   Local variables
   ---------------
   lcol   is the length of the column file.  It points to the last
          nonzero in the column list.
   lrow   is the analogous quantity for the row file.
   lfile  is the file length (lcol or lrow) after the most recent
          compression of the column list or row list.
   nrowd  and  ncold  are the number of rows and columns in the
          matrix defined by the pivot column and row.  They are the
          dimensions of the submatrix D being altered at this stage.
   melim  and  nelim  are the number of rows and columns in the
          same matrix D, excluding the pivot column and row.
   mleft  and  nleft  are the number of rows and columns
          still left to be factored.
   nzchng is the increase in nonzeros in the matrix that remains
          to be factored after the current elimination
          (usually negative).
   nzleft is the number of nonzeros still left to be factored.
   nspare is the space we leave at the end of the last row or
          column whenever a row or column is being moved to the }
          of its file.  nspare = 1 or 2 might help reduce the
          number of file compressions when storage is tight.
   The row and column ordering permutes A into the form
                      ------------------------
                       \                     |
                        \         U1         |
                         \                   |
                          --------------------
                          |\
                          | \
                          |  \
          P A Q   =       |   \
                          |    \
                          |     --------------
                          |     |            |
                          |     |            |
                          | L1  |     A2     |
                          |     |            |
                          |     |            |
                          --------------------
   where the block A2 is factored as  A2 = L2 U2.
   The phases of the factorization are as follows.
   Utri   is true when U1 is being determined.
          Any column of length 1 is accepted immediately (if TPP).
   Ltri   is true when L1 is being determined.
          lu1mar exits as soon as an acceptable pivot is found
          in a row of length 1.
   spars1 is true while the density of the (modified) A2 is less
          than the parameter dens1 = parmlu(7) = 0.3 say.
          lu1mar searches maxcol columns and maxrow rows,
          where  maxcol = luparm(3),  maxrow = maxcol - 1.
          lu1mxc is used to keep the biggest element at the top
          of all remaining columns.
   spars2 is true while the density of the modified A2 is less
          than the parameter dens2 = parmlu(8) = 0.6 say.
          lu1mar searches maxcol columns and no rows.
          lu1mxc could fix up only the first maxcol cols (with TPP).
          22 Sep 2000:  For simplicity, lu1mxc fixes all
                        modified cols.
   dense  is true once the density of A2 reaches dens2.
          lu1mar searches only 1 column (the shortest).
          lu1mxc could fix up only the first column (with TPP).
   ------------------------------------------------------------------
   00 Jan 1986  Version documented in LUSOL paper:
                Gill, Murray, Saunders and Wright (1987),
                Maintaining LU factors of a general sparse matrix,
                Linear algebra and its applications 88/89, 239-270.
   02 Feb 1989  Following Suhl and Aittoniemi (1987), the largest
                element in each column is now kept at the start of
                the column, i.e. in position locc(j) of a and indc.
                This should speed up the Markowitz searches.
                To save time on highly triangular matrices, we wait
                until there are no further columns of length 1
                before setting and maintaining that property.
   12 Apr 1989  ipinv and iqinv added (inverses of ip and iq)
                to save searching ip and iq for rows and columns
                altered in each elimination step.  (Used in lu1pq2)
   19 Apr 1989  Code segmented to reduce its size.
                lu1gau does most of the Gaussian elimination work.
                lu1mar does just the Markowitz search.
                lu1mxc moves biggest elements to top of columns.
                lu1pen deals with pending fill-in in the row list.
                lu1pq2 updates the row and column permutations.
   26 Apr 1989  maxtie replaced by maxcol, maxrow in the Markowitz
                search.  maxcol, maxrow change as density increases.
   25 Oct 1993  keepLU implemented.
   07 Feb 1994  Exit main loop early to finish off with a dense LU.
                densLU tells lu1fad whether to do it.
   21 Dec 1994  Bug fixed.  nrank was wrong after the call to lu1ful.
   12 Nov 1999  A parallel version of dcopy gave trouble in lu1ful
                during left-shift of dense matrix D within a(*).
                Fixed this unexpected problem here in lu1fad
                by making sure the first and second D don't overlap.
   13 Sep 2000  TCP (Threshold Complete Pivoting) implemented.
                lu2max added
                (finds aijmax from biggest elems in each col).
                Utri, Ltri and Spars1 phases apply.
                No switch to Dense CP yet.  (Only TPP switches.)
   14 Sep 2000  imax needed to remember row containing aijmax.
   22 Sep 2000  For simplicity, lu1mxc always fixes all modified cols.
                (TPP spars2 used to fix just the first maxcol cols.)
   08 Nov 2000: Speed up search for aijmax.
                Don't need to search all columns if the elimination
                didn't alter the col containing the current aijmax.
   21 Nov 2000: lu1slk implemented for Utri phase with TCP
                to guard against deceptive triangular matrices.
                (Utri used to have aijtol >= 0.9999 to include
                slacks, but this allows other 1s to be accepted.)
                Utri now accepts slacks, but applies normal aijtol
                test to other pivots.
   28 Nov 2000: TCP with empty cols must call lu1mxc and lu2max
                with ( lq1, n, ... ), not just ( 1, n, ... ).
   23 Mar 2001: lu1fad bug with TCP.
                A col of length 1 might not be accepted as a pivot.
                Later it appears in a pivot row and temporarily
                has length 0 (when pivot row is removed
                but before the column is filled in).  If it is the
                last column in storage, the preceding col also thinks
                it is "last".  Trouble arises when the preceding col
                needs fill-in -- it overlaps the real "last" column.
                (Very rarely, same trouble might have happened if
                the drop tolerance caused columns to have length 0.)
                Introduced ilast to record the last row in row file,
                           jlast to record the last col in col file.
                lu1rec returns ilast = indr(lrow + 1)
                            or jlast = indc(lcol + 1).
                (Should be an output parameter, but didn't want to
                alter lu1rec's parameter list.)
                lu1rec also treats empty rows or cols safely.
                (Doesn't eliminate them!)
   26 Apr 2002: Heap routines added for TCP.
                lu2max no longer needed.
                imax, jmax used only for printing.
   01 May 2002: lu1DCP implemented (dense complete pivoting).
                Both TPP and TCP now switch to dense LU
                when density exceeds dens2.
   06 May 2002: In dense mode, store diag(U) in natural order.
   09 May 2002: lu1mCP implemented (Markowitz TCP via heap).
   11 Jun 2002: lu1mRP implemented (Markowitz TRP).
   28 Jun 2002: Fixed call to lu1mxr.
   14 Dec 2002: lu1mSP implemented (Markowitz TSP).
   15 Dec 2002: Both TPP and TSP can grab cols of length 1
                during Utri.
   ================================================================== */
void LU1FAD(LUSOLrec *LUSOL,
#ifdef ClassicHamaxR
            int LENA2, int LENH, REAL HA[], int HJ[], int HK[], REAL AMAXR[],
#endif
            int *INFORM, int *LENL, int *LENU, int *MINLEN,
            int *MERSUM, int *NUTRI, int *NLTRI,
            int *NDENS1, int *NDENS2, int *NRANK,
            REAL *LMAX, REAL *UMAX, REAL *DUMAX, REAL *DUMIN, REAL *AKMAX)
{
  MYBOOL UTRI, LTRI, SPARS1, SPARS2, DENSE, DENSLU, KEEPLU, TCP, TPP, TRP,TSP;
  int    HLEN, HOPS, H, LPIV, LPRINT, MAXCOL, MAXROW, ILAST, JLAST, LFILE, LROW, LCOL,
         MINMN, MAXMN, NZLEFT, NSPARE, LU1, KK, J, LC, MLEFT, NLEFT, NROWU,
         LQ1, LQ2, JBEST, LQ, I, IBEST, MBEST, LEND, NFREE, LD, NCOLD, NROWD,
         MELIM, NELIM, JMAX, IMAX, LL1, LSAVE, LFREE, LIMIT, MINFRE, LPIVR, LPIVR1, LPIVR2,
         L, LPIVC, LPIVC1, LPIVC2, KBEST, LU, LR, LENJ, LC1, LAST, LL, LS,
         LENI, LR1, LFIRST, NFILL, NZCHNG, K, MRANK, NSING;
  REAL   LIJ, LTOL, SMALL, USPACE, DENS1, DENS2, AIJMAX, AIJTOL, AMAX, ABEST, DIAG, V;
#ifdef ClassicHamaxR
  int    LDIAGU;
#else
  int    LENA2 = LUSOL->lena;
#endif

#ifdef UseTimer
  int    eltime, mktime, ntime;
  timer ( "start", 3 );
  ntime = LUSOL->n / 4;
#endif

#ifdef ForceInitialization
  AIJMAX = 0;
  AIJTOL = 0;
  HLEN   = 0;
  JBEST  = 0;
  IBEST  = 0;
  MBEST  = 0;
  LEND   = 0;
  LD     = 0;
#endif

  LPRINT = LUSOL->luparm[LUSOL_IP_PRINTLEVEL];
  MAXCOL = LUSOL->luparm[LUSOL_IP_MARKOWITZ_MAXCOL];
  LPIV   = LUSOL->luparm[LUSOL_IP_PIVOTTYPE];
  KEEPLU = (MYBOOL) (LUSOL->luparm[LUSOL_IP_KEEPLU]!=FALSE);
/*      Threshold Partial   Pivoting (normal). */
  TPP = (MYBOOL) (LPIV==LUSOL_PIVMOD_TPP);
/*      Threshold Rook      Pivoting */
  TRP = (MYBOOL) (LPIV==LUSOL_PIVMOD_TRP);
/*      Threshold Complete  Pivoting. */
  TCP = (MYBOOL) (LPIV==LUSOL_PIVMOD_TCP);
/*      Threshold Symmetric Pivoting. */
  TSP = (MYBOOL) (LPIV==LUSOL_PIVMOD_TSP);
  DENSLU = FALSE;
  MAXROW = MAXCOL-1;
/*      Assume row m is last in the row file. */
  ILAST = LUSOL->m;
/*      Assume col n is last in the col file. */
  JLAST = LUSOL->n;
  LFILE = LUSOL->nelem;
  LROW = LUSOL->nelem;
  LCOL = LUSOL->nelem;
  MINMN = MIN(LUSOL->m,LUSOL->n);
  MAXMN = MAX(LUSOL->m,LUSOL->n);
  NZLEFT = LUSOL->nelem;
  NSPARE = 1;

  if(KEEPLU)
    LU1 = LENA2+1;
  else {
/*         Store only the diagonals of U in the top of memory. */
#ifdef ClassicdiagU
    LDIAGU = LENA2-LUSOL->n;
    LU1 = LDIAGU+1;
    LUSOL->diagU = LUSOL->a+LDIAGU;
#else
    LU1 = LENA2+1;
#endif
  }

  LTOL = LUSOL->parmlu[LUSOL_RP_FACTORMAX_Lij];
  SMALL = LUSOL->parmlu[LUSOL_RP_ZEROTOLERANCE];
  USPACE = LUSOL->parmlu[LUSOL_RP_COMPSPACE_U];
  DENS1 = LUSOL->parmlu[LUSOL_RP_MARKOWITZ_CONLY];
  DENS2 = LUSOL->parmlu[LUSOL_RP_MARKOWITZ_DENSE];
  UTRI = TRUE;
  LTRI = FALSE;
  SPARS1 = FALSE;
  SPARS2 = FALSE;
  DENSE = FALSE;
/*      Check parameters. */
  LTOL = MAX(LTOL,1.0001E+0);
  DENS1 = MIN(DENS1,DENS2);
/*      Initialize output parameters.
        lenL, lenU, minlen, mersum, nUtri, nLtri, ndens1, ndens2, nrank
        are already initialized by lu1fac. */
  *LMAX  = ZERO;
  *UMAX  = ZERO;
  *DUMAX = ZERO;
  *DUMIN = LUSOL_BIGNUM;
  if(LUSOL->nelem==0)
    *DUMIN = ZERO;
  *AKMAX = ZERO;
  HOPS = 0;
/*      More initialization.
        Don't worry yet about lu1mxc. */
  if(TPP || TSP) {
    AIJMAX = ZERO;
    AIJTOL = ZERO;
    HLEN = 1;
/*      TRP or TCP */
  }
  else {
/*      Move biggest element to top of each column.
        Set w(*) to mark slack columns (unit vectors). */
    LU1MXC(LUSOL, 1,LUSOL->n,LUSOL->iq);
    LU1SLK(LUSOL);
  }
  if(TRP)
/*      Find biggest element in each row. */
#ifdef ClassicHamaxR
    LU1MXR(LUSOL, 1,LUSOL->m,LUSOL->ip,AMAXR);
#else
    LU1MXR(LUSOL, 1,LUSOL->m,LUSOL->ip,LUSOL->amaxr);
#endif

  if(TCP) {
/*      Set Ha(1:Hlen) = biggest element in each column,
            Hj(1:Hlen) = corresponding column indices. */
    HLEN = 0;
    for(KK = 1; KK <= LUSOL->n; KK++) {
      HLEN++;
      J = LUSOL->iq[KK];
      LC = LUSOL->locc[J];
#ifdef ClassicHamaxR
      HA[HLEN] = fabs(LUSOL->a[LC]);
      HJ[HLEN] = J;
      HK[J] = HLEN;
#else
      LUSOL->Ha[HLEN] = fabs(LUSOL->a[LC]);
      LUSOL->Hj[HLEN] = J;
      LUSOL->Hk[J] = HLEN;
#endif
    }
/*      Build the heap, creating new Ha, Hj and setting Hk(1:Hlen). */
#ifdef ClassicHamaxR
    HBUILD(HA,HJ,HK,HLEN,&HOPS);
#else
    HBUILD(LUSOL->Ha,LUSOL->Hj,LUSOL->Hk,HLEN,&HOPS);
#endif
  }
/*      ------------------------------------------------------------------
        Start of main loop.
        ------------------------------------------------------------------ */
  MLEFT = LUSOL->m+1;
  NLEFT = LUSOL->n+1;
  for(NROWU = 1; NROWU <= MINMN; NROWU++) {
#ifdef UseTimer
    mktime = (nrowu / ntime) + 4;
    eltime = (nrowu / ntime) + 9;
#endif
    MLEFT--;
    NLEFT--;
/*         Bail out if there are no nonzero rows left. */
    if(LUSOL->iploc[1]>LUSOL->m)
      goto x900;
/*      For TCP, the largest Aij is at the top of the heap. */
   if(TCP) {
/*
              Marvelously easy */
#ifdef ClassicHamaxR
      AIJMAX = HA[1];
#else
      AIJMAX = LUSOL->Ha[1];
#endif
      *AKMAX = MAX(*AKMAX,AIJMAX);
      AIJTOL = AIJMAX/LTOL;
    }
/*         ===============================================================
           Find a suitable pivot element.
           =============================================================== */
    if(UTRI) {
/*            ------------------------------------------------------------
              So far all columns have had length 1.
              We are still looking for the (backward) triangular part of A
              that forms the first rows and columns of U.
              ------------------------------------------------------------ */
      LQ1 = LUSOL->iqloc[1];
      LQ2 = LUSOL->n;
      if(LUSOL->m>1)
        LQ2 = LUSOL->iqloc[2]-1;
/*      There are more cols of length 1. */
      if(LQ1<=LQ2) {
        if(TPP || TSP) {
/*      Grab the first one. */
          JBEST = LUSOL->iq[LQ1];
/*      Scan all columns of length 1 ... TRP or TCP */
        }
        else {
          JBEST = 0;
          for(LQ = LQ1; LQ <= LQ2; LQ++) {
            J = LUSOL->iq[LQ];
/*      Accept a slack */
            if(LUSOL->w[J]>ZERO) {
              JBEST = J;
              goto x250;
            }
            LC = LUSOL->locc[J];
            AMAX = fabs(LUSOL->a[LC]);
            if(TRP) {
              I = LUSOL->indc[LC];
#ifdef ClassicHamaxR
              AIJTOL = AMAXR[I]/LTOL;
#else
              AIJTOL = LUSOL->amaxr[I]/LTOL;
#endif
            }
            if(AMAX>=AIJTOL) {
              JBEST = J;
              goto x250;
            }
          }
        }
x250:
        if(JBEST>0) {
          LC = LUSOL->locc[JBEST];
          IBEST = LUSOL->indc[LC];
          MBEST = 0;
          goto x300;
        }
      }
/*            This is the end of the U triangle.
              We will not return to this part of the code.
              TPP and TSP call lu1mxc for the first time
              (to move biggest element to top of each column). */
      if(LPRINT>=LUSOL_MSG_PIVOT)
        LUSOL_report(LUSOL, 0, "Utri ended.  spars1 = TRUE\n");
      UTRI = FALSE;
      LTRI = TRUE;
      SPARS1 = TRUE;
      *NUTRI = NROWU-1;
      if(TPP || TSP)
        LU1MXC(LUSOL, LQ1,LUSOL->n,LUSOL->iq);
    }
    if(SPARS1) {
/*            ------------------------------------------------------------
              Perform a Markowitz search.
              Search cols of length 1, then rows of length 1,
              then   cols of length 2, then rows of length 2, etc.
              ------------------------------------------------------------ */
#ifdef UseTimer
        timer ( "start", mktime );
#endif
/*      12 Jun 2002: Next line disables lu1mCP below
              if (TPP) then */
      if(TPP || TCP) {
        LU1MAR(LUSOL, MAXMN,TCP,AIJTOL,LTOL,MAXCOL,MAXROW,&IBEST,&JBEST,&MBEST);
      }
      else if(TRP) {
#ifdef ClassicHamaxR
        LU1MRP(LUSOL, MAXMN,LTOL,MAXCOL,MAXROW,&IBEST,&JBEST,&MBEST,AMAXR);
#else
        LU1MRP(LUSOL, MAXMN,LTOL,MAXCOL,MAXROW,&IBEST,&JBEST,&MBEST,LUSOL->amaxr);
#endif
/*      else if (TCP) {
        lu1mCP( m    , n     , lena  , aijtol,
                      ibest, jbest , mbest ,
                      a    , indc  , indr  ,
                      lenc , lenr  , locc  ,
                      Hlen , Ha    , Hj    ) */
      }
      else if(TSP) {
        LU1MSP(LUSOL, MAXMN,LTOL,MAXCOL,&IBEST,&JBEST,&MBEST);
        if(IBEST==0)
          goto x990;
      }
#ifdef UseTimer
      timer ( "finish", mktime );
#endif
      if(LTRI) {
/*               So far all rows have had length 1.
                 We are still looking for the (forward) triangle of A
                 that forms the first rows and columns of L. */
        if(MBEST>0) {
          LTRI = FALSE;
          *NLTRI = NROWU-1-*NUTRI;
          if(LPRINT>=LUSOL_MSG_PIVOT)
            LUSOL_report(LUSOL, 0, "Ltri ended.\n");
        }
      }
      else {
/*               See if what's left is as dense as dens1. */
        if(NZLEFT>=(DENS1*MLEFT)*NLEFT) {
          SPARS1 = FALSE;
          SPARS2 = TRUE;
          *NDENS1 = NLEFT;
          MAXROW = 0;
          if(LPRINT>=LUSOL_MSG_PIVOT)
            LUSOL_report(LUSOL, 0, "spars1 ended.  spars2 = TRUE\n");
        }
      }
    }
    else if(SPARS2 || DENSE) {
/*            ------------------------------------------------------------
              Perform a restricted Markowitz search,
              looking at only the first maxcol columns.  (maxrow = 0.)
              ------------------------------------------------------------ */
#ifdef UseTimer
      timer ( "start", mktime );
#endif
/*      12 Jun 2002: Next line disables lu1mCP below
              if (TPP) then */
      if(TPP || TCP) {
        LU1MAR(LUSOL, MAXMN,TCP,AIJTOL,LTOL,MAXCOL,MAXROW,&IBEST,&JBEST,&MBEST);
      }
      else if(TRP) {
#ifdef ClassicHamaxR
        LU1MRP(LUSOL, MAXMN,LTOL,MAXCOL,MAXROW,&IBEST,&JBEST,&MBEST,AMAXR);
#else
        LU1MRP(LUSOL, MAXMN,LTOL,MAXCOL,MAXROW,&IBEST,&JBEST,&MBEST,LUSOL->amaxr);
#endif
/*      else if (TCP) {
        lu1mCP( m    , n     , lena  , aijtol,
                      ibest, jbest , mbest ,
                      a    , indc  , indr  ,
                      lenc , lenr  , locc  ,
                      Hlen , Ha    , Hj    ) */
      }
      else if(TSP) {
        LU1MSP(LUSOL, MAXMN,LTOL,MAXCOL,&IBEST,&JBEST,&MBEST);
        if(IBEST==0)
          goto x985;
      }
#ifdef UseTimer
      timer ( "finish", mktime );
#endif
/*            See if what's left is as dense as dens2. */
      if(SPARS2) {
        if(NZLEFT>=(DENS2*MLEFT)*NLEFT) {
          SPARS2 = FALSE;
          DENSE = TRUE;
          *NDENS2 = NLEFT;
          MAXCOL = 1;
          if(LPRINT>=LUSOL_MSG_PIVOT)
            LUSOL_report(LUSOL, 0, "spars2 ended.  dense = TRUE\n");
        }
      }
    }
/*         ---------------------------------------------------------------
           See if we can finish quickly.
           --------------------------------------------------------------- */
    if(DENSE) {
      LEND = MLEFT*NLEFT;
      NFREE = LU1-1;
      if(NFREE>=2*LEND) {
/*               There is room to treat the remaining matrix as
                 a dense matrix D.
                 We may have to compress the column file first.
                 12 Nov 1999: D used to be put at the
                              beginning of free storage (lD = lcol + 1).
                              Now put it at the end     (lD = lu1 - lenD)
                              so the left-shift in lu1ful will not
                              involve overlapping storage
                              (fatal with parallel dcopy).
   */
        DENSLU = TRUE;
        *NDENS2 = NLEFT;
        LD = LU1-LEND;
        if(LCOL>=LD) {
          LU1REC(LUSOL, LUSOL->n,TRUE,&LCOL,
                        LUSOL->indc,LUSOL->lenc,LUSOL->locc);
          LFILE = LCOL;
          JLAST = LUSOL->indc[LCOL+1];
        }
        goto x900;
      }
    }
/*         ===============================================================
           The best  aij  has been found.
           The pivot row  ibest  and the pivot column  jbest
           Define a dense matrix  D  of size  nrowd  by  ncold.
           =============================================================== */
x300:
    NCOLD = LUSOL->lenr[IBEST];
    NROWD = LUSOL->lenc[JBEST];
    MELIM = NROWD-1;
    NELIM = NCOLD-1;
    (*MERSUM) += MBEST;
    (*LENL) += MELIM;
    (*LENU) += NCOLD;
    if(LPRINT>=LUSOL_MSG_PIVOT) {
      if(NROWU==1)
        LUSOL_report(LUSOL, 0, "lu1fad debug:\n");
      if(TPP || TRP || TSP) {
        LUSOL_report(LUSOL, 0, "nrowu:%7d   i,jbest:%7d,%7d   nrowd,ncold:%6d,%6d\n",
                            NROWU, IBEST,JBEST, NROWD,NCOLD);
/*      TCP */
      }
      else {
#ifdef ClassicHamaxR
        JMAX = HJ[1];
#else
        JMAX = LUSOL->Hj[1];
#endif
        IMAX = LUSOL->indc[LUSOL->locc[JMAX]];
        LUSOL_report(LUSOL, 0, "nrowu:%7d   i,jbest:%7d,%7d   nrowd,ncold:%6d,%6d   i,jmax:%7d,%7d   aijmax:%g\n",
                            NROWU, IBEST,JBEST, NROWD,NCOLD, IMAX,JMAX, AIJMAX);
      }
    }
/*         ===============================================================
           Allocate storage for the next column of  L  and next row of  U.
           Initially the top of a, indc, indr are used as follows:
                      ncold       melim       ncold        melim
           a      |...........|...........|ujbest..ujn|li1......lim|
           indc   |...........|  lenr(i)  |  lenc(j)  |  markl(i)  |
           indr   |...........| iqloc(i)  |  jfill(j) |  ifill(i)  |
                 ^           ^             ^           ^            ^
                 lfree   lsave             lu1         ll1          oldlu1
           Later the correct indices are inserted:
           indc   |           |           |           |i1........im|
           indr   |           |           |jbest....jn|ibest..ibest|
           =============================================================== */
    if(!KEEPLU) {
/*            Always point to the top spot.
              Only the current column of L and row of U will
              take up space, overwriting the previous ones. */
#ifdef ClassicHamaxR
      LU1 = LDIAGU+1;
#else
      LU1 = LENA2+1;
#endif
    }
    /* Update (left-shift) pointers to make room for the new data */
    LL1 = LU1-MELIM;
    LU1 = LL1-NCOLD;
    LSAVE = LU1-NROWD;
    LFREE = LSAVE-NCOLD;

    /* Check if we need to allocate more memory, and allocate if necessary */
#if 0  /* Proposal by Michael A. Saunders (logic based on Markowitz' rule) */
    L = NROWD*NCOLD;
    
    /* Try to avoid future expansions by anticipating further updates - KE extension */
    if(LUSOL->luparm[LUSOL_IP_UPDATELIMIT] > 0)
#if 1
      L *= (int) (log(LUSOL->luparm[LUSOL_IP_UPDATELIMIT]-LUSOL->luparm[LUSOL_IP_UPDATECOUNT]+2.0) + 1);
#else
      L *= (LUSOL->luparm[LUSOL_IP_UPDATELIMIT]-LUSOL->luparm[LUSOL_IP_UPDATECOUNT]) / 2 + 1;
#endif

#else  /* Version by Kjell Eikland (from luparm[LUSOL_IP_MINIMUMLENA] and safety margin) */
    L  = MAX(LROW, LCOL) + 2*(LUSOL->m+LUSOL->n);
    L *= LUSOL_MULT_nz_a;
    L = MAX(L, NROWD*NCOLD);
#endif

    /* Do the memory expansion */  
    if((L > LFREE-LCOL) && LUSOL_expand_a(LUSOL, &L, &LFREE)) {
      LL1   += L;
      LU1   += L;
      LSAVE += L;
#ifdef ClassicdiagU
      LUSOL->diagU += L;
#endif
#ifdef ClassicHamaxR
      HA    += L;
      HJ    += L;
      HK    += L; 
      AMAXR += L;
#endif
    }
    LIMIT = (int) (USPACE*LFILE)+LUSOL->m+LUSOL->n+1000;

/*         Make sure the column file has room.
           Also force a compression if its length exceeds a certain limit. */
#ifdef StaticMemAlloc
    MINFRE = NCOLD+MELIM;
#else
    MINFRE = NROWD*NCOLD;
#endif
    NFREE = LFREE-LCOL;
    if(NFREE<MINFRE || LCOL>LIMIT) {
      LU1REC(LUSOL, LUSOL->n,TRUE,&LCOL,
                    LUSOL->indc,LUSOL->lenc,LUSOL->locc);
      LFILE = LCOL;
      JLAST = LUSOL->indc[LCOL+1];
      NFREE = LFREE-LCOL;
      if(NFREE<MINFRE)
        goto x970;
    }
/*         Make sure the row file has room. */
#ifdef StaticMemAlloc
    MINFRE = NCOLD+MELIM;
#else
    MINFRE = NROWD*NCOLD;
#endif
    NFREE = LFREE-LROW;
    if(NFREE<MINFRE || LROW>LIMIT) {
      LU1REC(LUSOL, LUSOL->m,FALSE,&LROW,
                    LUSOL->indr,LUSOL->lenr,LUSOL->locr);
      LFILE = LROW;
      ILAST = LUSOL->indr[LROW+1];
      NFREE = LFREE-LROW;
      if(NFREE<MINFRE)
        goto x970;
    }
/*         ===============================================================
           Move the pivot element to the front of its row
           and to the top of its column.
           =============================================================== */
    LPIVR = LUSOL->locr[IBEST];
    LPIVR1 = LPIVR+1;
    LPIVR2 = LPIVR+NELIM;
    for(L = LPIVR; L <= LPIVR2; L++) {
      if(LUSOL->indr[L]==JBEST)
        break;
    }

    LUSOL->indr[L] = LUSOL->indr[LPIVR];
    LUSOL->indr[LPIVR] = JBEST;
    LPIVC = LUSOL->locc[JBEST];
    LPIVC1 = LPIVC+1;
    LPIVC2 = LPIVC+MELIM;
    for(L = LPIVC; L <= LPIVC2; L++) {
      if(LUSOL->indc[L]==IBEST)
        break;
    }
    LUSOL->indc[L] = LUSOL->indc[LPIVC];
    LUSOL->indc[LPIVC] = IBEST;
    ABEST = LUSOL->a[L];
    LUSOL->a[L] = LUSOL->a[LPIVC];
    LUSOL->a[LPIVC] = ABEST;
    if(!KEEPLU)
/*            Store just the diagonal of U, in natural order.
   !!         a[ldiagU + nrowu] = abest ! This was in pivot order. */
      LUSOL->diagU[JBEST] = ABEST;

/*     ==============================================================
        Delete pivot col from heap.
        Hk tells us where it is in the heap.
       ============================================================== */
    if(TCP) {
#ifdef ClassicHamaxR
      KBEST = HK[JBEST];
      HDELETE(HA,HJ,HK,&HLEN,KBEST,&H);
#else
      KBEST = LUSOL->Hk[JBEST];
      HDELETE(LUSOL->Ha,LUSOL->Hj,LUSOL->Hk,&HLEN,KBEST,&H);
#endif
      HOPS += H;
    }
/*         ===============================================================
           Delete the pivot row from the column file
           and store it as the next row of  U.
           set  indr(lu) = 0     to initialize jfill ptrs on columns of D,
                indc(lu) = lenj  to save the original column lengths.
           =============================================================== */
    LUSOL->a[LU1] = ABEST;
    LUSOL->indr[LU1] = JBEST;
    LUSOL->indc[LU1] = NROWD;
    LU = LU1;
    DIAG = fabs(ABEST);
    *UMAX = MAX(*UMAX,DIAG);
    *DUMAX = MAX(*DUMAX,DIAG);
    *DUMIN = MIN(*DUMIN,DIAG);
    for(LR = LPIVR1; LR <= LPIVR2; LR++) {
      LU++;
      J = LUSOL->indr[LR];
      LENJ = LUSOL->lenc[J];
      LUSOL->lenc[J] = LENJ-1;
      LC1 = LUSOL->locc[J];
      LAST = LC1+LUSOL->lenc[J];
      for(L = LC1; L <= LAST; L++) {
        if(LUSOL->indc[L]==IBEST)
          break;
      }
      LUSOL->a[LU] = LUSOL->a[L];
      LUSOL->indr[LU] = 0;
      LUSOL->indc[LU] = LENJ;
      *UMAX = MAX(*UMAX,fabs(LUSOL->a[LU]));
      LUSOL->a[L] = LUSOL->a[LAST];
      LUSOL->indc[L] = LUSOL->indc[LAST];
/*      Free entry */
      LUSOL->indc[LAST] = 0;
/* ???        if (j .eq. jlast) lcol = lcol - 1 */
    }
/*         ===============================================================
           Delete the pivot column from the row file
           and store the nonzeros of the next column of  L.
           Set  indc(ll) = 0     to initialize markl(*) markers,
                indr(ll) = 0     to initialize ifill(*) row fill-in cntrs,
                indc(ls) = leni  to save the original row lengths,
                indr(ls) = iqloc(i)    to save parts of  iqloc(*),
                iqloc(i) = lsave - ls  to point to the nonzeros of  L
                         = -1, -2, -3, ... in mark(*).
           =============================================================== */
    LUSOL->indc[LSAVE] = NCOLD;
    if(MELIM==0)
      goto x700;
    LL = LL1-1;
    LS = LSAVE;
    ABEST = ONE/ABEST;
    for(LC = LPIVC1; LC <= LPIVC2; LC++) {
      LL++;
      LS++;
      I = LUSOL->indc[LC];
      LENI = LUSOL->lenr[I];
      LUSOL->lenr[I] = LENI-1;
      LR1 = LUSOL->locr[I];
      LAST = LR1+LUSOL->lenr[I];
      for(L = LR1; L <= LAST; L++) {
        if(LUSOL->indr[L]==JBEST)
          break;
      }
      LUSOL->indr[L] = LUSOL->indr[LAST];
/*      Free entry */
      LUSOL->indr[LAST] = 0;
      LUSOL->a[LL] = -LUSOL->a[LC]*ABEST;
      LIJ = fabs(LUSOL->a[LL]);
      *LMAX = MAX(*LMAX,LIJ);
      LUSOL->indc[LL] = 0;
      LUSOL->indr[LL] = 0;
      LUSOL->indc[LS] = LENI;
      LUSOL->indr[LS] = LUSOL->iqloc[I];
      LUSOL->iqloc[I] = LSAVE-LS;
    }
/*         ===============================================================
           Do the Gaussian elimination.
           This involves adding a multiple of the pivot column
           to all other columns in the pivot row.
           Sometimes more than one call to lu1gau is needed to allow
           compression of the column file.
           lfirst  says which column the elimination should start with.
           minfre  is a bound on the storage needed for any one column.
           lu      points to off-diagonals of u.
           nfill   keeps track of pending fill-in in the row file.
           =============================================================== */
    if(NELIM==0)
      goto x700;
    LFIRST = LPIVR1;
    MINFRE = MLEFT+NSPARE;
    LU = 1;
    NFILL = 0;

x400:
#ifdef UseTimer
    timer ( "start", eltime );
#endif
    LU1GAU(LUSOL, MELIM,NSPARE,SMALL,LPIVC1,LPIVC2,&LFIRST,LPIVR2,
           LFREE,MINFRE,ILAST,&JLAST,&LROW,&LCOL,&LU,&NFILL,
           LUSOL->iqloc, LUSOL->a+LL1-LUSOL_ARRAYOFFSET,
           LUSOL->indc+LL1-LUSOL_ARRAYOFFSET, LUSOL->a+LU1-LUSOL_ARRAYOFFSET,
           LUSOL->indr+LL1-LUSOL_ARRAYOFFSET, LUSOL->indr+LU1-LUSOL_ARRAYOFFSET);
#ifdef UseTimer
    timer ( "finish", eltime );
#endif
    if(LFIRST>0) {
/*            The elimination was interrupted.
              Compress the column file and try again.
              lfirst, lu and nfill have appropriate new values. */
      LU1REC(LUSOL, LUSOL->n,TRUE,&LCOL,
                    LUSOL->indc,LUSOL->lenc,LUSOL->locc);
      LFILE = LCOL;
      JLAST = LUSOL->indc[LCOL+1];
      LPIVC = LUSOL->locc[JBEST];
      LPIVC1 = LPIVC+1;
      LPIVC2 = LPIVC+MELIM;
      NFREE = LFREE-LCOL;
      if(NFREE<MINFRE) {
          goto x970;
      }
      goto x400;
    }
/*         ===============================================================
           The column file has been fully updated.
           Deal with any pending fill-in in the row file.
           =============================================================== */
    if(NFILL>0) {
/*            Compress the row file if necessary.
              lu1gau has set nfill to be the number of pending fill-ins
              plus the current length of any rows that need to be moved. */
      MINFRE = NFILL;
      NFREE = LFREE-LROW;
      if(NFREE<MINFRE) {
        LU1REC(LUSOL, LUSOL->m,FALSE,&LROW,
                      LUSOL->indr,LUSOL->lenr,LUSOL->locr);
        LFILE = LROW;
        ILAST = LUSOL->indr[LROW+1];
        LPIVR = LUSOL->locr[IBEST];
        LPIVR1 = LPIVR+1;
        LPIVR2 = LPIVR+NELIM;
        NFREE = LFREE-LROW;
        if(NFREE<MINFRE) {
            goto x970;
        }
      }
/*            Move rows that have pending fill-in to end of the row file.
              Then insert the fill-in. */
      LU1PEN(LUSOL, NSPARE,&ILAST,
             LPIVC1,LPIVC2,LPIVR1,LPIVR2,
             &LROW,LUSOL->indr+LL1-LUSOL_ARRAYOFFSET,LUSOL->indr+LU1-LUSOL_ARRAYOFFSET);
    }
/*         ===============================================================
           Restore the saved values of  iqloc.
           Insert the correct indices for the col of L and the row of U.
           =============================================================== */
x700:
    LUSOL->lenr[IBEST] = 0;
    LUSOL->lenc[JBEST] = 0;
    LL = LL1-1;
    LS = LSAVE;
    for(LC = LPIVC1; LC <= LPIVC2; LC++) {
      LL++;
      LS++;
      I = LUSOL->indc[LC];
      LUSOL->iqloc[I] = LUSOL->indr[LS];
      LUSOL->indc[LL] = I;
      LUSOL->indr[LL] = IBEST;
    }
    LU = LU1-1;
    for(LR = LPIVR; LR <= LPIVR2; LR++) {
      LU++;
      LUSOL->indr[LU] = LUSOL->indr[LR];
    }
/*         ===============================================================
           Free the space occupied by the pivot row
           and update the column permutation.
           Then free the space occupied by the pivot column
           and update the row permutation.
           nzchng is found in both calls to lu1pq2, but we use it only
           after the second.
           =============================================================== */
    LU1PQ2(LUSOL, NCOLD, &NZCHNG,
           LUSOL->indr+LPIVR-LUSOL_ARRAYOFFSET,
           LUSOL->indc+LU1-LUSOL_ARRAYOFFSET, LUSOL->lenc,
           LUSOL->iqloc, LUSOL->iq, LUSOL->iqinv);
    LU1PQ2(LUSOL, NROWD, &NZCHNG,
           LUSOL->indc+LPIVC-LUSOL_ARRAYOFFSET,
           LUSOL->indc+LSAVE-LUSOL_ARRAYOFFSET, LUSOL->lenr,
           LUSOL->iploc, LUSOL->ip, LUSOL->ipinv);
    NZLEFT += NZCHNG;

/*         ===============================================================
           lu1mxr resets Amaxr(i) in each modified row i.
           lu1mxc moves the largest aij to the top of each modified col j.
           28 Jun 2002: Note that cols of L have an implicit diag of 1.0,
                        so lu1mxr is called with ll1, not ll1+1, whereas
                           lu1mxc is called with          lu1+1.
           =============================================================== */
    if(UTRI && TPP) {
/*      Relax -- we're not keeping big elements at the top yet. */
    }
    else {
      if(TRP && MELIM>0)
#ifdef ClassicHamaxR
        LU1MXR(LUSOL, LL1,LL,LUSOL->indc,AMAXR);
#else
        LU1MXR(LUSOL, LL1,LL,LUSOL->indc,LUSOL->amaxr);
#endif

      if(NELIM>0) {
        LU1MXC(LUSOL, LU1+1,LU,LUSOL->indr);
/*      Update modified columns in heap */
        if(TCP) {
          for(KK = LU1+1; KK <= LU; KK++) {
            J = LUSOL->indr[KK];
#ifdef ClassicHamaxR
            K = HK[J];
#else
            K = LUSOL->Hk[J];
#endif
/*      Biggest aij in column j */
            V = fabs(LUSOL->a[LUSOL->locc[J]]);
#ifdef ClassicHamaxR
            HCHANGE(HA,HJ,HK,HLEN,K,V,J,&H);
#else
            HCHANGE(LUSOL->Ha,LUSOL->Hj,LUSOL->Hk,HLEN,K,V,J,&H);
#endif
            HOPS += H;
          }
        }
      }
    }
/*         ===============================================================
           Negate lengths of pivot row and column so they will be
           eliminated during compressions.
           =============================================================== */
    LUSOL->lenr[IBEST] = -NCOLD;
    LUSOL->lenc[JBEST] = -NROWD;

/*         Test for fatal bug: row or column lists overwriting L and U. */
    if(LROW>LSAVE || LCOL>LSAVE)
      goto x980;

/*         Reset the file lengths if pivot row or col was at the end. */
    if(IBEST==ILAST)
      LROW = LUSOL->locr[IBEST];

    if(JBEST==JLAST)
      LCOL = LUSOL->locc[JBEST];

  }
/*      ------------------------------------------------------------------
        End of main loop.
        ------------------------------------------------------------------
        ------------------------------------------------------------------
        Normal exit.
        Move empty rows and cols to the end of ip, iq.
        Then finish with a dense LU if necessary.
        ------------------------------------------------------------------ */
x900:
  *INFORM = LUSOL_INFORM_LUSUCCESS;
  LU1PQ3(LUSOL, LUSOL->m,LUSOL->lenr,LUSOL->ip,LUSOL->ipinv,&MRANK);
  LU1PQ3(LUSOL, LUSOL->n,LUSOL->lenc,LUSOL->iq,LUSOL->iqinv,NRANK);
  *NRANK = MIN(MRANK,*NRANK);
  if(DENSLU) {
#ifdef UseTimer
    timer ( "start", 17 );
#endif
    LU1FUL(LUSOL, LEND,LU1,TPP,MLEFT,NLEFT,*NRANK,NROWU,LENL,LENU,
           &NSING,KEEPLU,SMALL,LUSOL->a+LD-LUSOL_ARRAYOFFSET,LUSOL->locr);
/* ***     21 Dec 1994: Bug in next line.
   ***     nrank  = nrank - nsing */
    *NRANK = MINMN-NSING;
#ifdef UseTimer
    timer ( "finish", 17 );
#endif
  }
  *MINLEN = (*LENL)+(*LENU)+2*(LUSOL->m+LUSOL->n);
  goto x990;
/*      Not enough space free after a compress.
        Set  minlen  to an estimate of the necessary value of  lena. */
x970:
  *INFORM = LUSOL_INFORM_ANEEDMEM;
  *MINLEN = LENA2+LFILE+2*(LUSOL->m+LUSOL->n);
  goto x990;
/*      Fatal error.  This will never happen!
       (Famous last words.) */
x980:
  *INFORM = LUSOL_INFORM_FATALERR;
  goto x990;
/*      Fatal error with TSP.  Diagonal pivot not found. */
x985:
  *INFORM = LUSOL_INFORM_NOPIVOT;
/*      Exit. */
x990:
#ifdef UseTimer
  timer ( "finish", 3 );
#endif
;
}


/* ==================================================================
   lu1fac computes a factorization A = L*U, where A is a sparse
   matrix with m rows and n columns, P*L*P' is lower triangular
   and P*U*Q is upper triangular for certain permutations P, Q
   (which are returned in the arrays ip, iq).
   Stability is ensured by limiting the size of the elements of L.
   The nonzeros of A are input via the parallel arrays a, indc, indr,
   which should contain nelem entries of the form    aij,    i,    j
   in any order.  There should be no duplicate pairs         i,    j.

   ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   +        Beware !!!   The row indices i must be in indc,         +
   +              and the column indices j must be in indr.         +
   +              (Not the other way round!)                        +
   ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   It does not matter if some of the entries in a(*) are zero.
   Entries satisfying  abs( a(i) ) .le. parmlu(3)  are ignored.
   Other parameters in luparm and parmlu are described below.
   The matrix A may be singular.  On exit, nsing = luparm(11) gives
   the number of apparent singularities.  This is the number of
   "small" diagonals of the permuted factor U, as judged by
   the input tolerances Utol1 = parmlu(4) and  Utol2 = parmlu(5).
   The diagonal element diagj associated with column j of A is
   "small" if
               abs( diagj ) .le. Utol1
   or
               abs( diagj ) .le. Utol2 * max( uj ),
   where max( uj ) is the maximum element in the j-th column of U.
   The position of such elements is returned in w(*).  In general,
   w(j) = + max( uj ),  but if column j is a singularity,
   w(j) = - max( uj ).  Thus, w(j) .le. 0 if column j appears to be
   dependent on the other columns of A.
   NOTE: lu1fac (like certain other sparse LU packages) does not
   treat dense columns efficiently.  This means it will be slow
   on "arrow matrices" of the form
                A = (x       a)
                    (  x     b)
                    (    x   c)
                    (      x d)
                    (x x x x e)
   if the numerical values in the dense column allow it to be
   chosen LATE in the pivot order.
   With TPP (Threshold Partial Pivoting), the dense column is
   likely to be chosen late.
   With TCP (Threshold Complete Pivoting), if any of a,b,c,d
   is significantly larger than other elements of A, it will
   be chosen as the first pivot and the dense column will be
   eliminated, giving reasonably sparse factors.
   However, if element e is so big that TCP chooses it, the factors
   will become dense.  (It's hard to win on these examples!)
   ------------------------------------------------------------------

   Notes on the array names
   ------------------------
   During the LU factorization, the sparsity pattern of the matrix
   being factored is stored twice: in a column list and a row list.
   The column list is ( a, indc, locc, lenc )
   where
         a(*)    holds the nonzeros,
         indc(*) holds the indices for the column list,
         locc(j) points to the start of column j in a(*) and indc(*),
         lenc(j) is the number of nonzeros in column j.
   The row list is    (    indr, locr, lenr )
   where
         indr(*) holds the indices for the row list,
         locr(i) points to the start of row i in indr(*),
         lenr(i) is the number of nonzeros in row i.
   At all stages of the LU factorization, ip contains a complete
   row permutation.  At the start of stage k,  ip(1), ..., ip(k-1)
   are the first k-1 rows of the final row permutation P.
   The remaining rows are stored in an ordered list
                        ( ip, iploc, ipinv )
   where
         iploc(nz) points to the start in ip(*) of the set of rows
                   that currently contain nz nonzeros,
         ipinv(i)  points to the position of row i in ip(*).
   For example,
         iploc(1) = k   (and this is where rows of length 1 {),
         iploc(2) = k+p  if there are p rows of length 1
                        (and this is where rows of length 2 {).
   Similarly for iq, iqloc, iqinv.
   ---------------------------------------------------------------------
   INPUT PARAMETERS
   m      (not altered) is the number of rows in A.
   n      (not altered) is the number of columns in A.
   nelem  (not altered) is the number of matrix entries given in
          the arrays a, indc, indr.
   lena   (not altered) is the dimension of  a, indc, indr.
          This should be significantly larger than nelem.
          Typically one should have
             lena > max( 2*nelem, 10*m, 10*n, 10000 )
          but some applications may need more.
          On machines with virtual memory it is safe to have
          lena "far bigger than necessary", since not all of the
          arrays will be used.
   a      (overwritten) contains entries   Aij  in   a(1:nelem).
   indc   (overwritten) contains the indices i in indc(1:nelem).
   indr   (overwritten) contains the indices j in indr(1:nelem).
   luparm input parameters:                                Typical value
   luparm( 1) = nout     File number for printed messages.         6
   luparm( 2) = lprint   Print level.                              0
                    <  0 suppresses output.
                    =  0 gives error messages.
                   >= 10 gives statistics about the LU factors.
                   >= 50 gives debug output from lu1fac
                         (the pivot row and column and the
                         no. of rows and columns involved at
                         each elimination step).
   luparm( 3) = maxcol   lu1fac: maximum number of columns         5
                         searched allowed in a Markowitz-type
                         search for the next pivot element.
                         For some of the factorization, the
                         number of rows searched is
                         maxrow = maxcol - 1.
   luparm( 6) = 0    =>  TPP: Threshold Partial   Pivoting.        0
              = 1    =>  TRP: Threshold Rook      Pivoting.
              = 2    =>  TCP: Threshold Complete  Pivoting.
              = 3    =>  TSP: Threshold Symmetric Pivoting.
              = 4    =>  TDP: Threshold Diagonal  Pivoting.
                              (TDP not yet implemented).
                         TRP and TCP are more expensive than TPP but
                         more stable and better at revealing rank.
                         Take care with setting parmlu(1), especially
                         with TCP.
                         NOTE: TSP and TDP are for symmetric matrices
                         that are either definite or quasi-definite.
                         TSP is effectively TRP for symmetric matrices.
                         TDP is effectively TCP for symmetric matrices.
   luparm( 8) = keepLU   lu1fac: keepLU = 1 means the numerical    1
                         factors will be computed if possible.
                         keepLU = 0 means L and U will be discarded
                         but other information such as the row and
                         column permutations will be returned.
                         The latter option requires less storage.
   parmlu input parameters:                                Typical value
   parmlu( 1) = Ltol1    Max Lij allowed during Factor.
                                                   TPP     10.0 or 100.0
                                                   TRP      4.0 or  10.0
                                                   TCP      5.0 or  10.0
                                                   TSP      4.0 or  10.0
                         With TRP and TCP (Rook and Complete Pivoting),
                         values less than 25.0 may be expensive
                         on badly scaled data.  However,
                         values less than 10.0 may be needed
                         to obtain a reliable rank-revealing
                         factorization.
   parmlu( 2) = Ltol2    Max Lij allowed during Updates.            10.0
                         during updates.
   parmlu( 3) = small    Absolute tolerance for       eps**0.8 = 3.0d-13
                         treating reals as zero.
   parmlu( 4) = Utol1    Absolute tol for flagging    eps**0.67= 3.7d-11
                         small diagonals of U.
   parmlu( 5) = Utol2    Relative tol for flagging    eps**0.67= 3.7d-11
                         small diagonals of U.
                         (eps = machine precision)
   parmlu( 6) = Uspace   Factor limiting waste space in  U.      3.0
                         In lu1fac, the row or column lists
                         are compressed if their length
                         exceeds Uspace times the length of
                         either file after the last compression.
   parmlu( 7) = dens1    The density at which the Markowitz      0.3
                         pivot strategy should search maxcol
                         columns and no rows.
                         (Use 0.3 unless you are experimenting
                         with the pivot strategy.)
   parmlu( 8) = dens2    the density at which the Markowitz      0.5
                         strategy should search only 1 column,
                         or (if storage is available)
                         the density at which all remaining
                         rows and columns will be processed
                         by a dense LU code.
                         For example, if dens2 = 0.1 and lena is
                         large enough, a dense LU will be used
                         once more than 10 per cent of the
                         remaining matrix is nonzero.

   OUTPUT PARAMETERS
   a, indc, indr     contain the nonzero entries in the LU factors of A.
          If keepLU = 1, they are in a form suitable for use
          by other parts of the LUSOL package, such as lu6sol.
          U is stored by rows at the start of a, indr.
          L is stored by cols at the end   of a, indc.
          If keepLU = 0, only the diagonals of U are stored, at the
          end of a.
   ip, iq    are the row and column permutations defining the
          pivot order.  For example, row ip(1) and column iq(1)
          defines the first diagonal of U.
   lenc(1:numl0) contains the number of entries in nontrivial
          columns of L (in pivot order).
   lenr(1:m) contains the number of entries in each row of U
          (in original order).
   locc(1:n) = 0 (ready for the LU update routines).
   locr(1:m) points to the beginning of the rows of U in a, indr.
   iploc, iqloc, ipinv, iqinv  are undefined.
   w      indicates singularity as described above.
   inform = 0 if the LU factors were obtained successfully.
          = 1 if U appears to be singular, as judged by lu6chk.
          = 3 if some index pair indc(l), indr(l) lies outside
              the matrix dimensions 1:m , 1:n.
          = 4 if some index pair indc(l), indr(l) duplicates
              another such pair.
          = 7 if the arrays a, indc, indr were not large enough.
              Their length "lena" should be increase to at least
              the value "minlen" given in luparm(13).
          = 8 if there was some other fatal error.  (Shouldn't happen!)
          = 9 if no diagonal pivot could be found with TSP or TDP.
              The matrix must not be sufficiently definite
              or quasi-definite.
   luparm output parameters:
   luparm(10) = inform   Return code from last call to any LU routine.
   luparm(11) = nsing    No. of singularities marked in the
                         output array w(*).
   luparm(12) = jsing    Column index of last singularity.
   luparm(13) = minlen   Minimum recommended value for  lena.
   luparm(14) = maxlen   ?
   luparm(15) = nupdat   No. of updates performed by the lu8 routines.
   luparm(16) = nrank    No. of nonempty rows of U.
   luparm(17) = ndens1   No. of columns remaining when the density of
                         the matrix being factorized reached dens1.
   luparm(18) = ndens2   No. of columns remaining when the density of
                         the matrix being factorized reached dens2.
   luparm(19) = jumin    The column index associated with DUmin.
   luparm(20) = numL0    No. of columns in initial  L.
   luparm(21) = lenL0    Size of initial  L  (no. of nonzeros).
   luparm(22) = lenU0    Size of initial  U.
   luparm(23) = lenL     Size of current  L.
   luparm(24) = lenU     Size of current  U.
   luparm(25) = lrow     Length of row file.
   luparm(26) = ncp      No. of compressions of LU data structures.
   luparm(27) = mersum   lu1fac: sum of Markowitz merit counts.
   luparm(28) = nUtri    lu1fac: triangular rows in U.
   luparm(29) = nLtri    lu1fac: triangular rows in L.
   luparm(30) =
   parmlu output parameters:
   parmlu(10) = Amax     Maximum element in  A.
   parmlu(11) = Lmax     Maximum multiplier in current  L.
   parmlu(12) = Umax     Maximum element in current  U.
   parmlu(13) = DUmax    Maximum diagonal in  U.
   parmlu(14) = DUmin    Minimum diagonal in  U.
   parmlu(15) = Akmax    Maximum element generated at any stage
                         during TCP factorization.
   parmlu(16) = growth   TPP: Umax/Amax    TRP, TCP, TSP: Akmax/Amax
   parmlu(17) =
   parmlu(18) =
   parmlu(19) =
   parmlu(20) = resid    lu6sol: residual after solve with U or U'.
   ...
   parmlu(30) =
   ------------------------------------------------------------------
   00 Jun 1983  Original version.
   00 Jul 1987  nrank  saved in luparm(16).
   12 Apr 1989  ipinv, iqinv added as workspace.
   26 Apr 1989  maxtie replaced by maxcol in Markowitz search.
   16 Mar 1992  jumin  saved in luparm(19).
   10 Jun 1992  lu1fad has to move empty rows and cols to the bottom
                (via lu1pq3) before doing the dense LU.
   12 Jun 1992  Deleted dense LU (lu1ful, lu1vlu).
   25 Oct 1993  keepLU implemented.
   07 Feb 1994  Added new dense LU (lu1ful, lu1den).
   21 Dec 1994  Bugs fixed in lu1fad (nrank) and lu1ful (ipvt).
   08 Aug 1995  Use ip instead of w as parameter to lu1or3 (for F90).
   13 Sep 2000  TPP and TCP options implemented.
   17 Oct 2000  Fixed troubles due to A = empty matrix (Todd Munson).
   01 Dec 2000  Save Lmax, Umax, etc. after both lu1fad and lu6chk.
                lu1fad sets them when keepLU = false.
                lu6chk sets them otherwise, and includes items
                from the dense LU.
   11 Mar 2001  lu6chk now looks at diag(U) when keepLU = false.
   26 Apr 2002  New TCP implementation using heap routines to
                store largest element in each column.
                New workspace arrays Ha, Hj, Hk required.
                For compatibility, borrow space from a, indc, indr
                rather than adding new input parameters.
   01 May 2002  lu1den changed to lu1DPP (dense partial  pivoting).
                lu1DCP implemented       (dense complete pivoting).
                Both TPP and TCP now switch to dense mode and end.
   ================================================================== */
void LU1FAC(LUSOLrec *LUSOL, int *INFORM)
{
  MYBOOL  KEEPLU, TCP, TPP, TRP, TSP;
  int     LPIV, NELEM0, LPRINT, MINLEN, NUML0, LENL, LENU, LROW, MERSUM,
          NUTRI, NLTRI, NDENS1, NDENS2, NRANK, NSING, JSING, JUMIN, NUMNZ, LERR,
          LU, LL, LM, LTOPL, K, I, LENUK, J, LENLK, IDUMMY, LLSAVE, NMOVE, L2, L, NCP, NBUMP;
#ifdef ClassicHamaxR
  int     LENH, LENA2, LOCH, LMAXR;
#endif

  REAL    LMAX, LTOL, SMALL, AMAX, UMAX, DUMAX, DUMIN, AKMAX, DM, DN, DELEM, DENSTY,
          AGRWTH, UGRWTH, GROWTH, CONDU, DINCR, AVGMER;

/*      Free row-based version of L0 (regenerated by LUSOL_btran). */
  if(LUSOL->L0 != NULL)
    LUSOL_matfree(&(LUSOL->L0));

/*      Grab relevant input parameters. */
  NELEM0 = LUSOL->nelem;
  LPRINT = LUSOL->luparm[LUSOL_IP_PRINTLEVEL];
  LPIV   = LUSOL->luparm[LUSOL_IP_PIVOTTYPE];
  KEEPLU = (MYBOOL) (LUSOL->luparm[LUSOL_IP_KEEPLU]!=FALSE);
/*      Limit on size of Lij */
  LTOL   = LUSOL->parmlu[LUSOL_RP_FACTORMAX_Lij];
/*      Drop tolerance */
  SMALL  = LUSOL->parmlu[LUSOL_RP_ZEROTOLERANCE];
  TPP = (MYBOOL) (LPIV==LUSOL_PIVMOD_TPP);
  TRP = (MYBOOL) (LPIV==LUSOL_PIVMOD_TRP);
  TCP = (MYBOOL) (LPIV==LUSOL_PIVMOD_TCP);
  TSP = (MYBOOL) (LPIV==LUSOL_PIVMOD_TSP);
/*      Initialize output parameters. */
  *INFORM = LUSOL_INFORM_LUSUCCESS;
  LERR   = 0;
  MINLEN = LUSOL->nelem + 2*(LUSOL->m+LUSOL->n);
  NUML0  = 0;
  LENL   = 0;
  LENU   = 0;
  LROW   = 0;
  MERSUM = 0;
  NUTRI  = LUSOL->m;
  NLTRI  = 0;
  NDENS1 = 0;
  NDENS2 = 0;
  NRANK  = 0;
  NSING  = 0;
  JSING  = 0;
  JUMIN  = 0;
  AMAX   = ZERO;
  LMAX   = ZERO;
  UMAX   = ZERO;
  DUMAX  = ZERO;
  DUMIN  = ZERO;
  AKMAX  = ZERO;

/*      Float version of dimensions. */
  DM = LUSOL->m;
  DN = LUSOL->n;
  DELEM = LUSOL->nelem;

/*      Initialize workspace parameters. */
  LUSOL->luparm[LUSOL_IP_COMPRESSIONS_LU] = 0;
  if(LUSOL->lena < MINLEN) {
    if(!LUSOL_realloc_a(LUSOL, MINLEN))
      goto x970;
  }

/*      ------------------------------------------------------------------
        Organize the  aij's  in  a, indc, indr.
        lu1or1  deletes small entries, tests for illegal  i,j's,
                and counts the nonzeros in each row and column.
        lu1or2  reorders the elements of  A  by columns.
        lu1or3  uses the column list to test for duplicate entries
                (same indices  i,j).
        lu1or4  constructs a row list from the column list.
        ------------------------------------------------------------------ */
  LU1OR1(LUSOL, SMALL,&AMAX,&NUMNZ,&LERR,INFORM);
  if(LPRINT>=LUSOL_MSG_STATISTICS) {
    DENSTY = (100*DELEM)/(DM*DN);
    LUSOL_report(LUSOL, 0, "m:%6d %c n:%6d  nzcount:%9d  Amax:%g  Density:%g\n",
                           LUSOL->m, relationChar(LUSOL->m, LUSOL->n), LUSOL->n, 
                           LUSOL->nelem, AMAX, DENSTY);
  }
  if(*INFORM!=LUSOL_INFORM_LUSUCCESS)
    goto x930;
  LUSOL->nelem = NUMNZ;
  LU1OR2(LUSOL);
  LU1OR3(LUSOL, &LERR,INFORM);
  if(*INFORM!=LUSOL_INFORM_LUSUCCESS)
    goto x940;
  LU1OR4(LUSOL);
/*      ------------------------------------------------------------------
        Set up lists of rows and columns with equal numbers of nonzeros,
        using  indc(*)  as workspace.
        ------------------------------------------------------------------ */
  LU1PQ1(LUSOL, LUSOL->m,LUSOL->n,LUSOL->lenr,
         LUSOL->ip,LUSOL->iploc,LUSOL->ipinv,
         LUSOL->indc+LUSOL->nelem); /* LUSOL_ARRAYOFFSET implied */
  LU1PQ1(LUSOL, LUSOL->n,LUSOL->m,LUSOL->lenc,
         LUSOL->iq,LUSOL->iqloc,LUSOL->iqinv,
         LUSOL->indc+LUSOL->nelem); /* LUSOL_ARRAYOFFSET implied */
/*      ------------------------------------------------------------------
        For TCP, Ha, Hj, Hk are allocated separately, similarly amaxr
        for TRP. Then compute the factorization  A = L*U.
        ------------------------------------------------------------------ */
#ifdef ClassicHamaxR
  if(TPP || TSP) {
    LENH  = 1;
    LENA2 = LUSOL->lena;
    LOCH  = LUSOL->lena;
    LMAXR = 1;
  }
  else if(TRP) {
    LENH  = 1;                     /* Dummy                                */
    LENA2 = LUSOL->lena-LUSOL->m;  /* Reduced length of      a             */
    LOCH  = LUSOL->lena;           /* Dummy                                */
    LMAXR = LENA2+1;               /* Start of Amaxr      in a             */
  }
  else if(TCP) {
    LENH  = LUSOL->n;              /* Length of heap                       */
    LENA2 = LUSOL->lena-LENH;      /* Reduced length of      a, indc, indr */
    LOCH  = LENA2+1;               /* Start of Ha, Hj, Hk in a, indc, indr */
    LMAXR = 1;                     /* Dummy                                */
  }
  LU1FAD(LUSOL, 
         LENA2,LENH,
         LUSOL->a+LOCH-LUSOL_ARRAYOFFSET,
         LUSOL->indc+LOCH-LUSOL_ARRAYOFFSET,
         LUSOL->indr+LOCH-LUSOL_ARRAYOFFSET,
         LUSOL->a+LMAXR-LUSOL_ARRAYOFFSET,
         INFORM,&LENL,&LENU,
         &MINLEN,&MERSUM,&NUTRI,&NLTRI,&NDENS1,&NDENS2,
         &NRANK,&LMAX,&UMAX,&DUMAX,&DUMIN,&AKMAX);
#else
  LU1FAD(LUSOL,
         INFORM,&LENL,&LENU,
         &MINLEN,&MERSUM,&NUTRI,&NLTRI,&NDENS1,&NDENS2,
         &NRANK,&LMAX,&UMAX,&DUMAX,&DUMIN,&AKMAX);
#endif
  LUSOL->luparm[LUSOL_IP_RANK_U]     = NRANK;
  LUSOL->luparm[LUSOL_IP_NONZEROS_L] = LENL;
  if(*INFORM==LUSOL_INFORM_ANEEDMEM)
    goto x970;
  if(*INFORM==LUSOL_INFORM_NOPIVOT)
    goto x985;
  if(*INFORM>LUSOL_INFORM_LUSUCCESS)
    goto x980;
  if(KEEPLU) {
/*         ---------------------------------------------------------------
           The LU factors are at the top of  a, indc, indr,
           with the columns of  L  and the rows of  U  in the order
           ( free )   ... ( u3 ) ( l3 ) ( u2 ) ( l2 ) ( u1 ) ( l1 ).
           Starting with ( l1 ) and ( u1 ), move the rows of  U  to the
           left and the columns of  L  to the right, giving
           ( u1 ) ( u2 ) ( u3 ) ...   ( free )   ... ( l3 ) ( l2 ) ( l1 ).
           Also, set  numl0 = the number of nonempty columns of  U.
           --------------------------------------------------------------- */
    LU = 0;
    LL = LUSOL->lena+1;
#ifdef ClassicHamaxR
    LM = LENA2+1;
#else
    LM = LL;
#endif
    LTOPL = LL-LENL-LENU;
    LROW = LENU;
    for(K = 1; K <= NRANK; K++) {
      I = LUSOL->ip[K];
      LENUK = -LUSOL->lenr[I];
      LUSOL->lenr[I] = LENUK;
      J = LUSOL->iq[K];
      LENLK = -LUSOL->lenc[J]-1;
      if(LENLK>0) {
        NUML0++;
        LUSOL->iqloc[NUML0] = LENLK;
      }
      if(LU+LENUK<LTOPL) {
/*               =========================================================
                 There is room to move ( uk ).  Just right-shift ( lk ).
                 ========================================================= */
        for(IDUMMY = 1; IDUMMY <= LENLK; IDUMMY++) {
          LL--;
          LM--;
          LUSOL->a[LL] = LUSOL->a[LM];
          LUSOL->indc[LL] = LUSOL->indc[LM];
          LUSOL->indr[LL] = LUSOL->indr[LM];
        }
      }
      else {
/*               =========================================================
                 There is no room for ( uk ) yet.  We have to
                 right-shift the whole of the remaining LU file.
                 Note that ( lk ) ends up in the correct place.
                 ========================================================= */
        LLSAVE = LL-LENLK;
        NMOVE = LM-LTOPL;
        for(IDUMMY = 1; IDUMMY <= NMOVE; IDUMMY++) {
          LL--;
          LM--;
          LUSOL->a[LL] = LUSOL->a[LM];
          LUSOL->indc[LL] = LUSOL->indc[LM];
          LUSOL->indr[LL] = LUSOL->indr[LM];
        }
        LTOPL = LL;
        LL = LLSAVE;
        LM = LL;
      }
/*            ======================================================
              Left-shift ( uk ).
              ====================================================== */
      LUSOL->locr[I] = LU+1;
      L2 = LM-1;
      LM = LM-LENUK;
      for(L = LM; L <= L2; L++) {
        LU = LU+1;
        LUSOL->a[LU] = LUSOL->a[L];
        LUSOL->indr[LU] = LUSOL->indr[L];
      }
    }
/*         ---------------------------------------------------------------
           Save the lengths of the nonempty columns of  L,
           and initialize  locc(j)  for the LU update routines.
           --------------------------------------------------------------- */
    for(K = 1; K <= NUML0; K++) {
      LUSOL->lenc[K] = LUSOL->iqloc[K];
    }
    for(J = 1; J <= LUSOL->n; J++) {
      LUSOL->locc[J] = 0;
    }
/*         ---------------------------------------------------------------
           Test for singularity.
           lu6chk  sets  nsing, jsing, jumin, Lmax, Umax, DUmax, DUmin
           (including entries from the dense LU).
           inform = 1  if there are singularities (nsing gt 0).
           --------------------------------------------------------------- */
    LU6CHK(LUSOL, 1,LUSOL->lena,INFORM);
    NSING = LUSOL->luparm[LUSOL_IP_SINGULARITIES];
    JSING = LUSOL->luparm[LUSOL_IP_SINGULARINDEX];
    JUMIN = LUSOL->luparm[LUSOL_IP_COLINDEX_DUMIN];
    LMAX  = LUSOL->parmlu[LUSOL_RP_MAXMULT_L];
    UMAX  = LUSOL->parmlu[LUSOL_RP_MAXELEM_U];
    DUMAX = LUSOL->parmlu[LUSOL_RP_MAXELEM_DIAGU];
    DUMIN = LUSOL->parmlu[LUSOL_RP_MINELEM_DIAGU];
  }
  else {
/*         ---------------------------------------------------------------
           keepLU = 0.  L and U were not kept, just the diagonals of U.
           lu1fac will probably be called again soon with keepLU = .true.
           11 Mar 2001: lu6chk revised.  We can call it with keepLU = 0,
                        but we want to keep Lmax, Umax from lu1fad.
           05 May 2002: Allow for TCP with new lu1DCP.  Diag(U) starts
                        below lena2, not lena.  Need lena2 in next line.
           --------------------------------------------------------------- */
#ifdef ClassicHamaxR
    LU6CHK(LUSOL, 1,LENA2,INFORM);
#else
    LU6CHK(LUSOL, 1,LUSOL->lena,INFORM);
#endif
    NSING = LUSOL->luparm[LUSOL_IP_SINGULARITIES];
    JSING = LUSOL->luparm[LUSOL_IP_SINGULARINDEX];
    JUMIN = LUSOL->luparm[LUSOL_IP_COLINDEX_DUMIN];
    DUMAX = LUSOL->parmlu[LUSOL_RP_MAXELEM_DIAGU];
    DUMIN = LUSOL->parmlu[LUSOL_RP_MINELEM_DIAGU];
  }
  goto x990;
/*      ------------
        Error exits.
        ------------ */
x930:
  *INFORM = LUSOL_INFORM_ADIMERR;
  if(LPRINT>=LUSOL_MSG_SINGULARITY)
    LUSOL_report(LUSOL, 0, "lu1fac  error...\nentry  a[%d]  has an illegal row (%d) or column (%d) index\n",
                        LERR,LUSOL->indc[LERR],LUSOL->indr[LERR]);
  goto x990;
x940:
  *INFORM = LUSOL_INFORM_ADUPLICATE;
  if(LPRINT>=LUSOL_MSG_SINGULARITY)
    LUSOL_report(LUSOL, 0, "lu1fac  error...\nentry  a[%d]  is a duplicate with indeces indc=%d, indr=%d\n",
                        LERR,LUSOL->indc[LERR],LUSOL->indr[LERR]);
  goto x990;
x970:
  *INFORM = LUSOL_INFORM_ANEEDMEM;
  if(LPRINT>=LUSOL_MSG_SINGULARITY)
    LUSOL_report(LUSOL, 0, "lu1fac  error...\ninsufficient storage; increase  lena  from %d to at least %d\n",
                        LUSOL->lena, MINLEN);
  goto x990;
x980:
  *INFORM = LUSOL_INFORM_FATALERR;
  if(LPRINT>=LUSOL_MSG_SINGULARITY)
    LUSOL_report(LUSOL, 0, "lu1fac  error...\nfatal bug   (sorry --- this should never happen)\n");
  goto x990;
x985:
  *INFORM = LUSOL_INFORM_NOPIVOT;
  if(LPRINT>=LUSOL_MSG_SINGULARITY)
    LUSOL_report(LUSOL, 0, "lu1fac  error...\nTSP used but diagonal pivot could not be found\n");

/*      Finalize and store output parameters. */
x990:
  LUSOL->nelem = NELEM0;
  LUSOL->luparm[LUSOL_IP_SINGULARITIES]   = NSING;
  LUSOL->luparm[LUSOL_IP_SINGULARINDEX]   = JSING;
  LUSOL->luparm[LUSOL_IP_MINIMUMLENA]     = MINLEN;
  LUSOL->luparm[LUSOL_IP_UPDATECOUNT]     = 0;
  LUSOL->luparm[LUSOL_IP_RANK_U]          = NRANK;
  LUSOL->luparm[LUSOL_IP_COLCOUNT_DENSE1] = NDENS1;
  LUSOL->luparm[LUSOL_IP_COLCOUNT_DENSE2] = NDENS2;
  LUSOL->luparm[LUSOL_IP_COLINDEX_DUMIN]  = JUMIN;
  LUSOL->luparm[LUSOL_IP_COLCOUNT_L0]     = NUML0;
  LUSOL->luparm[LUSOL_IP_ROWCOUNT_L0]     = 0;
  LUSOL->luparm[LUSOL_IP_NONZEROS_L0]     = LENL;
  LUSOL->luparm[LUSOL_IP_NONZEROS_U0]     = LENU;
  LUSOL->luparm[LUSOL_IP_NONZEROS_L]      = LENL;
  LUSOL->luparm[LUSOL_IP_NONZEROS_U]      = LENU;
  LUSOL->luparm[LUSOL_IP_NONZEROS_ROW]    = LROW;
  LUSOL->luparm[LUSOL_IP_MARKOWITZ_MERIT] = MERSUM;
  LUSOL->luparm[LUSOL_IP_TRIANGROWS_U]    = NUTRI;
  LUSOL->luparm[LUSOL_IP_TRIANGROWS_L]    = NLTRI;
  LUSOL->parmlu[LUSOL_RP_MAXELEM_A]       = AMAX;
  LUSOL->parmlu[LUSOL_RP_MAXMULT_L]       = LMAX;
  LUSOL->parmlu[LUSOL_RP_MAXELEM_U]       = UMAX;
  LUSOL->parmlu[LUSOL_RP_MAXELEM_DIAGU]   = DUMAX;
  LUSOL->parmlu[LUSOL_RP_MINELEM_DIAGU]   = DUMIN;
  LUSOL->parmlu[LUSOL_RP_MAXELEM_TCP]     = AKMAX;
  AGRWTH = AKMAX/(AMAX+LUSOL_SMALLNUM);
  UGRWTH = UMAX/(AMAX+LUSOL_SMALLNUM);
  if(TPP)
    GROWTH = UGRWTH;
/*      TRP or TCP or TSP */
  else
    GROWTH = AGRWTH;
  LUSOL->parmlu[LUSOL_RP_GROWTHRATE]      = GROWTH;

  LUSOL->luparm[LUSOL_IP_FTRANCOUNT]      = 0;
  LUSOL->luparm[LUSOL_IP_BTRANCOUNT]      = 0;

/*      ------------------------------------------------------------------
        Set overall status variable.
        ------------------------------------------------------------------ */
  LUSOL->luparm[LUSOL_IP_INFORM]          = *INFORM;
  if(*INFORM == LUSOL_INFORM_NOMEMLEFT)
    LUSOL_report(LUSOL, 0, "lu1fac  error...\ninsufficient memory available\n");
  
/*      ------------------------------------------------------------------
        Print statistics for the LU factors.
        ------------------------------------------------------------------ */
  NCP   = LUSOL->luparm[LUSOL_IP_COMPRESSIONS_LU];
  CONDU = DUMAX/MAX(DUMIN,LUSOL_SMALLNUM);
  DINCR = (LENL+LENU)-LUSOL->nelem;
  DINCR = (DINCR*100)/MAX(DELEM,ONE);
  AVGMER = MERSUM;
  AVGMER = AVGMER/DM;
  NBUMP = LUSOL->m-NUTRI-NLTRI;
  if(LPRINT>=LUSOL_MSG_STATISTICS) {
    if(TPP) {
      LUSOL_report(LUSOL, 0, "Merit %g %d %d %d %g %d %d %g %g %d %d %d\n",
                          AVGMER,LENL,LENL+LENU,NCP,DINCR,NUTRI,LENU,
                          LTOL,UMAX,UGRWTH,NLTRI,NDENS1,LMAX);
    }
    else {
      LUSOL_report(LUSOL, 0, "Merit %s %g %d %d %d %g %d %d %g %g %d %d %d %g %g\n",
                          LUSOL_pivotLabel(LUSOL),
                          AVGMER,LENL,LENL+LENU,NCP,DINCR,NUTRI,LENU,
                          LTOL,UMAX,UGRWTH,NLTRI,NDENS1,LMAX,AKMAX,AGRWTH);
    }
    LUSOL_report(LUSOL, 0, "bump%9d  dense2%7d  DUmax%g DUmin%g  conDU%g\n",
                          NBUMP,NDENS2,DUMAX,DUMIN,CONDU);
  }
}



//#include "lusol7a.c"     /* Utility routines for updates */

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   File  lusol7a
      lu7add   lu7cyc   lu7elm   lu7for   lu7rnk   lu7zap
      Utilities for LUSOL's update routines.
      lu7for is the most important -- the forward sweep.
  01 May 2002: Derived from LUSOL's original lu7a.f file.
  01 May 2002: Current version of lusol7a.f.
   ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

/* ==================================================================
   lu7add  inserts the first nrank elements of the vector v(*)
   as column  jadd  of  U.  We assume that  U  does not yet have any
   entries in this column.
   Elements no larger than  parmlu(3)  are treated as zero.
   klast  will be set so that the last row to be affected
   (in pivotal order) is row  ip(klast).
   ------------------------------------------------------------------
   09 May 1988: First f77 version.
   ================================================================== */
void LU7ADD(LUSOLrec *LUSOL, int JADD, REAL V[], int LENL, int *LENU,
  int *LROW, int NRANK, int *INFORM, int *KLAST, REAL *VNORM)
{
  REAL SMALL;
  int  K, I, LENI, MINFRE, NFREE, LR1, LR2, L;
#ifndef LUSOLFastMove
  int J;
#endif

  SMALL = LUSOL->parmlu[LUSOL_RP_ZEROTOLERANCE];
  *VNORM = ZERO;
  *KLAST = 0;
  for(K = 1; K <= NRANK; K++) {
    I = LUSOL->ip[K];
    if(fabs(V[I])<=SMALL)
      continue;
    *KLAST = K;
    (*VNORM) += fabs(V[I]);
    LENI = LUSOL->lenr[I];
/*         Compress row file if necessary. */
    MINFRE = LENI+1;
    NFREE = LUSOL->lena - LENL - *LROW;
    if(NFREE<MINFRE) {
      LU1REC(LUSOL, LUSOL->m, TRUE,LROW,LUSOL->indr,LUSOL->lenr,LUSOL->locr);
      NFREE = LUSOL->lena - LENL - *LROW;
      if(NFREE<MINFRE)
        goto x970;
    }
/*         Move row  i  to the end of the row file,
           unless it is already there.
           No need to move if there is a gap already. */
    if(LENI==0)
      LUSOL->locr[I] = (*LROW) + 1;
    LR1 = LUSOL->locr[I];
    LR2 = (LR1+LENI)-1;
    if(LR2==*LROW)
      goto x150;
    if(LUSOL->indr[LR2+1]==0)
      goto x180;
    LUSOL->locr[I] = (*LROW) + 1;
#ifdef LUSOLFastMove
    L = LR2-LR1+1;
    if(L > 0) {
      LR2 = (*LROW)+1;
      MEMMOVE(LUSOL->a+LR2,    LUSOL->a+LR1, L);
      MEMMOVE(LUSOL->indr+LR2, LUSOL->indr+LR1, L);
      MEMCLEAR(LUSOL->indr+LR1, L);
      *LROW += L;
    }
#else
    for(L = LR1; L <= LR2; L++) {
      (*LROW)++;
      LUSOL->a[*LROW] = LUSOL->a[L];
      J = LUSOL->indr[L];
      LUSOL->indr[L] = 0;
      LUSOL->indr[*LROW] = J;
    }
#endif
x150:
    LR2 = *LROW;
    (*LROW)++;
/*         Add the element of  v. */
x180:
    LR2++;
    LUSOL->a[LR2] = V[I];
    LUSOL->indr[LR2] = JADD;
    LUSOL->lenr[I] = LENI+1;
    (*LENU)++;
  }
/*      Normal exit. */
  *INFORM = LUSOL_INFORM_LUSUCCESS;
  goto x990;
/*      Not enough storage. */
x970:
  *INFORM = LUSOL_INFORM_ANEEDMEM;
x990:
;
}

/* ==================================================================
   lu7cyc performs a cyclic permutation on the row or column ordering
   stored in ip, moving entry kfirst down to klast.
   If kfirst .ge. klast, lu7cyc should not be called.
   Sometimes klast = 0 and nothing should happen.
   ------------------------------------------------------------------
   09 May 1988: First f77 version.
   ================================================================== */
void LU7CYC(LUSOLrec *LUSOL, int KFIRST, int KLAST, int IX[])
{
  if(KFIRST<KLAST) {
    int IFIRST, K;
#ifdef LUSOLFastMove
#if 1
    IFIRST = IX[KFIRST];
    K = KLAST-KFIRST;
    MEMMOVE(IX+KFIRST, IX+KFIRST+1, K);
    IX[KLAST] = IFIRST;
#else
    int *IXK, *IXK1;
    IXK = IX+KFIRST;
    IFIRST = *IXK;
    for(K = KFIRST, IXK1 = IXK+1; K <= KLAST-1; K++, IXK++, IXK1++) {
      *IXK = *IXK1;
    }
    *IXK = IFIRST;
#endif
#else
    IFIRST = IX[KFIRST];
    for(K = KFIRST; K <= KLAST-1; K++) {
      IX[K] = IX[K+1];
    }
    IX[KLAST] = IFIRST;
#endif
  }
}

/* ==================================================================
   lu7elm  eliminates the subdiagonal elements of a vector  v(*),
   where  L*v = y  for some vector y.
   If  jelm > 0,  y  has just become column  jelm  of the matrix  A.
   lu7elm  should not be called unless  m  is greater than  nrank.
   inform = 0 if y contained no subdiagonal nonzeros to eliminate.
   inform = 1 if y contained at least one nontrivial subdiagonal.
   inform = 7 if there is insufficient storage.
   ------------------------------------------------------------------
   09 May 1988: First f77 version.
                No longer calls lu7for at end.  lu8rpc, lu8mod do so.
   ================================================================== */
void LU7ELM(LUSOLrec *LUSOL, int JELM, REAL V[], int *LENL,
            int *LROW, int NRANK, int *INFORM, REAL *DIAG)
{
  REAL VI, VMAX, SMALL;
  int  NRANK1, MINFRE, NFREE, KMAX, L, K, I, LMAX, IMAX, L1, L2;

#ifdef ForceInitialization
  LMAX = 0;
#endif

  SMALL = LUSOL->parmlu[LUSOL_RP_ZEROTOLERANCE];
  NRANK1 = NRANK+1;
  *DIAG = ZERO;
/*      Compress row file if necessary. */
  MINFRE = LUSOL->m-NRANK;
  NFREE = LUSOL->lena-(*LENL)-(*LROW);
  if(NFREE>=MINFRE)
    goto x100;
  LU1REC(LUSOL, LUSOL->m,TRUE,LROW,LUSOL->indr,LUSOL->lenr,LUSOL->locr);
  NFREE = LUSOL->lena-(*LENL)-(*LROW);
  if(NFREE<MINFRE)
    goto x970;
    
/*      Pack the subdiagonals of  v  into  L,  and find the largest. */
x100:
  VMAX = ZERO;
  KMAX = 0;
  L = (LUSOL->lena-(*LENL))+1;
  for(K = NRANK1; K <= LUSOL->m; K++) {
    I = LUSOL->ip[K];
    VI = fabs(V[I]);
    if(VI<=SMALL)
      continue;
    L--;
    LUSOL->a[L] = V[I];
    LUSOL->indc[L] = I;
    if(VMAX>=VI)
      continue;
    VMAX = VI;
    KMAX = K;
    LMAX = L;
  }
  if(KMAX==0)
    goto x900;
/*      ------------------------------------------------------------------
        Remove  vmax  by overwriting it with the last packed  v(i).
        Then set the multipliers in  L  for the other elements.
        ------------------------------------------------------------------ */
  IMAX = LUSOL->ip[KMAX];
  VMAX = LUSOL->a[LMAX];
  LUSOL->a[LMAX] = LUSOL->a[L];
  LUSOL->indc[LMAX] = LUSOL->indc[L];
  L1 = L+1;
  L2 = LUSOL->lena-(*LENL);
  *LENL = ((*LENL)+L2)-L;
  for(L = L1; L <= L2; L++) {
    LUSOL->a[L] /= -VMAX;
    LUSOL->indr[L] = IMAX;
  }
/*      Move the row containing vmax to pivotal position nrank + 1. */
  LUSOL->ip[KMAX] = LUSOL->ip[NRANK1];
  LUSOL->ip[NRANK1] = IMAX;
  *DIAG = VMAX;
/*      ------------------------------------------------------------------
        If jelm is positive, insert  vmax  into a new row of  U.
        This is now the only subdiagonal element.
        ------------------------------------------------------------------ */
  if(JELM>0) {
    (*LROW)++;
    LUSOL->locr[IMAX] = *LROW;
    LUSOL->lenr[IMAX] = 1;
    LUSOL->a[*LROW] = VMAX;
    LUSOL->indr[*LROW] = JELM;
  }
  *INFORM = LUSOL_INFORM_LUSINGULAR;
  goto x990;
/*      No elements to eliminate. */
x900:
  *INFORM = LUSOL_INFORM_LUSUCCESS;
  goto x990;
/*      Not enough storage. */
x970:
  *INFORM = LUSOL_INFORM_ANEEDMEM;
x990:
;
}

/* ==================================================================
   lu7for  (forward sweep) updates the LU factorization  A = L*U
   when row  iw = ip(klast)  of  U  is eliminated by a forward
   sweep of stabilized row operations, leaving  ip * U * iq  upper
   triangular.
   The row permutation  ip  is updated to preserve stability and/or
   sparsity.  The column permutation  iq  is not altered.
   kfirst  is such that row  ip(kfirst)  is the first row involved
   in eliminating row  iw.  (Hence,  kfirst  marks the first nonzero
   in row  iw  in pivotal order.)  If  kfirst  is unknown it may be
   input as  1.
   klast   is such that row  ip(klast)  is the row being eliminated.
   klast   is not altered.
   lu7for  should be called only if  kfirst .le. klast.
   If  kfirst = klast,  there are no nonzeros to eliminate, but the
   diagonal element of row  ip(klast)  may need to be moved to the
   front of the row.
   ------------------------------------------------------------------
   On entry,  locc(*)  must be zero.

   On exit:
   inform = 0  if row iw has a nonzero diagonal (could be small).
   inform = 1  if row iw has no diagonal.
   inform = 7  if there is not enough storage to finish the update.

   On a successful exit (inform le 1),  locc(*)  will again be zero.
   ------------------------------------------------------------------
      Jan 1985: Final f66 version.
   09 May 1988: First f77 version.
   ================================================================== */
void LU7FOR(LUSOLrec *LUSOL, int KFIRST, int KLAST, int *LENL, int *LENU,
                     int *LROW, int *INFORM, REAL *DIAG)
{
  MYBOOL SWAPPD;
  int    KBEGIN, IW, LENW, LW1, LW2, JFIRST, MINFRE, NFREE, L, J, KSTART, KSTOP, K,
         LFIRST, IV, LENV, LV1, JLAST, LV2, LV3, LV, JV, LW, LDIAG, LIMIT;
  REAL   AMULT, LTOL, USPACE, SMALL, VJ, WJ;

  LTOL   = LUSOL->parmlu[LUSOL_RP_UPDATEMAX_Lij];
  SMALL  = LUSOL->parmlu[LUSOL_RP_ZEROTOLERANCE];
  USPACE = LUSOL->parmlu[LUSOL_RP_COMPSPACE_U];
  KBEGIN = KFIRST;
  SWAPPD = FALSE;

/*      We come back here from below if a row interchange is performed. */
x100:
  IW = LUSOL->ip[KLAST];
  LENW = LUSOL->lenr[IW];
  if(LENW==0)
    goto x910;
  LW1 = LUSOL->locr[IW];
  LW2 = (LW1+LENW)-1;
  JFIRST = LUSOL->iq[KBEGIN];
  if(KBEGIN>=KLAST)
    goto x700;
/*      Make sure there is room at the end of the row file
        in case row  iw  is moved there and fills in completely. */
  MINFRE = LUSOL->n+1;
  NFREE = LUSOL->lena-(*LENL)-(*LROW);
  if(NFREE<MINFRE) {
    LU1REC(LUSOL, LUSOL->m,TRUE,LROW,LUSOL->indr,LUSOL->lenr,LUSOL->locr);
    LW1 = LUSOL->locr[IW];
    LW2 = (LW1+LENW)-1;
    NFREE = LUSOL->lena-(*LENL)-(*LROW);
    if(NFREE<MINFRE)
      goto x970;

  }
/*      Set markers on row  iw. */
  for(L = LW1; L <= LW2; L++) {
    J = LUSOL->indr[L];
    LUSOL->locc[J] = L;
  }
/*      ==================================================================
        Main elimination loop.
        ================================================================== */
  KSTART = KBEGIN;
  KSTOP = MIN(KLAST,LUSOL->n);
  for(K = KSTART; K <= KSTOP; K++) {
    JFIRST = LUSOL->iq[K];
    LFIRST = LUSOL->locc[JFIRST];
    if(LFIRST==0)
      goto x490;
/*         Row  iw  has its first element in column  jfirst. */
    WJ = LUSOL->a[LFIRST];
    if(K==KLAST)
      goto x490;
/*         ---------------------------------------------------------------
           We are about to use the first element of row  iv
                  to eliminate the first element of row  iw.
           However, we may wish to interchange the rows instead,
           to preserve stability and/or sparsity.
           --------------------------------------------------------------- */
    IV = LUSOL->ip[K];
    LENV = LUSOL->lenr[IV];
    LV1 = LUSOL->locr[IV];
    VJ = ZERO;
    if(LENV==0)
      goto x150;
    if(LUSOL->indr[LV1]!=JFIRST)
      goto x150;
    VJ = LUSOL->a[LV1];
    if(SWAPPD)
      goto x200;
    if(LTOL*fabs(WJ)<fabs(VJ))
      goto x200;
    if(LTOL*fabs(VJ)<fabs(WJ))
      goto x150;
    if(LENV<=LENW)
      goto x200;
/*         ---------------------------------------------------------------
           Interchange rows  iv  and  iw.
           --------------------------------------------------------------- */
x150:
    LUSOL->ip[KLAST] = IV;
    LUSOL->ip[K] = IW;
    KBEGIN = K;
    SWAPPD = TRUE;
    goto x600;
/*         ---------------------------------------------------------------
           Delete the eliminated element from row  iw
           by overwriting it with the last element.
           --------------------------------------------------------------- */
x200:
    LUSOL->a[LFIRST] = LUSOL->a[LW2];
    JLAST = LUSOL->indr[LW2];
    LUSOL->indr[LFIRST] = JLAST;
    LUSOL->indr[LW2] = 0;
    LUSOL->locc[JLAST] = LFIRST;
    LUSOL->locc[JFIRST] = 0;
    LENW--;
    (*LENU)--;
    if(*LROW==LW2)
      (*LROW)--;
    LW2 = LW2-1;
/*         ---------------------------------------------------------------
           Form the multiplier and store it in the  L  file.
           --------------------------------------------------------------- */
    if(fabs(WJ)<=SMALL)
      goto x490;
    AMULT = -WJ/VJ;
    L = LUSOL->lena-(*LENL);
    LUSOL->a[L] = AMULT;
    LUSOL->indr[L] = IV;
    LUSOL->indc[L] = IW;
    (*LENL)++;
/*         ---------------------------------------------------------------
           Add the appropriate multiple of row  iv  to row  iw.
           We use two different inner loops.  The first one is for the
           case where row  iw  is not at the end of storage.
           --------------------------------------------------------------- */
    if(LENV==1)
      goto x490;
    LV2 = LV1+1;
    LV3 = (LV1+LENV)-1;
    if(LW2==*LROW)
      goto x400;
/*         ...............................................................
           This inner loop will be interrupted only if
           fill-in occurs enough to bump into the next row.
           ............................................................... */
    for(LV = LV2; LV <= LV3; LV++) {
      JV = LUSOL->indr[LV];
      LW = LUSOL->locc[JV];
      if(LW>0) {
/*               No fill-in. */
        LUSOL->a[LW] += AMULT*LUSOL->a[LV];
        if(fabs(LUSOL->a[LW])<=SMALL) {
/*                  Delete small computed element. */
          LUSOL->a[LW] = LUSOL->a[LW2];
          J = LUSOL->indr[LW2];
          LUSOL->indr[LW] = J;
          LUSOL->indr[LW2] = 0;
          LUSOL->locc[J] = LW;
          LUSOL->locc[JV] = 0;
          (*LENU)--;
          LENW--;
          LW2--;
        }
      }
      else {
/*               Row  iw  doesn't have an element in column  jv  yet
                 so there is a fill-in. */
        if(LUSOL->indr[LW2+1]!=0)
          goto x360;
        (*LENU)++;
        LENW++;
        LW2++;
        LUSOL->a[LW2] = AMULT*LUSOL->a[LV];
        LUSOL->indr[LW2] = JV;
        LUSOL->locc[JV] = LW2;
      }
    }
    goto x490;
/*         Fill-in interrupted the previous loop.
           Move row  iw  to the end of the row file. */
x360:
    LV2 = LV;
    LUSOL->locr[IW] = (*LROW)+1;

#ifdef LUSOLFastMove
    L = LW2-LW1+1;
    if(L > 0) {
      int loci, *locp;
      for(loci = LW1, locp = LUSOL->indr+LW1;
          loci <= LW2; loci++, locp++) {
        (*LROW)++;
        LUSOL->locc[*locp] = *LROW;
      }
      LW2 = (*LROW)-L+1;
      MEMMOVE(LUSOL->a+LW2,    LUSOL->a+LW1, L);
      MEMMOVE(LUSOL->indr+LW2, LUSOL->indr+LW1, L);
      MEMCLEAR(LUSOL->indr+LW1, L);
    }
#else
    for(L = LW1; L <= LW2; L++) {
      (*LROW)++;
      LUSOL->a[*LROW] = LUSOL->a[L];
      J = LUSOL->indr[L];
      LUSOL->indr[L] = 0;
      LUSOL->indr[*LROW] = J;
      LUSOL->locc[J] = *LROW;
    }
#endif
    LW1 = LUSOL->locr[IW];
    LW2 = *LROW;
/*         ...............................................................
           Inner loop with row  iw  at the end of storage.
           ............................................................... */
x400:
    for(LV = LV2; LV <= LV3; LV++) {
      JV = LUSOL->indr[LV];
      LW = LUSOL->locc[JV];
      if(LW>0) {
/*               No fill-in. */
        LUSOL->a[LW] += AMULT*LUSOL->a[LV];
        if(fabs(LUSOL->a[LW])<=SMALL) {
/*                  Delete small computed element. */
          LUSOL->a[LW] = LUSOL->a[LW2];
          J = LUSOL->indr[LW2];
          LUSOL->indr[LW] = J;
          LUSOL->indr[LW2] = 0;
          LUSOL->locc[J] = LW;
          LUSOL->locc[JV] = 0;
          (*LENU)--;
          LENW--;
          LW2--;
        }
      }
      else {
/*               Row  iw  doesn't have an element in column  jv  yet
                 so there is a fill-in. */
        (*LENU)++;
        LENW++;
        LW2++;
        LUSOL->a[LW2] = AMULT*LUSOL->a[LV];
        LUSOL->indr[LW2] = JV;
        LUSOL->locc[JV] = LW2;
      }
    }
    *LROW = LW2;
/*         The  k-th  element of row  iw  has been processed.
           Reset  swappd  before looking at the next element. */
x490:
    SWAPPD = FALSE;
  }
/*      ==================================================================
        End of main elimination loop.
        ==================================================================

        Cancel markers on row  iw. */
x600:
  LUSOL->lenr[IW] = LENW;
  if(LENW==0)
    goto x910;
  for(L = LW1; L <= LW2; L++) {
    J = LUSOL->indr[L];
    LUSOL->locc[J] = 0;
  }
/*      Move the diagonal element to the front of row  iw.
        At this stage,  lenw gt 0  and  klast le n. */
x700:
  for(L = LW1; L <= LW2; L++) {
    LDIAG = L;
    if(LUSOL->indr[L]==JFIRST)
      goto x730;
  }
  goto x910;

x730:
  *DIAG = LUSOL->a[LDIAG];
  LUSOL->a[LDIAG] = LUSOL->a[LW1];
  LUSOL->a[LW1] = *DIAG;
  LUSOL->indr[LDIAG] = LUSOL->indr[LW1];
  LUSOL->indr[LW1] = JFIRST;
/*      If an interchange is needed, repeat from the beginning with the
        new row  iw,  knowing that the opposite interchange cannot occur. */
  if(SWAPPD)
    goto x100;
  *INFORM = LUSOL_INFORM_LUSUCCESS;
  goto x950;
/*      Singular. */
x910:
  *DIAG = ZERO;
  *INFORM = LUSOL_INFORM_LUSINGULAR;
/*      Force a compression if the file for  U  is much longer than the
        no. of nonzeros in  U  (i.e. if  lrow  is much bigger than  lenU).
        This should prevent memory fragmentation when there is far more
        memory than necessary  (i.e. when  lena  is huge). */
x950:
  LIMIT = (int) (USPACE*(*LENU))+LUSOL->m+LUSOL->n+1000;
  if(*LROW>LIMIT)
    LU1REC(LUSOL, LUSOL->m,TRUE,LROW,LUSOL->indr,LUSOL->lenr,LUSOL->locr);
  goto x990;
/*      Not enough storage. */
x970:
  *INFORM = LUSOL_INFORM_ANEEDMEM;
/*      Exit. */
x990:
;
}

/* ==================================================================
   lu7rnk (check rank) assumes U is currently nrank by n
   and determines if row nrank contains an acceptable pivot.
   If not, the row is deleted and nrank is decreased by 1.
   jsing is an input parameter (not altered).  If jsing is positive,
   column jsing has already been judged dependent.  A substitute
   (if any) must be some other column.
   ------------------------------------------------------------------
   -- Jul 1987: First version.
   09 May 1988: First f77 version.
   ================================================================== */
void LU7RNK(LUSOLrec *LUSOL, int JSING, int *LENU,
            int *LROW, int *NRANK, int *INFORM, REAL *DIAG)
{
  REAL UTOL1, UMAX;
  int  IW, LENW, L1, L2, LMAX, L, JMAX, KMAX;

#ifdef ForceInitialization
  L1 = 0;
  L2 = 0;
#endif

  UTOL1 = LUSOL->parmlu[LUSOL_RP_SMALLDIAG_U];
  *DIAG = ZERO;
/*      Find Umax, the largest element in row nrank. */
  IW = LUSOL->ip[*NRANK];
  LENW = LUSOL->lenr[IW];
  if(LENW==0)
    goto x400;
  L1 = LUSOL->locr[IW];
  L2 = (L1+LENW)-1;
  UMAX = ZERO;
  LMAX = L1;
  for(L = L1; L <= L2; L++) {
    if(UMAX<fabs(LUSOL->a[L])) {
      UMAX = fabs(LUSOL->a[L]);
      LMAX = L;
    }
  }
/*      Find which column that guy is in (in pivotal order).
        Interchange him with column nrank, then move him to be
        the new diagonal at the front of row nrank. */
  *DIAG = LUSOL->a[LMAX];
  JMAX = LUSOL->indr[LMAX];
  for(KMAX = *NRANK; KMAX <= LUSOL->n; KMAX++) {
    if(LUSOL->iq[KMAX]==JMAX)
      break;
  }
  LUSOL->iq[KMAX] = LUSOL->iq[*NRANK];
  LUSOL->iq[*NRANK] = JMAX;
  LUSOL->a[LMAX] = LUSOL->a[L1];
  LUSOL->a[L1] = *DIAG;
  LUSOL->indr[LMAX] = LUSOL->indr[L1];
  LUSOL->indr[L1] = JMAX;
/*      See if the new diagonal is big enough. */
  if(UMAX<=UTOL1)
    goto x400;
  if(JMAX==JSING)
    goto x400;
/*      ------------------------------------------------------------------
        The rank stays the same.
        ------------------------------------------------------------------ */
  *INFORM = LUSOL_INFORM_LUSUCCESS;
  return;
/*      ------------------------------------------------------------------
        The rank decreases by one.
        ------------------------------------------------------------------ */
x400:
  *INFORM = LUSOL_INFORM_RANKLOSS;
  (*NRANK)--;
  if(LENW>0) {
/*         Delete row nrank from U. */
    LENU = LENU-LENW;
    LUSOL->lenr[IW] = 0;
    for(L = L1; L <= L2; L++) {
      LUSOL->indr[L] = 0;
    }
    if(L2==*LROW) {
/*            This row was at the end of the data structure.
              We have to reset lrow.
              Preceding rows might already have been deleted, so we
              have to be prepared to go all the way back to 1. */
      for(L = 1; L <= L2; L++) {
        if(LUSOL->indr[*LROW]>0)
          goto x900;
        (*LROW)--;
      }
    }
  }
x900:
;
}

/* ==================================================================
   lu7zap  eliminates all nonzeros in column  jzap  of  U.
   It also sets  kzap  to the position of  jzap  in pivotal order.
   Thus, on exit we have  iq(kzap) = jzap.
   ------------------------------------------------------------------
   -- Jul 1987: nrank added.
   10 May 1988: First f77 version.
   ================================================================== */
void LU7ZAP(LUSOLrec *LUSOL, int JZAP, int *KZAP, int *LENU, int *LROW,
            int NRANK)
{
  int K, I, LENI, LR1, LR2, L;

  for(K = 1; K <= NRANK; K++) {
    I = LUSOL->ip[K];
    LENI = LUSOL->lenr[I];
    if(LENI==0)
      goto x90;
    LR1 = LUSOL->locr[I];
    LR2 = (LR1+LENI)-1;
    for(L = LR1; L <= LR2; L++) {
      if(LUSOL->indr[L]==JZAP)
        goto x60;
    }
    goto x90;
/*         Delete the old element. */
x60:
    LUSOL->a[L] = LUSOL->a[LR2];
    LUSOL->indr[L] = LUSOL->indr[LR2];
    LUSOL->indr[LR2] = 0;
    LUSOL->lenr[I] = LENI-1;
    (*LENU)--;
/*         Stop if we know there are no more rows containing  jzap. */
x90:
    *KZAP = K;
    if(LUSOL->iq[K]==JZAP)
      goto x800;
  }
/*      nrank must be smaller than n because we haven't found kzap yet. */
  L = LUSOL->n;
  for(K = NRANK+1; K <= L; K++) {
    *KZAP = K;
    if(LUSOL->iq[K]==JZAP)
      break;
  }
/*      See if we zapped the last element in the file. */
x800:
  if(*LROW>0) {
    if(LUSOL->indr[*LROW]==0)
      (*LROW)--;
  }

}




//#include "lusol8a.c"     /* Column update */

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   File  lusol8a
      lu8rpc
      Sparse LU update: Replace Column
      LUSOL's sparse implementation of the Bartels-Golub update.

   01 May 2002: Derived from LUSOL's original lu8a.f file.
   01 May 2002: Current version of lusol8a.f.
   ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

/* ==================================================================
   lu8rpc  updates the LU factorization  A = L*U  when column  jrep
   is replaced by some vector  a(new).
   lu8rpc  is an implementation of the Bartels-Golub update,
   designed for the case where A is rectangular and/or singular.
   L is a product of stabilized eliminations (m x m, nonsingular).
   P U Q is upper trapezoidal (m x n, rank nrank).

   If  mode1 = 0,  the old column is taken to be zero
                   (so it does not have to be removed from  U).
   If  mode1 = 1,  the old column need not have been zero.
   If  mode2 = 0,  the new column is taken to be zero.
                   v(*)  is not used or altered.
   If  mode2 = 1,  v(*)  must contain the new column  a(new).
                   On exit,  v(*)  will satisfy  L*v = a(new).
   If  mode2 = 2,  v(*)  must satisfy  L*v = a(new).

   The array  w(*)  is not used or altered.
   On entry, all elements of  locc  are assumed to be zero.
   On a successful exit (inform != 7), this will again be true.
   On exit:

   inform = -1  if the rank of U decreased by 1.
   inform =  0  if the rank of U stayed the same.
   inform =  1  if the rank of U increased by 1.
   inform =  2  if the update seemed to be unstable
                (diag much bigger than vnorm).
   inform =  7  if the update was not completed (lack of storage).
   inform =  8  if jrep is not between 1 and n.
   ------------------------------------------------------------------
   -- Jan 1985: Original F66 version.
   -- Jul 1987: Modified to maintain U in trapezoidal form.
   10 May 1988: First f77 version.
   16 Oct 2000: Added test for instability (inform = 2).
   ================================================================== */
void LU8RPC(LUSOLrec *LUSOL, int MODE1, int MODE2,
            int JREP, REAL V[], REAL W[],
            int *INFORM, REAL *DIAG, REAL *VNORM)
{
  MYBOOL SINGLR;
  int    LPRINT, NRANK, LENL, LENU, LROW, NRANK0, KREP, KLAST, IW, L1, J1, JSING;
  REAL   UTOL1, UTOL2;

  LPRINT = LUSOL->luparm[LUSOL_IP_PRINTLEVEL];
  NRANK  = LUSOL->luparm[LUSOL_IP_RANK_U];
  LENL   = LUSOL->luparm[LUSOL_IP_NONZEROS_L];
  LENU   = LUSOL->luparm[LUSOL_IP_NONZEROS_U];
  LROW   = LUSOL->luparm[LUSOL_IP_NONZEROS_ROW];
  UTOL1  = LUSOL->parmlu[LUSOL_RP_SMALLDIAG_U];
  UTOL2  = LUSOL->parmlu[LUSOL_RP_EPSDIAG_U];
  NRANK0 = NRANK;
  *DIAG  = ZERO;
  *VNORM = ZERO;
  if(JREP<1)
    goto x980;
  if(JREP>LUSOL->n)
    goto x980;

/*      ------------------------------------------------------------------
        If mode1 = 0, there are no elements to be removed from  U
        but we still have to set  krep  (using a backward loop).
        Otherwise, use lu7zap to remove column  jrep  from  U
        and set  krep  at the same time.
        ------------------------------------------------------------------ */
  if(MODE1==LUSOL_UPDATE_OLDEMPTY) {
    KREP = LUSOL->n+1;
x10:
    KREP--;
    if(LUSOL->iq[KREP]!=JREP)
      goto x10;
  }
  else
    LU7ZAP(LUSOL, JREP,&KREP,&LENU,&LROW,NRANK);

/*      ------------------------------------------------------------------
        Insert a new column of u and find klast.
        ------------------------------------------------------------------ */
  if(MODE2==LUSOL_UPDATE_NEWEMPTY) {
    KLAST = 0;
  }
  else {
    if(MODE2==LUSOL_UPDATE_NEWNONEMPTY) {
/*            Transform v = a(new) to satisfy  L*v = a(new). */
      LU6SOL(LUSOL, LUSOL_SOLVE_Lv_v, V,W, NULL, INFORM);
    }
    else if(V==NULL)
/* Otherwise, the V vector is taken to satisfy this already, or stored earlier. */
      V=LUSOL->vLU6L;
      

/*         Insert into  U  any nonzeros in the top of  v.
           row  ip(klast)  will contain the last nonzero in pivotal order.
           Note that  klast  will be in the range  ( 0, nrank ). */
    LU7ADD(LUSOL, JREP,V,LENL,&LENU,&LROW,NRANK,INFORM,&KLAST,VNORM);
    if(*INFORM==LUSOL_INFORM_ANEEDMEM)
      goto x970;
  }
/*      ------------------------------------------------------------------
        In general, the new column causes U to look like this:
                    krep        n                 krep  n
                   ....a.........          ..........a...
                    .  a        .           .        a  .
                     . a        .            .       a  .
                      .a        .             .      a  .
           P U Q =     a        .    or        .     a  .
                       b.       .               .    a  .
                       b .      .                .   a  .
                       b  .     .                 .  a  .
                       b   ......                  ..a...  nrank
                       c                             c
                       c                             c
                       c                             c     m
        klast points to the last nonzero "a" or "b".
        klast = 0 means all "a" and "b" entries are zero.
        ------------------------------------------------------------------ */
  if(MODE2==LUSOL_UPDATE_NEWEMPTY) {
    if(KREP>NRANK)
      goto x900;
  }
  else if(NRANK<LUSOL->m) {
/*         Eliminate any "c"s (in either case).
           Row nrank + 1 may end up containing one nonzero. */
    LU7ELM(LUSOL, JREP,V,&LENL,&LROW,NRANK,INFORM,DIAG);
    if(*INFORM==LUSOL_INFORM_ANEEDMEM)
      goto x970;
    if(*INFORM==LUSOL_INFORM_LUSINGULAR) {
/*            The nonzero is apparently significant.
              Increase nrank by 1 and make klast point to the bottom. */
      NRANK++;
      KLAST = NRANK;
    }
  }
  if(NRANK<LUSOL->n) {
/*         The column rank is low.
           In the first case, we want the new column to end up in
           position nrank, so the trapezoidal columns will have a chance
           later on (in lu7rnk) to pivot in that position.
           Otherwise the new column is not part of the triangle.  We
           swap it into position nrank so we can judge it for singularity.
           lu7rnk might choose some other trapezoidal column later. */
    if(KREP<NRANK)
      KLAST = NRANK;
    else {
      LUSOL->iq[KREP] = LUSOL->iq[NRANK];
      LUSOL->iq[NRANK] = JREP;
      KREP = NRANK;
    }
  }
/*      ------------------------------------------------------------------
        If krep .lt. klast, there are some "b"s to eliminate:
                     krep
                   ....a.........
                    .  a        .
                     . a        .
                      .a        .
           P U Q =     a        .  krep
                       b.       .
                       b .      .
                       b  .     .
                       b   ......  nrank
        If krep .eq. klast, there are no "b"s, but the last "a" still
        has to be moved to the front of row krep (by lu7for).
        ------------------------------------------------------------------ */
  if(KREP<=KLAST) {
/*         Perform a cyclic permutation on the current pivotal order,
           and eliminate the resulting row spike.  krep becomes klast.
           The final diagonal (if any) will be correctly positioned at
           the front of the new krep-th row.  nrank stays the same. */
    LU7CYC(LUSOL, KREP,KLAST,LUSOL->ip);
    LU7CYC(LUSOL, KREP,KLAST,LUSOL->iq);
    LU7FOR(LUSOL, KREP,KLAST,&LENL,&LENU,&LROW,INFORM,DIAG);
    if(*INFORM==LUSOL_INFORM_ANEEDMEM)
      goto x970;
    KREP = KLAST;
/*         Test for instability (diag much bigger than vnorm). */
    SINGLR = (MYBOOL) ((*VNORM)<UTOL2*fabs(*DIAG));
    if(SINGLR)
      goto x920;
  }
/*      ------------------------------------------------------------------
        Test for singularity in column krep (where krep .le. nrank).
        ------------------------------------------------------------------ */
  *DIAG = ZERO;
  IW = LUSOL->ip[KREP];
  SINGLR = (MYBOOL) (LUSOL->lenr[IW]==0);
  if(!SINGLR) {
    L1 = LUSOL->locr[IW];
    J1 = LUSOL->indr[L1];
    SINGLR = (MYBOOL) (J1!=JREP);
    if(!SINGLR) {
      *DIAG = LUSOL->a[L1];
      SINGLR = (MYBOOL) (fabs(*DIAG)<=UTOL1 || fabs(*DIAG)<=UTOL2*(*VNORM));
    }
  }
  if(SINGLR && KREP<NRANK) {
/*         Perform cyclic permutations to move column jrep to the end.
           Move the corresponding row to position nrank
           then eliminate the resulting row spike. */
    LU7CYC(LUSOL, KREP,NRANK,LUSOL->ip);
    LU7CYC(LUSOL, KREP,LUSOL->n,LUSOL->iq);
    LU7FOR(LUSOL, KREP,NRANK,&LENL,&LENU,&LROW,INFORM,DIAG);
    if(*INFORM==LUSOL_INFORM_ANEEDMEM)
      goto x970;
  }
/*      Find the best column to be in position nrank.
        If singlr, it can't be the new column, jrep.
        If nothing satisfactory exists, nrank will be decreased. */
  if(SINGLR || NRANK<LUSOL->n) {
    JSING = 0;
    if(SINGLR)
      JSING = JREP;
    LU7RNK(LUSOL, JSING,&LENU,&LROW,&NRANK,INFORM,DIAG);
  }

/*      ------------------------------------------------------------------
        Update indeces of optional row-based version of L0.
        ------------------------------------------------------------------ */
#if 0
  if(LUSOL->L0 != NULL)
    LU1L0UPD(LUSOL, INFORM);
#endif

/*      ------------------------------------------------------------------
        Set inform for exit.
        ------------------------------------------------------------------ */
x900:
  if(NRANK==NRANK0)
    *INFORM = LUSOL_INFORM_LUSUCCESS;
  else if(NRANK<NRANK0) {
    *INFORM = LUSOL_INFORM_RANKLOSS;
    if(NRANK0==LUSOL->n) {
      if(LPRINT>=LUSOL_MSG_SINGULARITY)
        LUSOL_report(LUSOL, 0, "lu8rpc  warning...\nSingularity after replacing column.    jrep=%8d    diag=%g\n",
                            JREP,DIAG);
    }
  }
  else
    *INFORM = LUSOL_INFORM_LUSINGULAR;
  goto x990;
/*      Instability. */
x920:
  *INFORM = LUSOL_INFORM_LUUNSTABLE;
  if(LPRINT>=LUSOL_MSG_SINGULARITY)
    LUSOL_report(LUSOL, 0, "lu8rpc  warning...\nInstability after replacing column.    jrep=%8d    diag=%g\n",
                        JREP,DIAG);
  goto x990;
/*      Not enough storage. */
x970:
  *INFORM = LUSOL_INFORM_ANEEDMEM;
  if(LPRINT>=LUSOL_MSG_SINGULARITY)
    LUSOL_report(LUSOL, 0, "lu8rpc  error...\nInsufficient memory.    lena=%8d\n",
                        LUSOL->lena);
  goto x990;
/*      jrep  is out of range. */
x980:
  *INFORM = LUSOL_INFORM_FATALERR;
  if(LPRINT>=LUSOL_MSG_SINGULARITY)
    LUSOL_report(LUSOL, 0, "lu8rpc  error...\njrep  is out of range.    m=%8d    n=%8d    jrep=%8d\n",
                        LUSOL->m,LUSOL->n,JREP);
/*      Exit. */
x990:
  LUSOL->luparm[LUSOL_IP_INFORM]       = *INFORM;
  LUSOL->luparm[LUSOL_IP_UPDATECOUNT]++;
  LUSOL->luparm[LUSOL_IP_RANK_U]       = NRANK;
  LUSOL->luparm[LUSOL_IP_NONZEROS_L]   = LENL;
  LUSOL->luparm[LUSOL_IP_NONZEROS_U]   = LENU;
  LUSOL->luparm[LUSOL_IP_NONZEROS_ROW] = LROW;
}


void LUSOL_dump(FILE *output, LUSOLrec *LUSOL)
{
  MYBOOL userfile = (MYBOOL) (output != NULL);

  if(!userfile)
    output = fopen("LUSOL.dbg", "w");

  blockWriteREAL(output, "a", LUSOL->a, 1, LUSOL->lena);
  blockWriteINT(output, "indc", LUSOL->indc, 1, LUSOL->lena);
  blockWriteINT(output, "indr", LUSOL->indr, 1, LUSOL->lena);

  blockWriteINT(output, "ip", LUSOL->ip, 1, LUSOL->m);
  blockWriteINT(output, "iq", LUSOL->iq, 1, LUSOL->n);
  blockWriteINT(output, "lenc", LUSOL->lenc, 1, LUSOL->n);
  blockWriteINT(output, "lenr", LUSOL->lenr, 1, LUSOL->m);

  blockWriteINT(output, "locc", LUSOL->locc, 1, LUSOL->n);
  blockWriteINT(output, "locr", LUSOL->locr, 1, LUSOL->m);
  blockWriteINT(output, "iploc", LUSOL->iploc, 1, LUSOL->n);
  blockWriteINT(output, "iqloc", LUSOL->iqloc, 1, LUSOL->m);

  blockWriteINT(output, "ipinv", LUSOL->ipinv, 1, LUSOL->m);
  blockWriteINT(output, "iqinv", LUSOL->iqinv, 1, LUSOL->n);

  if(!userfile)
    fclose(output);
}

LUSOLmat *LUSOL_matcreate(int dim, int nz)
{
  LUSOLmat *newm;

  newm = (LUSOLmat *) calloc(1, sizeof(*newm));
  if(newm != NULL) {
    newm->a    = (REAL *) malloc((nz+1)*sizeof(REAL));
    newm->lenx = (int *)  malloc((dim+1)*sizeof(int));
    newm->indx = (int *)  malloc((dim+1)*sizeof(int));
    newm->indr = (int *)  malloc((nz+1)*sizeof(int));
    newm->indc = (int *)  malloc((nz+1)*sizeof(int));
    if((newm->a == NULL) ||
       (newm->lenx == NULL) || (newm->indx == NULL) ||
       (newm->indr == NULL) || (newm->indc == NULL))
      LUSOL_matfree(&newm);
  }
  return(newm);
}
void LUSOL_matfree(LUSOLmat **mat)
{
  if((mat == NULL) || (*mat == NULL))
    return;
  FREE((*mat)->a);
  FREE((*mat)->indc);
  FREE((*mat)->indr);
  FREE((*mat)->lenx);
  FREE((*mat)->indx);
  FREE(*mat);
}

