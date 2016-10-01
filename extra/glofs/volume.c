#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>

#ifndef _UNDERSCORE
float cvolume(float *a, float *b, int *m, int *n, int *mmax, int *nmax);
void volume(float *a, float *b, int *m, int *n, int *mmax, int *nmax, float *volume);
float cvolumeb(float *a, float *b, int *m, int *n, int *mmax, int *nmax, int *opt, float *dvdb);
void volumeb(float *a, float *b, int *m, int *n, int *mmax, int *nmax, int *opt, float *volume, float *dvdb);
float cvolumebj(float *a, float *b, int *m, int *n, int *mmax, int *nmax, int con, float *dvdb);
void volumebj(float *a, float *b, int *m, int *n, int *mmax, int *nmax, int *con, float *volume, float *dvdb);
float cdvda(float *a, float *b, int *m, int *n, int *mmax, int *nmax, int *tdim, float *temp, int *jval, int *code);
void dvda(float *a, float *b, int *m, int *n, int *mmax, int *nmax, int *idim, float *dvda, int *code);
float cvolumef(float *a, float *b, int *m, int *n, int *mmax, int *nmax);
void volumef(float *a, float *b, int *m, int *n, int *mmax, int *nmax, float *volume);
float cdvdaf(float *a, float *b, int *m, int *n, int *mmax, int *nmax, int *tdim, float *temp, int *jval, int *code);
void dvdaf(float *a, float *b, int *m, int *n, int *mmax, int *nmax, int *idim, float *dvda, int *code);
float cdvda_debug(float *a, float *b, int *m, int *n, int *tdim, float *temp, int *jval, int *code); 
#else
float cvolume_ (float *a, float *b, int *m, int *n, int *mmax, int *nmax);
void volume_ (float *a, float *b, int *m, int *n, int *mmax, int *nmax, float *volume);
float cvolumeb_ (float *a, float *b, int *m, int *n, int *mmax, int *nmax, int *opt, float *dvdb);
void volumeb_ (float *a, float *b, int *m, int *n, int *mmax, int *nmax, int *opt, float *volume, float *dvdb);
float cvolumebj_ (float *a, float *b, int *m, int *n, int *mmax, int *nmax, int con, float *dvdb);
void volumebj_ (float *a, float *b, int *m, int *n, int *mmax, int *nmax, int *con, float *volume, float *dvdb);
float cdvda_ (float *a, float *b, int *m, int *n, int *mmax, int *nmax, int *tdim, float *temp, int *jval, int *code);
void dvda_ (float *a, float *b, int *m, int *n, int *mmax, int *nmax, int *idim, float *dvda, int *code);
float cvolumef_ (float *a, float *b, int *m, int *n, int *mmax, int *nmax);
void volumef_ (float *a, float *b, int *m, int *n, int *mmax, int *nmax, float *volume);
float cdvdaf_ (float *a, float *b, int *m, int *n, int *mmax, int *nmax, int *tdim, float *temp, int *jval, int *code);
void dvdaf_ (float *a, float *b, int *m, int *n, int *mmax, int *nmax, int *idim, float *dvda, int *code);
float cdvda_debug_ (float *a, float *b, int *m, int *n, int *tdim, float *temp, int *jval, int *code); 
#endif

/*--------------------------------------------------------------

    This file contains six recursive C routines for calculating
    the volume and derivatives of a convex polyhedron in n 
    dimensions. The routines volume, volumef, volumeb and dvda are 
    fortran callable interfaces to their C counterparts cvolume,
    cvolumeb and cdvda.

    ROUTINE      COMMENT			CALLS 

    cvolume	 calculates volume of region    cvolume 
                 which satisfies Ax <= b.

    cvolumef	 calculates volume of region    cvolumef 
                 which satisfies Ax <= b.
                 (faster version, see below)

    cvolumeb     calculates volume and the      cvolume
                 derivative d(vol)/db(0).

    cdvda        calculates derivatives		cdvda, cvolumeb, 
                 d(vol)/da(0,j) (j=1,n).        & cvolumebj

    cdvdaf        calculates derivatives	cdvdaf, cvolumeb, 
                 d(vol)/da(0,j) (j=1,n).        & cvolumebj
                 (faster version, see below)

    cvolumebj    calculates partial volumes
                 and derivatives d(vol)/d(b(j)
                 (j=1,n). (called by cdvda.)

    Here b(j) is the jth element of vector b, and a(i,j) is the
    (i,j)th coefficient of the matrix A corresponding to the
    ith constraint and the jth dimension.

    The two faster versions cvolumef and cdvdaf do not continue down
    the recursive tree if b(j) = 0 for any j at any time. This avoids
    extra work but restricts the applicability of the routines to
    cases where the origin does not pass through any inconsistent
    constraints.
     
				Malcolm Sambridge, April 1995.

  --------------------------------------------------------------*/

/*--------------------------------------------------------------

			ROUTINE: cvolume

    This routine calculates the `volume' in dimension n of the region
    bounded by a set of m linear inequality constraints of the form
    A x <= b, where a has m rows and n columns and is given by a(m,n),
    b is the n-vector and is contained in b(m). The recursive formula
    of Lasserre (1983) is used. Redundant constraints are allowed 
    If the inequality constraints are inconsistent then the volume 
    is returned as zero.  If the original polyhedron is unbounded 
    then a warning is issued and the volume is return as -1.0 
    (even though it is undefined). 

			Jean Braun and Malcolm Sambridge, Jan 1995.

    Lasserre, J. B., 1983. "An analytical expression and an Algorithm
    for the Volume of a Convex Polyhedron in R^n.", JOTA, vol. 39, 
    no. 3, 363-377.

  --------------------------------------------------------------*/

#ifndef _UNDERSCORE
float cvolume(a,b,m,n,mmax,nmax)
#else
float cvolume_ (a,b,m,n,mmax,nmax)
#endif
int   *n, *m, *mmax, *nmax;
float *a, *b;
{

float v,amax,pivot;
int   i,j,t,k,l;
int   jj,kk;
float   *ai, *aj, *ajt, *apjj, *bi;
int   kmm,tmm;
int   nn = *n,mm = *m, mm_max = *mmax;
int   nm1,mm1;
float  *ap, *bp;
int   firstmin,firstmax;
float bmin,bmax,bb;
float partialv;

/* one-dimensional case (full reduction) */

if (nn == 1) 
  {
  firstmin=0;
  firstmax=0;
  for (l=0;l<mm;l++)
    {
    if ( *(a+l) > 0.) 
      {
      bb= *(b+l)/ *(a+l);
      if (firstmin==0) {firstmin=1;bmin=bb;}
      else if (bb<bmin) bmin=bb;
      }
    else if ( *(a+l) < 0.)
      {
      bb= *(b+l)/ *(a+l);
      if (firstmax==0) {firstmax=1;bmax=bb;}
      else if (bb>bmax) bmax=bb;
      }
    else if ( *(b+l) < 0.)         /* Constraints are inconsistent */
      { 
       printf("Inconsistent constraints found after reduction to n = 1 \n");
       return(0.); 
      }
    }
  v=0.;
  if (firstmin*firstmax == 1) v=bmin-bmax;
  else 
    {
/*     printf("Volume is unbounded at 1-D; volume returned as -1\n");  */
       printf("Volume is unbounded at 1-D volume returned as -1\n");
     return(-1.0);
    }

  if (v<0.) {v=0.;} 

  return(v);
  }

nm1=nn-1;
mm1=mm-1;
v=0.;

for (i=0;i<mm;i++)
  {
   ai=a+i;
   bi=b+i;

/* find largest pivot */

  amax=0.;
  for (j=0;j<nn;j++) 
    if (fabs( *(ai+j*mm_max)) >= amax) {amax= fabs( *(ai+j*mm_max)); t=j;}
  tmm=t*mm_max;
  pivot=*(ai+tmm);

/* finds contribution to v from this pivot (if not nil) */

  if (amax == 0.)
  {
                                 /* Constraint is inconsistent */

    if ( *(bi) < 0.) 
      { printf("Constraint %d is inconsistent\n",i+1); return(0.); }

                                 /* otherwise constraint is redundant */

    printf("Constraint %d is redundant\n",i+1); 
  }
  else
  {

/* allocate memory */

  ap = (float *) malloc(4*nm1*mm1);
  bp = (float *) malloc(4*mm1);

/* reduce a and b into ap and bp eliminating variable t and constraint i */

  jj=-1;
    for (j=0;j<mm;j++)
      if (j != i)
      {
      jj=jj+1;
      aj=a+j;
      ajt=aj+tmm;
      *(bp+jj)= *(b+j) - *(bi) * *(ajt) / pivot;
      apjj=ap+jj;
      kk=-1;
      for (k=0;k<nn;k++)
        if (k != t)
        {
        kk=kk+1;
        kmm=k*mm_max;
        *(apjj+kk*mm1)= *(aj+kmm)- *(ajt) * *(ai+kmm)/ pivot;
        }
      }
  
/* add contribution to volume from volume calculated in smaller dimension */

#ifndef _UNDERSCORE
  partialv=cvolume(ap,bp, &mm1, &nm1, &mm1, &nm1);
#else
  partialv=cvolume_ (ap,bp, &mm1, &nm1, &mm1, &nm1);
#endif
  if(partialv == -1.0)return(-1.0);
  v=v+ *(bi)/amax*partialv/nn; 

  free(ap);
  free(bp);
  }
  }

return(v);

}

/*--------------------------------------------------------------

			ROUTINE: volume

   A dummy routine used to call cvolume from a fortran routine 


			Jean Braun and Malcolm Sambridge, Jan 1995.

----------------------------------------------------------------*/

#ifndef _UNDERSCORE
void volume(a,b,m,n,mmax,nmax,volume)
#else
void volume_ (a,b,m,n,mmax,nmax,volume)
#endif
int   *n, *m, *mmax, *nmax;
float *a, *b;
float *volume;
{
#ifndef _UNDERSCORE
*volume=cvolume(a,b,m,n,mmax,nmax);
#else
*volume = cvolume_ (a,b,m,n,mmax,nmax);
#endif
}

/*--------------------------------------------------------------

			ROUTINE: cvolumeb

    This routine calculates the `volume' in dimension n of the region
    bounded by a set of m linear inequality constraints of the form
    A x <= b, where a has m rows and n columns and is given by a(m,n),
    b is the n-vector and is contained in b(m). The recursive formula
    of Lasserre (1983) is used. Redundant constraints are allowed and
    a warning is issued if any are encountered. If the inequality 
    constraints are inconsistent then the volume is returned as zero.
    If the original polyhedron is unbounded then a warning is issued
    and the volume is return as zero (even though it is undefined). 

    This version also calculates the derivative of the volume with
    respect to the parameter b(0) using the simple formula of
    Lasserre (1983). If *opt == 2 then only the derivative is 
    calculated and not the volume.

    This routine is a variation from the routine `cvolume'.

    Calls are made to routine cvolume.

				Malcolm Sambridge, March 1995.


  --------------------------------------------------------------*/

#ifndef _UNDERSCORE
float cvolumeb(a,b,m,n,mmax,nmax,opt,dvdb)
#else
float cvolumeb_ (a,b,m,n,mmax,nmax,opt,dvdb)
#endif

int   *n, *m, *opt, *mmax, *nmax;
float *a, *b; 
float *dvdb; 

{

float v,amax,pivot;
int   i,j,t,k,l;
int   jj,kk;
float   *ai, *aj, *ajt, *apjj, *bi;
int   kmm,tmm;
int   nn = *n, mm = *m, mm_max = *mmax;
int   nm1,mm1;
int   lmin,lmax;
float  *ap, *bp;
int   firstmin,firstmax;
float bmin,bmax,bb;
float vol;

/* one-dimensional case (full reduction) */

if (nn == 1) 
  {
  firstmin=0;
  firstmax=0;
  lmax=0;
  lmin=0;
  for (l=0;l<mm;l++)
    {
    if ( *(a+l) > 0.) 
      {
      bb= *(b+l)/ *(a+l);
      if (firstmin==0) {firstmin=1;bmin=bb;lmin=l;}
      else if (bb<bmin) {bmin=bb;lmin=l;}
      }
    else if ( *(a+l) < 0.)
      {
      bb= *(b+l)/ *(a+l);
      if (firstmax==0) {firstmax=1;bmax=bb;lmax=l;}
      else if (bb>bmax) {bmax=bb;lmax=l;}
      }
    else if ( *(b+l) < 0.) 
      { 
                              /* Constraints are inconsistent.  
                                 Set volume and derivative to zero. */

       printf("Inconsistent constraints found after reduction to n = 1 \n");
      *dvdb = 0.; 
      return(0.);
      }
    }
  v=0.;
  *dvdb = 0.; 
  if (firstmin*firstmax == 1) v=bmin-bmax;
  else 
    {
     printf("Volume is unbounded; volume returned as -1\n derivatives returned as zero\n");
     return(-1.0);
    }

  if (v<0.) {v=0.;*dvdb=0.;} 
  else if (v>0. && lmin == 0) *dvdb = 1. / *a;
  else if (v>0. && lmax == 0) *dvdb = -1. / *a;
  return(v);
  }

nm1=nn-1;
mm1=mm-1;
v=0.;

/*  perform main loop over constraints */

for (i=0;i<mm;i++)
  {
   ai=a+i;
   bi=b+i;

/* find largest pivot */

  amax=0.;
  for (j=0;j<nn;j++) 
    if (fabs( *(ai+j*mm_max)) >= amax) {amax= fabs( *(ai+j*mm_max)); t=j;}
  tmm=t*mm_max;
  pivot=*(ai+tmm);

/* finds contribution to v from this pivot (if not nil) */

  if (amax == 0.)
  {
                                 /* Constraint is inconsistent */

    if ( *(bi) < 0.) 
       {
        printf("Constraint %d is inconsistent\n",i+1); 
        return(0.);
       } 

                                 /* otherwise constraint is redundant */

    printf("Constraint %d is redundant\n",i+1); 
  }
  else
  {

/* allocate memory */

  ap = (float *) malloc(4*nm1*mm1);
  bp = (float *) malloc(4*mm1);

/* reduce a and b into ap and bp eliminating variable t and constraint i */

  jj=-1;
    for (j=0;j<mm;j++)
      if (j != i)
      {
      jj=jj+1;
      aj=a+j;
      ajt=aj+tmm;
      *(bp+jj)= *(b+j) - *(bi) * *(ajt) / pivot;
      apjj=ap+jj;
      kk=-1;
      for (k=0;k<nn;k++)
        if (k != t)
        {
        kk=kk+1;
        kmm=k*mm_max;
        *(apjj+kk*mm1)= *(aj+kmm)- *(ajt) * *(ai+kmm)/ pivot;
        }
      }
  
/* add contribution to volume from volume calculated in smaller dimension */

#ifndef _UNDERSCORE
  vol=cvolume(ap,bp, &mm1, &nm1, &mm1, &nm1);
#else
  vol=cvolume_(ap,bp, &mm1, &nm1, &mm1, &nm1);
#endif
  if(vol == -1.0)
     {
      *dvdb = 0.;
      return(-1.0);
     }
  v=v+ *(bi)/amax*vol/nn;

  free(ap);
  free(bp);

/* calculate derivatives for first constraint only */

                                    /* calculate volume and derivative 
                                       with respect to b_0 */
  if (i == 0) *dvdb = vol/amax; 
                                    /* calculate derivative with respect 
                                       to b_0 but not volume */
  if (*opt == 2) return (0.); 
  }
  }

return(v);

}

/*--------------------------------------------------------------

			ROUTINE: volumeb

   A dummy routine used to call cvolumeb from a fortran routine 

----------------------------------------------------------------*/

#ifndef _UNDERSCORE
void volumeb(a,b,m,n,mmax,nmax,opt,volume,dvdb)
#else
void volumeb_ (a,b,m,n,mmax,nmax,opt,volume,dvdb)
#endif
int   *n, *m, *mmax, *nmax, *opt;
float *a, *b;
float *volume; 
float *dvdb;
{

if(*opt == 1)                      /* calculate volume and derivative with
                                       respect to b(0) */
  {
#ifndef _UNDERSCORE
   *volume=cvolumeb(a,b,m,n,mmax,nmax,opt,dvdb);
#else
   *volume=cvolumeb_(a,b,m,n,mmax,nmax,opt,dvdb);
#endif
  }
else if(*opt == 2)                /* calculate derivative with respect to b(0) 
                                       but not volume */
  {
#ifndef _UNDERSCORE
   *volume=cvolumeb(a,b,m,n,mmax,nmax,opt,dvdb);
#else
   *volume=cvolumeb_(a,b,m,n,mmax,nmax,opt,dvdb);
#endif
  }
else 
  {
   printf(" Warning: routine volumeb called with an invalid option parameter\n");
   printf("          valid options are 1 and 2, option given = %d \n",*opt);
  } 

}

/*--------------------------------------------------------------

			ROUTINE: cvolumebj

    This routine only calculates the j-th outer loop of the recursive
    routine cvolumeb. It returns the partial volume for projection 
    onto the j-th constraint and the derivative of the total volume 
    with respect to parameter b(j) using the simple formula of 
    Lasserre(1983). (See cvolumeb for more details.)
    
    The reason for this routine is so that the derivatives
    with respect to parameter b(j) can be calculated and passed back
    to a fortran routine without creating and an extra array of size
    m (because only one derivative is calculated per call). In
    this way it also avoids recalculating the total volume m times
    (which would be the case if we used a simple variation of routine
     cvolumeb).

    To calculate the total volume this routine must be called m times
    and the partial volume summed

    Constraint j is determined by the value of *con.

    Calls are made to routine cvolume.

				Malcolm Sambridge, March 1995.

  --------------------------------------------------------------*/

#ifndef _UNDERSCORE
float cvolumebj(a,b,m,n,mmax,nmax,con,dvdb)
#else
float cvolumebj_ (a,b,m,n,mmax,nmax,con,dvdb)
#endif
int   *n, *m, *mmax, *nmax, con;
float *a, *b; 
float *dvdb; 
{

float v,amax,pivot;
int   i,j,t,k,l;
int   jj,kk;
float   *ai, *aj, *ajt, *apjj, *bi;
int   kmm,tmm;
int   nn = *n,mm = *m, mm_max = *mmax;
int   nm1,mm1;
int   lmin,lmax;
float  *ap, *bp;
int   firstmin,firstmax;
float bmin,bmax,bb;
float vol;

/* one-dimensional case (full reduction) */

if (nn == 1) 
  {
  firstmin=0;
  firstmax=0;
  lmax=0;
  lmin=0;
  for (l=0;l<mm;l++) 
    {
    if ( *(a+l) > 0.) 
      {
      bb= *(b+l)/ *(a+l);
      if (firstmin==0) {firstmin=1;bmin=bb;lmin=l;}
      else if (bb<bmin) {bmin=bb;lmin=l;}
      }
    else if ( *(a+l) < 0.)
      {
      bb= *(b+l)/ *(a+l);
      if (firstmax==0) {firstmax=1;bmax=bb;lmax=l;}
      else if (bb>bmax) {bmax=bb;lmax=l;}
      }
    else if ( *(b+l) < 0.) 
      { 
                              /* Constraints are inconsistent.  
                                 Set volume and derivative to zero. */

       printf("Inconsistent constraints found after reduction to n = 1 \n");
      *dvdb = 0.; 
      return(0.);
      }
    }
  v=0.;
  *dvdb = 0.; 
  if (firstmin*firstmax == 1) v=bmin-bmax;
  else 
    {
     printf("Volume is unbounded; volume returned as -1\n derivatives returned as zero\n");
     return(-1.0);
    }

  if (v<0.) {v=0.;*dvdb=0.;} 
  else if (v>0. && lmin == con) *dvdb = 1. / *(a+lmin);
  else if (v>0. && lmax == con) *dvdb = -1. / *(a+lmax);
  return(v);
  }

nm1=nn-1;
mm1=mm-1;
v=0.;

/*  perform main loop over constraints */

/* for (i=0;i<mm;i++) */

  i = con;

   ai=a+i;
   bi=b+i;

/* find largest pivot */

  amax=0.;
  for (j=0;j<nn;j++) 
    if (fabs( *(ai+j*mm_max)) >= amax) {amax= fabs( *(ai+j*mm_max)); t=j;}
  tmm=t*mm_max;
  pivot=*(ai+tmm);

/* finds contribution to v from this pivot (if not nil) */

  if (amax == 0.)
  {
                                 /* Constraint is inconsistent */

    if ( *(bi) < 0.) 
       {
        printf("Constraint %d is inconsistent\n",i+1); 
        *dvdb = 0.; 
        return(0.);
       } 

                                 /* otherwise constraint is redundant */

    printf("Constraint %d is redundant\n",i+1); 
    *dvdb = 0.; 
  }
  else
  {

/* allocate memory */

  ap = (float *) malloc(4*nm1*mm1);
  bp = (float *) malloc(4*mm1);

/* reduce a and b into ap and bp eliminating variable t and constraint i */

  jj=-1;
    for (j=0;j<mm;j++)
      if (j != i)
      {
      jj=jj+1;
      aj=a+j;
      ajt=aj+tmm;
      *(bp+jj)= *(b+j) - *(bi) * *(ajt) / pivot;
      apjj=ap+jj;
      kk=-1;
      for (k=0;k<nn;k++)
        if (k != t)
        {
        kk=kk+1;
        kmm=k*mm_max;
        *(apjj+kk*mm1)= *(aj+kmm)- *(ajt) * *(ai+kmm)/ pivot;
        }
      }
  
                                    /* calculate partial volume */ 

#ifndef _UNDERSCORE
  vol=cvolume(ap,bp, &mm1, &nm1, &mm1, &nm1);
#else
  vol=cvolume_(ap,bp, &mm1, &nm1, &mm1, &nm1);
#endif
  if(vol == -1.0)
     {
      *dvdb = 0.;
      return(-1.0);
     }
  v=*(bi)/amax*vol/nn;

                                    /* calculate derivative of total volume
                                       with respect to current constraint */
  *dvdb = vol/amax; 

  free(ap);
  free(bp);
  }

return(v);

}

/*--------------------------------------------------------------

			ROUTINE: volumebj

   A dummy routine used to call cvolumebj from a fortran routine

----------------------------------------------------------------*/

#ifndef _UNDERSCORE
void volumebj(a,b,m,n,mmax,nmax,con,volume,dvdb)
#else
void volumebj_ (a,b,m,n,mmax,nmax,con,volume,dvdb)
#endif
int   *n, *m, *mmax, *nmax,*con;
float *a, *b;
float *volume; 
float *dvdb; 
{

int tcon;
tcon = *con - 1;
                                    /* calculate partial volume and 
                                       derivative with respect to b(con) */
#ifndef _UNDERSCORE
   *volume=cvolumebj(a,b,m,n,mmax,nmax,tcon,dvdb);
#else
   *volume=cvolumebj_(a,b,m,n,mmax,nmax,tcon,dvdb);
#endif

}
/*--------------------------------------------------------------

			ROUTINE: cdvda

    This routine calculates the derivative with respect to a(0,tdim)
    of the `volume' in dimension n of the region bounded by a set 
    of m linear inequality constraints of the form A x <= b, where 
    a has m rows and n columns and is given by a(m,n), b is the 
    n-vector and is contained in b(m). The derivative expression
    is recursive and derived from the formula of Lasserre (1983). 

    Redundant constraints are allowed and a warning is issued if any are 
    encountered. If the inequality constraints are inconsistent then the 
    derivative is returned as zero.  If any constraint is orthogonal to 
    the component a(0,idim) then the reduction can only take place onto 
    variable idim. A special case is used to handle this which involves
    no further recursive calls.

    If the original polyhedron is unbounded then a warning is issued
    and the derivative is return as zero. 

    Note: This code takes advantage of the fact that during recursive calls
    constraint 0 does not change its position in the list of remaining 
    constraints (if it has not been eliminated), i.e. it is always the 
    first constraint.  This would not be the case if the algorithm were 
    adapated to deal with other constraints, i.e. evaluate dvda_i,j 
    where i .ne. 0.

    Calls itself cvolumeb, and cvolumebj.
 
				Malcolm Sambridge, March 1995.

  --------------------------------------------------------------*/
#ifndef _UNDERSCORE
float cdvda(a,b,m,n,mmax,nmax,tdim,temp,jval,code)
#else
float cdvda_ (a,b,m,n,mmax,nmax,tdim,temp,jval,code)
#endif
int   *n, *m, *nmax, *mmax, *tdim, *jval, *code;
float *a, *b, *temp; 
{

float v,amax,pivot;
int   i,j,t,k,l;
int   jj,kk;
float   *ai, *aj, *ajt, *apjj, *bi;
int   kmm,tmm;
int   nn = *n, mm = *m, mm_max = *mmax, ttdim;
int   nm1,mm1;
int   lmin,lmax, jjval, kval;
float  *ap, *bp, *ttemp;
int   firstmin,firstmax;
float bmin,bmax,bb;
float deriv, junk, vol, dvdb, dbda;
int   special, opt;

/* one-dimensional case (full reduction) */

*code = 0;

if (nn == 1) 
  {
  firstmin=0;
  firstmax=0;
  lmax=0;
  lmin=0;
  for (l=0;l<mm;l++)
    {
    if ( *(a+l) > 0.) 
      {
      bb= *(b+l)/ *(a+l);
      if (firstmin==0) {firstmin=1;bmin=bb;lmin=l;}
      else if (bb<bmin) {bmin=bb;lmin=l;}
      }
    else if ( *(a+l) < 0.)
      {
      bb= *(b+l)/ *(a+l);
      if (firstmax==0) {firstmax=1;bmax=bb;lmax=l;}
      else if (bb>bmax) {bmax=bb;lmax=l;}
      }
    else if ( *(b+l) < 0.) 
      {
                                   /* Constraint is inconsistent.
                                      Set derivative to zero. */

       printf("Inconsistent constraints found after reduction to n = 1 \n");
       *code = 1;
       return(0.);
      }
    }
  v=0.;
  if (firstmin*firstmax == 1) v=bmin-bmax;
  else 
    {
     printf("Volume is unbounded; derivative returned is zero\n");
     *code = -1;
     return(0.);
    }

  if (v<0.) return(0.);       

  if(*jval == 1)           /* Constraint 0 has not yet been encountered */
    {
     if (lmin == 0) deriv = -bmin/ *a;
     else if (lmax == 0) deriv = bmax/ *a;
     else deriv = 0.;
     return(deriv); 
    }
  else if(*jval == 0)     /* Constraint 0  has already been encountered */
    {
     deriv =  ( *(temp+lmax) * (bmax/ *(a+lmax)) ) -
              ( *(temp+lmin) * (bmin/ *(a+lmin)) );
     return(deriv); 
    }
  }
nm1=nn-1;
mm1=mm-1;
v=0.;
 
                                 /*  perform main loop over constraints */

for (i=0;i<mm;i++)
  {
   ai=a+i;
   bi=b+i;
   ttdim = *tdim;
   special = 0;

/* find largest pivot */

  amax=0.;
  t = 0;
  for (j=0;j<nn;j++) 
    if (fabs( *(ai+j*mm_max)) >= amax && j != ttdim) 
        {amax= fabs( *(ai+j*mm_max)); t=j;}

                                 /* finds contribution to v from 
                                    this pivot (if not nil) */

  if (amax == 0.)
  {
   if(*(ai + ttdim * mm_max) == 0.0) 
     { 
                                 /* Constraint is inconsistent */
       if ( *(bi) < 0.)
         {
          printf("Constraint %d is inconsistent\n",i+1);
          *code = 1;
          return(0.);
         }

                                 /* otherwise constraint is redundant */

       printf("Constraint %d is redundant\n",i+1); 
     }
   else
     {
                                 /* if projection can only be peformed
                                    on dimension tdim then activate 
                                    special case */ 

     special = 1;
     t = ttdim;
     amax = fabs(*(ai+t * mm_max));

     }
  }

  tmm=t*mm_max;
  pivot=*(ai+tmm);

  if(t < ttdim) ttdim = ttdim -1;


  if (amax != 0)
  {

                 /* determine if constraint 0 has been encountered */
 
   kval = 0;
   if ( i == 0 && *jval == 1)
      {
                                     /* This is the first encounter of
                                        constraint 0 on this path so we 
                                        allocate memory and store parameters 
                                        to be used when n = 1 */

       if (special == 0) 
          {
           ttemp = (float *) malloc(4*mm1);
           for (j=0;j<mm1;j++) *(ttemp+j) = - *(a+j+1+tmm)/pivot;
           kval = 1;
          }
       jjval = 0;
      }
   else if (*jval == 0)             /* Constraint 0 has already been
                                        encountered */
      {
        jjval = 0;  
                                    /* perform recursive update of component
                                       derivative array temp. This eliminates 
                                       row i and copies into a new vector */

        if(special == 0)
          {
           ttemp = (float *) malloc(4*mm1);
           for (j=0;j<i;j++) *(ttemp+j) = *(temp+j) 
                                          -(*(temp+i) * *(a+j+tmm)/pivot);
           for (j=i;j<mm1;j++) *(ttemp+j) = *(temp+j+1)
                                          -(*(temp+i) * *(a+j+1+tmm)/pivot);
          }
      }
   else 
      {                             /* Constraint 0 has not yet been 
                                        encountered */
      jjval = 1;                
      }


/* allocate memory */

  ap = (float *) malloc(4*nm1*mm1);
  bp = (float *) malloc(4*mm1);

/* reduce a and b into ap and bp eliminating variable t and constraint i */

  jj=-1;
    for (j=0;j<mm;j++)
      if (j != i)
      {
      jj=jj+1;
      aj=a+j;
      ajt=aj+tmm;
      *(bp+jj)= *(b+j) - *(bi) * *(ajt) / pivot;
      apjj=ap+jj;
      kk=-1;
      for (k=0;k<nn;k++)
        if (k != t)
        {
        kk=kk+1;
        kmm=k*mm_max;
        *(apjj+kk*mm1)= *(aj+kmm)- *(ajt) * *(ai+kmm)/ pivot;
        }
      }
  
/* add contribution to derivative from that calculated in smaller dimension */

  					/* Normal case method */
  if(special == 0)
    {
#ifndef _UNDERSCORE
     deriv=cdvda(ap,bp, &mm1, &nm1, &mm1, &nm1, &ttdim, ttemp, &jjval,code);
#else
     deriv=cdvda_(ap,bp, &mm1, &nm1, &mm1, &nm1, &ttdim, ttemp, &jjval,code);
#endif
     v=v+ *(bi)/amax*deriv/nn;
     if (kval == 1 || *jval == 0) free(ttemp);
     if(*code != 0)return (0.);

    }
  else					/* Use special case method */
    {
     if( *jval == 1)                    
       {                     
        if(i == 0)                      /* This is constraint 0 */
          {
           deriv = 0.; 
           vol = 0.; 
           dvdb = 0.;
           for (j=1;j<mm;j++) 
              {
               k = j - 1;
#ifndef _UNDERSCORE
               junk=cvolumebj(ap,bp,&mm1,&nm1,&mm1,&nm1,k,&dvdb);
#else
               junk=cvolumebj_(ap,bp,&mm1,&nm1,&mm1,&nm1,k,&dvdb);
#endif
               if(junk == -1.)
                 {
                  *code = -1;
                  return(0.);
                 }
               deriv = deriv + dvdb * *(a + j + tmm) ;
               vol = vol + junk;
              }
           if(nm1 == 1)vol = junk;
           deriv = *(bi) * deriv/pivot;
           deriv = (deriv - vol) /pivot; 
           v=v+ *(bi)/amax*deriv/nn;

          }
        else				/* Constraint 0 not yet encountered */
          {
           opt = 2;
#ifndef _UNDERSCORE
           junk=cvolumeb(ap,bp,&mm1,&nm1,&mm1,&nm1,&opt,&deriv);
#else
           junk=cvolumeb_(ap,bp,&mm1,&nm1,&mm1,&nm1,&opt,&deriv);
#endif
           if(junk == -1.)
             {
              *code = -1;
              return(0.);
             }
           deriv= -(deriv *  *(bi)/pivot);
           v=v+ *(bi)/amax*deriv/nn;
          }
       }
     else if( *jval == 0)               /* Constraint 0 already encountered */
       {
        vol = 0.; 
        deriv = 0.; 
        for (j=0;j<i;j++) 
           {
#ifndef _UNDERSCORE
            junk=cvolumebj(ap,bp,&mm1,&nm1,&mm1,&nm1,j,&dvdb);
#else
            junk=cvolumebj_(ap,bp,&mm1,&nm1,&mm1,&nm1,j,&dvdb);
#endif
            if(junk == -1.)
              {
               *code = -1;
               return(0.);
              }
            dbda = (*(a+j+tmm) * *(temp+i)/pivot) - *(temp+j);
            deriv = deriv + dvdb*dbda;
            vol = vol + junk;
           }
        for (j=i+1;j<mm;j++) 
           {
            k = j - 1;
#ifndef _UNDERSCORE
            junk=cvolumebj(ap,bp,&mm1,&nm1,&mm1,&nm1,k,&dvdb);
#else
            junk=cvolumebj_(ap,bp,&mm1,&nm1,&mm1,&nm1,k,&dvdb);
#endif
            if(junk == -1.)
              {
               *code = -1;
               return(0.);
              }
            dbda = (*(a+j+tmm) * *(temp+i)/pivot) - *(temp+j);
            deriv = deriv + dvdb*dbda;
            vol = vol + junk;
           }
           if(nm1 == 1)vol = junk;
           deriv = (*(bi) * deriv - *(temp+i)*vol)/pivot;
           v=v+ *(bi)/amax*deriv/nn;
       }
    }

  free(ap);
  free(bp);
  }
  }

return(v);

}

/*--------------------------------------------------------------

			ROUTINE: dvda

   A dummy routine used to call cdvda from a fortran routine 

----------------------------------------------------------------*/
#ifndef _UNDERSCORE
void dvda(a,b,m,n,mmax,nmax,idim,dvda,code)
#else
void dvda_ (a,b,m,n,mmax,nmax,idim,dvda,code)
#endif
int   *n, *m, *mmax, *nmax, *idim, *code;
float *a, *b;
float *dvda; 
{

int    jval, tdim;
float *temp = 0;

jval = 1;
tdim = *idim - 1;
#ifndef _UNDERSCORE
*dvda=cdvda(a,b,m,n,mmax,nmax,&tdim,temp,&jval,code);
#else
*dvda=cdvda_(a,b,m,n,mmax,nmax,&tdim,temp,&jval,code);
#endif
free(temp);
}

/*--------------------------------------------------------------

			ROUTINE: cvolumef

    This routine calculates the `volume' in dimension n of the region
    bounded by a set of m linear inequality constraints of the form
    A x <= b, where a has m rows and n columns and is given by a(m,n),
    b is the n-vector and is contained in b(m). The recursive formula
    of Lasserre (1983) is used. Redundant constraints are allowed 
    If the inequality constraints are inconsistent then the volume 
    is returned as zero.  If the original polyhedron is unbounded 
    then a warning is issued and the volume is return as -1.0 
    (even though it is undefined). 

    This is a faster version which does not continue down the
    recursive tree if if encounters b(j) = 0 for any j. This
    restricts its use to cases where the origin does not
    pass through any inconsistent constraints. Also redundant
    constraints that pass through the origin will not be detected.
 
			Jean Braun and Malcolm Sambridge, Jan 1995.


    Lasserre, J. B., 1983. "An analytical expression and an Algorithm
    for the Volume of a Convex Polyhedron in R^n.", JOTA, vol. 39, 
    no. 3, 363-377.

  --------------------------------------------------------------*/

#ifndef _UNDERSCORE
float cvolumef(a,b,m,n,mmax,nmax)
#else
float cvolumef_ (a,b,m,n,mmax,nmax)
#endif
int   *n, *m, *mmax, *nmax;
float *a, *b;
{

float v,amax,pivot;
int   i,j,t,k,l;
int   jj,kk;
float   *ai, *aj, *ajt, *apjj, *bi;
int   kmm,tmm;
int   nn = *n, mm = *m, mm_max = *mmax;
int   nm1,mm1;
float  *ap, *bp;
int   firstmin,firstmax;
float bmin,bmax,bb;
float partialv;

/* one-dimensional case (full reduction) */

if (nn == 1) 
  {
  firstmin=0;
  firstmax=0;
  for (l=0;l<mm;l++)
    {
    if ( *(a+l) > 0.) 
      {
      bb= *(b+l)/ *(a+l);
      if (firstmin==0) {firstmin=1;bmin=bb;}
      else if (bb<bmin) bmin=bb;
      }
    else if ( *(a+l) < 0.)
      {
      bb= *(b+l)/ *(a+l);
      if (firstmax==0) {firstmax=1;bmax=bb;}
      else if (bb>bmax) bmax=bb;
      }
    else if ( *(b+l) < 0.)         /* Constraints are inconsistent */
      { 
       printf("Inconsistent constraints found after reduction to n = 1 \n");
       return(0.); 
      }
    }
  v=0.;
  if (firstmin*firstmax == 1) v=bmin-bmax;
  else 
    {
     printf("Volume is unbounded at 1-D; volume returned as -1\n");
     return(-1.0);
    }

  if (v<0.) {v=0.;} 

  return(v);
  }

nm1=nn-1;
mm1=mm-1;
v=0.;

for (i=0;i<mm;i++)
  {
   ai=a+i;
   bi=b+i;

/* find largest pivot */

  amax=0.;
  for (j=0;j<nn;j++) 
    if (fabs( *(ai+j*mm_max)) >= amax) {amax= fabs( *(ai+j*mm_max)); t=j;}
  tmm=t*mm_max;
  pivot=*(ai+tmm);

/* finds contribution to v from this pivot (if not nil) */

  if (amax == 0.)
  {
                                 /* Constraint is inconsistent */

    if ( *(bi) < 0.) 
      { printf("Constraint %d is inconsistent\n",i+1); return(0.); }

                                 /* otherwise constraint is redundant */

    printf("Constraint %d is redundant\n",i+1); 
  }
  else
  {

/* allocate memory */

  ap = (float *) malloc(4*nm1*mm1);
  bp = (float *) malloc(4*mm1);

/* reduce a and b into ap and bp eliminating variable t and constraint i */

  jj=-1;
    for (j=0;j<mm;j++)
      if (j != i)
      {
      jj=jj+1;
      aj=a+j;
      ajt=aj+tmm;
      *(bp+jj)= *(b+j) - *(bi) * *(ajt) / pivot;
      apjj=ap+jj;
      kk=-1;
      for (k=0;k<nn;k++)
        if (k != t)
        {
        kk=kk+1;
        kmm=k*mm_max;
        *(apjj+kk*mm1)= *(aj+kmm)- *(ajt) * *(ai+kmm)/ pivot;
        }
      }
  
/* add contribution to volume from volume calculated in smaller dimension */

  partialv = 0.;
  if (*(bi) != 0.)
#ifndef _UNDERSCORE
    partialv=cvolumef(ap,bp, &mm1, &nm1, &mm1, &nm1);
#else
    partialv=cvolumef_(ap,bp, &mm1, &nm1, &mm1, &nm1);
#endif
  if(partialv == -1.0)return(-1.0);
  v=v+ *(bi)/amax*partialv/nn; 

  free(ap);
  free(bp);
  }
  }

return(v);

}
/*--------------------------------------------------------------

			ROUTINE: volume

   A dummy routine used to call cvolume from a fortran routine 

----------------------------------------------------------------*/

#ifndef _UNDERSCORE
void volumef(a,b,m,n,mmax,nmax,volume)
#else
void volumef_ (a,b,m,n,mmax,nmax,volume)
#endif
int   *n, *m, *mmax, *nmax;
float *a, *b;
float *volume; 
{
#ifndef _UNDERSCORE
*volume=cvolumef(a,b,m,n,mmax,nmax);
#else
*volume=cvolumef_(a,b,m,n,mmax,nmax);
#endif

}
/*--------------------------------------------------------------

			ROUTINE: cdvdaf

    This routine calculates the derivative with respect to a(0,tdim)
    of the `volume' in dimension n of the region bounded by a set 
    of m linear inequality constraints of the form A x <= b, where 
    a has m rows and n columns and is given by a(m,n), b is the 
    n-vector and is contained in b(m). The derivative expression
    is recursive and derived from the formula of Lasserre (1983). 

    Redundant constraints are allowed and a warning is issued if any are 
    encountered. If the inequality constraints are inconsistent then the 
    derivative is returned as zero.  If any constraint is orthogonal to 
    the component a(0,idim) then the reduction can only take place onto 
    variable idim. A special case is used to handle this which involves
    no further recursive calls.

    If the original polyhedron is unbounded then a warning is issued
    and the derivative is return as zero. 

    Note: This code takes advantage of the fact that during recursive calls
    constraint 0 does not change its position in the list of remaining 
    constraints (if it has not been eliminated), i.e. it is always the 
    first constraint.  This would not be the case if the algorithm were 
    adapated to deal with other constraints, i.e. evaluate dvda_i,j 
    where i .ne. 0.

    This is a faster version which does not continue down the
    recursive tree if if encounters b(j) = 0 for any j. This
    restricts its use to cases where the origin does not
    pass through any inconsistent constraints. Also redundant
    constraints that pass through the origin will not be detected.
 
    Calls itself cvolumeb, and cvolumebj.
 
				Malcolm Sambridge, March 1995.

  --------------------------------------------------------------*/

#ifndef _UNDERSCORE
float cdvdaf(a,b,m,n,mmax,nmax,tdim,temp,jval,code)
#else
float cdvdaf_ (a,b,m,n,mmax,nmax,tdim,temp,jval,code)
#endif

int   *n, *m, *nmax, *mmax, *tdim, *jval, *code;
float *a, *b, *temp; 

{

float v,amax,pivot;
int   i,j,t,k,l;
int   jj,kk;
float   *ai, *aj, *ajt, *apjj, *bi;
int   kmm,tmm;
int   nn = *n, mm = *m, mm_max = *mmax, ttdim;
int   nm1,mm1;
int   lmin,lmax, jjval, kval;
float  *ap, *bp, *ttemp;
int   firstmin,firstmax;
float bmin,bmax,bb;
float deriv, junk, vol, dvdb, dbda;
int   special, opt;

/* one-dimensional case (full reduction) */

*code = 0;

if (nn == 1) 
  {
  firstmin=0;
  firstmax=0;
  lmax=0;
  lmin=0;
  for (l=0;l<mm;l++)
    {
    if ( *(a+l) > 0.) 
      {
      bb= *(b+l)/ *(a+l);
      if (firstmin==0) {firstmin=1;bmin=bb;lmin=l;}
      else if (bb<bmin) {bmin=bb;lmin=l;}
      }
    else if ( *(a+l) < 0.)
      {
      bb= *(b+l)/ *(a+l);
      if (firstmax==0) {firstmax=1;bmax=bb;lmax=l;}
      else if (bb>bmax) {bmax=bb;lmax=l;}
      }
    else if ( *(b+l) < 0.) 
      {
                                   /* Constraint is inconsistent.
                                      Set derivative to zero. */

       printf("Inconsistent constraints found after reduction to n = 1 \n");
       *code = 1;
       return(0.);
      }
    }
  v=0.;
  if (firstmin*firstmax == 1) v=bmin-bmax;
  else 
    {
     printf("Volume is unbounded; derivative returned is zero\n");
     *code = -1;
     return(0.);
    }

  if (v<0.) return(0.);       

  if(*jval == 1)           /* Constraint 0 has not yet been encountered */
    {
     if (lmin == 0) deriv = -bmin/ *a;
     else if (lmax == 0) deriv = bmax/ *a;
     else deriv = 0.;
     return(deriv); 
    }
  else if(*jval == 0)     /* Constraint 0  has already been encountered */
    {
     deriv =  ( *(temp+lmax) * (bmax/ *(a+lmax)) ) -
              ( *(temp+lmin) * (bmin/ *(a+lmin)) );
     return(deriv); 
    }
  }
nm1=nn-1;
mm1=mm-1;
v=0.;
 
                                 /*  perform main loop over constraints */

for (i=0;i<mm;i++)
  {
   ai=a+i;
   bi=b+i;
   ttdim = *tdim;
   special = 0;

/* find largest pivot */

  amax=0.;
  t = 0;
  for (j=0;j<nn;j++) 
    if (fabs( *(ai+j*mm_max)) >= amax && j != ttdim) 
        {amax= fabs( *(ai+j*mm_max)); t=j;}

                                 /* finds contribution to v from 
                                    this pivot (if not nil) */

  if (amax == 0.)
  {
   if(*(ai + ttdim * mm_max) == 0.0) 
     { 
                                 /* Constraint is inconsistent */
       if ( *(bi) < 0.)
         {
          printf("Constraint %d is inconsistent\n",i+1);
          *code = 1;
          return(0.);
         }

                                 /* otherwise constraint is redundant */

       printf("Constraint %d is redundant\n",i+1); 
     }
   else
     {
                                 /* if projection can only be peformed
                                    on dimension tdim then activate 
                                    special case */ 

     special = 1;
     t = ttdim;
     amax = fabs(*(ai+t * mm_max));

     }
  }

  tmm=t*mm_max;
  pivot=*(ai+tmm);

  if(t < ttdim) ttdim = ttdim -1;


  if (amax != 0)
  {

                 /* determine if constraint 0 has been encountered */
 
   kval = 0;
   if ( i == 0 && *jval == 1)
      {
                                     /* This is the first encounter of
                                        constraint 0 on this path so we 
                                        allocate memory and store parameters 
                                        to be used when n = 1 */

       if (special == 0) 
          {
           ttemp = (float *) malloc(4*mm1);
           for (j=0;j<mm1;j++) *(ttemp+j) = - *(a+j+1+tmm)/pivot;
           kval = 1;
          }
       jjval = 0;
      }
   else if (*jval == 0)             /* Constraint 0 has already been
                                        encountered */
      {
        jjval = 0;  
                                    /* perform recursive update of component
                                       derivative array temp. This eliminates 
                                       row i and copies into a new vector */

        if(special == 0)
          {
           ttemp = (float *) malloc(4*mm1);
           for (j=0;j<i;j++) *(ttemp+j) = *(temp+j) 
                                          -(*(temp+i) * *(a+j+tmm)/pivot);
           for (j=i;j<mm1;j++) *(ttemp+j) = *(temp+j+1)
                                          -(*(temp+i) * *(a+j+1+tmm)/pivot);
          }
      }
   else 
      {                             /* Constraint 0 has not yet been 
                                        encountered */
      jjval = 1;                
      }


/* allocate memory */

  ap = (float *) malloc(4*nm1*mm1);
  bp = (float *) malloc(4*mm1);

/* reduce a and b into ap and bp eliminating variable t and constraint i */

  jj=-1;
    for (j=0;j<mm;j++)
      if (j != i)
      {
      jj=jj+1;
      aj=a+j;
      ajt=aj+tmm;
      *(bp+jj)= *(b+j) - *(bi) * *(ajt) / pivot;
      apjj=ap+jj;
      kk=-1;
      for (k=0;k<nn;k++)
        if (k != t)
        {
        kk=kk+1;
        kmm=k*mm_max;
        *(apjj+kk*mm1)= *(aj+kmm)- *(ajt) * *(ai+kmm)/ pivot;
        }
      }
  
/* add contribution to derivative from that calculated in smaller dimension */

  					/* Normal case method */
  if(special == 0)
    {
     deriv = 0.;
     if (*(bi) != 0.)
#ifndef _UNDERSCORE
     {deriv=cdvdaf(ap,bp, &mm1, &nm1, &mm1, &nm1, &ttdim, ttemp, &jjval,code);}
#else
     {deriv=cdvdaf_(ap,bp, &mm1, &nm1, &mm1, &nm1, &ttdim, ttemp, &jjval,code);}
#endif
     v=v+ *(bi)/amax*deriv/nn;
     if (kval == 1 || *jval == 0) free(ttemp);
     if(*code != 0)return (0.);

    }
  else					/* Use special case method */
    {
     if( *jval == 1)                    
       {                     
        if(i == 0)                      /* This is constraint 0 */
          {
           deriv = 0.; 
           vol = 0.; 
           dvdb = 0.;
           for (j=1;j<mm;j++) 
              {
               k = j - 1;
#ifndef _UNDERSCORE
               junk=cvolumebj(ap,bp,&mm1,&nm1,&mm1,&nm1,k,&dvdb);
#else
               junk=cvolumebj_(ap,bp,&mm1,&nm1,&mm1,&nm1,k,&dvdb);
#endif
               if(junk == -1.)
                 {
                  *code = -1;
                  return(0.);
                 }
               deriv = deriv + dvdb * *(a + j + tmm) ;
               vol = vol + junk;
              }
           if(nm1 == 1)vol = junk;
           deriv = *(bi) * deriv/pivot;
           deriv = (deriv - vol) /pivot; 
           v=v+ *(bi)/amax*deriv/nn;

          }
        else				/* Constraint 0 not yet encountered */
          {
           opt = 2;
#ifndef _UNDERSCORE
           junk=cvolumeb(ap,bp,&mm1,&nm1,&mm1,&nm1,&opt,&deriv);
#else
           junk=cvolumeb_(ap,bp,&mm1,&nm1,&mm1,&nm1,&opt,&deriv);
#endif
           if(junk == -1.)
             {
              *code = -1;
              return(0.);
             }
           deriv= -(deriv *  *(bi)/pivot);
           v=v+ *(bi)/amax*deriv/nn;
          }
       }
     else if( *jval == 0)               /* Constraint 0 already encountered */
       {
        vol = 0.; 
        deriv = 0.; 
        for (j=0;j<i;j++) 
           {
#ifndef _UNDERSCORE
            junk=cvolumebj(ap,bp,&mm1,&nm1,&mm1,&nm1,j,&dvdb);
#else
            junk=cvolumebj_(ap,bp,&mm1,&nm1,&mm1,&nm1,j,&dvdb);
#endif
            if(junk == -1.)
              {
               *code = -1;
               return(0.);
              }
            dbda = (*(a+j+tmm) * *(temp+i)/pivot) - *(temp+j);
            deriv = deriv + dvdb*dbda;
            vol = vol + junk;
           }
        for (j=i+1;j<mm;j++) 
           {
            k = j - 1;
#ifndef _UNDERSCORE
            junk=cvolumebj(ap,bp,&mm1,&nm1,&mm1,&nm1,k,&dvdb);
#else
            junk=cvolumebj_(ap,bp,&mm1,&nm1,&mm1,&nm1,k,&dvdb);
#endif
            if(junk == -1.)
              {
               *code = -1;
               return(0.);
              }
            dbda = (*(a+j+tmm) * *(temp+i)/pivot) - *(temp+j);
            deriv = deriv + dvdb*dbda;
            vol = vol + junk;
           }
           if(nm1 == 1)vol = junk;
           deriv = (*(bi) * deriv - *(temp+i)*vol)/pivot;
           v=v+ *(bi)/amax*deriv/nn;
       }
    }

  free(ap);
  free(bp);
  }
  }

return(v);

}

/*--------------------------------------------------------------

			ROUTINE: dvda

   A dummy routine used to call cdvda from a fortran routine 

----------------------------------------------------------------*/

#ifndef _UNDERSCORE
void dvdaf(a,b,m,n,mmax,nmax,idim,dvda,code)
#else
void dvdaf_ (a,b,m,n,mmax,nmax,idim,dvda,code)
#endif
int   *n, *m, *mmax, *nmax, *idim, *code;
float *a, *b;
float *dvda;
{
int    jval, tdim;
float *temp = 0;

jval = 1;
tdim = *idim - 1;
#ifndef _UNDERSCORE
*dvda=cdvdaf(a,b,m,n,mmax,nmax,&tdim,temp,&jval,code);
#else
*dvda=cdvdaf_(a,b,m,n,mmax,nmax,&tdim,temp,&jval,code);
#endif
free(temp);
}


/*--------------------------------------------------------------

			ROUTINE: cdvda_debug

    This routine calculates the derivative with respect to a(0,tdim)
    of the `volume' in dimension n of the region bounded by a set 
    of m linear inequality constraints of the form A x <= b, where 
    a has m rows and n columns and is given by a(m,n), b is the 
    n-vector and is contained in b(m). The derivative expression
    is recursive and derived from the formula of Lasserre (1983). 

    Redundant constraints are allowed and a warning is issued if any are 
    encountered. If the inequality constraints are inconsistent then the 
    derivative is returned as zero.  If any constraint is orthogonal to 
    the component a(0,idim) then the reduction can only take place onto 
    variable idim. A special case is used to handle this which involves
    no further recursive calls.

    If the original polyhedron is unbounded then a warning is issued
    and the derivative is return as zero. 

    Note: This code takes advantage of the fact that during recursive calls
    constraint 0 does not change its position in the list of remaining 
    constraints (if it has not been eliminated), i.e. it is always the 
    first constraint.  This would not be the case if the algorithm were 
    adapated to deal with other constraints, i.e. evaluate dvda_i,j 
    where i .ne. 0.

    Calls itself cvolumeb, and cvolumebj.

    This is a version of cdvda used for debugging and contains extra
    code and print statements.
 
				Malcolm Sambridge, March 1995.

  --------------------------------------------------------------*/

#ifndef _UNDERSCORE
float cdvda_debug(a,b,m,n,tdim,temp,jval,code)
#else
float cdvda_debug_ (a,b,m,n,tdim,temp,jval,code)
#endif
int   *n, *m, *tdim, *jval, *code;
float *a, *b, *temp;
{

float v,amax,pivot;
int   i,j,t,k,l;
int   jj,kk;
float   *ai, *aj, *ajt, *apjj, *bi;
int   kmm,tmm;
int   nn,mm, ttdim;
int   nm1,mm1;
int   lmin,lmax, jjval, kval;
float  *ap, *bp, *ttemp;
int   firstmin,firstmax;
float bmin,bmax,bb;
float deriv, junk, vol, dvdb, dbda;
int   special, opt;
float volp, volm, da, FD, FDa;

nn= *n;
mm= *m;

/* write out a and b (debug) */

       printf("\n"); 
       printf(" Current matrix A and vector b \n");
   for (j=0;j<mm;j++)
      {
      aj=a+j;
      for (k=0;k<nn;k++) 
      {ajt=aj+k * mm; printf(" a %d %d = %f",j,k,*(ajt));}
/*      printf(" : b = %f \n", *(b+j));*/
      }

/* one-dimensional case (full reduction) */

*code = 0;

if (nn == 1) 
  {
  printf("One dimension left \n");
  firstmin=0;
  firstmax=0;
  lmax=0;
  lmin=0;
  for (l=0;l<mm;l++)
    {
    if ( *(a+l) > 0.) 
      {
      bb= *(b+l)/ *(a+l);
      if (firstmin==0) {firstmin=1;bmin=bb;lmin=l;}
      else if (bb<bmin) {bmin=bb;lmin=l;}
      }
    else if ( *(a+l) < 0.)
      {
      bb= *(b+l)/ *(a+l);
      if (firstmax==0) {firstmax=1;bmax=bb;lmax=l;}
      else if (bb>bmax) {bmax=bb;lmax=l;}
      }
    else if ( *(b+l) < 0.) 
      {
                                   /* Constraint is inconsistent.
                                      Set derivative to zero. */

       printf("Inconsistent constraints found after reduction to n = 1 \n");
       *code = 1;
       return(0.);
      }
    }
  v=0.;
  if (firstmin*firstmax == 1) v=bmin-bmax;
  else 
    {
     printf("Volume is unbounded; derivative returned is zero\n");
     *code = -1;
     return(0.);
    }

  if (v<0.) return(0.);       

  if(*jval == 1)           /* Constraint 0 has not yet been encountered */
    {
     if (lmin == 0) deriv = -bmin/ *a;
     else if (lmax == 0) deriv = bmax/ *a;
     else deriv = 0.;
     printf(" No projection onto constraint 0 \n");
     printf(" lmin = %d lmax = %d bmin = %f bmax = %f alpha = %f \n",
            lmin,lmax,bmin,bmax,*a);
     printf(" deriv contribution = %f \n",deriv);
     return(deriv); 
    }
  else if(*jval == 0)     /* Constraint 0  has already been encountered */
    {
     deriv =  ( *(temp+lmax) * (bmax/ *(a+lmax)) ) -
              ( *(temp+lmin) * (bmin/ *(a+lmin)) );
     printf(" Projection onto constraint 0 has occurred \n");
     printf(" lmin = %d bmin = %f alpha = %f temp = %f\n",
            lmin,bmin,*(a+lmin),*(temp+lmin));
     printf(" lmax = %d bmax = %f alpha = %f temp = %f\n",
            lmax,bmax,*(a+lmax),*(temp+lmax));
     printf(" deriv contribution = %f \n",deriv);
     return(deriv); 
    }
  }
nm1=nn-1;
mm1=mm-1;
v=0.;
 
printf(" About to start main loop n = %d m = %d \n",nn,mm);

                                 /*  perform main loop over constraints */

for (i=0;i<mm;i++)
  {
   ai=a+i;
   bi=b+i;
   ttdim = *tdim;
   special = 0;

/* find largest pivot */

  amax=0.;
  t = 0;
  for (j=0;j<nn;j++) 
    if (fabs( *(ai+j*mm)) >= amax && j != ttdim) 
        {amax= fabs( *(ai+j*mm)); t=j;}

                                 /* finds contribution to v from 
                                    this pivot (if not nil) */

  if (amax == 0.)
  {
   if(*(ai + ttdim * mm) == 0.0) 
     { 
                                 /* Constraint is inconsistent */
       if ( *(bi) < 0.)
         {
          printf("Constraint %d is inconsistent\n",i+1);
          *code = 1;
          return(0.);
         }

                                 /* otherwise constraint is redundant */

       printf("Constraint %d is redundant\n",i+1); 
     }
   else
     {
                                 /* if projection can only be peformed
                                    on dimension tdim then activate 
                                    special case */ 

     printf("Constraint %d can only be projected onto dimension %d\n",i+1,*tdim);
     printf("using special case method \n");
     special = 1;
     t = ttdim;
     amax = fabs(*(ai+t * mm));

     }
  }

  tmm=t*mm;
  pivot=*(ai+tmm);

  printf("\n Projection onto constraint %d variable %d k = %d\n",i,t,ttdim);

  if(t < ttdim) ttdim = ttdim -1;


  if (amax != 0)
  {

                 /* determine if constraint 0 has been encountered */
 
   kval = 0;
   if ( i == 0 && *jval == 1)
      {
                                     /* This is the first encounter of
                                        constraint 0 on this path so we 
                                        allocate memory and store parameters 
                                        to be used when n = 1 */

       if (special == 0) 
          {
           ttemp = (float *) malloc(4*mm);
           for (j=0;j<mm1;j++) *(ttemp+j) = - *(a+j+1+tmm)/pivot;
           kval = 1;
           printf("\n Generating temp n = %d m = %d i = %d \n",nn,mm,i); 
           printf(" ttemp : \n"); 
           for (j=0;j<mm1;j++) printf(" %f ",*(ttemp+j));
           printf("\n"); 
          }
       jjval = 0;
      }
   else if (*jval == 0)             /* Constraint 0 has already been
                                        encountered */
      {
        jjval = 0;  
                                    /* perform recursive update of component
                                       derivative array temp. This eliminates 
                                       row i and copies into a new vector */

        if(special == 0)
          {
           ttemp = (float *) malloc(4*mm1);
           for (j=0;j<i;j++) *(ttemp+j) = *(temp+j) 
                                          -(*(temp+i) * *(a+j+tmm)/pivot);
           for (j=i;j<mm1;j++) *(ttemp+j) = *(temp+j+1)
                                          -(*(temp+i) * *(a+j+1+tmm)/pivot);
           printf(" Reduced ttemp : \n"); 
           for (j=0;j<mm1;j++) printf(" %f ",*(ttemp+j));
           printf("\n"); 
          }
      }
   else 
      {                             /* Constraint 0 has not yet been 
                                        encountered */
      jjval = 1;                
      }


/* allocate memory */

  ap = (float *) malloc(4*nm1*mm1);
  bp = (float *) malloc(4*mm1);

/* reduce a and b into ap and bp eliminating variable t and constraint i */

  jj=-1;
    for (j=0;j<mm;j++)
      if (j != i)
      {
      jj=jj+1;
      aj=a+j;
      ajt=aj+tmm;
      *(bp+jj)= *(b+j) - *(bi) * *(ajt) / pivot;
      apjj=ap+jj;
      kk=-1;
      for (k=0;k<nn;k++)
        if (k != t)
        {
        kk=kk+1;
        kmm=k*mm;
        *(apjj+kk*mm1)= *(aj+kmm)- *(ajt) * *(ai+kmm)/ pivot;
        }
      }
  
/* add contribution to derivative from that calculated in smaller dimension */

  					/* Normal case method */
  if(special == 0)
    {
/*
     deriv = 0.;
     if (*(bi) != 0.) 
     {deriv=cdvda(ap,bp, &mm1, &nm1, &mm1, &nm1, &ttdim, ttemp, &jjval,code);}
*/
#ifndef _UNDERSCORE
     deriv=cdvda(ap,bp, &mm1, &nm1, &mm1, &nm1, &ttdim, ttemp, &jjval,code);
#else
     deriv=cdvda_(ap,bp, &mm1, &nm1, &mm1, &nm1, &ttdim, ttemp, &jjval,code);
#endif
     v=v+ *(bi)/amax*deriv/nn;
     if (kval == 1 || *jval == 0) free(ttemp);
     if(*code != 0)return (0.);

    }
  else					/* Use special case method */
    {
     if( *jval == 1)                    
       {                     
        if(i == 0)                      /* This is constraint 0 */
          {
           deriv = 0.; 
           vol = 0.; 
           dvdb = 0.;
           printf("\n special case 2: reduction by 0 and variable k\n"); 
           for (j=1;j<mm;j++) 
              {
               k = j - 1;
#ifndef _UNDERSCORE
               junk=cvolumebj(ap,bp,&mm1,&nm1,&mm1,&nm1,k,&dvdb);
#else
               junk=cvolumebj_(ap,bp,&mm1,&nm1,&mm1,&nm1,k,&dvdb);
#endif
               if(junk == -1.)
                 {
                  *code = -1;
                  return(0.);
                 }
               printf("\n j %d k %d dvdb %f vol %f",j,k,dvdb,junk); 
               deriv = deriv + dvdb * *(a + j + tmm) ;
               vol = vol + junk;
              }
           if(nm1 == 1)vol = junk;
           deriv = *(bi) * deriv/pivot;
           deriv = (deriv - vol) /pivot; 
           v=v+ *(bi)/amax*deriv/nn;
           printf("dcont %f",deriv); 
           printf("\n total volume of face = %f",vol); 

/* start debug section */

           FDa = *(bi)/amax*deriv/nn;
           da = 0.01* *(a + *tdim * mm);
           if(da == 0.0)da = 0.01;
           *(a + *tdim * mm) = *(a + *tdim * mm) + da;

  jj=-1;
    for (j=0;j<mm;j++)
      if (j != i)
      {
      jj=jj+1;
      aj=a+j;
      ajt=aj+tmm;
      *(bp+jj)= *(b+j) - *(bi) * *(ajt) / pivot;
      apjj=ap+jj;
      kk=-1;
      for (k=0;k<nn;k++)
        if (k != t)
        {
        kk=kk+1;
        kmm=k*mm;
        *(apjj+kk*mm1)= *(aj+kmm)- *(ajt) * *(ai+kmm)/ pivot;
        }
      }
#ifndef _UNDERSCORE
          volp = cvolume(ap,bp,&mm1,&nm1,&mm1,&nm1); 
#else
          volp = cvolume_(ap,bp,&mm1,&nm1,&mm1,&nm1); 
#endif
          amax = fabs(pivot+da);
          volp = (*(bi)/(amax*nn)) * volp;
           *(a + *tdim * mm) = *(a + *tdim * mm) - 2. * da;
  jj=-1;
    for (j=0;j<mm;j++)
      if (j != i)
      {
      jj=jj+1;
      aj=a+j;
      ajt=aj+tmm;
      *(bp+jj)= *(b+j) - *(bi) * *(ajt) / pivot;
      apjj=ap+jj;
      kk=-1;
      for (k=0;k<nn;k++)
        if (k != t)
        {
        kk=kk+1;
        kmm=k*mm;
        *(apjj+kk*mm1)= *(aj+kmm)- *(ajt) * *(ai+kmm)/ pivot;
        }
      }
#ifndef _UNDERSCORE
          volm = cvolume(ap,bp,&mm1,&nm1,&mm1,&nm1); 
#else
          volm = cvolume_(ap,bp,&mm1,&nm1,&mm1,&nm1); 
#endif
          amax = fabs(pivot-da);
           *(a + *tdim * mm) = *(a + *tdim * mm) + da;
          volm = (*(bi)/(amax*nn)) * volm;
          FD = (volp-volm)/(2. * da);
          printf("\n k  %d a = %f da = %f",*tdim,*(a + *tdim * mm),da); 
          printf("\n pivot = %f",pivot); 
          printf("\n volp = %f",volp); 
          printf("\n volm = %f",volm); 
          printf("\n FD estimate = %f",FD); 
          printf("\n analytical estimate = %f\n",FDa); 

/* end debug section */

          }
        else				/* Constraint 0 not yet encountered */
          {
           printf("\n special case 1: constraint 0 not yet reached\n"); 
           opt = 2;
#ifndef _UNDERSCORE
           junk=cvolumeb(ap,bp,&mm1,&nm1,&mm1,&nm1,&opt,&deriv);
#else
           junk=cvolumeb_(ap,bp,&mm1,&nm1,&mm1,&nm1,&opt,&deriv);
#endif
           if(junk == -1.)
             {
              *code = -1;
              return(0.);
             }
           deriv= -(deriv *  *(bi)/pivot);
           v=v+ *(bi)/amax*deriv/nn;
           printf("\n spec: 0 not yet reached i %d dcont %f",i,deriv); 
 
/* debug section: reduce system and repeat calculation with finite difference */

           FDa = *(bi)/amax*deriv/nn;
           da = 0.01* *(a + *tdim * mm);
           if(da == 0.0)da = 0.01;
           *(a + *tdim * mm) = *(a + *tdim * mm) + da;

/* reduce system */

  jj=-1;
    for (j=0;j<mm;j++)
      if (j != i)
      {
      jj=jj+1;
      aj=a+j;
      ajt=aj+tmm;
      *(bp+jj)= *(b+j) - *(bi) * *(ajt) / pivot;
      apjj=ap+jj;
      kk=-1;
      for (k=0;k<nn;k++)
        if (k != t)
        {
        kk=kk+1;
        kmm=k*mm;
        *(apjj+kk*mm1)= *(aj+kmm)- *(ajt) * *(ai+kmm)/ pivot;
        }
      }
#ifndef _UNDERSCORE
          volp = cvolume(ap,bp,&mm1,&nm1,&mm1,&nm1); 
#else
          volp = cvolume_(ap,bp,&mm1,&nm1,&mm1,&nm1); 
#endif
           *(a + *tdim * mm) = *(a + *tdim * mm) - 2. * da;

/* reduce system */

  jj=-1;
    for (j=0;j<mm;j++)
      if (j != i)
      {
      jj=jj+1;
      aj=a+j;
      ajt=aj+tmm;
      *(bp+jj)= *(b+j) - *(bi) * *(ajt) / pivot;
      apjj=ap+jj;
      kk=-1;
      for (k=0;k<nn;k++)
        if (k != t)
        {
        kk=kk+1;
        kmm=k*mm;
        *(apjj+kk*mm1)= *(aj+kmm)- *(ajt) * *(ai+kmm)/ pivot;
        }
      }
#ifndef _UNDERSCORE
          volm = cvolume(ap,bp,&mm1,&nm1,&mm1,&nm1); 
#else
          volm = cvolume_(ap,bp,&mm1,&nm1,&mm1,&nm1); 
#endif
           *(a + *tdim * mm) = *(a + *tdim * mm) + da;
          FD = (*(bi)/(amax*nn)) * (volp-volm)/(2. * da);
          printf("\n k  %d a = %f da = %f",*tdim,*(a + *tdim * mm),da); 
          printf("\n volp = %f",volp); 
          printf("\n volm = %f",volm); 
          printf("\n FD estimate = %f",FD); 
          printf("\n analytical estimate = %f\n",FDa); 

/* end debug section */

          }
       }
     else if( *jval == 0)               /* Constraint 0 already encountered */
       {
        printf("\n special case 3: constraint 0 already reduced\n"); 
        vol = 0.; 
        deriv = 0.; 
        for (j=0;j<i;j++) 
           {
#ifndef _UNDERSCORE
            junk=cvolumebj(ap,bp,&mm1,&nm1,&mm1,&nm1,j,&dvdb);
#else
            junk=cvolumebj_(ap,bp,&mm1,&nm1,&mm1,&nm1,j,&dvdb);
#endif
            if(junk == -1.)
              {
               *code = -1;
               return(0.);
              }
            dbda = (*(a+j+tmm) * *(temp+i)/pivot) - *(temp+j);
            deriv = deriv + dvdb*dbda;
            vol = vol + junk;
           }
        for (j=i+1;j<mm;j++) 
           {
            k = j - 1;
#ifndef _UNDERSCORE
            junk=cvolumebj(ap,bp,&mm1,&nm1,&mm1,&nm1,k,&dvdb);
#else
            junk=cvolumebj_(ap,bp,&mm1,&nm1,&mm1,&nm1,k,&dvdb);
#endif
            if(junk == -1.)
              {
               *code = -1;
               return(0.);
              }
            dbda = (*(a+j+tmm) * *(temp+i)/pivot) - *(temp+j);
            deriv = deriv + dvdb*dbda;
            vol = vol + junk;
           }
           if(nm1 == 1)vol = junk;
           deriv = (*(bi) * deriv - *(temp+i)*vol)/pivot;
           v=v+ *(bi)/amax*deriv/nn;
           printf("\n spec: 0 already encountered i = %d dcont %f",i,deriv); 
       }
    }

  free(ap);
  free(bp);
  }
  }

return(v);
}

