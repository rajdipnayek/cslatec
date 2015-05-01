/* scoef.f -- translated by f2c (version 12.02.01).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/

#include <stdlib.h> /* For exit() */
#include <f2c.h>

/* Common Block Declarations */

struct {
    real uro, sru, eps, sqovfl, twou, fouru;
    integer lpar;
} ml5mco_;

#define ml5mco_1 ml5mco_

/* Table of constant values */

static integer c__1 = 1;
static integer c__0 = 0;

/* DECK SCOEF */
/* Subroutine */ int scoef_(real *yh, real *yp, integer *ncomp, integer *
	nrowb, integer *nfc, integer *nic, real *b, real *beta, real *coef, 
	integer *inhomo, real *re, real *ae, real *by, real *cvec, real *work,
	 integer *iwork, integer *iflag, integer *nfcc)
{
    /* System generated locals */
    integer yh_dim1, yh_offset, b_dim1, b_offset, by_dim1, by_offset, i__1, 
	    i__2;
    real r__1, r__2, r__3;

    /* Local variables */
    static integer i__, j, k, l;
    static real bn;
    static integer nf, ki;
    static real un, bbn, gam, brn, bys, ypn, bykl, cons;
    extern doublereal sdot_(integer *, real *, integer *, real *, integer *);
    static integer mlso;
    extern /* Subroutine */ int suds_(real *, real *, real *, integer *, 
	    integer *, integer *, integer *, integer *, real *, integer *);
    static integer kflag;
    extern /* Subroutine */ int xgetf_(integer *), xsetf_(integer *);
    static integer nfccm1, ncomp2;

/* ***BEGIN PROLOGUE  SCOEF */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to BVSUP */
/* ***LIBRARY   SLATEC */
/* ***TYPE      SINGLE PRECISION (SCOEF-S, DCOEF-D) */
/* ***AUTHOR  Watts, H. A., (SNLA) */
/* ***DESCRIPTION */

/* ********************************************************************** */
/* INPUT TO SCOEF */
/* ********************************************************************** */

/*     YH = Matrix of homogeneous solutions. */
/*     YP = Vector containing particular solution. */
/*     NCOMP = Number of components per solution vector. */
/*     NROWB = First dimension of B in calling program. */
/*     NFC = Number of base solution vectors. */
/*     NFCC = 2*NFC for the special treatment of complex valued */
/*            equations. Otherwise, NFCC=NFC. */
/*     NIC = Number of specified initial conditions. */
/*     B = Boundary condition matrix at X = Xfinal. */
/*     BETA = Vector of nonhomogeneous boundary conditions at X = Xfinal. */
/*              1 - Nonzero particular solution */
/*     INHOMO = 2 - Zero particular solution */
/*              3 - Eigenvalue problem */
/*     RE = Relative error tolerance */
/*     AE = Absolute error tolerance */
/*     BY = Storage space for the matrix  B*YH */
/*     CVEC = Storage space for the vector  BETA-B*YP */
/*     WORK = Real array of internal storage. Dimension must be .GE. */
/*            NFCC*(NFCC+4) */
/*     IWORK = Integer array of internal storage. Dimension must be .GE. */
/*             3+NFCC */

/* ********************************************************************** */
/* OUTPUT FROM SCOEF */
/* ********************************************************************** */

/*     COEF = Array containing superposition constants. */
/*     IFLAG = Indicator of success from SUDS in solving the */
/*             boundary equations */
/*           = 0 Boundary equations are solved */
/*           = 1 Boundary equations appear to have many solutions */
/*           = 2 Boundary equations appear to be inconsistent */
/*           = 3 For this value of an eigenparameter, the boundary */
/*               equations have only the zero solution. */

/* ********************************************************************** */

/*     Subroutine SCOEF solves for the superposition constants from the */
/*     linear equations defined by the boundary conditions at X = Xfinal. */

/*                          B*YP + B*YH*COEF = BETA */

/* ********************************************************************** */

/* ***SEE ALSO  BVSUP */
/* ***ROUTINES CALLED  SDOT, SUDS, XGETF, XSETF */
/* ***COMMON BLOCKS    ML5MCO */
/* ***REVISION HISTORY  (YYMMDD) */
/*   750601  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890831  Modified array declarations.  (WRB) */
/*   890921  Realigned order of variables in certain COMMON blocks. */
/*           (WRB) */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900328  Added TYPE section.  (WRB) */
/*   910722  Updated AUTHOR section.  (ALS) */
/* ***END PROLOGUE  SCOEF */



/*     SET UP MATRIX  B*YH  AND VECTOR  BETA - B*YP */

/* ***FIRST EXECUTABLE STATEMENT  SCOEF */
    /* Parameter adjustments */
    --yp;
    yh_dim1 = *ncomp;
    yh_offset = 1 + yh_dim1;
    yh -= yh_offset;
    b_dim1 = *nrowb;
    b_offset = 1 + b_dim1;
    b -= b_offset;
    --beta;
    --coef;
    --cvec;
    --work;
    --iwork;
    by_dim1 = *nfcc;
    by_offset = 1 + by_dim1;
    by -= by_offset;

    /* Function Body */
    ncomp2 = *ncomp / 2;
    i__1 = *nfcc;
    for (k = 1; k <= i__1; ++k) {
	i__2 = *nfc;
	for (j = 1; j <= i__2; ++j) {
	    l = j;
	    if (*nfc != *nfcc) {
		l = (j << 1) - 1;
	    }
/* L1: */
	    by[k + l * by_dim1] = sdot_(ncomp, &b[k + b_dim1], nrowb, &yh[j * 
		    yh_dim1 + 1], &c__1);
	}
	if (*nfc == *nfcc) {
	    goto L3;
	}
	i__2 = *nfc;
	for (j = 1; j <= i__2; ++j) {
	    l = j << 1;
	    bykl = sdot_(&ncomp2, &b[k + b_dim1], nrowb, &yh[ncomp2 + 1 + j * 
		    yh_dim1], &c__1);
	    by[k + l * by_dim1] = sdot_(&ncomp2, &b[k + (ncomp2 + 1) * b_dim1]
		    , nrowb, &yh[j * yh_dim1 + 1], &c__1) - bykl;
/* L2: */
	}
L3:
	switch (*inhomo) {
	    case 1:  goto L4;
	    case 2:  goto L5;
	    case 3:  goto L6;
	}
/*     CASE 1 */
L4:
	cvec[k] = beta[k] - sdot_(ncomp, &b[k + b_dim1], nrowb, &yp[1], &c__1)
		;
	goto L7;
/*     CASE 2 */
L5:
	cvec[k] = beta[k];
	goto L7;
/*     CASE 3 */
L6:
	cvec[k] = 0.f;
L7:
	;
    }
    cons = dabs(cvec[1]);
    bys = (r__1 = by[by_dim1 + 1], dabs(r__1));

/* ********************************************************************** */
/*     SOLVE LINEAR SYSTEM */

    *iflag = 0;
    mlso = 0;
    if (*inhomo == 3) {
	mlso = 1;
    }
    kflag = r_lg10(&ml5mco_1.eps) * .5f;
    xgetf_(&nf);
    xsetf_(&c__0);
L10:
    suds_(&by[by_offset], &coef[1], &cvec[1], nfcc, nfcc, nfcc, &kflag, &mlso,
	     &work[1], &iwork[1]);
    if (kflag != 3) {
	goto L13;
    }
    kflag = 1;
    *iflag = 1;
    goto L10;
L13:
    if (kflag == 4) {
	*iflag = 2;
    }
    xsetf_(&nf);
    if (*nfcc == 1) {
	goto L25;
    }
    if (*inhomo != 3) {
	return 0;
    }
    if (iwork[1] < *nfcc) {
	goto L17;
    }
    *iflag = 3;
    i__1 = *nfcc;
    for (k = 1; k <= i__1; ++k) {
/* L14: */
	coef[k] = 0.f;
    }
    coef[*nfcc] = 1.f;
    nfccm1 = *nfcc - 1;
    i__1 = nfccm1;
    for (k = 1; k <= i__1; ++k) {
	j = *nfcc - k;
	l = *nfcc - j + 1;
	gam = sdot_(&l, &by[j + j * by_dim1], nfcc, &coef[j], &c__1) / (work[
		j] * by[j + j * by_dim1]);
	i__2 = *nfcc;
	for (i__ = j; i__ <= i__2; ++i__) {
/* L15: */
	    coef[i__] += gam * by[j + i__ * by_dim1];
	}
    }
    return 0;
L17:
    i__2 = *nfcc;
    for (k = 1; k <= i__2; ++k) {
	ki = (*nfcc << 2) + k;
/* L20: */
	coef[k] = work[ki];
    }
    return 0;

/* ********************************************************************** */
/*     TESTING FOR EXISTENCE AND UNIQUENESS OF BOUNDARY-VALUE PROBLEM */
/*     SOLUTION IN A SCALAR CASE */

L25:
    bn = 0.f;
    un = 0.f;
    ypn = 0.f;
    i__2 = *ncomp;
    for (k = 1; k <= i__2; ++k) {
/* Computing MAX */
	r__2 = un, r__3 = (r__1 = yh[k + yh_dim1], dabs(r__1));
	un = dmax(r__2,r__3);
/* Computing MAX */
	r__2 = ypn, r__3 = (r__1 = yp[k], dabs(r__1));
	ypn = dmax(r__2,r__3);
/* L30: */
/* Computing MAX */
	r__2 = bn, r__3 = (r__1 = b[k * b_dim1 + 1], dabs(r__1));
	bn = dmax(r__2,r__3);
    }
/* Computing MAX */
    r__1 = bn, r__2 = dabs(beta[1]);
    bbn = dmax(r__1,r__2);
    if (bys > (*re * un + *ae) * 10.f * bn) {
	goto L35;
    }
    brn = bbn / bn * bys;
    if (cons >= brn * .1f && cons <= brn * 10.f) {
	*iflag = 1;
    }
    if (cons > brn * 10.f) {
	*iflag = 2;
    }
    if (cons <= *re * dabs(beta[1]) + *ae + (*re * ypn + *ae) * bn) {
	*iflag = 1;
    }
    if (*inhomo == 3) {
	coef[1] = 1.f;
    }
    return 0;
L35:
    if (*inhomo != 3) {
	return 0;
    }
    *iflag = 3;
    coef[1] = 1.f;
    return 0;
} /* scoef_ */

