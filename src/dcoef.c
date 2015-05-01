/* dcoef.f -- translated by f2c (version 12.02.01).
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
    doublereal uro, sru, eps, sqovfl, twou, fouru;
    integer lpar;
} dml5mc_;

#define dml5mc_1 dml5mc_

/* Table of constant values */

static integer c__1 = 1;
static integer c__0 = 0;

/* DECK DCOEF */
/* Subroutine */ int dcoef_(doublereal *yh, doublereal *yp, integer *ncomp, 
	integer *nrowb, integer *nfc, integer *nic, doublereal *b, doublereal 
	*beta, doublereal *coef, integer *inhomo, doublereal *re, doublereal *
	ae, doublereal *by, doublereal *cvec, doublereal *work, integer *
	iwork, integer *iflag, integer *nfcc)
{
    /* System generated locals */
    integer b_dim1, b_offset, by_dim1, by_offset, yh_dim1, yh_offset, i__1, 
	    i__2;
    doublereal d__1, d__2, d__3;

    /* Local variables */
    static integer i__, j, k, l;
    static doublereal bn;
    static integer ki, nf;
    static doublereal un, bbn, gam, brn, bys, ypn;
    extern doublereal ddot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    static doublereal bykl, cons;
    static integer mlso, kflag;
    extern /* Subroutine */ int xgetf_(integer *), dsuds_(doublereal *, 
	    doublereal *, doublereal *, integer *, integer *, integer *, 
	    integer *, integer *, doublereal *, integer *), xsetf_(integer *);
    static integer nfccm1, ncomp2;

/* ***BEGIN PROLOGUE  DCOEF */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to DBVSUP */
/* ***LIBRARY   SLATEC */
/* ***TYPE      DOUBLE PRECISION (SCOEF-S, DCOEF-D) */
/* ***AUTHOR  Watts, H. A., (SNLA) */
/* ***DESCRIPTION */

/* ********************************************************************** */
/* INPUT to DCOEF */
/* ********************************************************************** */

/*     YH = matrix of homogeneous solutions. */
/*     YP = vector containing particular solution. */
/*     NCOMP = number of components per solution vector. */
/*     NROWB = first dimension of B in calling program. */
/*     NFC = number of base solution vectors. */
/*     NFCC = 2*NFC for the special treatment of COMPLEX*16 valued */
/*            equations. Otherwise, NFCC=NFC. */
/*     NIC = number of specified initial conditions. */
/*     B = boundary condition matrix at X = XFINAL. */
/*     BETA = vector of nonhomogeneous boundary conditions at X = XFINAL. */
/*              1 - nonzero particular solution */
/*     INHOMO = 2 - zero particular solution */
/*              3 - eigenvalue problem */
/*     RE = relative error tolerance. */
/*     AE = absolute error tolerance. */
/*     BY = storage space for the matrix  B*YH */
/*     CVEC = storage space for the vector  BETA-B*YP */
/*     WORK = double precision array of internal storage. Dimension must */
/*     be GE */
/*            NFCC*(NFCC+4) */
/*     IWORK = integer array of internal storage. Dimension must be GE */
/*             3+NFCC */

/* ********************************************************************** */
/* OUTPUT from DCOEF */
/* ********************************************************************** */

/*     COEF = array containing superposition constants. */
/*     IFLAG = indicator of success from DSUDS in solving the */
/*             boundary equations. */
/*           = 0 boundary equations are solved. */
/*           = 1 boundary equations appear to have many solutions. */
/*           = 2 boundary equations appear to be inconsistent. */
/*           = 3 for this value of an eigenparameter, the boundary */
/*               equations have only the zero solution. */

/* ********************************************************************** */

/*     Subroutine DCOEF solves for the superposition constants from the */
/*     linear equations defined by the boundary conditions at X = XFINAL. */

/*                          B*YP + B*YH*COEF = BETA */

/* ********************************************************************** */

/* ***SEE ALSO  DBVSUP */
/* ***ROUTINES CALLED  DDOT, DSUDS, XGETF, XSETF */
/* ***COMMON BLOCKS    DML5MC */
/* ***REVISION HISTORY  (YYMMDD) */
/*   750601  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890831  Modified array declarations.  (WRB) */
/*   890921  Realigned order of variables in certain COMMON blocks. */
/*           (WRB) */
/*   890921  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900328  Added TYPE section.  (WRB) */
/*   910722  Updated AUTHOR section.  (ALS) */
/* ***END PROLOGUE  DCOEF */


/* ***FIRST EXECUTABLE STATEMENT  DCOEF */

/*     SET UP MATRIX  B*YH  AND VECTOR  BETA - B*YP */

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
	    by[k + l * by_dim1] = ddot_(ncomp, &b[k + b_dim1], nrowb, &yh[j * 
		    yh_dim1 + 1], &c__1);
/* L10: */
	}
	if (*nfc == *nfcc) {
	    goto L30;
	}
	i__2 = *nfc;
	for (j = 1; j <= i__2; ++j) {
	    l = j << 1;
	    bykl = ddot_(&ncomp2, &b[k + b_dim1], nrowb, &yh[ncomp2 + 1 + j * 
		    yh_dim1], &c__1);
	    by[k + l * by_dim1] = ddot_(&ncomp2, &b[k + (ncomp2 + 1) * b_dim1]
		    , nrowb, &yh[j * yh_dim1 + 1], &c__1) - bykl;
/* L20: */
	}
L30:
	switch (*inhomo) {
	    case 1:  goto L40;
	    case 2:  goto L50;
	    case 3:  goto L60;
	}
/*        CASE 1 */
L40:
	cvec[k] = beta[k] - ddot_(ncomp, &b[k + b_dim1], nrowb, &yp[1], &c__1)
		;
	goto L70;
/*        CASE 2 */
L50:
	cvec[k] = beta[k];
	goto L70;
/*        CASE 3 */
L60:
	cvec[k] = 0.;
L70:
/* L80: */
	;
    }
    cons = abs(cvec[1]);
    bys = (d__1 = by[by_dim1 + 1], abs(d__1));

/*     ****************************************************************** */
/*         SOLVE LINEAR SYSTEM */

    *iflag = 0;
    mlso = 0;
    if (*inhomo == 3) {
	mlso = 1;
    }
    kflag = (integer) (d_lg10(&dml5mc_1.eps) * .5);
    xgetf_(&nf);
    xsetf_(&c__0);
L90:
    dsuds_(&by[by_offset], &coef[1], &cvec[1], nfcc, nfcc, nfcc, &kflag, &
	    mlso, &work[1], &iwork[1]);
    if (kflag != 3) {
	goto L100;
    }
    kflag = 1;
    *iflag = 1;
    goto L90;
L100:
    if (kflag == 4) {
	*iflag = 2;
    }
    xsetf_(&nf);
    if (*nfcc == 1) {
	goto L180;
    }
    if (*inhomo != 3) {
	goto L170;
    }
    if (iwork[1] < *nfcc) {
	goto L140;
    }
    *iflag = 3;
    i__1 = *nfcc;
    for (k = 1; k <= i__1; ++k) {
	coef[k] = 0.;
/* L110: */
    }
    coef[*nfcc] = 1.;
    nfccm1 = *nfcc - 1;
    i__1 = nfccm1;
    for (k = 1; k <= i__1; ++k) {
	j = *nfcc - k;
	l = *nfcc - j + 1;
	gam = ddot_(&l, &by[j + j * by_dim1], nfcc, &coef[j], &c__1) / (work[
		j] * by[j + j * by_dim1]);
	i__2 = *nfcc;
	for (i__ = j; i__ <= i__2; ++i__) {
	    coef[i__] += gam * by[j + i__ * by_dim1];
/* L120: */
	}
/* L130: */
    }
    goto L160;
L140:
    i__1 = *nfcc;
    for (k = 1; k <= i__1; ++k) {
	ki = (*nfcc << 2) + k;
	coef[k] = work[ki];
/* L150: */
    }
L160:
L170:
    goto L220;
L180:

/*        *************************************************************** */
/*            TESTING FOR EXISTENCE AND UNIQUENESS OF BOUNDARY-VALUE */
/*            PROBLEM SOLUTION IN A SCALAR CASE */

    bn = 0.;
    un = 0.;
    ypn = 0.;
    i__1 = *ncomp;
    for (k = 1; k <= i__1; ++k) {
/* Computing MAX */
	d__2 = un, d__3 = (d__1 = yh[k + yh_dim1], abs(d__1));
	un = max(d__2,d__3);
/* Computing MAX */
	d__2 = ypn, d__3 = (d__1 = yp[k], abs(d__1));
	ypn = max(d__2,d__3);
/* Computing MAX */
	d__2 = bn, d__3 = (d__1 = b[k * b_dim1 + 1], abs(d__1));
	bn = max(d__2,d__3);
/* L190: */
    }
/* Computing MAX */
    d__1 = bn, d__2 = abs(beta[1]);
    bbn = max(d__1,d__2);
    if (bys > (*re * un + *ae) * 10. * bn) {
	goto L200;
    }
    brn = bbn / bn * bys;
    if (cons >= brn * .1 && cons <= brn * 10.) {
	*iflag = 1;
    }
    if (cons > brn * 10.) {
	*iflag = 2;
    }
    if (cons <= *re * abs(beta[1]) + *ae + (*re * ypn + *ae) * bn) {
	*iflag = 1;
    }
    if (*inhomo == 3) {
	coef[1] = 1.;
    }
    goto L210;
L200:
    if (*inhomo != 3) {
	goto L210;
    }
    *iflag = 3;
    coef[1] = 1.;
L210:
L220:
    return 0;
} /* dcoef_ */

