/* slvs.f -- translated by f2c (version 12.02.01).
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
    real rownd, rowns[210], el0, h__, hmin, hmxi, hu, tn, uround;
    integer iownd[14], iowns[6], ier, jstart, kflag, l, meth, miter, maxord, 
	    n, nq, nst, nfe, nje, nqu;
} debdf1_;

#define debdf1_1 debdf1_

/* Table of constant values */

static integer c__0 = 0;

/* DECK SLVS */
/* Subroutine */ int slvs_(real *wm, integer *iwm, real *x, real *tem)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__;
    static real r__, di;
    static integer ml, mu;
    static real hl0, phl0;
    extern /* Subroutine */ int sgbsl_(real *, integer *, integer *, integer *
	    , integer *, integer *, real *, integer *), sgesl_(real *, 
	    integer *, integer *, integer *, real *, integer *);
    static integer meband;

/* ***BEGIN PROLOGUE  SLVS */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to DEBDF */
/* ***LIBRARY   SLATEC */
/* ***TYPE      SINGLE PRECISION (SLVS-S, DSLVS-D) */
/* ***AUTHOR  Watts, H. A., (SNLA) */
/* ***DESCRIPTION */

/*   SLVS solves the linear system in the iteration scheme for the */
/*   integrator package DEBDF. */

/* ***SEE ALSO  DEBDF */
/* ***ROUTINES CALLED  SGBSL, SGESL */
/* ***COMMON BLOCKS    DEBDF1 */
/* ***REVISION HISTORY  (YYMMDD) */
/*   800901  DATE WRITTEN */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900328  Added TYPE section.  (WRB) */
/*   910722  Updated AUTHOR section.  (ALS) */
/*   920422  Changed DIMENSION statement.  (WRB) */
/* ***END PROLOGUE  SLVS */

/* LLL. OPTIMIZE */
/* ----------------------------------------------------------------------- */
/* THIS ROUTINE MANAGES THE SOLUTION OF THE LINEAR SYSTEM ARISING FROM */
/* A CHORD ITERATION.  IT IS CALLED BY STOD  IF MITER .NE. 0. */
/* IF MITER IS 1 OR 2, IT CALLS SGESL TO ACCOMPLISH THIS. */
/* IF MITER = 3 IT UPDATES THE COEFFICIENT H*EL0 IN THE DIAGONAL */
/* MATRIX, AND THEN COMPUTES THE SOLUTION. */
/* IF MITER IS 4 OR 5, IT CALLS SGBSL. */
/* COMMUNICATION WITH SLVS USES THE FOLLOWING VARIABLES.. */
/* WM  = REAL WORK SPACE CONTAINING THE INVERSE DIAGONAL MATRIX IF MITER */
/*       IS 3 AND THE LU DECOMPOSITION OF THE MATRIX OTHERWISE. */
/*       STORAGE OF MATRIX ELEMENTS STARTS AT WM(3). */
/*       WM ALSO CONTAINS THE FOLLOWING MATRIX-RELATED DATA.. */
/*       WM(1) = SQRT(UROUND) (NOT USED HERE), */
/*       WM(2) = HL0, THE PREVIOUS VALUE OF H*EL0, USED IF MITER = 3. */
/* IWM = INTEGER WORK SPACE CONTAINING PIVOT INFORMATION, STARTING AT */
/*       IWM(21), IF MITER IS 1, 2, 4, OR 5.  IWM ALSO CONTAINS THE */
/*       BAND PARAMETERS ML = IWM(1) AND MU = IWM(2) IF MITER IS 4 OR 5. */
/* X   = THE RIGHT-HAND SIDE VECTOR ON INPUT, AND THE SOLUTION VECTOR */
/*       ON OUTPUT, OF LENGTH N. */
/* TEM = VECTOR OF WORK SPACE OF LENGTH N, NOT USED IN THIS VERSION. */
/* IER = OUTPUT FLAG (IN COMMON).  IER = 0 IF NO TROUBLE OCCURRED. */
/*       IER = -1 IF A SINGULAR MATRIX AROSE WITH MITER = 3. */
/* THIS ROUTINE ALSO USES THE COMMON VARIABLES EL0, H, MITER, AND N. */
/* ----------------------------------------------------------------------- */
/* ***FIRST EXECUTABLE STATEMENT  SLVS */
    /* Parameter adjustments */
    --tem;
    --x;
    --iwm;
    --wm;

    /* Function Body */
    debdf1_1.ier = 0;
    switch (debdf1_1.miter) {
	case 1:  goto L100;
	case 2:  goto L100;
	case 3:  goto L300;
	case 4:  goto L400;
	case 5:  goto L400;
    }
L100:
    sgesl_(&wm[3], &debdf1_1.n, &debdf1_1.n, &iwm[21], &x[1], &c__0);
    return 0;

L300:
    phl0 = wm[2];
    hl0 = debdf1_1.h__ * debdf1_1.el0;
    wm[2] = hl0;
    if (hl0 == phl0) {
	goto L330;
    }
    r__ = hl0 / phl0;
    i__1 = debdf1_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	di = 1.f - r__ * (1.f - 1.f / wm[i__ + 2]);
	if (dabs(di) == 0.f) {
	    goto L390;
	}
/* L320: */
	wm[i__ + 2] = 1.f / di;
    }
L330:
    i__1 = debdf1_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L340: */
	x[i__] = wm[i__ + 2] * x[i__];
    }
    return 0;
L390:
    debdf1_1.ier = -1;
    return 0;

L400:
    ml = iwm[1];
    mu = iwm[2];
    meband = (ml << 1) + mu + 1;
    sgbsl_(&wm[3], &meband, &debdf1_1.n, &ml, &mu, &iwm[21], &x[1], &c__0);
    return 0;
/* ----------------------- END OF SUBROUTINE SLVS ----------------------- */
} /* slvs_ */

