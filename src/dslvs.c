/* dslvs.f -- translated by f2c (version 12.02.01).
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
    doublereal rownd, rowns[210], el0, h__, hmin, hmxi, hu, tn, uround;
    integer iownd[14], iowns[6], ier, jstart, kflag, l, meth, miter, maxord, 
	    n, nq, nst, nfe, nje, nqu;
} ddebd1_;

#define ddebd1_1 ddebd1_

/* Table of constant values */

static integer c__0 = 0;

/* DECK DSLVS */
/* Subroutine */ int dslvs_(doublereal *wm, integer *iwm, doublereal *x, 
	doublereal *tem)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__;
    static doublereal r__, di;
    static integer ml, mu;
    static doublereal hl0, phl0;
    extern /* Subroutine */ int dgbsl_(doublereal *, integer *, integer *, 
	    integer *, integer *, integer *, doublereal *, integer *), dgesl_(
	    doublereal *, integer *, integer *, integer *, doublereal *, 
	    integer *);
    static integer meband;

/* ***BEGIN PROLOGUE  DSLVS */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to DDEBDF */
/* ***LIBRARY   SLATEC */
/* ***TYPE      DOUBLE PRECISION (SLVS-S, DSLVS-D) */
/* ***AUTHOR  Watts, H. A., (SNLA) */
/* ***DESCRIPTION */

/*   DSLVS solves the linear system in the iteration scheme for the */
/*   integrator package DDEBDF. */

/* ***SEE ALSO  DDEBDF */
/* ***ROUTINES CALLED  DGBSL, DGESL */
/* ***COMMON BLOCKS    DDEBD1 */
/* ***REVISION HISTORY  (YYMMDD) */
/*   820301  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900328  Added TYPE section.  (WRB) */
/*   910722  Updated AUTHOR section.  (ALS) */
/*   920422  Changed DIMENSION statement.  (WRB) */
/* ***END PROLOGUE  DSLVS */

/*     ------------------------------------------------------------------ */
/*      THIS ROUTINE MANAGES THE SOLUTION OF THE LINEAR SYSTEM ARISING */
/*      FROM A CHORD ITERATION.  IT IS CALLED BY DSTOD  IF MITER .NE. 0. */
/*      IF MITER IS 1 OR 2, IT CALLS DGESL TO ACCOMPLISH THIS. */
/*      IF MITER = 3 IT UPDATES THE COEFFICIENT H*EL0 IN THE DIAGONAL */
/*      MATRIX, AND THEN COMPUTES THE SOLUTION. */
/*      IF MITER IS 4 OR 5, IT CALLS DGBSL. */
/*      COMMUNICATION WITH DSLVS USES THE FOLLOWING VARIABLES.. */
/*      WM  = DOUBLE PRECISION WORK SPACE CONTAINING THE INVERSE DIAGONAL */
/*      MATRIX IF MITER */
/*            IS 3 AND THE LU DECOMPOSITION OF THE MATRIX OTHERWISE. */
/*            STORAGE OF MATRIX ELEMENTS STARTS AT WM(3). */
/*            WM ALSO CONTAINS THE FOLLOWING MATRIX-RELATED DATA.. */
/*            WM(1) = SQRT(UROUND) (NOT USED HERE), */
/*            WM(2) = HL0, THE PREVIOUS VALUE OF H*EL0, USED IF MITER = */
/*            3. */
/*      IWM = INTEGER WORK SPACE CONTAINING PIVOT INFORMATION, STARTING */
/*            AT IWM(21), IF MITER IS 1, 2, 4, OR 5.  IWM ALSO CONTAINS */
/*            THE BAND PARAMETERS ML = IWM(1) AND MU = IWM(2) IF MITER IS */
/*            4 OR 5. */
/*      X   = THE RIGHT-HAND SIDE VECTOR ON INPUT, AND THE SOLUTION */
/*            VECTOR ON OUTPUT, OF LENGTH N. */
/*      TEM = VECTOR OF WORK SPACE OF LENGTH N, NOT USED IN THIS VERSION. */
/*      IER = OUTPUT FLAG (IN COMMON).  IER = 0 IF NO TROUBLE OCCURRED. */
/*            IER = -1 IF A SINGULAR MATRIX AROSE WITH MITER = 3. */
/*      THIS ROUTINE ALSO USES THE COMMON VARIABLES EL0, H, MITER, AND N. */
/* ----------------------------------------------------------------------- */
/*     BEGIN BLOCK PERMITTING ...EXITS TO 80 */
/*        BEGIN BLOCK PERMITTING ...EXITS TO 60 */
/* ***FIRST EXECUTABLE STATEMENT  DSLVS */
    /* Parameter adjustments */
    --tem;
    --x;
    --iwm;
    --wm;

    /* Function Body */
    ddebd1_1.ier = 0;
    switch (ddebd1_1.miter) {
	case 1:  goto L10;
	case 2:  goto L10;
	case 3:  goto L20;
	case 4:  goto L70;
	case 5:  goto L70;
    }
L10:
    dgesl_(&wm[3], &ddebd1_1.n, &ddebd1_1.n, &iwm[21], &x[1], &c__0);
/*     ......EXIT */
    goto L80;

L20:
    phl0 = wm[2];
    hl0 = ddebd1_1.h__ * ddebd1_1.el0;
    wm[2] = hl0;
    if (hl0 == phl0) {
	goto L40;
    }
    r__ = hl0 / phl0;
    i__1 = ddebd1_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	di = 1. - r__ * (1. - 1. / wm[i__ + 2]);
/*        .........EXIT */
	if (abs(di) == 0.) {
	    goto L60;
	}
	wm[i__ + 2] = 1. / di;
/* L30: */
    }
L40:
    i__1 = ddebd1_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	x[i__] = wm[i__ + 2] * x[i__];
/* L50: */
    }
/*     ......EXIT */
    goto L80;
L60:
    ddebd1_1.ier = -1;
/*     ...EXIT */
    goto L80;

L70:
    ml = iwm[1];
    mu = iwm[2];
    meband = (ml << 1) + mu + 1;
    dgbsl_(&wm[3], &meband, &ddebd1_1.n, &ml, &mu, &iwm[21], &x[1], &c__0);
L80:
    return 0;
/*     ----------------------- END OF SUBROUTINE DSLVS */
/*     ----------------------- */
} /* dslvs_ */

