/* chfev.f -- translated by f2c (version 12.02.01).
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

/* Table of constant values */

static integer c__1 = 1;

/* DECK CHFEV */
/* Subroutine */ int chfev_(real *x1, real *x2, real *f1, real *f2, real *d1, 
	real *d2, integer *ne, real *xe, real *fe, integer *next, integer *
	ierr)
{
    /* Initialized data */

    static real zero = 0.f;

    /* System generated locals */
    integer i__1;

    /* Local variables */
    static real h__;
    static integer i__;
    static real x, c2, c3, xma, xmi, del1, del2, delta;
    extern /* Subroutine */ int xermsg_(char *, char *, char *, integer *, 
	    integer *, ftnlen, ftnlen, ftnlen);

/* ***BEGIN PROLOGUE  CHFEV */
/* ***PURPOSE  Evaluate a cubic polynomial given in Hermite form at an */
/*            array of points.  While designed for use by PCHFE, it may */
/*            be useful directly as an evaluator for a piecewise cubic */
/*            Hermite function in applications, such as graphing, where */
/*            the interval is known in advance. */
/* ***LIBRARY   SLATEC (PCHIP) */
/* ***CATEGORY  E3 */
/* ***TYPE      SINGLE PRECISION (CHFEV-S, DCHFEV-D) */
/* ***KEYWORDS  CUBIC HERMITE EVALUATION, CUBIC POLYNOMIAL EVALUATION, */
/*             PCHIP */
/* ***AUTHOR  Fritsch, F. N., (LLNL) */
/*             Lawrence Livermore National Laboratory */
/*             P.O. Box 808  (L-316) */
/*             Livermore, CA  94550 */
/*             FTS 532-4275, (510) 422-4275 */
/* ***DESCRIPTION */

/*          CHFEV:  Cubic Hermite Function EValuator */

/*     Evaluates the cubic polynomial determined by function values */
/*     F1,F2 and derivatives D1,D2 on interval (X1,X2) at the points */
/*     XE(J), J=1(1)NE. */

/* ---------------------------------------------------------------------- */

/*  Calling sequence: */

/*        INTEGER  NE, NEXT(2), IERR */
/*        REAL  X1, X2, F1, F2, D1, D2, XE(NE), FE(NE) */

/*        CALL  CHFEV (X1,X2, F1,F2, D1,D2, NE, XE, FE, NEXT, IERR) */

/*   Parameters: */

/*     X1,X2 -- (input) endpoints of interval of definition of cubic. */
/*           (Error return if  X1.EQ.X2 .) */

/*     F1,F2 -- (input) values of function at X1 and X2, respectively. */

/*     D1,D2 -- (input) values of derivative at X1 and X2, respectively. */

/*     NE -- (input) number of evaluation points.  (Error return if */
/*           NE.LT.1 .) */

/*     XE -- (input) real array of points at which the function is to be */
/*           evaluated.  If any of the XE are outside the interval */
/*           [X1,X2], a warning error is returned in NEXT. */

/*     FE -- (output) real array of values of the cubic function defined */
/*           by  X1,X2, F1,F2, D1,D2  at the points  XE. */

/*     NEXT -- (output) integer array indicating number of extrapolation */
/*           points: */
/*            NEXT(1) = number of evaluation points to left of interval. */
/*            NEXT(2) = number of evaluation points to right of interval. */

/*     IERR -- (output) error flag. */
/*           Normal return: */
/*              IERR = 0  (no errors). */
/*           "Recoverable" errors: */
/*              IERR = -1  if NE.LT.1 . */
/*              IERR = -2  if X1.EQ.X2 . */
/*                (The FE-array has not been changed in either case.) */

/* ***REFERENCES  (NONE) */
/* ***ROUTINES CALLED  XERMSG */
/* ***REVISION HISTORY  (YYMMDD) */
/*   811019  DATE WRITTEN */
/*   820803  Minor cosmetic changes for release 1. */
/*   890411  Added SAVE statements (Vers. 3.2). */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890703  Corrected category record.  (WRB) */
/*   890703  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ) */
/* ***END PROLOGUE  CHFEV */
/*  Programming notes: */

/*     To produce a double precision version, simply: */
/*        a. Change CHFEV to DCHFEV wherever it occurs, */
/*        b. Change the real declaration to double precision, and */
/*        c. Change the constant ZERO to double precision. */

/*  DECLARE ARGUMENTS. */


/*  DECLARE LOCAL VARIABLES. */

    /* Parameter adjustments */
    --next;
    --fe;
    --xe;

    /* Function Body */

/*  VALIDITY-CHECK ARGUMENTS. */

/* ***FIRST EXECUTABLE STATEMENT  CHFEV */
    if (*ne < 1) {
	goto L5001;
    }
    h__ = *x2 - *x1;
    if (h__ == zero) {
	goto L5002;
    }

/*  INITIALIZE. */

    *ierr = 0;
    next[1] = 0;
    next[2] = 0;
    xmi = dmin(zero,h__);
    xma = dmax(zero,h__);

/*  COMPUTE CUBIC COEFFICIENTS (EXPANDED ABOUT X1). */

    delta = (*f2 - *f1) / h__;
    del1 = (*d1 - delta) / h__;
    del2 = (*d2 - delta) / h__;
/*                                           (DELTA IS NO LONGER NEEDED.) */
    c2 = -(del1 + del1 + del2);
    c3 = (del1 + del2) / h__;
/*                               (H, DEL1 AND DEL2 ARE NO LONGER NEEDED.) */

/*  EVALUATION LOOP. */

    i__1 = *ne;
    for (i__ = 1; i__ <= i__1; ++i__) {
	x = xe[i__] - *x1;
	fe[i__] = *f1 + x * (*d1 + x * (c2 + x * c3));
/*          COUNT EXTRAPOLATION POINTS. */
	if (x < xmi) {
	    ++next[1];
	}
	if (x > xma) {
	    ++next[2];
	}
/*        (NOTE REDUNDANCY--IF EITHER CONDITION IS TRUE, OTHER IS FALSE.) */
/* L500: */
    }

/*  NORMAL RETURN. */

    return 0;

/*  ERROR RETURNS. */

L5001:
/*     NE.LT.1 RETURN. */
    *ierr = -1;
    xermsg_("SLATEC", "CHFEV", "NUMBER OF EVALUATION POINTS LESS THAN ONE", 
	    ierr, &c__1, (ftnlen)6, (ftnlen)5, (ftnlen)41);
    return 0;

L5002:
/*     X1.EQ.X2 RETURN. */
    *ierr = -2;
    xermsg_("SLATEC", "CHFEV", "INTERVAL ENDPOINTS EQUAL", ierr, &c__1, (
	    ftnlen)6, (ftnlen)5, (ftnlen)24);
    return 0;
/* ------------- LAST LINE OF CHFEV FOLLOWS ------------------------------ */
} /* chfev_ */

