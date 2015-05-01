/* dintrv.f -- translated by f2c (version 12.02.01).
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

/* DECK DINTRV */
/* Subroutine */ int dintrv_(doublereal *xt, integer *lxt, doublereal *x, 
	integer *ilo, integer *ileft, integer *mflag)
{
    static integer ihi, istep, middle;

/* ***BEGIN PROLOGUE  DINTRV */
/* ***PURPOSE  Compute the largest integer ILEFT in 1 .LE. ILEFT .LE. LXT */
/*            such that XT(ILEFT) .LE. X where XT(*) is a subdivision of */
/*            the X interval. */
/* ***LIBRARY   SLATEC */
/* ***CATEGORY  E3, K6 */
/* ***TYPE      DOUBLE PRECISION (INTRV-S, DINTRV-D) */
/* ***KEYWORDS  B-SPLINE, DATA FITTING, INTERPOLATION, SPLINES */
/* ***AUTHOR  Amos, D. E., (SNLA) */
/* ***DESCRIPTION */

/*     Written by Carl de Boor and modified by D. E. Amos */

/*     Abstract    **** a double precision routine **** */
/*         DINTRV is the INTERV routine of the reference. */

/*         DINTRV computes the largest integer ILEFT in 1 .LE. ILEFT .LE. */
/*         LXT such that XT(ILEFT) .LE. X where XT(*) is a subdivision of */
/*         the X interval.  Precisely, */

/*                      X .LT. XT(1)                1         -1 */
/*         if  XT(I) .LE. X .LT. XT(I+1)  then  ILEFT=I  , MFLAG=0 */
/*           XT(LXT) .LE. X                         LXT        1, */

/*         That is, when multiplicities are present in the break point */
/*         to the left of X, the largest index is taken for ILEFT. */

/*     Description of Arguments */

/*         Input      XT,X are double precision */
/*          XT      - XT is a knot or break point vector of length LXT */
/*          LXT     - length of the XT vector */
/*          X       - argument */
/*          ILO     - an initialization parameter which must be set */
/*                    to 1 the first time the spline array XT is */
/*                    processed by DINTRV. */

/*         Output */
/*          ILO     - ILO contains information for efficient process- */
/*                    ing after the initial call and ILO must not be */
/*                    changed by the user.  Distinct splines require */
/*                    distinct ILO parameters. */
/*          ILEFT   - largest integer satisfying XT(ILEFT) .LE. X */
/*          MFLAG   - signals when X lies out of bounds */

/*     Error Conditions */
/*         None */

/* ***REFERENCES  Carl de Boor, Package for calculating with B-splines, */
/*                 SIAM Journal on Numerical Analysis 14, 3 (June 1977), */
/*                 pp. 441-472. */
/* ***ROUTINES CALLED  (NONE) */
/* ***REVISION HISTORY  (YYMMDD) */
/*   800901  DATE WRITTEN */
/*   890831  Modified array declarations.  (WRB) */
/*   890831  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   920501  Reformatted the REFERENCES section.  (WRB) */
/* ***END PROLOGUE  DINTRV */

/* ***FIRST EXECUTABLE STATEMENT  DINTRV */
    /* Parameter adjustments */
    --xt;

    /* Function Body */
    ihi = *ilo + 1;
    if (ihi < *lxt) {
	goto L10;
    }
    if (*x >= xt[*lxt]) {
	goto L110;
    }
    if (*lxt <= 1) {
	goto L90;
    }
    *ilo = *lxt - 1;
    ihi = *lxt;

L10:
    if (*x >= xt[ihi]) {
	goto L40;
    }
    if (*x >= xt[*ilo]) {
	goto L100;
    }

/* *** NOW X .LT. XT(IHI) . FIND LOWER BOUND */
    istep = 1;
L20:
    ihi = *ilo;
    *ilo = ihi - istep;
    if (*ilo <= 1) {
	goto L30;
    }
    if (*x >= xt[*ilo]) {
	goto L70;
    }
    istep <<= 1;
    goto L20;
L30:
    *ilo = 1;
    if (*x < xt[1]) {
	goto L90;
    }
    goto L70;
/* *** NOW X .GE. XT(ILO) . FIND UPPER BOUND */
L40:
    istep = 1;
L50:
    *ilo = ihi;
    ihi = *ilo + istep;
    if (ihi >= *lxt) {
	goto L60;
    }
    if (*x < xt[ihi]) {
	goto L70;
    }
    istep <<= 1;
    goto L50;
L60:
    if (*x >= xt[*lxt]) {
	goto L110;
    }
    ihi = *lxt;

/* *** NOW XT(ILO) .LE. X .LT. XT(IHI) . NARROW THE INTERVAL */
L70:
    middle = (*ilo + ihi) / 2;
    if (middle == *ilo) {
	goto L100;
    }
/*     NOTE. IT IS ASSUMED THAT MIDDLE = ILO IN CASE IHI = ILO+1 */
    if (*x < xt[middle]) {
	goto L80;
    }
    *ilo = middle;
    goto L70;
L80:
    ihi = middle;
    goto L70;
/* *** SET OUTPUT AND RETURN */
L90:
    *mflag = -1;
    *ileft = 1;
    return 0;
L100:
    *mflag = 0;
    *ileft = *ilo;
    return 0;
L110:
    *mflag = 1;
    *ileft = *lxt;
    return 0;
} /* dintrv_ */

