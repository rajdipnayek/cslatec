/* ddanrm.f -- translated by f2c (version 12.02.01).
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

/* DECK DDANRM */
doublereal ddanrm_(integer *neq, doublereal *v, doublereal *wt, doublereal *
	rpar, integer *ipar)
{
    /* System generated locals */
    integer i__1;
    doublereal ret_val, d__1, d__2;

    /* Local variables */
    static integer i__;
    static doublereal sum, vmax;

/* ***BEGIN PROLOGUE  DDANRM */
/* ***SUBSIDIARY */
/* ***PURPOSE  Compute vector norm for DDASSL. */
/* ***LIBRARY   SLATEC (DASSL) */
/* ***TYPE      DOUBLE PRECISION (SDANRM-S, DDANRM-D) */
/* ***AUTHOR  Petzold, Linda R., (LLNL) */
/* ***DESCRIPTION */
/* ----------------------------------------------------------------------- */
/*     THIS FUNCTION ROUTINE COMPUTES THE WEIGHTED */
/*     ROOT-MEAN-SQUARE NORM OF THE VECTOR OF LENGTH */
/*     NEQ CONTAINED IN THE ARRAY V,WITH WEIGHTS */
/*     CONTAINED IN THE ARRAY WT OF LENGTH NEQ. */
/*        DDANRM=SQRT((1/NEQ)*SUM(V(I)/WT(I))**2) */
/* ----------------------------------------------------------------------- */
/* ***ROUTINES CALLED  (NONE) */
/* ***REVISION HISTORY  (YYMMDD) */
/*   830315  DATE WRITTEN */
/*   901009  Finished conversion to SLATEC 4.0 format (F.N.Fritsch) */
/*   901019  Merged changes made by C. Ulrich with SLATEC 4.0 format. */
/*   901026  Added explicit declarations for all variables and minor */
/*           cosmetic changes to prologue.  (FNF) */
/* ***END PROLOGUE  DDANRM */



/* ***FIRST EXECUTABLE STATEMENT  DDANRM */
    /* Parameter adjustments */
    --wt;
    --v;
    --rpar;
    --ipar;

    /* Function Body */
    ret_val = 0.;
    vmax = 0.;
    i__1 = *neq;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if ((d__1 = v[i__] / wt[i__], abs(d__1)) > vmax) {
	    vmax = (d__2 = v[i__] / wt[i__], abs(d__2));
	}
/* L10: */
    }
    if (vmax <= 0.) {
	goto L30;
    }
    sum = 0.;
    i__1 = *neq;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L20: */
/* Computing 2nd power */
	d__1 = v[i__] / wt[i__] / vmax;
	sum += d__1 * d__1;
    }
    ret_val = vmax * sqrt(sum / *neq);
L30:
    return ret_val;
/* ------END OF FUNCTION DDANRM------ */
} /* ddanrm_ */

