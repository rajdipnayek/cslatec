/* vnwrms.f -- translated by f2c (version 12.02.01).
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

/* DECK VNWRMS */
doublereal vnwrms_(integer *n, real *v, real *w)
{
    /* System generated locals */
    integer i__1;
    real ret_val, r__1;

    /* Local variables */
    static integer i__;
    static real sum;

/* ***BEGIN PROLOGUE  VNWRMS */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to DEBDF */
/* ***LIBRARY   SLATEC */
/* ***TYPE      SINGLE PRECISION (VNWRMS-S, DVNRMS-D) */
/* ***AUTHOR  (UNKNOWN) */
/* ***DESCRIPTION */

/*   VNWRMS computes a weighted root-mean-square vector norm for the */
/*   integrator package DEBDF. */

/* ***SEE ALSO  DEBDF */
/* ***ROUTINES CALLED  (NONE) */
/* ***REVISION HISTORY  (YYMMDD) */
/*   800901  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890831  Modified array declarations.  (WRB) */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900328  Added TYPE section.  (WRB) */
/* ***END PROLOGUE  VNWRMS */


/* LLL. OPTIMIZE */
/* ----------------------------------------------------------------------- */
/* THIS FUNCTION ROUTINE COMPUTES THE WEIGHTED ROOT-MEAN-SQUARE NORM */
/* OF THE VECTOR OF LENGTH N CONTAINED IN THE ARRAY V, WITH WEIGHTS */
/* CONTAINED IN THE ARRAY W OF LENGTH N.. */
/*   VNWRMS = SQRT( (1/N) * SUM( V(I)/W(I) )**2 ) */
/* ----------------------------------------------------------------------- */
/* ***FIRST EXECUTABLE STATEMENT  VNWRMS */
    /* Parameter adjustments */
    --w;
    --v;

    /* Function Body */
    sum = 0.f;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L10: */
/* Computing 2nd power */
	r__1 = v[i__] / w[i__];
	sum += r__1 * r__1;
    }
    ret_val = sqrt(sum / *n);
    return ret_val;
/* ----------------------- END OF FUNCTION VNWRMS ------------------------ */
} /* vnwrms_ */

