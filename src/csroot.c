/* csroot.f -- translated by f2c (version 12.02.01).
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

/* DECK CSROOT */
/* Subroutine */ int csroot_(real *xr, real *xi, real *yr, real *yi)
{
    /* Local variables */
    static real s, ti, tr;
    extern doublereal pythag_(real *, real *);

/* ***BEGIN PROLOGUE  CSROOT */
/* ***SUBSIDIARY */
/* ***PURPOSE  Compute the complex square root of a complex number. */
/* ***LIBRARY   SLATEC */
/* ***TYPE      SINGLE PRECISION (CSROOT-S) */
/* ***AUTHOR  (UNKNOWN) */
/* ***DESCRIPTION */

/*     (YR,YI) = complex sqrt(XR,XI) */

/* ***SEE ALSO  EISDOC */
/* ***ROUTINES CALLED  PYTHAG */
/* ***REVISION HISTORY  (YYMMDD) */
/*   811101  DATE WRITTEN */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900402  Added TYPE section.  (WRB) */
/* ***END PROLOGUE  CSROOT */

/*     BRANCH CHOSEN SO THAT YR .GE. 0.0 AND SIGN(YI) .EQ. SIGN(XI) */
/* ***FIRST EXECUTABLE STATEMENT  CSROOT */
    tr = *xr;
    ti = *xi;
    s = sqrt((pythag_(&tr, &ti) + dabs(tr)) * .5f);
    if (tr >= 0.f) {
	*yr = s;
    }
    if (ti < 0.f) {
	s = -s;
    }
    if (tr <= 0.f) {
	*yi = s;
    }
    if (tr < 0.f) {
	*yr = ti / *yi * .5f;
    }
    if (tr > 0.f) {
	*yi = ti / *yr * .5f;
    }
    return 0;
} /* csroot_ */

