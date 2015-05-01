/* sindg.f -- translated by f2c (version 12.02.01).
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

static real c_b2 = 90.f;
static real c_b3 = 1.f;

/* DECK SINDG */
doublereal sindg_(real *x)
{
    /* Initialized data */

    static real raddeg = .017453292519943296f;

    /* System generated locals */
    real ret_val;

    /* Local variables */
    static integer n;

/* ***BEGIN PROLOGUE  SINDG */
/* ***PURPOSE  Compute the sine of an argument in degrees. */
/* ***LIBRARY   SLATEC (FNLIB) */
/* ***CATEGORY  C4A */
/* ***TYPE      SINGLE PRECISION (SINDG-S, DSINDG-D) */
/* ***KEYWORDS  DEGREES, ELEMENTARY FUNCTIONS, FNLIB, SINE, TRIGONOMETRIC */
/* ***AUTHOR  Fullerton, W., (LANL) */
/* ***DESCRIPTION */

/* SINDG(X) evaluates the single precision sine of X where */
/* X is in degrees. */

/* ***REFERENCES  (NONE) */
/* ***ROUTINES CALLED  (NONE) */
/* ***REVISION HISTORY  (YYMMDD) */
/*   770601  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890531  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/* ***END PROLOGUE  SINDG */
/* JUNE 1977 EDITION.   W. FULLERTON, C3, LOS ALAMOS SCIENTIFIC LAB. */

/* ***FIRST EXECUTABLE STATEMENT  SINDG */
    ret_val = sin(raddeg * *x);

    if (r_mod(x, &c_b2) != 0.f) {
	return ret_val;
    }
    n = dabs(*x) / 90.f + .5f;
    n %= 2;
    if (n == 0) {
	ret_val = 0.f;
    }
    if (n == 1) {
	ret_val = r_sign(&c_b3, &ret_val);
    }

    return ret_val;
} /* sindg_ */

