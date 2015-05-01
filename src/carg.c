/* carg.f -- translated by f2c (version 12.02.01).
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

/* DECK CARG */
doublereal carg_(complex *z__)
{
    /* System generated locals */
    real ret_val;

/* ***BEGIN PROLOGUE  CARG */
/* ***PURPOSE  Compute the argument of a complex number. */
/* ***LIBRARY   SLATEC (FNLIB) */
/* ***CATEGORY  A4A */
/* ***TYPE      COMPLEX (CARG-C) */
/* ***KEYWORDS  ARGUMENT OF A COMPLEX NUMBER, ELEMENTARY FUNCTIONS, FNLIB */
/* ***AUTHOR  Fullerton, W., (LANL) */
/* ***DESCRIPTION */

/* CARG(Z) calculates the argument of the complex number Z.  Note */
/* that CARG returns a real result.  If Z = X+iY, then CARG is ATAN(Y/X), */
/* except when both X and Y are zero, in which case the result */
/* will be zero. */

/* ***REFERENCES  (NONE) */
/* ***ROUTINES CALLED  (NONE) */
/* ***REVISION HISTORY  (YYMMDD) */
/*   770401  DATE WRITTEN */
/*   861211  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/* ***END PROLOGUE  CARG */
/* ***FIRST EXECUTABLE STATEMENT  CARG */
    ret_val = 0.f;
    if (z__->r != 0.f || r_imag(z__) != 0.f) {
	ret_val = atan2(r_imag(z__), (real) z__->r);
    }

    return ret_val;
} /* carg_ */

