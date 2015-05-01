/* clog10.f -- translated by f2c (version 12.02.01).
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

/* DECK CLOG10 */
/* Complex */ void clog10_(complex * ret_val, complex *z__)
{
    /* Initialized data */

    static real aloge = .43429448190325182765f;

    /* System generated locals */
    complex q__1, q__2;

/* ***BEGIN PROLOGUE  CLOG10 */
/* ***PURPOSE  Compute the principal value of the complex base 10 */
/*            logarithm. */
/* ***LIBRARY   SLATEC (FNLIB) */
/* ***CATEGORY  C4B */
/* ***TYPE      COMPLEX (CLOG10-C) */
/* ***KEYWORDS  BASE TEN LOGARITHM, ELEMENTARY FUNCTIONS, FNLIB */
/* ***AUTHOR  Fullerton, W., (LANL) */
/* ***DESCRIPTION */

/* CLOG10(Z) calculates the principal value of the complex common */
/* or base 10 logarithm of Z for -PI .LT. arg(Z) .LE. +PI. */

/* ***REFERENCES  (NONE) */
/* ***ROUTINES CALLED  (NONE) */
/* ***REVISION HISTORY  (YYMMDD) */
/*   770401  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890531  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/* ***END PROLOGUE  CLOG10 */
/* ***FIRST EXECUTABLE STATEMENT  CLOG10 */
    c_log(&q__2, z__);
    q__1.r = aloge * q__2.r, q__1.i = aloge * q__2.i;
     ret_val->r = q__1.r,  ret_val->i = q__1.i;

    return ;
} /* clog10_ */

