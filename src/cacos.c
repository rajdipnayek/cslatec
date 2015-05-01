/* cacos.f -- translated by f2c (version 12.02.01).
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

/* DECK CACOS */
/* Complex */ void cacos_(complex * ret_val, complex *z__)
{
    /* Initialized data */

    static real pi2 = 1.57079632679489661923f;

    /* System generated locals */
    complex q__1, q__2;

    /* Local variables */
    extern /* Complex */ void casin_(complex *, complex *);

/* ***BEGIN PROLOGUE  CACOS */
/* ***PURPOSE  Compute the complex arc cosine. */
/* ***LIBRARY   SLATEC (FNLIB) */
/* ***CATEGORY  C4A */
/* ***TYPE      COMPLEX (CACOS-C) */
/* ***KEYWORDS  ARC COSINE, ELEMENTARY FUNCTIONS, FNLIB, TRIGONOMETRIC */
/* ***AUTHOR  Fullerton, W., (LANL) */
/* ***DESCRIPTION */

/* CACOS(Z) calculates the complex trigonometric arc cosine of Z. */
/* The result is in units of radians, and the real part is in the */
/* first or second quadrant. */

/* ***REFERENCES  (NONE) */
/* ***ROUTINES CALLED  CASIN */
/* ***REVISION HISTORY  (YYMMDD) */
/*   770401  DATE WRITTEN */
/*   861211  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/* ***END PROLOGUE  CACOS */
/* ***FIRST EXECUTABLE STATEMENT  CACOS */
    casin_(&q__2, z__);
    q__1.r = pi2 - q__2.r, q__1.i = -q__2.i;
     ret_val->r = q__1.r,  ret_val->i = q__1.i;

    return ;
} /* cacos_ */

