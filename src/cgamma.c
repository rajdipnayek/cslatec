/* cgamma.f -- translated by f2c (version 12.02.01).
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

/* DECK CGAMMA */
/* Complex */ void cgamma_(complex * ret_val, complex *z__)
{
    /* System generated locals */
    complex q__1, q__2;

    /* Local variables */
    extern /* Complex */ void clngam_(complex *, complex *);

/* ***BEGIN PROLOGUE  CGAMMA */
/* ***PURPOSE  Compute the complete Gamma function. */
/* ***LIBRARY   SLATEC (FNLIB) */
/* ***CATEGORY  C7A */
/* ***TYPE      COMPLEX (GAMMA-S, DGAMMA-D, CGAMMA-C) */
/* ***KEYWORDS  COMPLETE GAMMA FUNCTION, FNLIB, SPECIAL FUNCTIONS */
/* ***AUTHOR  Fullerton, W., (LANL) */
/* ***DESCRIPTION */

/* CGAMMA(Z) calculates the complete gamma function for COMPLEX */
/* argument Z.  This is a preliminary version that is portable */
/* but not accurate. */

/* ***REFERENCES  (NONE) */
/* ***ROUTINES CALLED  CLNGAM */
/* ***REVISION HISTORY  (YYMMDD) */
/*   770701  DATE WRITTEN */
/*   861211  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/* ***END PROLOGUE  CGAMMA */
/* ***FIRST EXECUTABLE STATEMENT  CGAMMA */
    clngam_(&q__2, z__);
    c_exp(&q__1, &q__2);
     ret_val->r = q__1.r,  ret_val->i = q__1.i;

    return ;
} /* cgamma_ */

