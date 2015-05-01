/* ccosh.f -- translated by f2c (version 12.02.01).
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

/* DECK CCOSH */
/* Complex */ void ccosh_(complex * ret_val, complex *z__)
{
    /* Initialized data */

    static complex ci = {0.f,1.f};

    /* System generated locals */
    complex q__1, q__2;

/* ***BEGIN PROLOGUE  CCOSH */
/* ***PURPOSE  Compute the complex hyperbolic cosine. */
/* ***LIBRARY   SLATEC (FNLIB) */
/* ***CATEGORY  C4C */
/* ***TYPE      COMPLEX (CCOSH-C) */
/* ***KEYWORDS  ELEMENTARY FUNCTIONS, FNLIB, HYPERBOLIC COSINE */
/* ***AUTHOR  Fullerton, W., (LANL) */
/* ***DESCRIPTION */

/* CCOSH(Z) calculates the complex hyperbolic cosine of Z. */

/* ***REFERENCES  (NONE) */
/* ***ROUTINES CALLED  (NONE) */
/* ***REVISION HISTORY  (YYMMDD) */
/*   770401  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890531  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/* ***END PROLOGUE  CCOSH */
/* ***FIRST EXECUTABLE STATEMENT  CCOSH */
    q__2.r = ci.r * z__->r - ci.i * z__->i, q__2.i = ci.r * z__->i + ci.i * 
	    z__->r;
    c_cos(&q__1, &q__2);
     ret_val->r = q__1.r,  ret_val->i = q__1.i;

    return ;
} /* ccosh_ */

