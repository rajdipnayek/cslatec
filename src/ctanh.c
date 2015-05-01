/* ctanh.f -- translated by f2c (version 12.02.01).
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

/* DECK CTANH */
/* Complex */ void ctanh_(complex * ret_val, complex *z__)
{
    /* Initialized data */

    static complex ci = {0.f,1.f};

    /* System generated locals */
    complex q__1, q__2, q__3, q__4;

    /* Local variables */
    extern /* Complex */ void ctan_(complex *, complex *);

/* ***BEGIN PROLOGUE  CTANH */
/* ***PURPOSE  Compute the complex hyperbolic tangent. */
/* ***LIBRARY   SLATEC (FNLIB) */
/* ***CATEGORY  C4C */
/* ***TYPE      COMPLEX (CTANH-C) */
/* ***KEYWORDS  ELEMENTARY FUNCTIONS, FNLIB, HYPERBOLIC TANGENT */
/* ***AUTHOR  Fullerton, W., (LANL) */
/* ***DESCRIPTION */

/* CTANH(Z) calculates the complex hyperbolic tangent of complex */
/* argument Z.  Z is in units of radians. */

/* ***REFERENCES  (NONE) */
/* ***ROUTINES CALLED  CTAN */
/* ***REVISION HISTORY  (YYMMDD) */
/*   770401  DATE WRITTEN */
/*   861211  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/* ***END PROLOGUE  CTANH */
/* ***FIRST EXECUTABLE STATEMENT  CTANH */
    q__2.r = -ci.r, q__2.i = -ci.i;
    q__4.r = ci.r * z__->r - ci.i * z__->i, q__4.i = ci.r * z__->i + ci.i * 
	    z__->r;
    ctan_(&q__3, &q__4);
    q__1.r = q__2.r * q__3.r - q__2.i * q__3.i, q__1.i = q__2.r * q__3.i + 
	    q__2.i * q__3.r;
     ret_val->r = q__1.r,  ret_val->i = q__1.i;

    return ;
} /* ctanh_ */

