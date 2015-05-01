/* casinh.f -- translated by f2c (version 12.02.01).
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

/* DECK CASINH */
/* Complex */ void casinh_(complex * ret_val, complex *z__)
{
    /* Initialized data */

    static complex ci = {0.f,1.f};

    /* System generated locals */
    complex q__1, q__2, q__3, q__4;

    /* Local variables */
    extern /* Complex */ void casin_(complex *, complex *);

/* ***BEGIN PROLOGUE  CASINH */
/* ***PURPOSE  Compute the arc hyperbolic sine. */
/* ***LIBRARY   SLATEC (FNLIB) */
/* ***CATEGORY  C4C */
/* ***TYPE      COMPLEX (ASINH-S, DASINH-D, CASINH-C) */
/* ***KEYWORDS  ARC HYPERBOLIC SINE, ASINH, ELEMENTARY FUNCTIONS, FNLIB, */
/*             INVERSE HYPERBOLIC SINE */
/* ***AUTHOR  Fullerton, W., (LANL) */
/* ***DESCRIPTION */

/* CASINH(Z) calculates the complex arc hyperbolic sine of Z. */

/* ***REFERENCES  (NONE) */
/* ***ROUTINES CALLED  CASIN */
/* ***REVISION HISTORY  (YYMMDD) */
/*   770401  DATE WRITTEN */
/*   861211  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/* ***END PROLOGUE  CASINH */
/* ***FIRST EXECUTABLE STATEMENT  CASINH */
    q__2.r = -ci.r, q__2.i = -ci.i;
    q__4.r = ci.r * z__->r - ci.i * z__->i, q__4.i = ci.r * z__->i + ci.i * 
	    z__->r;
    casin_(&q__3, &q__4);
    q__1.r = q__2.r * q__3.r - q__2.i * q__3.i, q__1.i = q__2.r * q__3.i + 
	    q__2.i * q__3.r;
     ret_val->r = q__1.r,  ret_val->i = q__1.i;

    return ;
} /* casinh_ */

