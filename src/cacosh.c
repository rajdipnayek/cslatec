/* cacosh.f -- translated by f2c (version 12.02.01).
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

/* DECK CACOSH */
/* Complex */ void cacosh_(complex * ret_val, complex *z__)
{
    /* Initialized data */

    static complex ci = {0.f,1.f};

    /* System generated locals */
    complex q__1, q__2;

    /* Local variables */
    extern /* Complex */ void cacos_(complex *, complex *);

/* ***BEGIN PROLOGUE  CACOSH */
/* ***PURPOSE  Compute the arc hyperbolic cosine. */
/* ***LIBRARY   SLATEC (FNLIB) */
/* ***CATEGORY  C4C */
/* ***TYPE      COMPLEX (ACOSH-S, DACOSH-D, CACOSH-C) */
/* ***KEYWORDS  ACOSH, ARC HYPERBOLIC COSINE, ELEMENTARY FUNCTIONS, FNLIB, */
/*             INVERSE HYPERBOLIC COSINE */
/* ***AUTHOR  Fullerton, W., (LANL) */
/* ***DESCRIPTION */

/* CACOSH(Z) calculates the complex arc hyperbolic cosine of Z. */

/* ***REFERENCES  (NONE) */
/* ***ROUTINES CALLED  CACOS */
/* ***REVISION HISTORY  (YYMMDD) */
/*   770401  DATE WRITTEN */
/*   861211  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/* ***END PROLOGUE  CACOSH */
/* ***FIRST EXECUTABLE STATEMENT  CACOSH */
    cacos_(&q__2, z__);
    q__1.r = ci.r * q__2.r - ci.i * q__2.i, q__1.i = ci.r * q__2.i + ci.i * 
	    q__2.r;
     ret_val->r = q__1.r,  ret_val->i = q__1.i;

    return ;
} /* cacosh_ */

