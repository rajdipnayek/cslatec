/* ccbrt.f -- translated by f2c (version 12.02.01).
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

/* DECK CCBRT */
/* Complex */ void ccbrt_(complex * ret_val, complex *z__)
{
    /* System generated locals */
    real r__1, r__2;
    complex q__1;

    /* Local variables */
    static real r__;
    extern doublereal carg_(complex *), cbrt_(real *);
    static real theta;

/* ***BEGIN PROLOGUE  CCBRT */
/* ***PURPOSE  Compute the cube root. */
/* ***LIBRARY   SLATEC (FNLIB) */
/* ***CATEGORY  C2 */
/* ***TYPE      COMPLEX (CBRT-S, DCBRT-D, CCBRT-C) */
/* ***KEYWORDS  CUBE ROOT, ELEMENTARY FUNCTIONS, FNLIB, ROOTS */
/* ***AUTHOR  Fullerton, W., (LANL) */
/* ***DESCRIPTION */

/* CCBRT(Z) calculates the complex cube root of Z.  The principal root */
/* for which -PI .LT. arg(Z) .LE. +PI is returned. */

/* ***REFERENCES  (NONE) */
/* ***ROUTINES CALLED  CARG, CBRT */
/* ***REVISION HISTORY  (YYMMDD) */
/*   770401  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890531  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/* ***END PROLOGUE  CCBRT */
/* ***FIRST EXECUTABLE STATEMENT  CCBRT */
    theta = carg_(z__) / 3.f;
    r__1 = c_abs(z__);
    r__ = cbrt_(&r__1);

    r__1 = r__ * cos(theta);
    r__2 = r__ * sin(theta);
    q__1.r = r__1, q__1.i = r__2;
     ret_val->r = q__1.r,  ret_val->i = q__1.i;

    return ;
} /* ccbrt_ */

