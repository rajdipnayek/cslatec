/* cexprl.f -- translated by f2c (version 12.02.01).
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

static integer c__3 = 3;

/* DECK CEXPRL */
/* Complex */ void cexprl_(complex * ret_val, complex *z__)
{
    /* Initialized data */

    static logical first = TRUE_;

    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1;
    complex q__1, q__2, q__3;

    /* Local variables */
    static integer i__;
    static real r__, xn, xln, rbnd;
    extern doublereal r1mach_(integer *);
    static real alneps;
    static integer nterms;

/* ***BEGIN PROLOGUE  CEXPRL */
/* ***PURPOSE  Calculate the relative error exponential (EXP(X)-1)/X. */
/* ***LIBRARY   SLATEC (FNLIB) */
/* ***CATEGORY  C4B */
/* ***TYPE      COMPLEX (EXPREL-S, DEXPRL-D, CEXPRL-C) */
/* ***KEYWORDS  ELEMENTARY FUNCTIONS, EXPONENTIAL, FIRST ORDER, FNLIB */
/* ***AUTHOR  Fullerton, W., (LANL) */
/* ***DESCRIPTION */

/* Evaluate  (EXP(Z)-1)/Z .  For small ABS(Z), we use the Taylor */
/* series.  We could instead use the expression */
/*        CEXPRL(Z) = (EXP(X)*EXP(I*Y)-1)/Z */
/*                  = (X*EXPREL(X) * (1 - 2*SIN(Y/2)**2) - 2*SIN(Y/2)**2 */
/*                                    + I*SIN(Y)*(1+X*EXPREL(X))) / Z */

/* ***REFERENCES  (NONE) */
/* ***ROUTINES CALLED  R1MACH */
/* ***REVISION HISTORY  (YYMMDD) */
/*   770801  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890531  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/* ***END PROLOGUE  CEXPRL */
/* ***FIRST EXECUTABLE STATEMENT  CEXPRL */
    if (first) {
	alneps = log(r1mach_(&c__3));
	xn = 3.72f - alneps * .3f;
	xln = log((xn + 1.f) / 1.36f);
	nterms = xn - (xn * xln + alneps) / (xln + 1.36f) + 1.5f;
	rbnd = r1mach_(&c__3);
    }
    first = FALSE_;

    r__ = c_abs(z__);
    if (r__ > .5f) {
	c_exp(&q__3, z__);
	q__2.r = q__3.r - 1.f, q__2.i = q__3.i;
	c_div(&q__1, &q__2, z__);
	 ret_val->r = q__1.r,  ret_val->i = q__1.i;
    }
    if (r__ > .5f) {
	return ;
    }

     ret_val->r = 1.f,  ret_val->i = 0.f;
    if (r__ < rbnd) {
	return ;
    }

     ret_val->r = 0.f,  ret_val->i = 0.f;
    i__1 = nterms;
    for (i__ = 1; i__ <= i__1; ++i__) {
	q__3.r =  ret_val->r * z__->r -  ret_val->i * z__->i, q__3.i =  
		ret_val->r * z__->i +  ret_val->i * z__->r;
	i__2 = nterms + 2 - i__;
	d__1 = (doublereal) i__2;
	q__2.r = q__3.r / d__1, q__2.i = q__3.i / d__1;
	q__1.r = q__2.r + 1.f, q__1.i = q__2.i;
	 ret_val->r = q__1.r,  ret_val->i = q__1.i;
/* L20: */
    }

    return ;
} /* cexprl_ */

