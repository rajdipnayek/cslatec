/* c9ln2r.f -- translated by f2c (version 12.02.01).
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

/* DECK C9LN2R */
/* Complex */ void c9ln2r_(complex * ret_val, complex *z__)
{
    /* System generated locals */
    real r__1, r__2, r__3;
    complex q__1, q__2, q__3, q__4, q__5, q__6, q__7, q__8;

    /* Local variables */
    static real x, y, xz, yz, y1x, arg, cabsz;
    extern doublereal r9atn1_(real *);
    static real rpart;
    extern doublereal r9ln2r_(real *);
    static real aipart;

/* ***BEGIN PROLOGUE  C9LN2R */
/* ***SUBSIDIARY */
/* ***PURPOSE  Evaluate LOG(1+Z) from second order relative accuracy so */
/*            that  LOG(1+Z) = Z - Z**2/2 + Z**3*C9LN2R(Z). */
/* ***LIBRARY   SLATEC (FNLIB) */
/* ***CATEGORY  C4B */
/* ***TYPE      COMPLEX (R9LN2R-S, D9LN2R-D, C9LN2R-C) */
/* ***KEYWORDS  ELEMENTARY FUNCTIONS, FNLIB, LOGARITHM, SECOND ORDER */
/* ***AUTHOR  Fullerton, W., (LANL) */
/* ***DESCRIPTION */

/* Evaluate  LOG(1+Z)  from 2-nd order with relative error accuracy so */
/* that     LOG(1+Z) = Z - Z**2/2 + Z**3*C9LN2R(Z). */

/* Now  LOG(1+Z) = 0.5*LOG(1+2*X+ABS(Z)**2) + I*CARG(1+Z), */
/* where X = REAL(Z)  and  Y = AIMAG(Z). */
/* We find */
/*     Z**3 * C9LN2R(Z) = -X*ABS(Z)**2 - 0.25*ABS(Z)**4 */
/*        + (2*X+ABS(Z)**2)**3 * R9LN2R(2*X+ABS(Z)**2) */
/*        + I * (CARG(1+Z) + (X-1)*Y) */
/* The imaginary part must be evaluated carefully as */
/*     (ATAN(Y/(1+X)) - Y/(1+X)) + Y/(1+X) - (1-X)*Y */
/*       = (Y/(1+X))**3 * R9ATN1(Y/(1+X)) + X**2*Y/(1+X) */

/* Now we divide through by Z**3 carefully.  Write */
/*     1/Z**3 = (X-I*Y)/ABS(Z)**3 * (1/ABS(Z)**3) */
/* then   C9LN2R(Z) = ((X-I*Y)/ABS(Z))**3 * (-X/ABS(Z) - ABS(Z)/4 */
/*        + 0.5*((2*X+ABS(Z)**2)/ABS(Z))**3 * R9LN2R(2*X+ABS(Z)**2) */
/*        + I*Y/(ABS(Z)*(1+X)) * ((X/ABS(Z))**2 + */
/*          + (Y/(ABS(Z)*(1+X)))**2 * R9ATN1(Y/(1+X)) ) ) */

/* If we let  XZ = X/ABS(Z)  and  YZ = Y/ABS(Z)  we may write */
/*     C9LN2R(Z) = (XZ-I*YZ)**3 * (-XZ - ABS(Z)/4 */
/*        + 0.5*(2*XZ+ABS(Z))**3 * R9LN2R(2*X+ABS(Z)**2) */
/*        + I*YZ/(1+X) * (XZ**2 + (YZ/(1+X))**2*R9ATN1(Y/(1+X)) )) */

/* ***REFERENCES  (NONE) */
/* ***ROUTINES CALLED  R9ATN1, R9LN2R */
/* ***REVISION HISTORY  (YYMMDD) */
/*   780401  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890531  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900720  Routine changed from user-callable to subsidiary.  (WRB) */
/* ***END PROLOGUE  C9LN2R */
/* ***FIRST EXECUTABLE STATEMENT  C9LN2R */
    x = z__->r;
    y = r_imag(z__);

    cabsz = c_abs(z__);
    if (cabsz > .8125f) {
	goto L20;
    }

     ret_val->r = .33333333333333331f,  ret_val->i = 0.f;
    if (cabsz == 0.f) {
	return ;
    }

    xz = x / cabsz;
    yz = y / cabsz;

    arg = xz * 2.f + cabsz;
/* Computing 3rd power */
    r__1 = arg;
    r__2 = cabsz * arg;
    rpart = r__1 * (r__1 * r__1) * .5f * r9ln2r_(&r__2) - xz - cabsz * .25f;
    y1x = yz / (x + 1.f);
/* Computing 2nd power */
    r__1 = xz;
/* Computing 2nd power */
    r__2 = y1x;
    r__3 = cabsz * y1x;
    aipart = y1x * (r__1 * r__1 + r__2 * r__2 * r9atn1_(&r__3));

    r__1 = -yz;
    q__3.r = xz, q__3.i = r__1;
    pow_ci(&q__2, &q__3, &c__3);
    q__4.r = rpart, q__4.i = aipart;
    q__1.r = q__2.r * q__4.r - q__2.i * q__4.i, q__1.i = q__2.r * q__4.i + 
	    q__2.i * q__4.r;
     ret_val->r = q__1.r,  ret_val->i = q__1.i;
    return ;

L20:
    q__4.r = z__->r + 1.f, q__4.i = z__->i;
    c_log(&q__3, &q__4);
    q__7.r = z__->r * .5f, q__7.i = z__->i * .5f;
    q__6.r = 1.f - q__7.r, q__6.i = -q__7.i;
    q__5.r = z__->r * q__6.r - z__->i * q__6.i, q__5.i = z__->r * q__6.i + 
	    z__->i * q__6.r;
    q__2.r = q__3.r - q__5.r, q__2.i = q__3.i - q__5.i;
    pow_ci(&q__8, z__, &c__3);
    c_div(&q__1, &q__2, &q__8);
     ret_val->r = q__1.r,  ret_val->i = q__1.i;
    return ;

} /* c9ln2r_ */

