/* catan.f -- translated by f2c (version 12.02.01).
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
static integer c__4 = 4;
static integer c__2 = 2;
static integer c__1 = 1;

/* DECK CATAN */
/* Complex */ void catan_(complex * ret_val, complex *z__)
{
    /* Initialized data */

    static real pi2 = 1.57079632679489661923f;
    static logical first = TRUE_;

    /* System generated locals */
    integer i__1;
    real r__1;
    complex q__1, q__2;

    /* Local variables */
    static integer i__;
    static real r__, x, y, r2;
    static complex z2;
    static real rmin, rmax, xans, yans, twoi, sqeps;
    extern doublereal r1mach_(integer *);
    extern /* Subroutine */ int xermsg_(char *, char *, char *, integer *, 
	    integer *, ftnlen, ftnlen, ftnlen);
    static integer nterms;

/* ***BEGIN PROLOGUE  CATAN */
/* ***PURPOSE  Compute the complex arc tangent. */
/* ***LIBRARY   SLATEC (FNLIB) */
/* ***CATEGORY  C4A */
/* ***TYPE      COMPLEX (CATAN-C) */
/* ***KEYWORDS  ARC TANGENT, ELEMENTARY FUNCTIONS, FNLIB, TRIGONOMETRIC */
/* ***AUTHOR  Fullerton, W., (LANL) */
/* ***DESCRIPTION */

/* CATAN(Z) calculates the complex trigonometric arc tangent of Z. */
/* The result is in units of radians, and the real part is in the first */
/* or fourth quadrant. */

/* ***REFERENCES  (NONE) */
/* ***ROUTINES CALLED  R1MACH, XERMSG */
/* ***REVISION HISTORY  (YYMMDD) */
/*   770801  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890531  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ) */
/*   900326  Removed duplicate information from DESCRIPTION section. */
/*           (WRB) */
/* ***END PROLOGUE  CATAN */
/* ***FIRST EXECUTABLE STATEMENT  CATAN */
    if (first) {
/* NTERMS = LOG(EPS)/LOG(RBND) WHERE RBND = 0.1 */
	nterms = log(r1mach_(&c__3)) * -.4343f + 1.f;
	sqeps = sqrt(r1mach_(&c__4));
	rmin = sqrt(r1mach_(&c__3) * 3.f);
	rmax = 1.f / r1mach_(&c__3);
    }
    first = FALSE_;

    r__ = c_abs(z__);
    if (r__ > .1f) {
	goto L30;
    }

     ret_val->r = z__->r,  ret_val->i = z__->i;
    if (r__ < rmin) {
	return ;
    }

     ret_val->r = 0.f,  ret_val->i = 0.f;
    q__1.r = z__->r * z__->r - z__->i * z__->i, q__1.i = z__->r * z__->i + 
	    z__->i * z__->r;
    z2.r = q__1.r, z2.i = q__1.i;
    i__1 = nterms;
    for (i__ = 1; i__ <= i__1; ++i__) {
	twoi = (real) ((nterms - i__ << 1) + 1);
	r__1 = 1.f / twoi;
	q__2.r = z2.r *  ret_val->r - z2.i *  ret_val->i, q__2.i = z2.r *  
		ret_val->i + z2.i *  ret_val->r;
	q__1.r = r__1 - q__2.r, q__1.i = -q__2.i;
	 ret_val->r = q__1.r,  ret_val->i = q__1.i;
/* L20: */
    }
    q__1.r = z__->r *  ret_val->r - z__->i *  ret_val->i, q__1.i = z__->r *  
	    ret_val->i + z__->i *  ret_val->r;
     ret_val->r = q__1.r,  ret_val->i = q__1.i;
    return ;

L30:
    if (r__ > rmax) {
	goto L50;
    }
    x = z__->r;
    y = r_imag(z__);
    r2 = r__ * r__;
    if (r2 == 1.f && x == 0.f) {
	xermsg_("SLATEC", "CATAN", "Z IS +I OR -I", &c__2, &c__2, (ftnlen)6, (
		ftnlen)5, (ftnlen)13);
    }
    if ((r__1 = r2 - 1.f, dabs(r__1)) > sqeps) {
	goto L40;
    }
    q__2.r = z__->r * z__->r - z__->i * z__->i, q__2.i = z__->r * z__->i + 
	    z__->i * z__->r;
    q__1.r = q__2.r + 1.f, q__1.i = q__2.i + 0.f;
    if (c_abs(&q__1) < sqeps) {
	xermsg_("SLATEC", "CATAN", "ANSWER LT HALF PRECISION, Z**2 CLOSE TO "
		"-1", &c__1, &c__1, (ftnlen)6, (ftnlen)5, (ftnlen)42);
    }

L40:
    xans = atan2(x * 2.f, 1.f - r2) * .5f;
    yans = log((r2 + y * 2.f + 1.f) / (r2 - y * 2.f + 1.f)) * .25f;
    q__1.r = xans, q__1.i = yans;
     ret_val->r = q__1.r,  ret_val->i = q__1.i;
    return ;

L50:
    q__1.r = pi2, q__1.i = 0.f;
     ret_val->r = q__1.r,  ret_val->i = q__1.i;
    if (z__->r < 0.f) {
	r__1 = -pi2;
	q__1.r = r__1, q__1.i = 0.f;
	 ret_val->r = q__1.r,  ret_val->i = q__1.i;
    }
    return ;

} /* catan_ */

