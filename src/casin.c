/* casin.f -- translated by f2c (version 12.02.01).
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

/* DECK CASIN */
/* Complex */ void casin_(complex * ret_val, complex *zinp)
{
    /* Initialized data */

    static real pi2 = 1.57079632679489661923f;
    static real pi = 3.14159265358979324f;
    static complex ci = {0.f,1.f};
    static logical first = TRUE_;

    /* System generated locals */
    integer i__1;
    real r__1, r__2;
    complex q__1, q__2, q__3, q__4, q__5, q__6, q__7;

    /* Local variables */
    static integer i__;
    static real r__;
    static complex z__, z2;
    static real rmin, twoi;
    static complex sqzp1;
    extern doublereal r1mach_(integer *);
    static integer nterms;

/* ***BEGIN PROLOGUE  CASIN */
/* ***PURPOSE  Compute the complex arc sine. */
/* ***LIBRARY   SLATEC (FNLIB) */
/* ***CATEGORY  C4A */
/* ***TYPE      COMPLEX (CASIN-C) */
/* ***KEYWORDS  ARC SINE, ELEMENTARY FUNCTIONS, FNLIB, TRIGONOMETRIC */
/* ***AUTHOR  Fullerton, W., (LANL) */
/* ***DESCRIPTION */

/* CASIN(ZINP) calculates the complex trigonometric arc sine of ZINP. */
/* The result is in units of radians, and the real part is in the first */
/* or fourth quadrant. */

/* ***REFERENCES  (NONE) */
/* ***ROUTINES CALLED  R1MACH */
/* ***REVISION HISTORY  (YYMMDD) */
/*   770701  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890531  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/* ***END PROLOGUE  CASIN */
/* ***FIRST EXECUTABLE STATEMENT  CASIN */
    if (first) {
/* NTERMS = LOG(EPS)/LOG(RMAX)  WHERE RMAX = 0.1 */
	nterms = log(r1mach_(&c__3)) * -.4343f;
	rmin = sqrt(r1mach_(&c__3) * 6.f);
    }
    first = FALSE_;

    z__.r = zinp->r, z__.i = zinp->i;
    r__ = c_abs(&z__);
    if (r__ > .1f) {
	goto L30;
    }

     ret_val->r = z__.r,  ret_val->i = z__.i;
    if (r__ < rmin) {
	return ;
    }

     ret_val->r = 0.f,  ret_val->i = 0.f;
    q__1.r = z__.r * z__.r - z__.i * z__.i, q__1.i = z__.r * z__.i + z__.i * 
	    z__.r;
    z2.r = q__1.r, z2.i = q__1.i;
    i__1 = nterms;
    for (i__ = 1; i__ <= i__1; ++i__) {
	twoi = (real) ((nterms - i__ << 1) + 1);
	r__1 = 1.f / twoi;
	q__4.r = twoi *  ret_val->r, q__4.i = twoi *  ret_val->i;
	q__3.r = q__4.r * z2.r - q__4.i * z2.i, q__3.i = q__4.r * z2.i + 
		q__4.i * z2.r;
	r__2 = twoi + 1.f;
	q__2.r = q__3.r / r__2, q__2.i = q__3.i / r__2;
	q__1.r = r__1 + q__2.r, q__1.i = q__2.i;
	 ret_val->r = q__1.r,  ret_val->i = q__1.i;
/* L20: */
    }
    q__1.r = z__.r *  ret_val->r - z__.i *  ret_val->i, q__1.i = z__.r *  
	    ret_val->i + z__.i *  ret_val->r;
     ret_val->r = q__1.r,  ret_val->i = q__1.i;
    return ;

L30:
    if (zinp->r < 0.f) {
	q__1.r = -zinp->r, q__1.i = -zinp->i;
	z__.r = q__1.r, z__.i = q__1.i;
    }

    q__2.r = z__.r + 1.f, q__2.i = z__.i;
    c_sqrt(&q__1, &q__2);
    sqzp1.r = q__1.r, sqzp1.i = q__1.i;
    if (r_imag(&sqzp1) < 0.f) {
	q__1.r = -sqzp1.r, q__1.i = -sqzp1.i;
	sqzp1.r = q__1.r, sqzp1.i = q__1.i;
    }
    q__7.r = z__.r - 1.f, q__7.i = z__.i;
    c_sqrt(&q__6, &q__7);
    q__5.r = sqzp1.r * q__6.r - sqzp1.i * q__6.i, q__5.i = sqzp1.r * q__6.i + 
	    sqzp1.i * q__6.r;
    q__4.r = z__.r + q__5.r, q__4.i = z__.i + q__5.i;
    c_log(&q__3, &q__4);
    q__2.r = ci.r * q__3.r - ci.i * q__3.i, q__2.i = ci.r * q__3.i + ci.i * 
	    q__3.r;
    q__1.r = pi2 - q__2.r, q__1.i = -q__2.i;
     ret_val->r = q__1.r,  ret_val->i = q__1.i;

    if ( ret_val->r > pi2) {
	q__1.r = pi -  ret_val->r, q__1.i = - ret_val->i;
	 ret_val->r = q__1.r,  ret_val->i = q__1.i;
    }
    if ( ret_val->r <= -pi2) {
	r__1 = -pi;
	q__1.r = r__1 -  ret_val->r, q__1.i = - ret_val->i;
	 ret_val->r = q__1.r,  ret_val->i = q__1.i;
    }
    if (zinp->r < 0.f) {
	q__1.r = - ret_val->r, q__1.i = - ret_val->i;
	 ret_val->r = q__1.r,  ret_val->i = q__1.i;
    }

    return ;
} /* casin_ */

