/* c9lgmc.f -- translated by f2c (version 12.02.01).
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
static integer c__1 = 1;
static integer c__2 = 2;
static complex c_b18 = {1.f,0.f};

/* DECK C9LGMC */
/* Complex */ void c9lgmc_(complex * ret_val, complex *zin)
{
    /* Initialized data */

    static real bern[11] = { .083333333333333333f,-.0027777777777777778f,
	    7.9365079365079365e-4f,-5.9523809523809524e-4f,
	    8.4175084175084175e-4f,-.0019175269175269175f,
	    .0064102564102564103f,-.029550653594771242f,.17964437236883057f,
	    -1.3924322169059011f,13.402864044168392f };
    static logical first = TRUE_;

    /* System generated locals */
    integer i__1, i__2;
    real r__1, r__2;
    doublereal d__1, d__2;
    complex q__1, q__2;

    /* Local variables */
    static integer i__;
    static real x, y;
    static complex z__;
    static integer ndx;
    static real xbig, xmax;
    static complex z2inv;
    static real cabsz, bound;
    static integer nterm;
    extern doublereal r1mach_(integer *);
    extern /* Subroutine */ int xermsg_(char *, char *, char *, integer *, 
	    integer *, ftnlen, ftnlen, ftnlen);

/* ***BEGIN PROLOGUE  C9LGMC */
/* ***SUBSIDIARY */
/* ***PURPOSE  Compute the log gamma correction factor so that */
/*            LOG(CGAMMA(Z)) = 0.5*LOG(2.*PI) + (Z-0.5)*LOG(Z) - Z */
/*            + C9LGMC(Z). */
/* ***LIBRARY   SLATEC (FNLIB) */
/* ***CATEGORY  C7A */
/* ***TYPE      COMPLEX (R9LGMC-S, D9LGMC-D, C9LGMC-C) */
/* ***KEYWORDS  COMPLETE GAMMA FUNCTION, CORRECTION TERM, FNLIB, */
/*             LOG GAMMA, LOGARITHM, SPECIAL FUNCTIONS */
/* ***AUTHOR  Fullerton, W., (LANL) */
/* ***DESCRIPTION */

/* Compute the LOG GAMMA correction term for large ABS(Z) when REAL(Z) */
/* .GE. 0.0 and for large ABS(AIMAG(Y)) when REAL(Z) .LT. 0.0.  We find */
/* C9LGMC so that */
/*   LOG(Z) = 0.5*LOG(2.*PI) + (Z-0.5)*LOG(Z) - Z + C9LGMC(Z) */

/* ***REFERENCES  (NONE) */
/* ***ROUTINES CALLED  R1MACH, XERMSG */
/* ***REVISION HISTORY  (YYMMDD) */
/*   780401  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890531  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ) */
/*   900326  Removed duplicate information from DESCRIPTION section. */
/*           (WRB) */
/*   900720  Routine changed from user-callable to subsidiary.  (WRB) */
/* ***END PROLOGUE  C9LGMC */
/* ***FIRST EXECUTABLE STATEMENT  C9LGMC */
    if (first) {
	nterm = log(r1mach_(&c__3)) * -.3f;
	d__1 = (doublereal) (r1mach_(&c__3) * .1f);
	d__2 = (doublereal) (-1.f / ((nterm << 1) - 1));
	bound = nterm * .117f * pow_dd(&d__1, &d__2);
	xbig = 1.f / sqrt(r1mach_(&c__3));
/* Computing MIN */
	r__1 = log(r1mach_(&c__2) / 12.f), r__2 = -log(r1mach_(&c__1) * 12.f);
	xmax = exp((dmin(r__1,r__2)));
    }
    first = FALSE_;

    z__.r = zin->r, z__.i = zin->i;
    x = z__.r;
    y = r_imag(&z__);
    cabsz = c_abs(&z__);

    if (x < 0.f && dabs(y) < bound) {
	xermsg_("SLATEC", "C9LGMC", "NOT VALID FOR NEGATIVE REAL(Z) AND SMAL"
		"L ABS(AIMAG(Z))", &c__2, &c__2, (ftnlen)6, (ftnlen)6, (ftnlen)
		54);
    }
    if (cabsz < bound) {
	xermsg_("SLATEC", "C9LGMC", "NOT VALID FOR SMALL ABS(Z)", &c__3, &
		c__2, (ftnlen)6, (ftnlen)6, (ftnlen)26);
    }

    if (cabsz >= xmax) {
	goto L50;
    }

    if (cabsz >= xbig) {
	q__2.r = z__.r * 12.f, q__2.i = z__.i * 12.f;
	c_div(&q__1, &c_b18, &q__2);
	 ret_val->r = q__1.r,  ret_val->i = q__1.i;
    }
    if (cabsz >= xbig) {
	return ;
    }

    pow_ci(&q__2, &z__, &c__2);
    c_div(&q__1, &c_b18, &q__2);
    z2inv.r = q__1.r, z2inv.i = q__1.i;
     ret_val->r = 0.f,  ret_val->i = 0.f;
    i__1 = nterm;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ndx = nterm + 1 - i__;
	i__2 = ndx - 1;
	q__2.r =  ret_val->r * z2inv.r -  ret_val->i * z2inv.i, q__2.i =  
		ret_val->r * z2inv.i +  ret_val->i * z2inv.r;
	q__1.r = bern[i__2] + q__2.r, q__1.i = q__2.i;
	 ret_val->r = q__1.r,  ret_val->i = q__1.i;
/* L40: */
    }

    c_div(&q__1,  ret_val, &z__);
     ret_val->r = q__1.r,  ret_val->i = q__1.i;
    return ;

L50:
     ret_val->r = 0.f,  ret_val->i = 0.f;
    xermsg_("SLATEC", "C9LGMC", "Z SO BIG C9LGMC UNDERFLOWS", &c__1, &c__1, (
	    ftnlen)6, (ftnlen)6, (ftnlen)26);
    return ;

} /* c9lgmc_ */

