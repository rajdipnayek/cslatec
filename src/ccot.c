/* ccot.f -- translated by f2c (version 12.02.01).
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

static integer c__4 = 4;
static integer c__2 = 2;
static integer c__1 = 1;

/* DECK CCOT */
/* Complex */ void ccot_(complex * ret_val, complex *z__)
{
    /* Initialized data */

    static real sqeps = 0.f;

    /* System generated locals */
    real r__1, r__2;
    complex q__1;

    /* Local variables */
    static real x2, y2, den, sn2x;
    extern doublereal r1mach_(integer *);
    extern /* Subroutine */ int xerclr_(void), xermsg_(char *, char *, char *,
	     integer *, integer *, ftnlen, ftnlen, ftnlen);

/* ***BEGIN PROLOGUE  CCOT */
/* ***PURPOSE  Compute the cotangent. */
/* ***LIBRARY   SLATEC (FNLIB) */
/* ***CATEGORY  C4A */
/* ***TYPE      COMPLEX (COT-S, DCOT-D, CCOT-C) */
/* ***KEYWORDS  COTANGENT, ELEMENTARY FUNCTIONS, FNLIB, TRIGONOMETRIC */
/* ***AUTHOR  Fullerton, W., (LANL) */
/* ***DESCRIPTION */

/* CCOT(Z) calculates the complex trigonometric cotangent of Z. */

/* ***REFERENCES  (NONE) */
/* ***ROUTINES CALLED  R1MACH, XERCLR, XERMSG */
/* ***REVISION HISTORY  (YYMMDD) */
/*   770401  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890531  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ) */
/*   900326  Removed duplicate information from DESCRIPTION section. */
/*           (WRB) */
/* ***END PROLOGUE  CCOT */
/* ***FIRST EXECUTABLE STATEMENT  CCOT */
    if (sqeps == 0.f) {
	sqeps = sqrt(r1mach_(&c__4));
    }

    x2 = z__->r * 2.f;
    y2 = r_imag(z__) * 2.f;

    sn2x = sin(x2);
    xerclr_();

    den = cosh(y2) - cos(x2);
    if (den == 0.f) {
	xermsg_("SLATEC", "CCOT", "COT IS SINGULAR FOR INPUT Z (X IS 0 OR PI"
		" AND Y IS 0)", &c__2, &c__2, (ftnlen)6, (ftnlen)4, (ftnlen)53)
		;
    }

/* Computing MAX */
    r__1 = dabs(x2);
    if (dabs(den) > dmax(r__1,1.f) * sqeps) {
	goto L10;
    }
    xerclr_();
    xermsg_("SLATEC", "CCOT", "ANSWER LT HALF PRECISION, ABS(X) TOO BIG OR X"
	    " TOO NEAR 0 OR PI", &c__1, &c__1, (ftnlen)6, (ftnlen)4, (ftnlen)
	    62);

L10:
    r__1 = sn2x / den;
    r__2 = -sinh(y2) / den;
    q__1.r = r__1, q__1.i = r__2;
     ret_val->r = q__1.r,  ret_val->i = q__1.i;

    return ;
} /* ccot_ */

