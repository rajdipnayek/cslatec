/* ctan.f -- translated by f2c (version 12.02.01).
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

/* DECK CTAN */
/* Complex */ void ctan_(complex * ret_val, complex *z__)
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

/* ***BEGIN PROLOGUE  CTAN */
/* ***PURPOSE  Compute the complex tangent. */
/* ***LIBRARY   SLATEC (FNLIB) */
/* ***CATEGORY  C4A */
/* ***TYPE      COMPLEX (CTAN-C) */
/* ***KEYWORDS  ELEMENTARY FUNCTIONS, FNLIB, TANGENT, TRIGONOMETRIC */
/* ***AUTHOR  Fullerton, W., (LANL) */
/* ***DESCRIPTION */

/* CTAN(Z) calculates the complex trigonometric tangent of complex */
/* argument Z.  Z is in units of radians. */

/* ***REFERENCES  (NONE) */
/* ***ROUTINES CALLED  R1MACH, XERCLR, XERMSG */
/* ***REVISION HISTORY  (YYMMDD) */
/*   770401  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890531  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ) */
/* ***END PROLOGUE  CTAN */
/* ***FIRST EXECUTABLE STATEMENT  CTAN */
    if (sqeps == 0.f) {
	sqeps = sqrt(r1mach_(&c__4));
    }

    x2 = z__->r * 2.f;
    y2 = r_imag(z__) * 2.f;

    sn2x = sin(x2);
    xerclr_();

    den = cos(x2) + cosh(y2);
    if (den == 0.f) {
	xermsg_("SLATEC", "CTAN", "TAN IS SINGULAR FOR INPUT Z (X IS PI/2 OR"
		" 3*PI/2 AND Y IS 0)", &c__2, &c__2, (ftnlen)6, (ftnlen)4, (
		ftnlen)60);
    }

/* Computing MAX */
    r__1 = dabs(x2);
    if (dabs(den) > dmax(r__1,1.f) * sqeps) {
	goto L10;
    }
    xerclr_();
    xermsg_("SLATEC", "CTAN", "ANSWER LT HALF PRECISION, ABS(X) TOO BIG OR X"
	    " TOO NEAR PI/2 OR 3*PI/2", &c__1, &c__1, (ftnlen)6, (ftnlen)4, (
	    ftnlen)69);

L10:
    r__1 = sn2x / den;
    r__2 = sinh(y2) / den;
    q__1.r = r__1, q__1.i = r__2;
     ret_val->r = q__1.r,  ret_val->i = q__1.i;

    return ;
} /* ctan_ */

