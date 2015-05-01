/* clnrel.f -- translated by f2c (version 12.02.01).
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
static integer c__1 = 1;

/* DECK CLNREL */
/* Complex */ void clnrel_(complex * ret_val, complex *z__)
{
    /* Initialized data */

    static real sqeps = 0.f;

    /* System generated locals */
    real r__1, r__2, r__3, r__4;
    complex q__1, q__2;

    /* Local variables */
    static real x, rho;
    extern doublereal carg_(complex *), r1mach_(integer *), alnrel_(real *);
    extern /* Subroutine */ int xermsg_(char *, char *, char *, integer *, 
	    integer *, ftnlen, ftnlen, ftnlen);

/* ***BEGIN PROLOGUE  CLNREL */
/* ***PURPOSE  Evaluate ln(1+X) accurate in the sense of relative error. */
/* ***LIBRARY   SLATEC (FNLIB) */
/* ***CATEGORY  C4B */
/* ***TYPE      COMPLEX (ALNREL-S, DLNREL-D, CLNREL-C) */
/* ***KEYWORDS  ELEMENTARY FUNCTIONS, FNLIB, LOGARITHM */
/* ***AUTHOR  Fullerton, W., (LANL) */
/* ***DESCRIPTION */

/* CLNREL(Z) = LOG(1+Z) with relative error accuracy near Z = 0. */
/* Let   RHO = ABS(Z)  and */
/*       R**2 = ABS(1+Z)**2 = (1+X)**2 + Y**2 = 1 + 2*X + RHO**2 . */
/* Now if RHO is small we may evaluate CLNREL(Z) accurately by */
/*       LOG(1+Z) = CMPLX  (LOG(R), CARG(1+Z)) */
/*                 = CMPLX  (0.5*LOG(R**2), CARG(1+Z)) */
/*                 = CMPLX  (0.5*ALNREL(2*X+RHO**2), CARG(1+Z)) */

/* ***REFERENCES  (NONE) */
/* ***ROUTINES CALLED  ALNREL, CARG, R1MACH, XERMSG */
/* ***REVISION HISTORY  (YYMMDD) */
/*   770401  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890531  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ) */
/* ***END PROLOGUE  CLNREL */
/* ***FIRST EXECUTABLE STATEMENT  CLNREL */
    if (sqeps == 0.f) {
	sqeps = sqrt(r1mach_(&c__4));
    }

    q__1.r = z__->r + 1.f, q__1.i = z__->i;
    if (c_abs(&q__1) < sqeps) {
	xermsg_("SLATEC", "CLNREL", "ANSWER LT HALF PRECISION BECAUSE Z TOO "
		"NEAR -1", &c__1, &c__1, (ftnlen)6, (ftnlen)6, (ftnlen)46);
    }

    rho = c_abs(z__);
    if (rho > .375f) {
	q__2.r = z__->r + 1.f, q__2.i = z__->i;
	c_log(&q__1, &q__2);
	 ret_val->r = q__1.r,  ret_val->i = q__1.i;
    }
    if (rho > .375f) {
	return ;
    }

    x = z__->r;
/* Computing 2nd power */
    r__3 = rho;
    r__2 = x * 2.f + r__3 * r__3;
    r__1 = alnrel_(&r__2) * .5f;
    q__2.r = z__->r + 1.f, q__2.i = z__->i;
    r__4 = carg_(&q__2);
    q__1.r = r__1, q__1.i = r__4;
     ret_val->r = q__1.r,  ret_val->i = q__1.i;

    return ;
} /* clnrel_ */

