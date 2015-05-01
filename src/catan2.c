/* catan2.f -- translated by f2c (version 12.02.01).
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

static integer c__1 = 1;
static integer c__2 = 2;

/* DECK CATAN2 */
/* Complex */ void catan2_(complex * ret_val, complex *csn, complex *ccs)
{
    /* Initialized data */

    static real pi = 3.14159265358979323846f;

    /* System generated locals */
    real r__1, r__2, r__3;
    complex q__1, q__2;

    /* Local variables */
    extern /* Complex */ void catan_(complex *, complex *);
    extern /* Subroutine */ int xermsg_(char *, char *, char *, integer *, 
	    integer *, ftnlen, ftnlen, ftnlen);

/* ***BEGIN PROLOGUE  CATAN2 */
/* ***PURPOSE  Compute the complex arc tangent in the proper quadrant. */
/* ***LIBRARY   SLATEC (FNLIB) */
/* ***CATEGORY  C4A */
/* ***TYPE      COMPLEX (CATAN2-C) */
/* ***KEYWORDS  ARC TANGENT, ELEMENTARY FUNCTIONS, FNLIB, POLAR ANGEL, */
/*             QUADRANT, TRIGONOMETRIC */
/* ***AUTHOR  Fullerton, W., (LANL) */
/* ***DESCRIPTION */

/* CATAN2(CSN,CCS) calculates the complex trigonometric arc */
/* tangent of the ratio CSN/CCS and returns a result whose real */
/* part is in the correct quadrant (within a multiple of 2*PI).  The */
/* result is in units of radians and the real part is between -PI */
/* and +PI. */

/* ***REFERENCES  (NONE) */
/* ***ROUTINES CALLED  CATAN, XERMSG */
/* ***REVISION HISTORY  (YYMMDD) */
/*   770401  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890531  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ) */
/*   900326  Removed duplicate information from DESCRIPTION section. */
/*           (WRB) */
/* ***END PROLOGUE  CATAN2 */
/* ***FIRST EXECUTABLE STATEMENT  CATAN2 */
    if (c_abs(ccs) == 0.f) {
	goto L10;
    }

    c_div(&q__2, csn, ccs);
    catan_(&q__1, &q__2);
     ret_val->r = q__1.r,  ret_val->i = q__1.i;
    if (ccs->r < 0.f) {
	q__1.r =  ret_val->r + pi, q__1.i =  ret_val->i;
	 ret_val->r = q__1.r,  ret_val->i = q__1.i;
    }
    if ( ret_val->r > pi) {
	r__1 = pi * 2.f;
	q__1.r =  ret_val->r - r__1, q__1.i =  ret_val->i;
	 ret_val->r = q__1.r,  ret_val->i = q__1.i;
    }
    return ;

L10:
    if (c_abs(csn) == 0.f) {
	xermsg_("SLATEC", "CATAN2", "CALLED WITH BOTH ARGUMENTS ZERO", &c__1, 
		&c__2, (ftnlen)6, (ftnlen)6, (ftnlen)31);
    }

    r__2 = pi * .5f;
    r__3 = csn->r;
    r__1 = r_sign(&r__2, &r__3);
    q__1.r = r__1, q__1.i = 0.f;
     ret_val->r = q__1.r,  ret_val->i = q__1.i;

    return ;
} /* catan2_ */

