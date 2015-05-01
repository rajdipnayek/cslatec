/* besy1.f -- translated by f2c (version 12.02.01).
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
static integer c__14 = 14;
static integer c__21 = 21;
static integer c__24 = 24;
static integer c__2 = 2;
static integer c__1 = 1;
static integer c__4 = 4;

/* DECK BESY1 */
doublereal besy1_(real *x)
{
    /* Initialized data */

    static real by1cs[14] = { .03208047100611908629f,1.26270789743350045f,
	    .006499961899923175f,-.08936164528860504117f,
	    .01325088122175709545f,-8.9790591196483523e-4f,
	    3.647361487958306e-5f,-1.001374381666e-6f,1.99453965739e-8f,
	    -3.0230656018e-10f,3.60987815e-12f,-3.487488e-14f,2.7838e-16f,
	    -1.86e-18f };
    static real bm1cs[21] = { .1047362510931285f,.00442443893702345f,
	    -5.661639504035e-5f,2.31349417339e-6f,-1.7377182007e-7f,
	    1.89320993e-8f,-2.65416023e-9f,4.4740209e-10f,-8.691795e-11f,
	    1.891492e-11f,-4.51884e-12f,1.16765e-12f,-3.2265e-13f,9.45e-14f,
	    -2.913e-14f,9.39e-15f,-3.15e-15f,1.09e-15f,-3.9e-16f,1.4e-16f,
	    -5e-17f };
    static real bth1cs[24] = { .7406014102631385f,-.00457175565963769f,
	    1.19818510964326e-4f,-6.964561891648e-6f,6.55495621447e-7f,
	    -8.4066228945e-8f,1.3376886564e-8f,-2.499565654e-9f,5.294951e-10f,
	    -1.24135944e-10f,3.1656485e-11f,-8.66864e-12f,2.523758e-12f,
	    -7.75085e-13f,2.49527e-13f,-8.3773e-14f,2.9205e-14f,-1.0534e-14f,
	    3.919e-15f,-1.5e-15f,5.89e-16f,-2.37e-16f,9.7e-17f,-4e-17f };
    static real twodpi = .63661977236758134f;
    static real pi4 = .78539816339744831f;
    static logical first = TRUE_;

    /* System generated locals */
    real ret_val, r__1, r__2;

    /* Local variables */
    static real y, z__;
    static integer ntm1, nty1;
    static real ampl, xmin, xmax, xsml;
    extern doublereal besj1_(real *);
    static integer ntth1;
    static real theta;
    extern doublereal csevl_(real *, real *, integer *);
    extern integer inits_(real *, integer *, real *);
    extern doublereal r1mach_(integer *);
    extern /* Subroutine */ int xermsg_(char *, char *, char *, integer *, 
	    integer *, ftnlen, ftnlen, ftnlen);

/* ***BEGIN PROLOGUE  BESY1 */
/* ***PURPOSE  Compute the Bessel function of the second kind of order */
/*            one. */
/* ***LIBRARY   SLATEC (FNLIB) */
/* ***CATEGORY  C10A1 */
/* ***TYPE      SINGLE PRECISION (BESY1-S, DBESY1-D) */
/* ***KEYWORDS  BESSEL FUNCTION, FNLIB, ORDER ONE, SECOND KIND, */
/*             SPECIAL FUNCTIONS */
/* ***AUTHOR  Fullerton, W., (LANL) */
/* ***DESCRIPTION */

/* BESY1(X) calculates the Bessel function of the second kind of */
/* order one for real argument X. */

/* Series for BY1        on the interval  0.          to  1.60000D+01 */
/*                                        with weighted error   1.87E-18 */
/*                                         log weighted error  17.73 */
/*                               significant figures required  17.83 */
/*                                    decimal places required  18.30 */

/* Series for BM1        on the interval  0.          to  6.25000D-02 */
/*                                        with weighted error   5.61E-17 */
/*                                         log weighted error  16.25 */
/*                               significant figures required  14.97 */
/*                                    decimal places required  16.91 */

/* Series for BTH1       on the interval  0.          to  6.25000D-02 */
/*                                        with weighted error   4.10E-17 */
/*                                         log weighted error  16.39 */
/*                               significant figures required  15.96 */
/*                                    decimal places required  17.08 */

/* ***REFERENCES  (NONE) */
/* ***ROUTINES CALLED  BESJ1, CSEVL, INITS, R1MACH, XERMSG */
/* ***REVISION HISTORY  (YYMMDD) */
/*   770401  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890531  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ) */
/*   900326  Removed duplicate information from DESCRIPTION section. */
/*           (WRB) */
/* ***END PROLOGUE  BESY1 */
/* ***FIRST EXECUTABLE STATEMENT  BESY1 */
    if (first) {
	r__1 = r1mach_(&c__3) * .1f;
	nty1 = inits_(by1cs, &c__14, &r__1);
	r__1 = r1mach_(&c__3) * .1f;
	ntm1 = inits_(bm1cs, &c__21, &r__1);
	r__1 = r1mach_(&c__3) * .1f;
	ntth1 = inits_(bth1cs, &c__24, &r__1);

/* Computing MAX */
	r__1 = log(r1mach_(&c__1)), r__2 = -log(r1mach_(&c__2));
	xmin = exp(dmax(r__1,r__2) + .01f) * 1.571f;
	xsml = sqrt(r1mach_(&c__3) * 4.f);
	xmax = 1.f / r1mach_(&c__4);
    }
    first = FALSE_;

    if (*x <= 0.f) {
	xermsg_("SLATEC", "BESY1", "X IS ZERO OR NEGATIVE", &c__1, &c__2, (
		ftnlen)6, (ftnlen)5, (ftnlen)21);
    }
    if (*x > 4.f) {
	goto L20;
    }

    if (*x < xmin) {
	xermsg_("SLATEC", "BESY1", "X SO SMALL Y1 OVERFLOWS", &c__3, &c__2, (
		ftnlen)6, (ftnlen)5, (ftnlen)23);
    }
    y = 0.f;
    if (*x > xsml) {
	y = *x * *x;
    }
    r__1 = y * .125f - 1.f;
    ret_val = twodpi * log(*x * .5f) * besj1_(x) + (csevl_(&r__1, by1cs, &
	    nty1) + .5f) / *x;
    return ret_val;

L20:
    if (*x > xmax) {
	xermsg_("SLATEC", "BESY1", "NO PRECISION BECAUSE X IS BIG", &c__2, &
		c__2, (ftnlen)6, (ftnlen)5, (ftnlen)29);
    }

/* Computing 2nd power */
    r__1 = *x;
    z__ = 32.f / (r__1 * r__1) - 1.f;
    ampl = (csevl_(&z__, bm1cs, &ntm1) + .75f) / sqrt(*x);
    theta = *x - pi4 * 3.f + csevl_(&z__, bth1cs, &ntth1) / *x;
    ret_val = ampl * sin(theta);

    return ret_val;
} /* besy1_ */

