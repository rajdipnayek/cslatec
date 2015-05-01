/* besj1.f -- translated by f2c (version 12.02.01).
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
static integer c__12 = 12;
static integer c__21 = 21;
static integer c__24 = 24;
static integer c__1 = 1;
static integer c__4 = 4;
static integer c__2 = 2;

/* DECK BESJ1 */
doublereal besj1_(real *x)
{
    /* Initialized data */

    static real bj1cs[12] = { -.11726141513332787f,-.2536152183079064f,
	    .050127080984469569f,-.004631514809625081f,2.47996229415914e-4f,
	    -8.678948686278e-6f,2.14293917143e-7f,-3.936093079e-9f,
	    5.5911823e-11f,-6.32761e-13f,5.84e-15f,-4.4e-17f };
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
    static real pi4 = .78539816339744831f;
    static logical first = TRUE_;

    /* System generated locals */
    real ret_val, r__1;

    /* Local variables */
    static real y, z__;
    static integer ntj1, ntm1;
    static real ampl, xmin, xmax, xsml;
    static integer ntth1;
    static real theta;
    extern doublereal csevl_(real *, real *, integer *);
    extern integer inits_(real *, integer *, real *);
    extern doublereal r1mach_(integer *);
    extern /* Subroutine */ int xermsg_(char *, char *, char *, integer *, 
	    integer *, ftnlen, ftnlen, ftnlen);

/* ***BEGIN PROLOGUE  BESJ1 */
/* ***PURPOSE  Compute the Bessel function of the first kind of order one. */
/* ***LIBRARY   SLATEC (FNLIB) */
/* ***CATEGORY  C10A1 */
/* ***TYPE      SINGLE PRECISION (BESJ1-S, DBESJ1-D) */
/* ***KEYWORDS  BESSEL FUNCTION, FIRST KIND, FNLIB, ORDER ONE, */
/*             SPECIAL FUNCTIONS */
/* ***AUTHOR  Fullerton, W., (LANL) */
/* ***DESCRIPTION */

/* BESJ1(X) calculates the Bessel function of the first kind of */
/* order one for real argument X. */

/* Series for BJ1        on the interval  0.          to  1.60000D+01 */
/*                                        with weighted error   4.48E-17 */
/*                                         log weighted error  16.35 */
/*                               significant figures required  15.77 */
/*                                    decimal places required  16.89 */

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
/* ***ROUTINES CALLED  CSEVL, INITS, R1MACH, XERMSG */
/* ***REVISION HISTORY  (YYMMDD) */
/*   780601  DATE WRITTEN */
/*   890210  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ) */
/*   900326  Removed duplicate information from DESCRIPTION section. */
/*           (WRB) */
/* ***END PROLOGUE  BESJ1 */
/* ***FIRST EXECUTABLE STATEMENT  BESJ1 */
    if (first) {
	r__1 = r1mach_(&c__3) * .1f;
	ntj1 = inits_(bj1cs, &c__12, &r__1);
	r__1 = r1mach_(&c__3) * .1f;
	ntm1 = inits_(bm1cs, &c__21, &r__1);
	r__1 = r1mach_(&c__3) * .1f;
	ntth1 = inits_(bth1cs, &c__24, &r__1);

	xsml = sqrt(r1mach_(&c__3) * 8.f);
	xmin = r1mach_(&c__1) * 2.f;
	xmax = 1.f / r1mach_(&c__4);
    }
    first = FALSE_;

    y = dabs(*x);
    if (y > 4.f) {
	goto L20;
    }

    ret_val = 0.f;
    if (y == 0.f) {
	return ret_val;
    }
    if (y <= xmin) {
	xermsg_("SLATEC", "BESJ1", "ABS(X) SO SMALL J1 UNDERFLOWS", &c__1, &
		c__1, (ftnlen)6, (ftnlen)5, (ftnlen)29);
    }
    if (y > xmin) {
	ret_val = *x * .5f;
    }
    if (y > xsml) {
	r__1 = y * .125f * y - 1.f;
	ret_val = *x * (csevl_(&r__1, bj1cs, &ntj1) + .25f);
    }
    return ret_val;

L20:
    if (y > xmax) {
	xermsg_("SLATEC", "BESJ1", "NO PRECISION BECAUSE ABS(X) IS TOO BIG", &
		c__2, &c__2, (ftnlen)6, (ftnlen)5, (ftnlen)38);
    }
/* Computing 2nd power */
    r__1 = y;
    z__ = 32.f / (r__1 * r__1) - 1.f;
    ampl = (csevl_(&z__, bm1cs, &ntm1) + .75f) / sqrt(y);
    theta = y - pi4 * 3.f + csevl_(&z__, bth1cs, &ntth1) / y;
    ret_val = r_sign(&ampl, x) * cos(theta);

    return ret_val;
} /* besj1_ */

