/* dbesj1.f -- translated by f2c (version 12.02.01).
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
static integer c__19 = 19;
static integer c__1 = 1;

/* DECK DBESJ1 */
doublereal dbesj1_(doublereal *x)
{
    /* Initialized data */

    static doublereal bj1cs[19] = { -.117261415133327865606240574524003,
	    -.253615218307906395623030884554698,
	    .0501270809844695685053656363203743,
	    -.00463151480962508191842619728789772,
	    2.47996229415914024539124064592364e-4,
	    -8.67894868627882584521246435176416e-6,
	    2.14293917143793691502766250991292e-7,
	    -3.93609307918317979229322764073061e-9,
	    5.59118231794688004018248059864032e-11,
	    -6.3276164046613930247769527401488e-13,
	    5.84099161085724700326945563268266e-15,
	    -4.48253381870125819039135059199999e-17,
	    2.90538449262502466306018688e-19,
	    -1.61173219784144165412118186666666e-21,
	    7.73947881939274637298346666666666e-24,
	    -3.24869378211199841143466666666666e-26,1.2022376772274102272e-28,
	    -3.95201221265134933333333333333333e-31,
	    1.16167808226645333333333333333333e-33 };
    static logical first = TRUE_;

    /* System generated locals */
    real r__1;
    doublereal ret_val, d__1;

    /* Local variables */
    static doublereal y;
    static integer ntj1;
    static doublereal ampl, xmin, xsml;
    extern /* Subroutine */ int d9b1mp_(doublereal *, doublereal *, 
	    doublereal *);
    static doublereal theta;
    extern doublereal d1mach_(integer *), dcsevl_(doublereal *, doublereal *, 
	    integer *);
    extern integer initds_(doublereal *, integer *, real *);
    extern /* Subroutine */ int xermsg_(char *, char *, char *, integer *, 
	    integer *, ftnlen, ftnlen, ftnlen);

/* ***BEGIN PROLOGUE  DBESJ1 */
/* ***PURPOSE  Compute the Bessel function of the first kind of order one. */
/* ***LIBRARY   SLATEC (FNLIB) */
/* ***CATEGORY  C10A1 */
/* ***TYPE      DOUBLE PRECISION (BESJ1-S, DBESJ1-D) */
/* ***KEYWORDS  BESSEL FUNCTION, FIRST KIND, FNLIB, ORDER ONE, */
/*             SPECIAL FUNCTIONS */
/* ***AUTHOR  Fullerton, W., (LANL) */
/* ***DESCRIPTION */

/* DBESJ1(X) calculates the double precision Bessel function of the */
/* first kind of order one for double precision argument X. */

/* Series for BJ1        on the interval  0.          to  1.60000E+01 */
/*                                        with weighted error   1.16E-33 */
/*                                         log weighted error  32.93 */
/*                               significant figures required  32.36 */
/*                                    decimal places required  33.57 */

/* ***REFERENCES  (NONE) */
/* ***ROUTINES CALLED  D1MACH, D9B1MP, DCSEVL, INITDS, XERMSG */
/* ***REVISION HISTORY  (YYMMDD) */
/*   780601  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890531  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ) */
/*   910401  Corrected error in code which caused values to have the */
/*           wrong sign for arguments less than 4.0.  (WRB) */
/* ***END PROLOGUE  DBESJ1 */
/* ***FIRST EXECUTABLE STATEMENT  DBESJ1 */
    if (first) {
	r__1 = (real) d1mach_(&c__3) * .1f;
	ntj1 = initds_(bj1cs, &c__19, &r__1);

	xsml = sqrt(d1mach_(&c__3) * 8.);
	xmin = d1mach_(&c__1) * 2.;
    }
    first = FALSE_;

    y = abs(*x);
    if (y > 4.) {
	goto L20;
    }

    ret_val = 0.;
    if (y == 0.) {
	return ret_val;
    }
    if (y <= xmin) {
	xermsg_("SLATEC", "DBESJ1", "ABS(X) SO SMALL J1 UNDERFLOWS", &c__1, &
		c__1, (ftnlen)6, (ftnlen)6, (ftnlen)29);
    }
    if (y > xmin) {
	ret_val = *x * .5;
    }
    if (y > xsml) {
	d__1 = y * .125 * y - 1.;
	ret_val = *x * (dcsevl_(&d__1, bj1cs, &ntj1) + .25);
    }
    return ret_val;

L20:
    d9b1mp_(&y, &ampl, &theta);
    ret_val = d_sign(&ampl, x) * cos(theta);

    return ret_val;
} /* dbesj1_ */

