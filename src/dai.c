/* dai.f -- translated by f2c (version 12.02.01).
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
static integer c__13 = 13;
static doublereal c_b7 = .3334;
static integer c__1 = 1;
static doublereal c_b9 = .6667;

/* DECK DAI */
doublereal dai_(doublereal *x)
{
    /* Initialized data */

    static doublereal aifcs[13] = { -.037971358496669997496197089469414,
	    .059191888537263638574319728013777,
	    9.862928057727997536560389104406e-4,
	    6.8488438190765667554854830182412e-6,
	    2.5942025962194713019489279081403e-8,
	    6.1766127740813750329445749697236e-11,
	    1.0092454172466117901429556224601e-13,
	    1.2014792511179938141288033225333e-16,
	    1.0882945588716991878525295466666e-19,
	    7.75137721966848870392384e-23,
	    4.4548112037175638391466666666666e-26,
	    2.1092845231692343466666666666666e-29,
	    8.3701735910741333333333333333333e-33 };
    static doublereal aigcs[13] = { .018152365581161273011556209957864,
	    .021572563166010755534030638819968,
	    2.5678356987483249659052428090133e-4,
	    1.4265214119792403898829496921721e-6,
	    4.5721149200180426070434097558191e-9,
	    9.5251708435647098607392278840592e-12,
	    1.392563460577139905115042068619e-14,
	    1.5070999142762379592306991138666e-17,
	    1.2559148312567778822703205333333e-20,
	    8.3063073770821340343829333333333e-24,
	    4.4657538493718567445333333333333e-27,
	    1.9900855034518869333333333333333e-30,
	    7.4702885256533333333333333333333e-34 };
    static logical first = TRUE_;

    /* System generated locals */
    real r__1;
    doublereal ret_val, d__1;

    /* Local variables */
    static doublereal z__, xm;
    extern doublereal daie_(doublereal *);
    static integer naif, naig;
    static doublereal xmax, x3sml, theta;
    extern doublereal d1mach_(integer *);
    static doublereal xmaxt;
    extern /* Subroutine */ int d9aimp_(doublereal *, doublereal *, 
	    doublereal *);
    extern doublereal dcsevl_(doublereal *, doublereal *, integer *);
    extern integer initds_(doublereal *, integer *, real *);
    extern /* Subroutine */ int xermsg_(char *, char *, char *, integer *, 
	    integer *, ftnlen, ftnlen, ftnlen);

/* ***BEGIN PROLOGUE  DAI */
/* ***PURPOSE  Evaluate the Airy function. */
/* ***LIBRARY   SLATEC (FNLIB) */
/* ***CATEGORY  C10D */
/* ***TYPE      DOUBLE PRECISION (AI-S, DAI-D) */
/* ***KEYWORDS  AIRY FUNCTION, FNLIB, SPECIAL FUNCTIONS */
/* ***AUTHOR  Fullerton, W., (LANL) */
/* ***DESCRIPTION */

/* DAI(X) calculates the double precision Airy function for double */
/* precision argument X. */

/* Series for AIF        on the interval -1.00000E+00 to  1.00000E+00 */
/*                                        with weighted error   8.37E-33 */
/*                                         log weighted error  32.08 */
/*                               significant figures required  30.87 */
/*                                    decimal places required  32.63 */

/* Series for AIG        on the interval -1.00000E+00 to  1.00000E+00 */
/*                                        with weighted error   7.47E-34 */
/*                                         log weighted error  33.13 */
/*                               significant figures required  31.50 */
/*                                    decimal places required  33.68 */

/* ***REFERENCES  (NONE) */
/* ***ROUTINES CALLED  D1MACH, D9AIMP, DAIE, DCSEVL, INITDS, XERMSG */
/* ***REVISION HISTORY  (YYMMDD) */
/*   770701  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890531  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ) */
/*   920618  Removed space from variable names.  (RWC, WRB) */
/* ***END PROLOGUE  DAI */
/* ***FIRST EXECUTABLE STATEMENT  DAI */
    if (first) {
	r__1 = (real) d1mach_(&c__3) * .1f;
	naif = initds_(aifcs, &c__13, &r__1);
	r__1 = (real) d1mach_(&c__3) * .1f;
	naig = initds_(aigcs, &c__13, &r__1);

	d__1 = d1mach_(&c__3);
	x3sml = pow_dd(&d__1, &c_b7);
	d__1 = log(d1mach_(&c__1)) * -1.5;
	xmaxt = pow_dd(&d__1, &c_b9);
	xmax = xmaxt - xmaxt * log(xmaxt) / (sqrt(xmaxt) * 4. + 1.) - .01;
    }
    first = FALSE_;

    if (*x >= -1.) {
	goto L20;
    }
    d9aimp_(x, &xm, &theta);
    ret_val = xm * cos(theta);
    return ret_val;

L20:
    if (*x > 1.) {
	goto L30;
    }
    z__ = 0.;
    if (abs(*x) > x3sml) {
/* Computing 3rd power */
	d__1 = *x;
	z__ = d__1 * (d__1 * d__1);
    }
    ret_val = dcsevl_(&z__, aifcs, &naif) - *x * (dcsevl_(&z__, aigcs, &naig) 
	    + .25) + .375;
    return ret_val;

L30:
    if (*x > xmax) {
	goto L40;
    }
    ret_val = daie_(x) * exp(*x * -2. * sqrt(*x) / 3.);
    return ret_val;

L40:
    ret_val = 0.;
    xermsg_("SLATEC", "DAI", "X SO BIG AI UNDERFLOWS", &c__1, &c__1, (ftnlen)
	    6, (ftnlen)3, (ftnlen)22);
    return ret_val;

} /* dai_ */

