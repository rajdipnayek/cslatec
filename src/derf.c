/* derf.f -- translated by f2c (version 12.02.01).
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
static integer c__21 = 21;
static doublereal c_b7 = 1.;

/* DECK DERF */
doublereal derf_(doublereal *x)
{
    /* Initialized data */

    static doublereal erfcs[21] = { -.049046121234691808039984544033376,
	    -.14226120510371364237824741899631,
	    .010035582187599795575754676712933,
	    -5.7687646997674847650827025509167e-4,
	    2.7419931252196061034422160791471e-5,
	    -1.1043175507344507604135381295905e-6,
	    3.8488755420345036949961311498174e-8,
	    -1.1808582533875466969631751801581e-9,
	    3.2334215826050909646402930953354e-11,
	    -7.9910159470045487581607374708595e-13,
	    1.7990725113961455611967245486634e-14,
	    -3.7186354878186926382316828209493e-16,
	    7.1035990037142529711689908394666e-18,
	    -1.2612455119155225832495424853333e-19,
	    2.0916406941769294369170500266666e-21,
	    -3.253973102931407298236416e-23,
	    4.7668672097976748332373333333333e-25,
	    -6.5980120782851343155199999999999e-27,
	    8.6550114699637626197333333333333e-29,
	    -1.0788925177498064213333333333333e-30,
	    1.2811883993017002666666666666666e-32 };
    static doublereal sqrtpi = 1.77245385090551602729816748334115;
    static logical first = TRUE_;

    /* System generated locals */
    real r__1;
    doublereal ret_val, d__1;

    /* Local variables */
    static doublereal y, xbig;
    extern doublereal derfc_(doublereal *);
    static integer nterf;
    static doublereal sqeps;
    extern doublereal d1mach_(integer *), dcsevl_(doublereal *, doublereal *, 
	    integer *);
    extern integer initds_(doublereal *, integer *, real *);

/* ***BEGIN PROLOGUE  DERF */
/* ***PURPOSE  Compute the error function. */
/* ***LIBRARY   SLATEC (FNLIB) */
/* ***CATEGORY  C8A, L5A1E */
/* ***TYPE      DOUBLE PRECISION (ERF-S, DERF-D) */
/* ***KEYWORDS  ERF, ERROR FUNCTION, FNLIB, SPECIAL FUNCTIONS */
/* ***AUTHOR  Fullerton, W., (LANL) */
/* ***DESCRIPTION */

/* DERF(X) calculates the double precision error function for double */
/* precision argument X. */

/* Series for ERF        on the interval  0.          to  1.00000E+00 */
/*                                        with weighted error   1.28E-32 */
/*                                         log weighted error  31.89 */
/*                               significant figures required  31.05 */
/*                                    decimal places required  32.55 */

/* ***REFERENCES  (NONE) */
/* ***ROUTINES CALLED  D1MACH, DCSEVL, DERFC, INITDS */
/* ***REVISION HISTORY  (YYMMDD) */
/*   770701  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890531  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900727  Added EXTERNAL statement.  (WRB) */
/*   920618  Removed space from variable name.  (RWC, WRB) */
/* ***END PROLOGUE  DERF */
/* ***FIRST EXECUTABLE STATEMENT  DERF */
    if (first) {
	r__1 = (real) d1mach_(&c__3) * .1f;
	nterf = initds_(erfcs, &c__21, &r__1);
	xbig = sqrt(-log(sqrtpi * d1mach_(&c__3)));
	sqeps = sqrt(d1mach_(&c__3) * 2.);
    }
    first = FALSE_;

    y = abs(*x);
    if (y > 1.) {
	goto L20;
    }

/* ERF(X) = 1.0 - ERFC(X)  FOR  -1.0 .LE. X .LE. 1.0 */

    if (y <= sqeps) {
	ret_val = *x * 2. * *x / sqrtpi;
    }
    if (y > sqeps) {
	d__1 = *x * 2. * *x - 1.;
	ret_val = *x * (dcsevl_(&d__1, erfcs, &nterf) + 1.);
    }
    return ret_val;

/* ERF(X) = 1.0 - ERFC(X) FOR ABS(X) .GT. 1.0 */

L20:
    if (y <= xbig) {
	d__1 = 1. - derfc_(&y);
	ret_val = d_sign(&d__1, x);
    }
    if (y > xbig) {
	ret_val = d_sign(&c_b7, x);
    }

    return ret_val;
} /* derf_ */

