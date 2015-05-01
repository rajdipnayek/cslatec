/* dspenc.f -- translated by f2c (version 12.02.01).
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
static integer c__38 = 38;

/* DECK DSPENC */
doublereal dspenc_(doublereal *x)
{
    /* Initialized data */

    static doublereal spencs[38] = { .1527365598892405872946684910028,
	    .08169658058051014403501838185271,
	    .005814157140778730872977350641182,
	    5.371619814541527542247889005319e-4,
	    5.724704675185826233210603054782e-5,
	    6.674546121649336343607835438589e-6,
	    8.276467339715676981584391689011e-7,
	    1.073315673030678951270005873354e-7,
	    1.440077294303239402334590331513e-8,
	    1.984442029965906367898877139608e-9,
	    2.794005822163638720201994821615e-10,
	    4.003991310883311823072580445908e-11,
	    5.823462892044638471368135835757e-12,
	    8.576708692638689278097914771224e-13,
	    1.276862586280193045989483033433e-13,
	    1.918826209042517081162380416062e-14,
	    2.907319206977138177795799719673e-15,
	    4.437112685276780462557473641745e-16,
	    6.815727787414599527867359135607e-17,
	    1.053017386015574429547019416644e-17,
	    1.63538980675237710005182173457e-18,
	    2.551852874940463932310901642581e-19,
	    3.999020621999360112770470379519e-20,
	    6.291501645216811876514149171199e-21,
	    9.933827435675677643803887752533e-22,
	    1.573679570749964816721763805866e-22,
	    2.500595316849476129369270954666e-23,
	    3.984740918383811139210663253333e-24,
	    6.366473210082843892691326293333e-25,
	    1.019674287239678367077061973333e-25,
	    1.636881058913518841111074133333e-26,
	    2.633310439417650117345279999999e-27,
	    4.244811560123976817224362666666e-28,
	    6.855411983680052916824746666666e-29,
	    1.109122433438056434018986666666e-29,
	    1.797431304999891457365333333333e-30,
	    2.917505845976095173290666666666e-31,
	    4.742646808928671061333333333333e-32 };
    static doublereal pi26 = 1.644934066848226436472415166646025189219;
    static logical first = TRUE_;

    /* System generated locals */
    real r__1;
    doublereal ret_val, d__1, d__2;

    /* Local variables */
    static doublereal aln, xbig;
    extern doublereal d1mach_(integer *), dcsevl_(doublereal *, doublereal *, 
	    integer *);
    static integer nspenc;
    extern integer initds_(doublereal *, integer *, real *);

/* ***BEGIN PROLOGUE  DSPENC */
/* ***PURPOSE  Compute a form of Spence's integral due to K. Mitchell. */
/* ***LIBRARY   SLATEC (FNLIB) */
/* ***CATEGORY  C5 */
/* ***TYPE      DOUBLE PRECISION (SPENC-S, DSPENC-D) */
/* ***KEYWORDS  FNLIB, SPECIAL FUNCTIONS, SPENCE'S INTEGRAL */
/* ***AUTHOR  Fullerton, W., (LANL) */
/* ***DESCRIPTION */

/* DSPENC(X) calculates the double precision Spence's integral */
/* for double precision argument X.  Spence's function defined by */
/*        integral from 0 to X of  -LOG(1-Y)/Y  DY. */
/* For ABS(X) .LE. 1, the uniformly convergent expansion */
/*        DSPENC = sum K=1,infinity  X**K / K**2     is valid. */
/* This is a form of Spence's integral due to K. Mitchell which differs */
/* from the definition in the NBS Handbook of Mathematical Functions. */

/* Spence's function can be used to evaluate much more general integral */
/* forms.  For example, */
/*        integral from 0 to Z of  LOG(A*X+B)/(C*X+D)  DX  = */
/*             LOG(ABS(B-A*D/C))*LOG(ABS(A*(C*X+D)/(A*D-B*C)))/C */
/*             - DSPENC (A*(C*Z+D)/(A*D-B*C)) / C. */

/* Ref -- K. Mitchell, Philosophical Magazine, 40, p.351 (1949). */
/*        Stegun and Abromowitz, AMS 55, p.1004. */


/* Series for SPEN       on the interval  0.          to  5.00000E-01 */
/*                                        with weighted error   4.74E-32 */
/*                                         log weighted error  31.32 */
/*                               significant figures required  30.37 */
/*                                    decimal places required  32.11 */

/* ***REFERENCES  (NONE) */
/* ***ROUTINES CALLED  D1MACH, DCSEVL, INITDS */
/* ***REVISION HISTORY  (YYMMDD) */
/*   780201  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   891115  Corrected third argument in reference to INITDS.  (WRB) */
/*   891115  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/* ***END PROLOGUE  DSPENC */
/* ***FIRST EXECUTABLE STATEMENT  DSPENC */
    if (first) {
	r__1 = (real) d1mach_(&c__3) * .1f;
	nspenc = initds_(spencs, &c__38, &r__1);
	xbig = 1. / d1mach_(&c__3);
    }
    first = FALSE_;

    if (*x > 2.) {
	goto L60;
    }
    if (*x > 1.) {
	goto L50;
    }
    if (*x > .5) {
	goto L40;
    }
    if (*x >= 0.) {
	goto L30;
    }
    if (*x > -1.) {
	goto L20;
    }

/* HERE IF X .LE. -1.0 */

    aln = log(1. - *x);
    ret_val = -pi26 - aln * .5 * (log(-(*x)) * 2. - aln);
    if (*x > -xbig) {
	d__1 = 4. / (1. - *x) - 1.;
	ret_val += (dcsevl_(&d__1, spencs, &nspenc) + 1.) / (1. - *x);
    }
    return ret_val;

/* -1.0 .LT. X .LT. 0.0 */

L20:
/* Computing 2nd power */
    d__1 = log(1. - *x);
    d__2 = *x * 4. / (*x - 1.) - 1.;
    ret_val = d__1 * d__1 * -.5 - *x * (dcsevl_(&d__2, spencs, &nspenc) + 1.) 
	    / (*x - 1.);
    return ret_val;

/* 0.0 .LE. X .LE. 0.5 */

L30:
    d__1 = *x * 4. - 1.;
    ret_val = *x * (dcsevl_(&d__1, spencs, &nspenc) + 1.);
    return ret_val;

/* 0.5 .LT. X .LE. 1.0 */

L40:
    ret_val = pi26;
    if (*x != 1.) {
	d__1 = (1. - *x) * 4. - 1.;
	ret_val = pi26 - log(*x) * log(1. - *x) - (1. - *x) * (dcsevl_(&d__1, 
		spencs, &nspenc) + 1.);
    }
    return ret_val;

/* 1.0 .LT. X .LE. 2.0 */

L50:
/* Computing 2nd power */
    d__1 = *x - 1.;
    d__2 = (*x - 1.) * 4. / *x - 1.;
    ret_val = pi26 - log(*x) * .5 * log(d__1 * d__1 / *x) + (*x - 1.) * (
	    dcsevl_(&d__2, spencs, &nspenc) + 1.) / *x;
    return ret_val;

/* X .GT. 2.0 */

L60:
/* Computing 2nd power */
    d__1 = log(*x);
    ret_val = pi26 * 2. - d__1 * d__1 * .5;
    if (*x < xbig) {
	d__1 = 4. / *x - 1.;
	ret_val -= (dcsevl_(&d__1, spencs, &nspenc) + 1.) / *x;
    }
    return ret_val;

} /* dspenc_ */

