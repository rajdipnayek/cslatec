/* e1.f -- translated by f2c (version 12.02.01).
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
static integer c__39 = 39;
static integer c__25 = 25;
static integer c__19 = 19;
static integer c__16 = 16;
static integer c__26 = 26;
static integer c__1 = 1;
static integer c__2 = 2;

/* DECK E1 */
doublereal e1_(real *x)
{
    /* Initialized data */

    static real ae11cs[39] = { .12150323971606579f,-.06508877851355015f,
	    .00489765135745967f,-6.49237843027216e-4f,9.3840434587471e-5f,
	    4.20236380882e-7f,-8.113374735904e-6f,2.804247688663e-6f,
	    5.6487164441e-8f,-3.4480917445e-7f,5.8209273578e-8f,
	    3.8711426349e-8f,-1.2453235014e-8f,-5.118504888e-9f,
	    2.148771527e-9f,8.68459898e-10f,-3.43650105e-10f,-1.79796603e-10f,
	    4.744206e-11f,4.0423282e-11f,-3.543928e-12f,-8.853444e-12f,
	    -9.60151e-13f,1.692921e-12f,6.0799e-13f,-2.24338e-13f,
	    -2.00327e-13f,-6.246e-15f,4.5571e-14f,1.6383e-14f,-5.561e-15f,
	    -6.074e-15f,-8.62e-16f,1.223e-15f,7.16e-16f,-2.4e-17f,-2.01e-16f,
	    -8.2e-17f,1.7e-17f };
    static real ae12cs[25] = { .58241749513472674f,-.15834885090578275f,
	    -.006764275590323141f,.005125843950185725f,4.35232492169391e-4f,
	    -1.43613366305483e-4f,-4.1801320556301e-5f,-2.71339575864e-6f,
	    1.151381913647e-6f,4.20650022012e-7f,6.6581901391e-8f,
	    6.62143777e-10f,-2.84410487e-9f,-9.40724197e-10f,-1.77476602e-10f,
	    -1.5830222e-11f,2.905732e-12f,1.769356e-12f,4.92735e-13f,
	    9.3709e-14f,1.0707e-14f,-5.37e-16f,-7.16e-16f,-2.44e-16f,
	    -5.8e-17f };
    static real e11cs[19] = { -16.113461655571494026f,7.7940727787426802769f,
	    -1.9554058188631419507f,.37337293866277945612f,
	    -.05692503191092901938f,.00721107776966009185f,
	    -7.8104901449841593e-4f,7.388093356262168e-5f,
	    -6.2028618758082e-6f,4.6816002303176e-7f,-3.209288853329e-8f,
	    2.01519974874e-9f,-1.1673686816e-10f,6.27627066e-12f,
	    -3.1481541e-13f,1.479904e-14f,-6.5457e-16f,2.733e-17f,-1.08e-18f }
	    ;
    static real e12cs[16] = { -.037390214792202795f,.042723986062209577f,
	    -.1303182079849700544f,.01441912402469889073f,
	    -.00134617078051068022f,1.073102925306378e-4f,
	    -7.42999951611943e-6f,4.5377325690753e-7f,-2.47641721139e-8f,
	    1.22076581374e-9f,-5.48514148e-11f,2.26362142e-12f,-8.635897e-14f,
	    3.06291e-15f,-1.0148e-16f,3.15e-18f };
    static real ae13cs[25] = { -.60577324664060346f,-.1125352434836609f,
	    .013432266247902779f,-.001926845187381145f,3.09118337720603e-4f,
	    -5.3564132129618e-5f,9.827812880247e-6f,-1.885368984916e-6f,
	    3.74943193568e-7f,-7.682345587e-8f,1.6143270567e-8f,
	    -3.466802211e-9f,7.58754209e-10f,-1.68864333e-10f,3.8145706e-11f,
	    -8.733026e-12f,2.023672e-12f,-4.74132e-13f,1.12211e-13f,
	    -2.6804e-14f,6.457e-15f,-1.568e-15f,3.83e-16f,-9.4e-17f,2.3e-17f }
	    ;
    static real ae14cs[26] = { -.1892918000753017f,-.08648117855259871f,
	    .00722410154374659f,-8.0975594575573e-4f,1.0999134432661e-4f,
	    -1.717332998937e-5f,2.98562751447e-6f,-5.6596491457e-7f,
	    1.1526808397e-7f,-2.49503044e-8f,5.6923242e-9f,-1.35995766e-9f,
	    3.3846628e-10f,-8.737853e-11f,2.331588e-11f,-6.41148e-12f,
	    1.81224e-12f,-5.2538e-13f,1.5592e-13f,-4.729e-14f,1.463e-14f,
	    -4.61e-15f,1.48e-15f,-4.8e-16f,1.6e-16f,-5e-17f };
    static logical first = TRUE_;

    /* System generated locals */
    real ret_val, r__1;

    /* Local variables */
    static real eta;
    static integer nte11, nte12;
    static real xmax;
    static integer ntae11, ntae12, ntae13, ntae14;
    extern doublereal csevl_(real *, real *, integer *);
    extern integer inits_(real *, integer *, real *);
    static real xmaxt;
    extern doublereal r1mach_(integer *);
    extern /* Subroutine */ int xermsg_(char *, char *, char *, integer *, 
	    integer *, ftnlen, ftnlen, ftnlen);

/* ***BEGIN PROLOGUE  E1 */
/* ***PURPOSE  Compute the exponential integral E1(X). */
/* ***LIBRARY   SLATEC (FNLIB) */
/* ***CATEGORY  C5 */
/* ***TYPE      SINGLE PRECISION (E1-S, DE1-D) */
/* ***KEYWORDS  E1 FUNCTION, EXPONENTIAL INTEGRAL, FNLIB, */
/*             SPECIAL FUNCTIONS */
/* ***AUTHOR  Fullerton, W., (LANL) */
/* ***DESCRIPTION */

/* E1 calculates the single precision exponential integral, E1(X), for */
/* positive single precision argument X and the Cauchy principal value */
/* for negative X.  If principal values are used everywhere, then, for */
/* all X, */

/*    E1(X) = -Ei(-X) */
/* or */
/*    Ei(X) = -E1(-X). */


/* Series for AE11       on the interval -1.00000D-01 to  0. */
/*                                        with weighted error   1.76E-17 */
/*                                         log weighted error  16.75 */
/*                               significant figures required  15.70 */
/*                                    decimal places required  17.55 */


/* Series for AE12       on the interval -2.50000D-01 to -1.00000D-01 */
/*                                        with weighted error   5.83E-17 */
/*                                         log weighted error  16.23 */
/*                               significant figures required  15.76 */
/*                                    decimal places required  16.93 */


/* Series for E11        on the interval -4.00000D+00 to -1.00000D+00 */
/*                                        with weighted error   1.08E-18 */
/*                                         log weighted error  17.97 */
/*                               significant figures required  19.02 */
/*                                    decimal places required  18.61 */


/* Series for E12        on the interval -1.00000D+00 to  1.00000D+00 */
/*                                        with weighted error   3.15E-18 */
/*                                         log weighted error  17.50 */
/*                        approx significant figures required  15.8 */
/*                                    decimal places required  18.10 */


/* Series for AE13       on the interval  2.50000D-01 to  1.00000D+00 */
/*                                        with weighted error   2.34E-17 */
/*                                         log weighted error  16.63 */
/*                               significant figures required  16.14 */
/*                                    decimal places required  17.33 */


/* Series for AE14       on the interval  0.          to  2.50000D-01 */
/*                                        with weighted error   5.41E-17 */
/*                                         log weighted error  16.27 */
/*                               significant figures required  15.38 */
/*                                    decimal places required  16.97 */

/* ***REFERENCES  (NONE) */
/* ***ROUTINES CALLED  CSEVL, INITS, R1MACH, XERMSG */
/* ***REVISION HISTORY  (YYMMDD) */
/*   770701  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   891115  Modified prologue description.  (WRB) */
/*   891115  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ) */
/*   920618  Removed space from variable names.  (RWC, WRB) */
/* ***END PROLOGUE  E1 */
/* ***FIRST EXECUTABLE STATEMENT  E1 */
    if (first) {
	eta = r1mach_(&c__3) * .1f;
	ntae11 = inits_(ae11cs, &c__39, &eta);
	ntae12 = inits_(ae12cs, &c__25, &eta);
	nte11 = inits_(e11cs, &c__19, &eta);
	nte12 = inits_(e12cs, &c__16, &eta);
	ntae13 = inits_(ae13cs, &c__25, &eta);
	ntae14 = inits_(ae14cs, &c__26, &eta);

	xmaxt = -log(r1mach_(&c__1));
	xmax = xmaxt - log(xmaxt);
    }
    first = FALSE_;

    if (*x > -10.f) {
	goto L20;
    }

/* E1(X) = -EI(-X) FOR X .LE. -10. */

    r__1 = 20.f / *x + 1.f;
    ret_val = exp(-(*x)) / *x * (csevl_(&r__1, ae11cs, &ntae11) + 1.f);
    return ret_val;

L20:
    if (*x > -4.f) {
	goto L30;
    }

/* E1(X) = -EI(-X) FOR -10. .LT. X .LE. -4. */

    r__1 = (40.f / *x + 7.f) / 3.f;
    ret_val = exp(-(*x)) / *x * (csevl_(&r__1, ae12cs, &ntae12) + 1.f);
    return ret_val;

L30:
    if (*x > -1.f) {
	goto L40;
    }

/* E1(X) = -EI(-X) FOR -4. .LT. X .LE. -1. */

    r__1 = (*x * 2.f + 5.f) / 3.f;
    ret_val = -log((dabs(*x))) + csevl_(&r__1, e11cs, &nte11);
    return ret_val;

L40:
    if (*x > 1.f) {
	goto L50;
    }
    if (*x == 0.f) {
	xermsg_("SLATEC", "E1", "X IS 0", &c__2, &c__2, (ftnlen)6, (ftnlen)2, 
		(ftnlen)6);
    }

/* E1(X) = -EI(-X) FOR -1. .LT. X .LE. 1.,  X .NE. 0. */

    ret_val = -log((dabs(*x))) - .6875f + *x + csevl_(x, e12cs, &nte12);
    return ret_val;

L50:
    if (*x > 4.f) {
	goto L60;
    }

/* E1(X) = -EI(-X) FOR 1. .LT. X .LE. 4. */

    r__1 = (8.f / *x - 5.f) / 3.f;
    ret_val = exp(-(*x)) / *x * (csevl_(&r__1, ae13cs, &ntae13) + 1.f);
    return ret_val;

L60:
    if (*x > xmax) {
	goto L70;
    }

/* E1(X) = -EI(-X) FOR 4. .LT. X .LE. XMAX */

    r__1 = 8.f / *x - 1.f;
    ret_val = exp(-(*x)) / *x * (csevl_(&r__1, ae14cs, &ntae14) + 1.f);
    return ret_val;

/* E1(X) = -EI(-X) FOR X .GT. XMAX */

L70:
    xermsg_("SLATEC", "E1", "X SO BIG E1 UNDERFLOWS", &c__1, &c__1, (ftnlen)6,
	     (ftnlen)2, (ftnlen)22);
    ret_val = 0.f;
    return ret_val;

} /* e1_ */

