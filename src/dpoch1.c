/* dpoch1.f -- translated by f2c (version 12.02.01).
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
static integer c__3 = 3;
static integer c__2 = 2;

/* DECK DPOCH1 */
doublereal dpoch1_(doublereal *a, doublereal *x)
{
    /* Initialized data */

    static doublereal bern[20] = { .0833333333333333333333333333333333,
	    -.00138888888888888888888888888888888,
	    3.3068783068783068783068783068783e-5,
	    -8.26719576719576719576719576719576e-7,
	    2.08767569878680989792100903212014e-8,
	    -5.28419013868749318484768220217955e-10,
	    1.33825365306846788328269809751291e-11,
	    -3.38968029632258286683019539124944e-13,
	    8.58606205627784456413590545042562e-15,
	    -2.17486869855806187304151642386591e-16,
	    5.50900282836022951520265260890225e-18,
	    -1.39544646858125233407076862640635e-19,
	    3.53470703962946747169322997780379e-21,
	    -8.95351742703754685040261131811274e-23,
	    2.26795245233768306031095073886816e-24,
	    -5.744724395202645238348479719434e-25,
	    1.45517247561486490186626486727132e-27,
	    -3.68599494066531017818178247990866e-29,
	    9.33673425709504467203255515278562e-31,
	    -2.36502241570062993455963519636983e-32 };
    static doublereal pi = 3.141592653589793238462643383279503;
    static logical first = TRUE_;

    /* System generated locals */
    integer i__1, i__2;
    doublereal ret_val, d__1;

    /* Local variables */
    static doublereal b;
    static integer i__, j, k;
    static doublereal q, bp;
    static integer ii;
    static doublereal gbk, rho, var;
    static integer ndx;
    static doublereal var2, absa;
    extern doublereal dcot_(doublereal *);
    static integer incr;
    static doublereal absx, binv;
    extern doublereal dpsi_(doublereal *);
    static doublereal trig, term, poly1, gbern[21];
    extern doublereal dpoch_(doublereal *, doublereal *), d1mach_(integer *);
    static doublereal sinpx2, alneps, alnvar, sqtbig;
    extern doublereal dexprl_(doublereal *);
    extern /* Subroutine */ int xermsg_(char *, char *, char *, integer *, 
	    integer *, ftnlen, ftnlen, ftnlen);
    static integer nterms;
    static doublereal sinpxx;

/* ***BEGIN PROLOGUE  DPOCH1 */
/* ***PURPOSE  Calculate a generalization of Pochhammer's symbol starting */
/*            from first order. */
/* ***LIBRARY   SLATEC (FNLIB) */
/* ***CATEGORY  C1, C7A */
/* ***TYPE      DOUBLE PRECISION (POCH1-S, DPOCH1-D) */
/* ***KEYWORDS  FIRST ORDER, FNLIB, POCHHAMMER, SPECIAL FUNCTIONS */
/* ***AUTHOR  Fullerton, W., (LANL) */
/* ***DESCRIPTION */

/* Evaluate a double precision generalization of Pochhammer's symbol */
/* for double precision A and X for special situations that require */
/* especially accurate values when X is small in */
/*        POCH1(A,X) = (POCH(A,X)-1)/X */
/*                   = (GAMMA(A+X)/GAMMA(A) - 1.0)/X . */
/* This specification is particularly suited for stably computing */
/* expressions such as */
/*        (GAMMA(A+X)/GAMMA(A) - GAMMA(B+X)/GAMMA(B))/X */
/*             = POCH1(A,X) - POCH1(B,X) */
/* Note that POCH1(A,0.0) = PSI(A) */

/* When ABS(X) is so small that substantial cancellation will occur if */
/* the straightforward formula is used, we use an expansion due */
/* to Fields and discussed by Y. L. Luke, The Special Functions and Their */
/* Approximations, Vol. 1, Academic Press, 1969, page 34. */

/* The ratio POCH(A,X) = GAMMA(A+X)/GAMMA(A) is written by Luke as */
/*        (A+(X-1)/2)**X * polynomial in (A+(X-1)/2)**(-2) . */
/* In order to maintain significance in POCH1, we write for positive a */
/*        (A+(X-1)/2)**X = EXP(X*LOG(A+(X-1)/2)) = EXP(Q) */
/*                       = 1.0 + Q*EXPREL(Q) . */
/* Likewise the polynomial is written */
/*        POLY = 1.0 + X*POLY1(A,X) . */
/* Thus, */
/*        POCH1(A,X) = (POCH(A,X) - 1) / X */
/*                   = EXPREL(Q)*(Q/X + Q*POLY1(A,X)) + POLY1(A,X) */

/* ***REFERENCES  (NONE) */
/* ***ROUTINES CALLED  D1MACH, DCOT, DEXPRL, DPOCH, DPSI, XERMSG */
/* ***REVISION HISTORY  (YYMMDD) */
/*   770801  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890911  Removed unnecessary intrinsics.  (WRB) */
/*   890911  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ) */
/*   900727  Added EXTERNAL statement.  (WRB) */
/* ***END PROLOGUE  DPOCH1 */
/* ***FIRST EXECUTABLE STATEMENT  DPOCH1 */
    if (first) {
	sqtbig = 1. / sqrt(d1mach_(&c__1) * 24.);
	alneps = log(d1mach_(&c__3));
    }
    first = FALSE_;

    if (*x == 0.) {
	ret_val = dpsi_(a);
    }
    if (*x == 0.) {
	return ret_val;
    }

    absx = abs(*x);
    absa = abs(*a);
    if (absx > absa * .1) {
	goto L70;
    }
    if (absx * log((max(absa,2.))) > .1) {
	goto L70;
    }

    bp = *a;
    if (*a < -.5) {
	bp = 1. - *a - *x;
    }
    incr = 0;
    if (bp < 10.) {
	incr = (integer) (11. - bp);
    }
    b = bp + incr;

    var = b + (*x - 1.) * .5;
    alnvar = log(var);
    q = *x * alnvar;

    poly1 = 0.;
    if (var >= sqtbig) {
	goto L40;
    }
/* Computing 2nd power */
    d__1 = 1. / var;
    var2 = d__1 * d__1;

    rho = (*x + 1.) * .5;
    gbern[0] = 1.;
    gbern[1] = -rho / 12.;
    term = var2;
    poly1 = gbern[1] * term;

    nterms = (integer) (alneps * -.5 / alnvar + 1.);
    if (nterms > 20) {
	xermsg_("SLATEC", "DPOCH1", "NTERMS IS TOO BIG, MAYBE D1MACH(3) IS B"
		"AD", &c__1, &c__2, (ftnlen)6, (ftnlen)6, (ftnlen)41);
    }
    if (nterms < 2) {
	goto L40;
    }

    i__1 = nterms;
    for (k = 2; k <= i__1; ++k) {
	gbk = 0.;
	i__2 = k;
	for (j = 1; j <= i__2; ++j) {
	    ndx = k - j + 1;
	    gbk += bern[ndx - 1] * gbern[j - 1];
/* L20: */
	}
	gbern[k] = -rho * gbk / k;

	term = term * ((k << 1) - 2 - *x) * ((k << 1) - 1 - *x) * var2;
	poly1 += gbern[k] * term;
/* L30: */
    }

L40:
    poly1 = (*x - 1.) * poly1;
    ret_val = dexprl_(&q) * (alnvar + q * poly1) + poly1;

    if (incr == 0) {
	goto L60;
    }

/* WE HAVE DPOCH1(B,X), BUT BP IS SMALL, SO WE USE BACKWARDS RECURSION */
/* TO OBTAIN DPOCH1(BP,X). */

    i__1 = incr;
    for (ii = 1; ii <= i__1; ++ii) {
	i__ = incr - ii;
	binv = 1. / (bp + i__);
	ret_val = (ret_val - binv) / (*x * binv + 1.);
/* L50: */
    }

L60:
    if (bp == *a) {
	return ret_val;
    }

/* WE HAVE DPOCH1(BP,X), BUT A IS LT -0.5.  WE THEREFORE USE A REFLECTION */
/* FORMULA TO OBTAIN DPOCH1(A,X). */

    sinpxx = sin(pi * *x) / *x;
    sinpx2 = sin(pi * .5 * *x);
    d__1 = pi * b;
    trig = sinpxx * dcot_(&d__1) - sinpx2 * 2. * (sinpx2 / *x);

    ret_val = trig + (*x * trig + 1.) * ret_val;
    return ret_val;

L70:
    ret_val = (dpoch_(a, x) - 1.) / *x;
    return ret_val;

} /* dpoch1_ */

