/* poch1.f -- translated by f2c (version 12.02.01).
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

/* DECK POCH1 */
doublereal poch1_(real *a, real *x)
{
    /* Initialized data */

    static real bern[9] = { .083333333333333333f,-.0013888888888888889f,
	    3.3068783068783069e-5f,-8.2671957671957672e-7f,
	    2.0876756987868099e-8f,-5.2841901386874932e-10f,
	    1.3382536530684679e-11f,-3.3896802963225829e-13f,
	    8.5860620562778446e-15f };
    static real pi = 3.14159265358979324f;
    static logical first = TRUE_;

    /* System generated locals */
    integer i__1, i__2;
    real ret_val, r__1;

    /* Local variables */
    static real b;
    static integer i__, j, k;
    static real q, bp;
    static integer ii;
    static real gbk;
    extern doublereal cot_(real *);
    static real var, rho;
    static integer ndx;
    extern doublereal psi_(real *);
    static real var2, absa;
    extern doublereal poch_(real *, real *);
    static integer incr;
    static real absx, binv, trig, term, poly1, gbern[10];
    extern doublereal r1mach_(integer *);
    static real sinpx2, alneps, alnvar, sqtbig;
    extern doublereal exprel_(real *);
    extern /* Subroutine */ int xermsg_(char *, char *, char *, integer *, 
	    integer *, ftnlen, ftnlen, ftnlen);
    static integer nterms;
    static real sinpxx;

/* ***BEGIN PROLOGUE  POCH1 */
/* ***PURPOSE  Calculate a generalization of Pochhammer's symbol starting */
/*            from first order. */
/* ***LIBRARY   SLATEC (FNLIB) */
/* ***CATEGORY  C1, C7A */
/* ***TYPE      SINGLE PRECISION (POCH1-S, DPOCH1-D) */
/* ***KEYWORDS  FIRST ORDER, FNLIB, POCHHAMMER, SPECIAL FUNCTIONS */
/* ***AUTHOR  Fullerton, W., (LANL) */
/* ***DESCRIPTION */

/* Evaluate a generalization of Pochhammer's symbol for special */
/* situations that require especially accurate values when X is small in */
/*        POCH1(A,X) = (POCH(A,X)-1)/X */
/*                   = (GAMMA(A+X)/GAMMA(A) - 1.0)/X . */
/* This specification is particularly suited for stably computing */
/* expressions such as */
/*        (GAMMA(A+X)/GAMMA(A) - GAMMA(B+X)/GAMMA(B))/X */
/*             = POCH1(A,X) - POCH1(B,X) */
/* Note that POCH1(A,0.0) = PSI(A) */

/* When ABS(X) is so small that substantial cancellation will occur if */
/* the straightforward formula is used, we  use an expansion due */
/* to Fields and discussed by Y. L. Luke, The Special Functions and Their */
/* Approximations, Vol. 1, Academic Press, 1969, page 34. */

/* The ratio POCH(A,X) = GAMMA(A+X)/GAMMA(A) is written by Luke as */
/*        (A+(X-1)/2)**X * polynomial in (A+(X-1)/2)**(-2) . */
/* In order to maintain significance in POCH1, we write for positive A */
/*        (A+(X-1)/2)**X = EXP(X*LOG(A+(X-1)/2)) = EXP(Q) */
/*                       = 1.0 + Q*EXPREL(Q) . */
/* Likewise the polynomial is written */
/*        POLY = 1.0 + X*POLY1(A,X) . */
/* Thus, */
/*        POCH1(A,X) = (POCH(A,X) - 1) / X */
/*                   = EXPREL(Q)*(Q/X + Q*POLY1(A,X)) + POLY1(A,X) */

/* ***REFERENCES  (NONE) */
/* ***ROUTINES CALLED  COT, EXPREL, POCH, PSI, R1MACH, XERMSG */
/* ***REVISION HISTORY  (YYMMDD) */
/*   770801  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890531  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ) */
/*   900727  Added EXTERNAL statement.  (WRB) */
/* ***END PROLOGUE  POCH1 */
/* ***FIRST EXECUTABLE STATEMENT  POCH1 */
    if (first) {
	sqtbig = 1.f / sqrt(r1mach_(&c__1) * 24.f);
	alneps = log(r1mach_(&c__3));
    }
    first = FALSE_;

    if (*x == 0.f) {
	ret_val = psi_(a);
    }
    if (*x == 0.f) {
	return ret_val;
    }

    absx = dabs(*x);
    absa = dabs(*a);
    if (absx > absa * .1f) {
	goto L70;
    }
    if (absx * log((dmax(absa,2.f))) > .1f) {
	goto L70;
    }

    bp = *a;
    if (*a < -.5f) {
	bp = 1.f - *a - *x;
    }
    incr = 0;
    if (bp < 10.f) {
	incr = 11.f - bp;
    }
    b = bp + incr;

    var = b + (*x - 1.f) * .5f;
    alnvar = log(var);
    q = *x * alnvar;

    poly1 = 0.f;
    if (var >= sqtbig) {
	goto L40;
    }
/* Computing 2nd power */
    r__1 = 1.f / var;
    var2 = r__1 * r__1;

    rho = (*x + 1.f) * .5f;
    gbern[0] = 1.f;
    gbern[1] = -rho / 12.f;
    term = var2;
    poly1 = gbern[1] * term;

    nterms = alneps * -.5f / alnvar + 1.f;
    if (nterms > 9) {
	xermsg_("SLATEC", "POCH1", "NTERMS IS TOO BIG, MAYBE R1MACH(3) IS BAD"
		, &c__1, &c__2, (ftnlen)6, (ftnlen)5, (ftnlen)41);
    }
    if (nterms < 2) {
	goto L40;
    }

    i__1 = nterms;
    for (k = 2; k <= i__1; ++k) {
	gbk = 0.f;
	i__2 = k;
	for (j = 1; j <= i__2; ++j) {
	    ndx = k - j + 1;
	    gbk += bern[ndx - 1] * gbern[j - 1];
/* L20: */
	}
	gbern[k] = -rho * gbk / k;

	term = term * ((k << 1) - 2.f - *x) * ((k << 1) - 1.f - *x) * var2;
	poly1 += gbern[k] * term;
/* L30: */
    }

L40:
    poly1 = (*x - 1.f) * poly1;
    ret_val = exprel_(&q) * (alnvar + q * poly1) + poly1;

    if (incr == 0) {
	goto L60;
    }

/* WE HAVE POCH1(B,X).  BUT BP IS SMALL, SO WE USE BACKWARDS RECURSION */
/* TO OBTAIN POCH1(BP,X). */

    i__1 = incr;
    for (ii = 1; ii <= i__1; ++ii) {
	i__ = incr - ii;
	binv = 1.f / (bp + i__);
	ret_val = (ret_val - binv) / (*x * binv + 1.f);
/* L50: */
    }

L60:
    if (bp == *a) {
	return ret_val;
    }

/* WE HAVE POCH1(BP,X), BUT A IS LT -0.5.  WE THEREFORE USE A REFLECTION */
/* FORMULA TO OBTAIN POCH1(A,X). */

    sinpxx = sin(pi * *x) / *x;
    sinpx2 = sin(pi * .5f * *x);
    r__1 = pi * b;
    trig = sinpxx * cot_(&r__1) - sinpx2 * 2.f * (sinpx2 / *x);

    ret_val = trig + (*x * trig + 1.f) * ret_val;
    return ret_val;

L70:
    ret_val = (poch_(a, x) - 1.f) / *x;
    return ret_val;

} /* poch1_ */

