/* r9knus.f -- translated by f2c (version 12.02.01).
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
static integer c__16 = 16;
static integer c__12 = 12;
static integer c__1 = 1;
static integer c__2 = 2;

/* DECK R9KNUS */
/* Subroutine */ int r9knus_(real *xnu, real *x, real *bknu, real *bknu1, 
	integer *iswtch)
{
    /* Initialized data */

    static real c0kcs[16] = { .060183057242626108f,-.15364871433017286f,
	    -.011751176008210492f,-8.52487888919795e-4f,-6.1329838767496e-5f,
	    -4.405228124551e-6f,-3.16312467283e-7f,-2.2710719382e-8f,
	    -1.63056446e-9f,-1.17069392e-10f,-8.405206e-12f,-6.03466e-13f,
	    -4.3326e-14f,-3.11e-15f,-2.23e-16f,-1.6e-17f };
    static real znu1cs[12] = { .20330675699419173f,.14007793341321977f,
	    .007916796961001613f,3.39801182532104e-4f,1.1741975688989e-5f,
	    3.39357570612e-7f,8.425941769e-9f,1.83336677e-10f,3.549698e-12f,
	    6.1903e-14f,9.81e-16f,1.4e-17f };
    static real euler = .57721566490153286f;
    static real sqpi2 = 1.2533141373155003f;
    static real aln2 = .69314718055994531f;
    static logical first = TRUE_;

    /* System generated locals */
    integer i__1, i__2;
    real r__1, r__2;

    /* Local variables */
    static real a[15];
    static integer i__, n;
    static real v, z__, a0, b0, c0, p1, p2, p3, an, bn;
    static integer ii;
    static real xi, qq, x2n;
    static integer inu;
    static real xmu, beta[15], alnz, xsml, expx, vlnz, ztov;
    static integer ntc0k;
    static real bknu0;
    extern doublereal gamma_(real *);
    static real x2tov, alpha[15], bknud;
    extern doublereal csevl_(real *, real *, integer *);
    extern integer inits_(real *, integer *, real *);
    extern doublereal r1mach_(integer *);
    static real sqrtx, alnbig;
    static integer ntznu1;
    static real alneps, alnsml;
    extern /* Subroutine */ int xermsg_(char *, char *, char *, integer *, 
	    integer *, ftnlen, ftnlen, ftnlen);
    static integer nterms;
    static real result, xnusml;

/* ***BEGIN PROLOGUE  R9KNUS */
/* ***SUBSIDIARY */
/* ***PURPOSE  Compute Bessel functions EXP(X)*K-SUB-XNU(X) and EXP(X)* */
/*            K-SUB-XNU+1(X) for 0.0 .LE. XNU .LT. 1.0. */
/* ***LIBRARY   SLATEC (FNLIB) */
/* ***CATEGORY  C10B3 */
/* ***TYPE      SINGLE PRECISION (R9KNUS-S, D9KNUS-D) */
/* ***KEYWORDS  BESSEL FUNCTION, FNLIB, SPECIAL FUNCTIONS */
/* ***AUTHOR  Fullerton, W., (LANL) */
/* ***DESCRIPTION */

/* Compute Bessel functions EXP(X) * K-sub-XNU (X)  and */
/* EXP(X) * K-sub-XNU+1 (X) for 0.0 .LE. XNU .LT. 1.0 . */

/* Series for C0K        on the interval  0.          to  2.50000D-01 */
/*                                        with weighted error   1.60E-17 */
/*                                         log weighted error  16.79 */
/*                               significant figures required  15.99 */
/*                                    decimal places required  17.40 */

/* Series for ZNU1       on the interval -7.00000D-01 to  0. */
/*                                        with weighted error   1.43E-17 */
/*                                         log weighted error  16.85 */
/*                               significant figures required  16.08 */
/*                                    decimal places required  17.38 */

/* ***REFERENCES  (NONE) */
/* ***ROUTINES CALLED  CSEVL, GAMMA, INITS, R1MACH, XERMSG */
/* ***REVISION HISTORY  (YYMMDD) */
/*   770601  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890531  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ) */
/*   900720  Routine changed from user-callable to subsidiary.  (WRB) */
/*   900727  Added EXTERNAL statement.  (WRB) */
/*   920618  Removed space from variable names.  (RWC, WRB) */
/* ***END PROLOGUE  R9KNUS */
/* ***FIRST EXECUTABLE STATEMENT  R9KNUS */
    if (first) {
	r__1 = r1mach_(&c__3) * .1f;
	ntc0k = inits_(c0kcs, &c__16, &r__1);
	r__1 = r1mach_(&c__3) * .1f;
	ntznu1 = inits_(znu1cs, &c__12, &r__1);

	xnusml = sqrt(r1mach_(&c__3) / 8.f);
	xsml = r1mach_(&c__3) * .1f;
	alnsml = log(r1mach_(&c__1));
	alnbig = log(r1mach_(&c__2));
	alneps = log(r1mach_(&c__3) * .1f);
    }
    first = FALSE_;

    if (*xnu < 0.f || *xnu >= 1.f) {
	xermsg_("SLATEC", "R9KNUS", "XNU MUST BE GE 0 AND LT 1", &c__1, &c__2,
		 (ftnlen)6, (ftnlen)6, (ftnlen)25);
    }
    if (*x <= 0.f) {
	xermsg_("SLATEC", "R9KNUS", "X MUST BE GT 0", &c__2, &c__2, (ftnlen)6,
		 (ftnlen)6, (ftnlen)14);
    }

    *iswtch = 0;
    if (*x > 2.f) {
	goto L50;
    }

/* X IS SMALL.  COMPUTE K-SUB-XNU (X) AND THE DERIVATIVE OF K-SUB-XNU (X) */
/* THEN FIND K-SUB-XNU+1 (X).  XNU IS REDUCED TO THE INTERVAL (-.5,+.5) */
/* THEN TO (0., .5), BECAUSE K OF NEGATIVE ORDER (-NU) = K OF POSITIVE */
/* ORDER (+NU). */

    v = *xnu;
    if (*xnu > .5f) {
	v = 1.f - *xnu;
    }

/* CAREFULLY FIND (X/2)**XNU AND Z**XNU WHERE Z = X*X/4. */
    alnz = (log(*x) - aln2) * 2.f;

    if (*x > *xnu) {
	goto L20;
    }
    if (*xnu * -.5f * alnz - aln2 - log(*xnu) > alnbig) {
	xermsg_("SLATEC", "R9KNUS", "X SO SMALL BESSEL K-SUB-XNU OVERFLOWS", &
		c__3, &c__2, (ftnlen)6, (ftnlen)6, (ftnlen)37);
    }

L20:
    vlnz = v * alnz;
    x2tov = exp(vlnz * .5f);
    ztov = 0.f;
    if (vlnz > alnsml) {
/* Computing 2nd power */
	r__1 = x2tov;
	ztov = r__1 * r__1;
    }

    r__1 = v + 1.f;
    a0 = gamma_(&r__1) * .5f;
    r__1 = 1.f - v;
    b0 = gamma_(&r__1) * .5f;
    c0 = -euler;
    if (ztov > .5f && v > xnusml) {
	r__1 = v * 8.f * v - 1.f;
	c0 = csevl_(&r__1, c0kcs, &ntc0k) - .75f;
    }

    if (ztov <= .5f) {
	alpha[0] = (a0 - ztov * b0) / v;
    }
    if (ztov > .5f) {
	r__1 = vlnz / .35f + 1.f;
	alpha[0] = c0 - alnz * (csevl_(&r__1, znu1cs, &ntznu1) + .75f) * b0;
    }
    beta[0] = (a0 + ztov * b0) * -.5f;

    z__ = 0.f;
    if (*x > xsml) {
	z__ = *x * .25f * *x;
    }
/* Computing MAX */
    r__1 = 2.f, r__2 = (alnz * 8.f - 25.19f - alneps) / (4.28f - alnz) + 11.f;
    nterms = dmax(r__1,r__2);
    i__1 = nterms;
    for (i__ = 2; i__ <= i__1; ++i__) {
	xi = (real) (i__ - 1);
	a0 /= xi * (xi - v);
	b0 /= xi * (xi + v);
	alpha[i__ - 1] = (alpha[i__ - 2] + xi * 2.f * a0) / (xi * (xi + v));
	beta[i__ - 1] = (xi - v * .5f) * alpha[i__ - 1] - ztov * b0;
/* L30: */
    }

    *bknu = alpha[nterms - 1];
    bknud = beta[nterms - 1];
    i__1 = nterms;
    for (ii = 2; ii <= i__1; ++ii) {
	i__ = nterms + 1 - ii;
	*bknu = alpha[i__ - 1] + *bknu * z__;
	bknud = beta[i__ - 1] + bknud * z__;
/* L40: */
    }

    expx = exp(*x);
    *bknu = expx * *bknu / x2tov;

    if ((*xnu + 1.f) * -.5f * alnz - aln2 * 2.f > alnbig) {
	*iswtch = 1;
    }
    if (*iswtch == 1) {
	return 0;
    }
    bknud = expx * bknud * 2.f / (x2tov * *x);

    if (*xnu <= .5f) {
	*bknu1 = v * *bknu / *x - bknud;
    }
    if (*xnu <= .5f) {
	return 0;
    }

    bknu0 = *bknu;
    *bknu = -v * *bknu / *x - bknud;
    *bknu1 = *xnu * 2.f * *bknu / *x + bknu0;
    return 0;

/* X IS LARGE.  FIND K-SUB-XNU (X) AND K-SUB-XNU+1 (X) WITH Y. L. LUKE-S */
/* RATIONAL EXPANSION. */

L50:
    sqrtx = sqrt(*x);
    if (*x > 1.f / xsml) {
	goto L90;
    }
    an = 4.f / *x - 1.56f;
    bn = -.29f - .22f / *x;
/* Computing MIN */
/* Computing MAX */
    r__1 = 3.f, r__2 = an + bn * alneps;
    i__1 = 15, i__2 = (integer) dmax(r__1,r__2);
    nterms = min(i__1,i__2);

    for (inu = 1; inu <= 2; ++inu) {
	xmu = 0.f;
	if (inu == 1 && *xnu > xnusml) {
	    xmu = *xnu * 4.f * *xnu;
	}
	if (inu == 2) {
/* Computing 2nd power */
	    r__1 = dabs(*xnu) + 1.f;
	    xmu = r__1 * r__1 * 4.f;
	}

	a[0] = 1.f - xmu;
	a[1] = 9.f - xmu;
	a[2] = 25.f - xmu;
	if (a[1] == 0.f) {
	    result = sqpi2 * (*x * 16.f + xmu + 7.f) / (*x * 16.f * sqrtx);
	}
	if (a[1] == 0.f) {
	    goto L70;
	}

	alpha[0] = 1.f;
	alpha[1] = (*x * 16.f + a[1]) / a[1];
	alpha[2] = ((*x * 768.f + a[2] * 48.f) * *x + a[1] * a[2]) / (a[1] * 
		a[2]);

	beta[0] = 1.f;
	beta[1] = (*x * 16.f + (xmu + 7.f)) / a[1];
	beta[2] = ((*x * 768.f + (xmu + 23.f) * 48.f) * *x + ((xmu + 62.f) * 
		xmu + 129.f)) / (a[1] * a[2]);

	if (nterms < 4) {
	    goto L65;
	}
	i__1 = nterms;
	for (i__ = 4; i__ <= i__1; ++i__) {
	    n = i__ - 1;
	    x2n = (real) ((n << 1) - 1);

/* Computing 2nd power */
	    r__1 = x2n + 2.f;
	    a[i__ - 1] = r__1 * r__1 - xmu;
	    qq = x2n * 16.f / a[i__ - 1];
	    p1 = -x2n * (n * 12 * n - n * 20 - a[0]) / ((x2n - 2.f) * a[i__ - 
		    1]) - qq * *x;
	    p2 = (n * 12 * n - n * 28 + 8 - a[0]) / a[i__ - 1] - qq * *x;
	    p3 = -x2n * a[i__ - 4] / ((x2n - 2.f) * a[i__ - 1]);

	    alpha[i__ - 1] = -p1 * alpha[i__ - 2] - p2 * alpha[i__ - 3] - p3 *
		     alpha[i__ - 4];
	    beta[i__ - 1] = -p1 * beta[i__ - 2] - p2 * beta[i__ - 3] - p3 * 
		    beta[i__ - 4];
/* L60: */
	}

L65:
	result = sqpi2 * beta[nterms - 1] / (sqrtx * alpha[nterms - 1]);

L70:
	if (inu == 1) {
	    *bknu = result;
	}
	if (inu == 2) {
	    *bknu1 = result;
	}
/* L80: */
    }
    return 0;

L90:
    *bknu = sqpi2 / sqrtx;
    *bknu1 = *bknu;
    return 0;

} /* r9knus_ */

