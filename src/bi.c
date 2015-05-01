/* bi.f -- translated by f2c (version 12.02.01).
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
static integer c__9 = 9;
static integer c__8 = 8;
static integer c__10 = 10;
static doublereal c_b7 = .3333;
static integer c__2 = 2;
static doublereal c_b9 = .6666;
static integer c__1 = 1;

/* DECK BI */
doublereal bi_(real *x)
{
    /* Initialized data */

    static real bifcs[9] = { -.01673021647198664948f,.1025233583424944561f,
	    .00170830925073815165f,1.186254546774468e-5f,4.493290701779e-8f,
	    1.0698207143e-10f,1.7480643e-13f,2.081e-16f,1.8e-19f };
    static real bigcs[8] = { .02246622324857452f,.03736477545301955f,
	    4.4476218957212e-4f,2.47080756363e-6f,7.91913533e-9f,
	    1.649807e-11f,2.411e-14f,2e-17f };
    static real bif2cs[10] = { .09984572693816041f,.478624977863005538f,
	    .0251552119604330118f,5.820693885232645e-4f,7.4997659644377e-6f,
	    6.13460287034e-8f,3.462753885e-10f,1.428891e-12f,4.4962e-15f,
	    1.11e-17f };
    static real big2cs[10] = { .03330566214551434f,.161309215123197068f,
	    .0063190073096134286f,1.187904568162517e-4f,1.30453458862e-6f,
	    9.3741259955e-9f,4.74580188e-11f,1.783107e-13f,5.167e-16f,
	    1.1e-18f };
    static logical first = TRUE_;

    /* System generated locals */
    real ret_val, r__1;
    doublereal d__1;

    /* Local variables */
    static real z__, xm;
    extern doublereal bie_(real *);
    static real eta;
    static integer nbif, nbig;
    static real xmax;
    static integer nbif2, nbig2;
    static real x3sml, theta;
    extern doublereal csevl_(real *, real *, integer *);
    extern integer inits_(real *, integer *, real *);
    extern doublereal r1mach_(integer *);
    extern /* Subroutine */ int r9aimp_(real *, real *, real *), xermsg_(char 
	    *, char *, char *, integer *, integer *, ftnlen, ftnlen, ftnlen);

/* ***BEGIN PROLOGUE  BI */
/* ***PURPOSE  Evaluate the Bairy function (the Airy function of the */
/*            second kind). */
/* ***LIBRARY   SLATEC (FNLIB) */
/* ***CATEGORY  C10D */
/* ***TYPE      SINGLE PRECISION (BI-S, DBI-D) */
/* ***KEYWORDS  BAIRY FUNCTION, FNLIB, SPECIAL FUNCTIONS */
/* ***AUTHOR  Fullerton, W., (LANL) */
/* ***DESCRIPTION */

/* BI(X) calculates the Airy function of the second kind for real */
/* argument X. */

/* Series for BIF        on the interval -1.00000D+00 to  1.00000D+00 */
/*                                        with weighted error   1.88E-19 */
/*                                         log weighted error  18.72 */
/*                               significant figures required  17.74 */
/*                                    decimal places required  19.20 */

/* Series for BIG        on the interval -1.00000D+00 to  1.00000D+00 */
/*                                        with weighted error   2.61E-17 */
/*                                         log weighted error  16.58 */
/*                               significant figures required  15.17 */
/*                                    decimal places required  17.03 */

/* Series for BIF2       on the interval  1.00000D+00 to  8.00000D+00 */
/*                                        with weighted error   1.11E-17 */
/*                                         log weighted error  16.95 */
/*                        approx significant figures required  16.5 */
/*                                    decimal places required  17.45 */

/* Series for BIG2       on the interval  1.00000D+00 to  8.00000D+00 */
/*                                        with weighted error   1.19E-18 */
/*                                         log weighted error  17.92 */
/*                        approx significant figures required  17.2 */
/*                                    decimal places required  18.42 */

/* ***REFERENCES  (NONE) */
/* ***ROUTINES CALLED  BIE, CSEVL, INITS, R1MACH, R9AIMP, XERMSG */
/* ***REVISION HISTORY  (YYMMDD) */
/*   770701  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890531  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ) */
/*   900326  Removed duplicate information from DESCRIPTION section. */
/*           (WRB) */
/* ***END PROLOGUE  BI */
/* ***FIRST EXECUTABLE STATEMENT  BI */
    if (first) {
	eta = r1mach_(&c__3) * .1f;
	nbif = inits_(bifcs, &c__9, &eta);
	nbig = inits_(bigcs, &c__8, &eta);
	nbif2 = inits_(bif2cs, &c__10, &eta);
	nbig2 = inits_(big2cs, &c__10, &eta);

	d__1 = (doublereal) eta;
	x3sml = pow_dd(&d__1, &c_b7);
	d__1 = (doublereal) (log(r1mach_(&c__2)) * 1.5f);
	xmax = pow_dd(&d__1, &c_b9);
    }
    first = FALSE_;

    if (*x >= -1.f) {
	goto L20;
    }
    r9aimp_(x, &xm, &theta);
    ret_val = xm * sin(theta);
    return ret_val;

L20:
    if (*x > 1.f) {
	goto L30;
    }
    z__ = 0.f;
    if (dabs(*x) > x3sml) {
/* Computing 3rd power */
	r__1 = *x;
	z__ = r__1 * (r__1 * r__1);
    }
    ret_val = csevl_(&z__, bifcs, &nbif) + .625f + *x * (csevl_(&z__, bigcs, &
	    nbig) + .4375f);
    return ret_val;

L30:
    if (*x > 2.f) {
	goto L40;
    }
/* Computing 3rd power */
    r__1 = *x;
    z__ = (r__1 * (r__1 * r__1) * 2.f - 9.f) / 7.f;
    ret_val = csevl_(&z__, bif2cs, &nbif2) + 1.125f + *x * (csevl_(&z__, 
	    big2cs, &nbig2) + .625f);
    return ret_val;

L40:
    if (*x > xmax) {
	xermsg_("SLATEC", "BI", "X SO BIG THAT BI OVERFLOWS", &c__1, &c__2, (
		ftnlen)6, (ftnlen)2, (ftnlen)26);
    }

    ret_val = bie_(x) * exp(*x * 2.f * sqrt(*x) / 3.f);
    return ret_val;

} /* bi_ */

