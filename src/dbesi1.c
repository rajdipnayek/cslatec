/* dbesi1.f -- translated by f2c (version 12.02.01).
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
static integer c__17 = 17;
static integer c__1 = 1;
static integer c__2 = 2;

/* DECK DBESI1 */
doublereal dbesi1_(doublereal *x)
{
    /* Initialized data */

    static doublereal bi1cs[17] = { -.0019717132610998597316138503218149,
	    .40734887667546480608155393652014,
	    .034838994299959455866245037783787,
	    .0015453945563001236038598401058489,
	    4.188852109837778412945883200412e-5,
	    7.6490267648362114741959703966069e-7,
	    1.0042493924741178689179808037238e-8,
	    9.9322077919238106481371298054863e-11,
	    7.6638017918447637275200171681349e-13,
	    4.741418923816739498038809194816e-15,
	    2.4041144040745181799863172032e-17,
	    1.0171505007093713649121100799999e-19,
	    3.6450935657866949458491733333333e-22,
	    1.1205749502562039344810666666666e-24,2.9875441934468088832e-27,
	    6.9732310939194709333333333333333e-30,1.43679482206208e-32 };
    static logical first = TRUE_;

    /* System generated locals */
    real r__1;
    doublereal ret_val, d__1;

    /* Local variables */
    static doublereal y;
    static integer nti1;
    static doublereal xmin, xmax, xsml;
    extern doublereal d1mach_(integer *), dbsi1e_(doublereal *), dcsevl_(
	    doublereal *, doublereal *, integer *);
    extern integer initds_(doublereal *, integer *, real *);
    extern /* Subroutine */ int xermsg_(char *, char *, char *, integer *, 
	    integer *, ftnlen, ftnlen, ftnlen);

/* ***BEGIN PROLOGUE  DBESI1 */
/* ***PURPOSE  Compute the modified (hyperbolic) Bessel function of the */
/*            first kind of order one. */
/* ***LIBRARY   SLATEC (FNLIB) */
/* ***CATEGORY  C10B1 */
/* ***TYPE      DOUBLE PRECISION (BESI1-S, DBESI1-D) */
/* ***KEYWORDS  FIRST KIND, FNLIB, HYPERBOLIC BESSEL FUNCTION, */
/*             MODIFIED BESSEL FUNCTION, ORDER ONE, SPECIAL FUNCTIONS */
/* ***AUTHOR  Fullerton, W., (LANL) */
/* ***DESCRIPTION */

/* DBESI1(X) calculates the double precision modified (hyperbolic) */
/* Bessel function of the first kind of order one and double precision */
/* argument X. */

/* Series for BI1        on the interval  0.          to  9.00000E+00 */
/*                                        with weighted error   1.44E-32 */
/*                                         log weighted error  31.84 */
/*                               significant figures required  31.45 */
/*                                    decimal places required  32.46 */

/* ***REFERENCES  (NONE) */
/* ***ROUTINES CALLED  D1MACH, DBSI1E, DCSEVL, INITDS, XERMSG */
/* ***REVISION HISTORY  (YYMMDD) */
/*   770701  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890531  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ) */
/* ***END PROLOGUE  DBESI1 */
/* ***FIRST EXECUTABLE STATEMENT  DBESI1 */
    if (first) {
	r__1 = (real) d1mach_(&c__3) * .1f;
	nti1 = initds_(bi1cs, &c__17, &r__1);
	xmin = d1mach_(&c__1) * 2.;
	xsml = sqrt(d1mach_(&c__3) * 4.5);
	xmax = log(d1mach_(&c__2));
    }
    first = FALSE_;

    y = abs(*x);
    if (y > 3.) {
	goto L20;
    }

    ret_val = 0.;
    if (y == 0.) {
	return ret_val;
    }

    if (y <= xmin) {
	xermsg_("SLATEC", "DBESI1", "ABS(X) SO SMALL I1 UNDERFLOWS", &c__1, &
		c__1, (ftnlen)6, (ftnlen)6, (ftnlen)29);
    }
    if (y > xmin) {
	ret_val = *x * .5;
    }
    if (y > xsml) {
	d__1 = y * y / 4.5 - 1.;
	ret_val = *x * (dcsevl_(&d__1, bi1cs, &nti1) + .875);
    }
    return ret_val;

L20:
    if (y > xmax) {
	xermsg_("SLATEC", "DBESI1", "ABS(X) SO BIG I1 OVERFLOWS", &c__2, &
		c__2, (ftnlen)6, (ftnlen)6, (ftnlen)26);
    }

    ret_val = exp(y) * dbsi1e_(x);

    return ret_val;
} /* dbesi1_ */

