/* daws.f -- translated by f2c (version 12.02.01).
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
static integer c__29 = 29;
static integer c__26 = 26;
static integer c__1 = 1;
static integer c__2 = 2;

/* DECK DAWS */
doublereal daws_(real *x)
{
    /* Initialized data */

    static real dawcs[13] = { -.006351734375145949f,-.22940714796773869f,
	    .022130500939084764f,-.001549265453892985f,8.4973277156849e-5f,
	    -3.828266270972e-6f,1.46285480625e-7f,-4.851982381e-9f,
	    1.42146357e-10f,-3.728836e-12f,8.8549e-14f,-1.92e-15f,3.8e-17f };
    static real daw2cs[29] = { -.056886544105215527f,-.31811346996168131f,
	    .20873845413642237f,-.12475409913779131f,.067869305186676777f,
	    -.03365914489527094f,.015260781271987972f,-.006348370962596214f,
	    .002432674092074852f,-8.6219541491065e-4f,2.83765733363216e-4f,
	    -8.705754987417e-5f,2.4986849985481e-5f,-6.731928676416e-6f,
	    1.707857878557e-6f,-4.09175512264e-7f,9.2828292216e-8f,
	    -1.999140361e-8f,4.096349064e-9f,-8.00324095e-10f,1.49385031e-10f,
	    -2.6687999e-11f,4.571221e-12f,-7.51873e-13f,1.18931e-13f,
	    -1.8116e-14f,2.661e-15f,-3.77e-16f,5.1e-17f };
    static real dawacs[26] = { .01690485637765704f,.00868325227840695f,
	    2.4248640424177e-4f,1.261182399572e-5f,1.06645331463e-6f,
	    1.3581597947e-7f,2.171042356e-8f,2.8670105e-9f,-1.9013363e-10f,
	    -3.0977804e-10f,-1.0294148e-10f,-6.26035e-12f,8.56313e-12f,
	    3.03304e-12f,-2.5236e-13f,-4.2106e-13f,-4.431e-14f,4.911e-14f,
	    1.235e-14f,-5.78e-15f,-2.28e-15f,7.6e-16f,3.8e-16f,-1.1e-16f,
	    -6e-17f,2e-17f };
    static logical first = TRUE_;

    /* System generated locals */
    real ret_val, r__1, r__2;

    /* Local variables */
    static real y, eps, xbig, xmax, xsml;
    extern doublereal csevl_(real *, real *, integer *);
    static integer ntdaw;
    extern integer inits_(real *, integer *, real *);
    extern doublereal r1mach_(integer *);
    static integer ntdaw2, ntdawa;
    extern /* Subroutine */ int xermsg_(char *, char *, char *, integer *, 
	    integer *, ftnlen, ftnlen, ftnlen);

/* ***BEGIN PROLOGUE  DAWS */
/* ***PURPOSE  Compute Dawson's function. */
/* ***LIBRARY   SLATEC (FNLIB) */
/* ***CATEGORY  C8C */
/* ***TYPE      SINGLE PRECISION (DAWS-S, DDAWS-D) */
/* ***KEYWORDS  DAWSON'S FUNCTION, FNLIB, SPECIAL FUNCTIONS */
/* ***AUTHOR  Fullerton, W., (LANL) */
/* ***DESCRIPTION */

/* DAWS(X) calculates Dawson's integral for real argument X. */

/* Series for DAW        on the interval  0.          to  1.00000D+00 */
/*                                        with weighted error   3.83E-17 */
/*                                         log weighted error  16.42 */
/*                               significant figures required  15.78 */
/*                                    decimal places required  16.97 */

/* Series for DAW2       on the interval  0.          to  1.60000D+01 */
/*                                        with weighted error   5.17E-17 */
/*                                         log weighted error  16.29 */
/*                               significant figures required  15.90 */
/*                                    decimal places required  17.02 */

/* Series for DAWA       on the interval  0.          to  6.25000D-02 */
/*                                        with weighted error   2.24E-17 */
/*                                         log weighted error  16.65 */
/*                               significant figures required  14.73 */
/*                                    decimal places required  17.36 */

/* ***REFERENCES  (NONE) */
/* ***ROUTINES CALLED  CSEVL, INITS, R1MACH, XERMSG */
/* ***REVISION HISTORY  (YYMMDD) */
/*   780401  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890531  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ) */
/*   920618  Removed space from variable names.  (RWC, WRB) */
/* ***END PROLOGUE  DAWS */
/* ***FIRST EXECUTABLE STATEMENT  DAWS */
    if (first) {
	eps = r1mach_(&c__3);
	r__1 = eps * .1f;
	ntdaw = inits_(dawcs, &c__13, &r__1);
	r__1 = eps * .1f;
	ntdaw2 = inits_(daw2cs, &c__29, &r__1);
	r__1 = eps * .1f;
	ntdawa = inits_(dawacs, &c__26, &r__1);

	xsml = sqrt(eps * 1.5f);
	xbig = sqrt(.5f / eps);
/* Computing MIN */
	r__1 = -log(r1mach_(&c__1) * 2.f), r__2 = log(r1mach_(&c__2));
	xmax = exp(dmin(r__1,r__2) - 1.f);
    }
    first = FALSE_;

    y = dabs(*x);
    if (y > 1.f) {
	goto L20;
    }

    ret_val = *x;
    if (y <= xsml) {
	return ret_val;
    }

    r__1 = y * 2.f * y - 1.f;
    ret_val = *x * (csevl_(&r__1, dawcs, &ntdaw) + .75f);
    return ret_val;

L20:
    if (y > 4.f) {
	goto L30;
    }
    r__1 = y * .125f * y - 1.f;
    ret_val = *x * (csevl_(&r__1, daw2cs, &ntdaw2) + .25f);
    return ret_val;

L30:
    if (y > xmax) {
	goto L40;
    }
    ret_val = .5f / *x;
    if (y > xbig) {
	return ret_val;
    }

/* Computing 2nd power */
    r__2 = y;
    r__1 = 32.f / (r__2 * r__2) - 1.f;
    ret_val = (csevl_(&r__1, dawacs, &ntdawa) + .5f) / *x;
    return ret_val;

L40:
    xermsg_("SLATEC", "DAWS", "ABS(X) SO LARGE DAWS UNDERFLOWS", &c__1, &c__1,
	     (ftnlen)6, (ftnlen)4, (ftnlen)31);
    ret_val = 0.f;
    return ret_val;

} /* daws_ */

