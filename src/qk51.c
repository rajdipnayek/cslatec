/* qk51.f -- translated by f2c (version 12.02.01).
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

static integer c__4 = 4;
static integer c__1 = 1;
static doublereal c_b7 = 1.5;

/* DECK QK51 */
/* Subroutine */ int qk51_(E_fp f, real *a, real *b, real *result, real *
	abserr, real *resabs, real *resasc)
{
    /* Initialized data */

    static real xgk[26] = { .9992621049926098f,.9955569697904981f,
	    .9880357945340772f,.9766639214595175f,.9616149864258425f,
	    .9429745712289743f,.9207471152817016f,.8949919978782754f,
	    .8658470652932756f,.833442628760834f,.7978737979985001f,
	    .7592592630373576f,.7177664068130844f,.6735663684734684f,
	    .6268100990103174f,.577662930241223f,.5263252843347192f,
	    .473002731445715f,.4178853821930377f,.3611723058093878f,
	    .3030895389311078f,.2438668837209884f,.1837189394210489f,
	    .1228646926107104f,.06154448300568508f,0.f };
    static real wgk[26] = { .001987383892330316f,.005561932135356714f,
	    .009473973386174152f,.01323622919557167f,.0168478177091283f,
	    .02043537114588284f,.02400994560695322f,.02747531758785174f,
	    .03079230016738749f,.03400213027432934f,.03711627148341554f,
	    .04008382550403238f,.04287284502017005f,.04550291304992179f,
	    .04798253713883671f,.05027767908071567f,.05236288580640748f,
	    .05425112988854549f,.05595081122041232f,.05743711636156783f,
	    .05868968002239421f,.05972034032417406f,.06053945537604586f,
	    .06112850971705305f,.06147118987142532f,.06158081806783294f };
    static real wg[13] = { .01139379850102629f,.02635498661503214f,
	    .04093915670130631f,.05490469597583519f,.06803833381235692f,
	    .08014070033500102f,.09102826198296365f,.1005359490670506f,
	    .1085196244742637f,.1148582591457116f,.1194557635357848f,
	    .12224244299031f,.1231760537267155f };

    /* System generated locals */
    real r__1, r__2;
    doublereal d__1;

    /* Local variables */
    static integer j;
    static real fc, fv1[25], fv2[25];
    static integer jtw;
    static real absc, resg, resk, fsum, fval1, fval2;
    static integer jtwm1;
    static real hlgth, centr, reskh, uflow;
    extern doublereal r1mach_(integer *);
    static real epmach, dhlgth;

/* ***BEGIN PROLOGUE  QK51 */
/* ***PURPOSE  To compute I = Integral of F over (A,B) with error */
/*                           estimate */
/*                       J = Integral of ABS(F) over (A,B) */
/* ***LIBRARY   SLATEC (QUADPACK) */
/* ***CATEGORY  H2A1A2 */
/* ***TYPE      SINGLE PRECISION (QK51-S, DQK51-D) */
/* ***KEYWORDS  51-POINT GAUSS-KRONROD RULES, QUADPACK, QUADRATURE */
/* ***AUTHOR  Piessens, Robert */
/*             Applied Mathematics and Programming Division */
/*             K. U. Leuven */
/*           de Doncker, Elise */
/*             Applied Mathematics and Programming Division */
/*             K. U. Leuven */
/* ***DESCRIPTION */

/*           Integration rules */
/*           Standard fortran subroutine */
/*           Real version */

/*           PARAMETERS */
/*            ON ENTRY */
/*              F      - Real */
/*                       Function subroutine defining the integrand */
/*                       function F(X). The actual name for F needs to be */
/*                       declared E X T E R N A L in the calling program. */

/*              A      - Real */
/*                       Lower limit of integration */

/*              B      - Real */
/*                       Upper limit of integration */

/*            ON RETURN */
/*              RESULT - Real */
/*                       Approximation to the integral I */
/*                       RESULT is computed by applying the 51-point */
/*                       Kronrod rule (RESK) obtained by optimal addition */
/*                       of abscissae to the 25-point Gauss rule (RESG). */

/*              ABSERR - Real */
/*                       Estimate of the modulus of the absolute error, */
/*                       which should not exceed ABS(I-RESULT) */

/*              RESABS - Real */
/*                       Approximation to the integral J */

/*              RESASC - Real */
/*                       Approximation to the integral of ABS(F-I/(B-A)) */
/*                       over (A,B) */

/* ***REFERENCES  (NONE) */
/* ***ROUTINES CALLED  R1MACH */
/* ***REVISION HISTORY  (YYMMDD) */
/*   800101  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890531  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/* ***END PROLOGUE  QK51 */



/*           THE ABSCISSAE AND WEIGHTS ARE GIVEN FOR THE INTERVAL (-1,1). */
/*           BECAUSE OF SYMMETRY ONLY THE POSITIVE ABSCISSAE AND THEIR */
/*           CORRESPONDING WEIGHTS ARE GIVEN. */

/*           XGK    - ABSCISSAE OF THE 51-POINT KRONROD RULE */
/*                    XGK(2), XGK(4), ...  ABSCISSAE OF THE 25-POINT */
/*                    GAUSS RULE */
/*                    XGK(1), XGK(3), ...  ABSCISSAE WHICH ARE OPTIMALLY */
/*                    ADDED TO THE 25-POINT GAUSS RULE */

/*           WGK    - WEIGHTS OF THE 51-POINT KRONROD RULE */

/*           WG     - WEIGHTS OF THE 25-POINT GAUSS RULE */



/*           LIST OF MAJOR VARIABLES */
/*           ----------------------- */

/*           CENTR  - MID POINT OF THE INTERVAL */
/*           HLGTH  - HALF-LENGTH OF THE INTERVAL */
/*           ABSC   - ABSCISSA */
/*           FVAL*  - FUNCTION VALUE */
/*           RESG   - RESULT OF THE 25-POINT GAUSS FORMULA */
/*           RESK   - RESULT OF THE 51-POINT KRONROD FORMULA */
/*           RESKH  - APPROXIMATION TO THE MEAN VALUE OF F OVER (A,B), */
/*                    I.E. TO I/(B-A) */

/*           MACHINE DEPENDENT CONSTANTS */
/*           --------------------------- */

/*           EPMACH IS THE LARGEST RELATIVE SPACING. */
/*           UFLOW IS THE SMALLEST POSITIVE MAGNITUDE. */

/* ***FIRST EXECUTABLE STATEMENT  QK51 */
    epmach = r1mach_(&c__4);
    uflow = r1mach_(&c__1);

    centr = (*a + *b) * .5f;
    hlgth = (*b - *a) * .5f;
    dhlgth = dabs(hlgth);

/*           COMPUTE THE 51-POINT KRONROD APPROXIMATION TO */
/*           THE INTEGRAL, AND ESTIMATE THE ABSOLUTE ERROR. */

    fc = (*f)(&centr);
    resg = wg[12] * fc;
    resk = wgk[25] * fc;
    *resabs = dabs(resk);
    for (j = 1; j <= 12; ++j) {
	jtw = j << 1;
	absc = hlgth * xgk[jtw - 1];
	r__1 = centr - absc;
	fval1 = (*f)(&r__1);
	r__1 = centr + absc;
	fval2 = (*f)(&r__1);
	fv1[jtw - 1] = fval1;
	fv2[jtw - 1] = fval2;
	fsum = fval1 + fval2;
	resg += wg[j - 1] * fsum;
	resk += wgk[jtw - 1] * fsum;
	*resabs += wgk[jtw - 1] * (dabs(fval1) + dabs(fval2));
/* L10: */
    }
    for (j = 1; j <= 13; ++j) {
	jtwm1 = (j << 1) - 1;
	absc = hlgth * xgk[jtwm1 - 1];
	r__1 = centr - absc;
	fval1 = (*f)(&r__1);
	r__1 = centr + absc;
	fval2 = (*f)(&r__1);
	fv1[jtwm1 - 1] = fval1;
	fv2[jtwm1 - 1] = fval2;
	fsum = fval1 + fval2;
	resk += wgk[jtwm1 - 1] * fsum;
	*resabs += wgk[jtwm1 - 1] * (dabs(fval1) + dabs(fval2));
/* L15: */
    }
    reskh = resk * .5f;
    *resasc = wgk[25] * (r__1 = fc - reskh, dabs(r__1));
    for (j = 1; j <= 25; ++j) {
	*resasc += wgk[j - 1] * ((r__1 = fv1[j - 1] - reskh, dabs(r__1)) + (
		r__2 = fv2[j - 1] - reskh, dabs(r__2)));
/* L20: */
    }
    *result = resk * hlgth;
    *resabs *= dhlgth;
    *resasc *= dhlgth;
    *abserr = (r__1 = (resk - resg) * hlgth, dabs(r__1));
    if (*resasc != 0.f && *abserr != 0.f) {
/* Computing MIN */
	d__1 = (doublereal) (*abserr * 200.f / *resasc);
	r__1 = 1.f, r__2 = pow_dd(&d__1, &c_b7);
	*abserr = *resasc * dmin(r__1,r__2);
    }
    if (*resabs > uflow / (epmach * 50.f)) {
/* Computing MAX */
	r__1 = epmach * 50.f * *resabs;
	*abserr = dmax(r__1,*abserr);
    }
    return 0;
} /* qk51_ */

