/* qk31.f -- translated by f2c (version 12.02.01).
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

/* DECK QK31 */
/* Subroutine */ int qk31_(E_fp f, real *a, real *b, real *result, real *
	abserr, real *resabs, real *resasc)
{
    /* Initialized data */

    static real xgk[16] = { .9980022986933971f,.9879925180204854f,
	    .9677390756791391f,.9372733924007059f,.8972645323440819f,
	    .8482065834104272f,.7904185014424659f,.72441773136017f,
	    .650996741297417f,.5709721726085388f,.4850818636402397f,
	    .3941513470775634f,.2991800071531688f,.2011940939974345f,
	    .1011420669187175f,0.f };
    static real wgk[16] = { .005377479872923349f,.01500794732931612f,
	    .02546084732671532f,.03534636079137585f,.04458975132476488f,
	    .05348152469092809f,.06200956780067064f,.06985412131872826f,
	    .07684968075772038f,.08308050282313302f,.08856444305621177f,
	    .09312659817082532f,.09664272698362368f,.09917359872179196f,
	    .1007698455238756f,.1013300070147915f };
    static real wg[8] = { .03075324199611727f,.07036604748810812f,
	    .1071592204671719f,.1395706779261543f,.1662692058169939f,
	    .1861610000155622f,.1984314853271116f,.2025782419255613f };

    /* System generated locals */
    real r__1, r__2;
    doublereal d__1;

    /* Local variables */
    static integer j;
    static real fc, fv1[15], fv2[15];
    static integer jtw;
    static real absc, resg, resk, fsum, fval1, fval2;
    static integer jtwm1;
    static real hlgth, centr, reskh, uflow;
    extern doublereal r1mach_(integer *);
    static real epmach, dhlgth;

/* ***BEGIN PROLOGUE  QK31 */
/* ***PURPOSE  To compute I = Integral of F over (A,B) with error */
/*                           estimate */
/*                       J = Integral of ABS(F) over (A,B) */
/* ***LIBRARY   SLATEC (QUADPACK) */
/* ***CATEGORY  H2A1A2 */
/* ***TYPE      SINGLE PRECISION (QK31-S, DQK31-D) */
/* ***KEYWORDS  31-POINT GAUSS-KRONROD RULES, QUADPACK, QUADRATURE */
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
/*                       Function subprogram defining the integrand */
/*                       FUNCTION F(X). The actual name for F needs to be */
/*                       Declared E X T E R N A L in the calling program. */

/*              A      - Real */
/*                       Lower limit of integration */

/*              B      - Real */
/*                       Upper limit of integration */

/*            ON RETURN */
/*              RESULT - Real */
/*                       Approximation to the integral I */
/*                       RESULT is computed by applying the 31-POINT */
/*                       GAUSS-KRONROD RULE (RESK), obtained by optimal */
/*                       addition of abscissae to the 15-POINT GAUSS */
/*                       RULE (RESG). */

/*              ABSERR - Real */
/*                       Estimate of the modulus of the modulus, */
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
/* ***END PROLOGUE  QK31 */


/*           THE ABSCISSAE AND WEIGHTS ARE GIVEN FOR THE INTERVAL (-1,1). */
/*           BECAUSE OF SYMMETRY ONLY THE POSITIVE ABSCISSAE AND THEIR */
/*           CORRESPONDING WEIGHTS ARE GIVEN. */

/*           XGK    - ABSCISSAE OF THE 31-POINT KRONROD RULE */
/*                    XGK(2), XGK(4), ...  ABSCISSAE OF THE 15-POINT */
/*                    GAUSS RULE */
/*                    XGK(1), XGK(3), ...  ABSCISSAE WHICH ARE OPTIMALLY */
/*                    ADDED TO THE 15-POINT GAUSS RULE */

/*           WGK    - WEIGHTS OF THE 31-POINT KRONROD RULE */

/*           WG     - WEIGHTS OF THE 15-POINT GAUSS RULE */



/*           LIST OF MAJOR VARIABLES */
/*           ----------------------- */
/*           CENTR  - MID POINT OF THE INTERVAL */
/*           HLGTH  - HALF-LENGTH OF THE INTERVAL */
/*           ABSC   - ABSCISSA */
/*           FVAL*  - FUNCTION VALUE */
/*           RESG   - RESULT OF THE 15-POINT GAUSS FORMULA */
/*           RESK   - RESULT OF THE 31-POINT KRONROD FORMULA */
/*           RESKH  - APPROXIMATION TO THE MEAN VALUE OF F OVER (A,B), */
/*                    I.E. TO I/(B-A) */

/*           MACHINE DEPENDENT CONSTANTS */
/*           --------------------------- */
/*           EPMACH IS THE LARGEST RELATIVE SPACING. */
/*           UFLOW IS THE SMALLEST POSITIVE MAGNITUDE. */

/* ***FIRST EXECUTABLE STATEMENT  QK31 */
    epmach = r1mach_(&c__4);
    uflow = r1mach_(&c__1);

    centr = (*a + *b) * .5f;
    hlgth = (*b - *a) * .5f;
    dhlgth = dabs(hlgth);

/*           COMPUTE THE 31-POINT KRONROD APPROXIMATION TO */
/*           THE INTEGRAL, AND ESTIMATE THE ABSOLUTE ERROR. */

    fc = (*f)(&centr);
    resg = wg[7] * fc;
    resk = wgk[15] * fc;
    *resabs = dabs(resk);
    for (j = 1; j <= 7; ++j) {
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
    for (j = 1; j <= 8; ++j) {
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
    *resasc = wgk[15] * (r__1 = fc - reskh, dabs(r__1));
    for (j = 1; j <= 15; ++j) {
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
} /* qk31_ */

