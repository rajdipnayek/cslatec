/* qk41.f -- translated by f2c (version 12.02.01).
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

/* DECK QK41 */
/* Subroutine */ int qk41_(E_fp f, real *a, real *b, real *result, real *
	abserr, real *resabs, real *resasc)
{
    /* Initialized data */

    static real xgk[21] = { .9988590315882777f,.9931285991850949f,
	    .9815078774502503f,.9639719272779138f,.9408226338317548f,
	    .9122344282513259f,.878276811252282f,.8391169718222188f,
	    .7950414288375512f,.7463319064601508f,.6932376563347514f,
	    .636053680726515f,.5751404468197103f,.5108670019508271f,
	    .4435931752387251f,.3737060887154196f,.301627868114913f,
	    .2277858511416451f,.1526054652409227f,.07652652113349733f,0.f };
    static real wgk[21] = { .003073583718520532f,.008600269855642942f,
	    .01462616925697125f,.02038837346126652f,.02588213360495116f,
	    .0312873067770328f,.0366001697582008f,.04166887332797369f,
	    .04643482186749767f,.05094457392372869f,.05519510534828599f,
	    .05911140088063957f,.06265323755478117f,.06583459713361842f,
	    .06864867292852162f,.07105442355344407f,.07303069033278667f,
	    .07458287540049919f,.07570449768455667f,.07637786767208074f,
	    .07660071191799966f };
    static real wg[10] = { .01761400713915212f,.04060142980038694f,
	    .06267204833410906f,.08327674157670475f,.1019301198172404f,
	    .1181945319615184f,.1316886384491766f,.1420961093183821f,
	    .1491729864726037f,.1527533871307259f };

    /* System generated locals */
    real r__1, r__2;
    doublereal d__1;

    /* Local variables */
    static integer j;
    static real fc, fv1[20], fv2[20];
    static integer jtw;
    static real absc, resg, resk, fsum, fval1, fval2;
    static integer jtwm1;
    static real hlgth, centr, reskh, uflow;
    extern doublereal r1mach_(integer *);
    static real epmach, dhlgth;

/* ***BEGIN PROLOGUE  QK41 */
/* ***PURPOSE  To compute I = Integral of F over (A,B), with error */
/*                           estimate */
/*                       J = Integral of ABS(F) over (A,B) */
/* ***LIBRARY   SLATEC (QUADPACK) */
/* ***CATEGORY  H2A1A2 */
/* ***TYPE      SINGLE PRECISION (QK41-S, DQK41-D) */
/* ***KEYWORDS  41-POINT GAUSS-KRONROD RULES, QUADPACK, QUADRATURE */
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
/*                       declared E X T E R N A L in the calling program. */

/*              A      - Real */
/*                       Lower limit of integration */

/*              B      - Real */
/*                       Upper limit of integration */

/*            ON RETURN */
/*              RESULT - Real */
/*                       Approximation to the integral I */
/*                       RESULT is computed by applying the 41-POINT */
/*                       GAUSS-KRONROD RULE (RESK) obtained by optimal */
/*                       addition of abscissae to the 20-POINT GAUSS */
/*                       RULE (RESG). */

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
/* ***END PROLOGUE  QK41 */



/*           THE ABSCISSAE AND WEIGHTS ARE GIVEN FOR THE INTERVAL (-1,1). */
/*           BECAUSE OF SYMMETRY ONLY THE POSITIVE ABSCISSAE AND THEIR */
/*           CORRESPONDING WEIGHTS ARE GIVEN. */

/*           XGK    - ABSCISSAE OF THE 41-POINT GAUSS-KRONROD RULE */
/*                    XGK(2), XGK(4), ...  ABSCISSAE OF THE 20-POINT */
/*                    GAUSS RULE */
/*                    XGK(1), XGK(3), ...  ABSCISSAE WHICH ARE OPTIMALLY */
/*                    ADDED TO THE 20-POINT GAUSS RULE */

/*           WGK    - WEIGHTS OF THE 41-POINT GAUSS-KRONROD RULE */

/*           WG     - WEIGHTS OF THE 20-POINT GAUSS RULE */



/*           LIST OF MAJOR VARIABLES */
/*           ----------------------- */

/*           CENTR  - MID POINT OF THE INTERVAL */
/*           HLGTH  - HALF-LENGTH OF THE INTERVAL */
/*           ABSC   - ABSCISSA */
/*           FVAL*  - FUNCTION VALUE */
/*           RESG   - RESULT OF THE 20-POINT GAUSS FORMULA */
/*           RESK   - RESULT OF THE 41-POINT KRONROD FORMULA */
/*           RESKH  - APPROXIMATION TO MEAN VALUE OF F OVER (A,B), I.E. */
/*                    TO I/(B-A) */

/*           MACHINE DEPENDENT CONSTANTS */
/*           --------------------------- */

/*           EPMACH IS THE LARGEST RELATIVE SPACING. */
/*           UFLOW IS THE SMALLEST POSITIVE MAGNITUDE. */

/* ***FIRST EXECUTABLE STATEMENT  QK41 */
    epmach = r1mach_(&c__4);
    uflow = r1mach_(&c__1);

    centr = (*a + *b) * .5f;
    hlgth = (*b - *a) * .5f;
    dhlgth = dabs(hlgth);

/*           COMPUTE THE 41-POINT GAUSS-KRONROD APPROXIMATION TO */
/*           THE INTEGRAL, AND ESTIMATE THE ABSOLUTE ERROR. */

    resg = 0.f;
    fc = (*f)(&centr);
    resk = wgk[20] * fc;
    *resabs = dabs(resk);
    for (j = 1; j <= 10; ++j) {
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
    for (j = 1; j <= 10; ++j) {
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
    *resasc = wgk[20] * (r__1 = fc - reskh, dabs(r__1));
    for (j = 1; j <= 20; ++j) {
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
} /* qk41_ */

