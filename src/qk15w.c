/* qk15w.f -- translated by f2c (version 12.02.01).
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

/* DECK QK15W */
/* Subroutine */ int qk15w_(E_fp f, E_fp w, real *p1, real *p2, real *p3, 
	real *p4, integer *kp, real *a, real *b, real *result, real *abserr, 
	real *resabs, real *resasc)
{
    /* Initialized data */

    static real xgk[8] = { .9914553711208126f,.9491079123427585f,
	    .8648644233597691f,.7415311855993944f,.5860872354676911f,
	    .4058451513773972f,.2077849550078985f,0.f };
    static real wgk[8] = { .02293532201052922f,.06309209262997855f,
	    .1047900103222502f,.1406532597155259f,.1690047266392679f,
	    .1903505780647854f,.2044329400752989f,.2094821410847278f };
    static real wg[4] = { .1294849661688697f,.2797053914892767f,
	    .3818300505051889f,.4179591836734694f };

    /* System generated locals */
    real r__1, r__2;
    doublereal d__1;

    /* Local variables */
    static integer j;
    static real fc, fv1[7], fv2[7];
    static integer jtw;
    static real absc, resg, resk, fsum, absc1, absc2, fval1, fval2;
    static integer jtwm1;
    static real hlgth, centr, reskh, uflow;
    extern doublereal r1mach_(integer *);
    static real epmach, dhlgth;

/* ***BEGIN PROLOGUE  QK15W */
/* ***PURPOSE  To compute I = Integral of F*W over (A,B), with error */
/*                           estimate */
/*                       J = Integral of ABS(F*W) over (A,B) */
/* ***LIBRARY   SLATEC (QUADPACK) */
/* ***CATEGORY  H2A2A2 */
/* ***TYPE      SINGLE PRECISION (QK15W-S, DQK15W-D) */
/* ***KEYWORDS  15-POINT GAUSS-KRONROD RULES, QUADPACK, QUADRATURE */
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
/*             ON ENTRY */
/*              F      - Real */
/*                       Function subprogram defining the integrand */
/*                       function F(X). The actual name for F needs to be */
/*                       declared E X T E R N A L in the driver program. */

/*              W      - Real */
/*                       Function subprogram defining the integrand */
/*                       WEIGHT function W(X). The actual name for W */
/*                       needs to be declared E X T E R N A L in the */
/*                       calling program. */

/*              P1, P2, P3, P4 - Real */
/*                       Parameters in the WEIGHT function */

/*              KP     - Integer */
/*                       Key for indicating the type of WEIGHT function */

/*              A      - Real */
/*                       Lower limit of integration */

/*              B      - Real */
/*                       Upper limit of integration */

/*            ON RETURN */
/*              RESULT - Real */
/*                       Approximation to the integral I */
/*                       RESULT is computed by applying the 15-point */
/*                       Kronrod rule (RESK) obtained by optimal addition */
/*                       of abscissae to the 7-point Gauss rule (RESG). */

/*              ABSERR - Real */
/*                       Estimate of the modulus of the absolute error, */
/*                       which should equal or exceed ABS(I-RESULT) */

/*              RESABS - Real */
/*                       Approximation to the integral of ABS(F) */

/*              RESASC - Real */
/*                       Approximation to the integral of ABS(F-I/(B-A)) */

/* ***REFERENCES  (NONE) */
/* ***ROUTINES CALLED  R1MACH */
/* ***REVISION HISTORY  (YYMMDD) */
/*   810101  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890531  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/* ***END PROLOGUE  QK15W */



/*           THE ABSCISSAE AND WEIGHTS ARE GIVEN FOR THE INTERVAL (-1,1). */
/*           BECAUSE OF SYMMETRY ONLY THE POSITIVE ABSCISSAE AND THEIR */
/*           CORRESPONDING WEIGHTS ARE GIVEN. */

/*           XGK    - ABSCISSAE OF THE 15-POINT GAUSS-KRONROD RULE */
/*                    XGK(2), XGK(4), ... ABSCISSAE OF THE 7-POINT */
/*                    GAUSS RULE */
/*                    XGK(1), XGK(3), ... ABSCISSAE WHICH ARE OPTIMALLY */
/*                    ADDED TO THE 7-POINT GAUSS RULE */

/*           WGK    - WEIGHTS OF THE 15-POINT GAUSS-KRONROD RULE */

/*           WG     - WEIGHTS OF THE 7-POINT GAUSS RULE */





/*           LIST OF MAJOR VARIABLES */
/*           ----------------------- */

/*           CENTR  - MID POINT OF THE INTERVAL */
/*           HLGTH  - HALF-LENGTH OF THE INTERVAL */
/*           ABSC*  - ABSCISSA */
/*           FVAL*  - FUNCTION VALUE */
/*           RESG   - RESULT OF THE 7-POINT GAUSS FORMULA */
/*           RESK   - RESULT OF THE 15-POINT KRONROD FORMULA */
/*           RESKH  - APPROXIMATION TO THE MEAN VALUE OF F*W OVER (A,B), */
/*                    I.E. TO I/(B-A) */

/*           MACHINE DEPENDENT CONSTANTS */
/*           --------------------------- */

/*           EPMACH IS THE LARGEST RELATIVE SPACING. */
/*           UFLOW IS THE SMALLEST POSITIVE MAGNITUDE. */

/* ***FIRST EXECUTABLE STATEMENT  QK15W */
    epmach = r1mach_(&c__4);
    uflow = r1mach_(&c__1);

    centr = (*a + *b) * .5f;
    hlgth = (*b - *a) * .5f;
    dhlgth = dabs(hlgth);

/*           COMPUTE THE 15-POINT KRONROD APPROXIMATION TO THE */
/*           INTEGRAL, AND ESTIMATE THE ERROR. */

    fc = (*f)(&centr) * (*w)(&centr, p1, p2, p3, p4, kp);
    resg = wg[3] * fc;
    resk = wgk[7] * fc;
    *resabs = dabs(resk);
    for (j = 1; j <= 3; ++j) {
	jtw = j << 1;
	absc = hlgth * xgk[jtw - 1];
	absc1 = centr - absc;
	absc2 = centr + absc;
	fval1 = (*f)(&absc1) * (*w)(&absc1, p1, p2, p3, p4, kp);
	fval2 = (*f)(&absc2) * (*w)(&absc2, p1, p2, p3, p4, kp);
	fv1[jtw - 1] = fval1;
	fv2[jtw - 1] = fval2;
	fsum = fval1 + fval2;
	resg += wg[j - 1] * fsum;
	resk += wgk[jtw - 1] * fsum;
	*resabs += wgk[jtw - 1] * (dabs(fval1) + dabs(fval2));
/* L10: */
    }
    for (j = 1; j <= 4; ++j) {
	jtwm1 = (j << 1) - 1;
	absc = hlgth * xgk[jtwm1 - 1];
	absc1 = centr - absc;
	absc2 = centr + absc;
	fval1 = (*f)(&absc1) * (*w)(&absc1, p1, p2, p3, p4, kp);
	fval2 = (*f)(&absc2) * (*w)(&absc2, p1, p2, p3, p4, kp);
	fv1[jtwm1 - 1] = fval1;
	fv2[jtwm1 - 1] = fval2;
	fsum = fval1 + fval2;
	resk += wgk[jtwm1 - 1] * fsum;
	*resabs += wgk[jtwm1 - 1] * (dabs(fval1) + dabs(fval2));
/* L15: */
    }
    reskh = resk * .5f;
    *resasc = wgk[7] * (r__1 = fc - reskh, dabs(r__1));
    for (j = 1; j <= 7; ++j) {
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
} /* qk15w_ */

