/* dqk15w.f -- translated by f2c (version 12.02.01).
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

/* DECK DQK15W */
/* Subroutine */ int dqk15w_(D_fp f, D_fp w, doublereal *p1, doublereal *p2, 
	doublereal *p3, doublereal *p4, integer *kp, doublereal *a, 
	doublereal *b, doublereal *result, doublereal *abserr, doublereal *
	resabs, doublereal *resasc)
{
    /* Initialized data */

    static doublereal xgk[8] = { .9914553711208126,.9491079123427585,
	    .8648644233597691,.7415311855993944,.5860872354676911,
	    .4058451513773972,.2077849550078985,0. };
    static doublereal wgk[8] = { .02293532201052922,.06309209262997855,
	    .1047900103222502,.1406532597155259,.1690047266392679,
	    .1903505780647854,.2044329400752989,.2094821410847278 };
    static doublereal wg[4] = { .1294849661688697,.2797053914892767,
	    .3818300505051889,.4179591836734694 };

    /* System generated locals */
    doublereal d__1, d__2, d__3;

    /* Local variables */
    static integer j;
    static doublereal fc, fv1[7], fv2[7];
    static integer jtw;
    static doublereal absc, resg, resk, fsum, absc1, absc2, fval1, fval2;
    static integer jtwm1;
    static doublereal hlgth, centr, reskh, uflow;
    extern doublereal d1mach_(integer *);
    static doublereal epmach, dhlgth;

/* ***BEGIN PROLOGUE  DQK15W */
/* ***PURPOSE  To compute I = Integral of F*W over (A,B), with error */
/*                           estimate */
/*                       J = Integral of ABS(F*W) over (A,B) */
/* ***LIBRARY   SLATEC (QUADPACK) */
/* ***CATEGORY  H2A2A2 */
/* ***TYPE      DOUBLE PRECISION (QK15W-S, DQK15W-D) */
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
/*           Double precision version */

/*           PARAMETERS */
/*             ON ENTRY */
/*              F      - Double precision */
/*                       Function subprogram defining the integrand */
/*                       function F(X). The actual name for F needs to be */
/*                       declared E X T E R N A L in the driver program. */

/*              W      - Double precision */
/*                       Function subprogram defining the integrand */
/*                       WEIGHT function W(X). The actual name for W */
/*                       needs to be declared E X T E R N A L in the */
/*                       calling program. */

/*              P1, P2, P3, P4 - Double precision */
/*                       Parameters in the WEIGHT function */

/*              KP     - Integer */
/*                       Key for indicating the type of WEIGHT function */

/*              A      - Double precision */
/*                       Lower limit of integration */

/*              B      - Double precision */
/*                       Upper limit of integration */

/*            ON RETURN */
/*              RESULT - Double precision */
/*                       Approximation to the integral I */
/*                       RESULT is computed by applying the 15-point */
/*                       Kronrod rule (RESK) obtained by optimal addition */
/*                       of abscissae to the 7-point Gauss rule (RESG). */

/*              ABSERR - Double precision */
/*                       Estimate of the modulus of the absolute error, */
/*                       which should equal or exceed ABS(I-RESULT) */

/*              RESABS - Double precision */
/*                       Approximation to the integral of ABS(F) */

/*              RESASC - Double precision */
/*                       Approximation to the integral of ABS(F-I/(B-A)) */

/* ***REFERENCES  (NONE) */
/* ***ROUTINES CALLED  D1MACH */
/* ***REVISION HISTORY  (YYMMDD) */
/*   810101  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890531  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/* ***END PROLOGUE  DQK15W */



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

/* ***FIRST EXECUTABLE STATEMENT  DQK15W */
    epmach = d1mach_(&c__4);
    uflow = d1mach_(&c__1);

    centr = (*a + *b) * .5;
    hlgth = (*b - *a) * .5;
    dhlgth = abs(hlgth);

/*           COMPUTE THE 15-POINT KRONROD APPROXIMATION TO THE */
/*           INTEGRAL, AND ESTIMATE THE ERROR. */

    fc = (*f)(&centr) * (*w)(&centr, p1, p2, p3, p4, kp);
    resg = wg[3] * fc;
    resk = wgk[7] * fc;
    *resabs = abs(resk);
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
	*resabs += wgk[jtw - 1] * (abs(fval1) + abs(fval2));
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
	*resabs += wgk[jtwm1 - 1] * (abs(fval1) + abs(fval2));
/* L15: */
    }
    reskh = resk * .5;
    *resasc = wgk[7] * (d__1 = fc - reskh, abs(d__1));
    for (j = 1; j <= 7; ++j) {
	*resasc += wgk[j - 1] * ((d__1 = fv1[j - 1] - reskh, abs(d__1)) + (
		d__2 = fv2[j - 1] - reskh, abs(d__2)));
/* L20: */
    }
    *result = resk * hlgth;
    *resabs *= dhlgth;
    *resasc *= dhlgth;
    *abserr = (d__1 = (resk - resg) * hlgth, abs(d__1));
    if (*resasc != 0. && *abserr != 0.) {
/* Computing MIN */
	d__3 = *abserr * 200. / *resasc;
	d__1 = 1., d__2 = pow_dd(&d__3, &c_b7);
	*abserr = *resasc * min(d__1,d__2);
    }
    if (*resabs > uflow / (epmach * 50.)) {
/* Computing MAX */
	d__1 = epmach * 50. * *resabs;
	*abserr = max(d__1,*abserr);
    }
    return 0;
} /* dqk15w_ */

