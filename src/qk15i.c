/* qk15i.f -- translated by f2c (version 12.02.01).
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
static doublereal c_b6 = 1.5;

/* DECK QK15I */
/* Subroutine */ int qk15i_(E_fp f, real *boun, integer *inf, real *a, real *
	b, real *result, real *abserr, real *resabs, real *resasc)
{
    /* Initialized data */

    static real xgk[8] = { .9914553711208126f,.9491079123427585f,
	    .8648644233597691f,.7415311855993944f,.5860872354676911f,
	    .4058451513773972f,.2077849550078985f,0.f };
    static real wgk[8] = { .02293532201052922f,.06309209262997855f,
	    .1047900103222502f,.1406532597155259f,.1690047266392679f,
	    .1903505780647854f,.2044329400752989f,.2094821410847278f };
    static real wg[8] = { 0.f,.1294849661688697f,0.f,.2797053914892767f,0.f,
	    .3818300505051189f,0.f,.4179591836734694f };

    /* System generated locals */
    real r__1, r__2;
    doublereal d__1;

    /* Local variables */
    static integer j;
    static real fc, fv1[7], fv2[7], absc, dinf, resg, resk, fsum, absc1, 
	    absc2, fval1, fval2, hlgth, centr, reskh, uflow;
    extern doublereal r1mach_(integer *);
    static real tabsc1, tabsc2, epmach;

/* ***BEGIN PROLOGUE  QK15I */
/* ***PURPOSE  The original (infinite integration range is mapped */
/*            onto the interval (0,1) and (A,B) is a part of (0,1). */
/*            it is the purpose to compute */
/*            I = Integral of transformed integrand over (A,B), */
/*            J = Integral of ABS(Transformed Integrand) over (A,B). */
/* ***LIBRARY   SLATEC (QUADPACK) */
/* ***CATEGORY  H2A3A2, H2A4A2 */
/* ***TYPE      SINGLE PRECISION (QK15I-S, DQK15I-D) */
/* ***KEYWORDS  15-POINT GAUSS-KRONROD RULES, QUADPACK, QUADRATURE */
/* ***AUTHOR  Piessens, Robert */
/*             Applied Mathematics and Programming Division */
/*             K. U. Leuven */
/*           de Doncker, Elise */
/*             Applied Mathematics and Programming Division */
/*             K. U. Leuven */
/* ***DESCRIPTION */

/*           Integration Rule */
/*           Standard Fortran subroutine */
/*           Real version */

/*           PARAMETERS */
/*            ON ENTRY */
/*              F      - Real */
/*                       Function subprogram defining the integrand */
/*                       FUNCTION F(X). The actual name for F needs to be */
/*                       Declared E X T E R N A L in the calling program. */

/*              BOUN   - Real */
/*                       Finite bound of original integration */
/*                       Range (SET TO ZERO IF INF = +2) */

/*              INF    - Integer */
/*                       If INF = -1, the original interval is */
/*                                   (-INFINITY,BOUND), */
/*                       If INF = +1, the original interval is */
/*                                   (BOUND,+INFINITY), */
/*                       If INF = +2, the original interval is */
/*                                   (-INFINITY,+INFINITY) AND */
/*                       The integral is computed as the sum of two */
/*                       integrals, one over (-INFINITY,0) and one over */
/*                       (0,+INFINITY). */

/*              A      - Real */
/*                       Lower limit for integration over subrange */
/*                       of (0,1) */

/*              B      - Real */
/*                       Upper limit for integration over subrange */
/*                       of (0,1) */

/*            ON RETURN */
/*              RESULT - Real */
/*                       Approximation to the integral I */
/*                       Result is computed by applying the 15-POINT */
/*                       KRONROD RULE(RESK) obtained by optimal addition */
/*                       of abscissae to the 7-POINT GAUSS RULE(RESG). */

/*              ABSERR - Real */
/*                       Estimate of the modulus of the absolute error, */
/*                       WHICH SHOULD EQUAL or EXCEED ABS(I-RESULT) */

/*              RESABS - Real */
/*                       Approximation to the integral J */

/*              RESASC - Real */
/*                       Approximation to the integral of */
/*                       ABS((TRANSFORMED INTEGRAND)-I/(B-A)) over (A,B) */

/* ***REFERENCES  (NONE) */
/* ***ROUTINES CALLED  R1MACH */
/* ***REVISION HISTORY  (YYMMDD) */
/*   800101  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890531  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/* ***END PROLOGUE  QK15I */



/*           THE ABSCISSAE AND WEIGHTS ARE SUPPLIED FOR THE INTERVAL */
/*           (-1,1).  BECAUSE OF SYMMETRY ONLY THE POSITIVE ABSCISSAE AND */
/*           THEIR CORRESPONDING WEIGHTS ARE GIVEN. */

/*           XGK    - ABSCISSAE OF THE 15-POINT KRONROD RULE */
/*                    XGK(2), XGK(4), ... ABSCISSAE OF THE 7-POINT */
/*                    GAUSS RULE */
/*                    XGK(1), XGK(3), ...  ABSCISSAE WHICH ARE OPTIMALLY */
/*                    ADDED TO THE 7-POINT GAUSS RULE */

/*           WGK    - WEIGHTS OF THE 15-POINT KRONROD RULE */

/*           WG     - WEIGHTS OF THE 7-POINT GAUSS RULE, CORRESPONDING */
/*                    TO THE ABSCISSAE XGK(2), XGK(4), ... */
/*                    WG(1), WG(3), ... ARE SET TO ZERO. */





/*           LIST OF MAJOR VARIABLES */
/*           ----------------------- */

/*           CENTR  - MID POINT OF THE INTERVAL */
/*           HLGTH  - HALF-LENGTH OF THE INTERVAL */
/*           ABSC*  - ABSCISSA */
/*           TABSC* - TRANSFORMED ABSCISSA */
/*           FVAL*  - FUNCTION VALUE */
/*           RESG   - RESULT OF THE 7-POINT GAUSS FORMULA */
/*           RESK   - RESULT OF THE 15-POINT KRONROD FORMULA */
/*           RESKH  - APPROXIMATION TO THE MEAN VALUE OF THE TRANSFORMED */
/*                    INTEGRAND OVER (A,B), I.E. TO I/(B-A) */

/*           MACHINE DEPENDENT CONSTANTS */
/*           --------------------------- */

/*           EPMACH IS THE LARGEST RELATIVE SPACING. */
/*           UFLOW IS THE SMALLEST POSITIVE MAGNITUDE. */

/* ***FIRST EXECUTABLE STATEMENT  QK15I */
    epmach = r1mach_(&c__4);
    uflow = r1mach_(&c__1);
    dinf = (real) min(1,*inf);

    centr = (*a + *b) * .5f;
    hlgth = (*b - *a) * .5f;
    tabsc1 = *boun + dinf * (1.f - centr) / centr;
    fval1 = (*f)(&tabsc1);
    if (*inf == 2) {
	r__1 = -tabsc1;
	fval1 += (*f)(&r__1);
    }
    fc = fval1 / centr / centr;

/*           COMPUTE THE 15-POINT KRONROD APPROXIMATION TO */
/*           THE INTEGRAL, AND ESTIMATE THE ERROR. */

    resg = wg[7] * fc;
    resk = wgk[7] * fc;
    *resabs = dabs(resk);
    for (j = 1; j <= 7; ++j) {
	absc = hlgth * xgk[j - 1];
	absc1 = centr - absc;
	absc2 = centr + absc;
	tabsc1 = *boun + dinf * (1.f - absc1) / absc1;
	tabsc2 = *boun + dinf * (1.f - absc2) / absc2;
	fval1 = (*f)(&tabsc1);
	fval2 = (*f)(&tabsc2);
	if (*inf == 2) {
	    r__1 = -tabsc1;
	    fval1 += (*f)(&r__1);
	}
	if (*inf == 2) {
	    r__1 = -tabsc2;
	    fval2 += (*f)(&r__1);
	}
	fval1 = fval1 / absc1 / absc1;
	fval2 = fval2 / absc2 / absc2;
	fv1[j - 1] = fval1;
	fv2[j - 1] = fval2;
	fsum = fval1 + fval2;
	resg += wg[j - 1] * fsum;
	resk += wgk[j - 1] * fsum;
	*resabs += wgk[j - 1] * (dabs(fval1) + dabs(fval2));
/* L10: */
    }
    reskh = resk * .5f;
    *resasc = wgk[7] * (r__1 = fc - reskh, dabs(r__1));
    for (j = 1; j <= 7; ++j) {
	*resasc += wgk[j - 1] * ((r__1 = fv1[j - 1] - reskh, dabs(r__1)) + (
		r__2 = fv2[j - 1] - reskh, dabs(r__2)));
/* L20: */
    }
    *result = resk * hlgth;
    *resasc *= hlgth;
    *resabs *= hlgth;
    *abserr = (r__1 = (resk - resg) * hlgth, dabs(r__1));
    if (*resasc != 0.f && *abserr != 0.f) {
/* Computing MIN */
	d__1 = (doublereal) (*abserr * 200.f / *resasc);
	r__1 = 1.f, r__2 = pow_dd(&d__1, &c_b6);
	*abserr = *resasc * dmin(r__1,r__2);
    }
    if (*resabs > uflow / (epmach * 50.f)) {
/* Computing MAX */
	r__1 = epmach * 50.f * *resabs;
	*abserr = dmax(r__1,*abserr);
    }
    return 0;
} /* qk15i_ */

