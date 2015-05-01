/* dqc25c.f -- translated by f2c (version 12.02.01).
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

/* DECK DQC25C */
/* Subroutine */ int dqc25c_(D_fp f, doublereal *a, doublereal *b, doublereal 
	*c__, doublereal *result, doublereal *abserr, integer *krul, integer *
	neval)
{
    /* Initialized data */

    static doublereal x[11] = { .9914448613738104,.9659258262890683,
	    .9238795325112868,.8660254037844386,.7933533402912352,
	    .7071067811865475,.6087614290087206,.5,.3826834323650898,
	    .2588190451025208,.1305261922200516 };

    /* System generated locals */
    doublereal d__1;

    /* Local variables */
    static integer i__, k;
    static doublereal u, p2, p3, p4, cc;
    static integer kp;
    static doublereal ak22, fval[25], res12, res24;
    static integer isym;
    static doublereal amom0, amom1, amom2, cheb12[13], cheb24[25], hlgth, 
	    centr;
    extern /* Subroutine */ int dqk15w_(D_fp, D_fp, doublereal *, doublereal *
	    , doublereal *, doublereal *, integer *, doublereal *, doublereal 
	    *, doublereal *, doublereal *, doublereal *, doublereal *), 
	    dqcheb_(doublereal *, doublereal *, doublereal *, doublereal *);
    static doublereal resabs, resasc;
    extern doublereal dqwgtc_();

/* ***BEGIN PROLOGUE  DQC25C */
/* ***PURPOSE  To compute I = Integral of F*W over (A,B) with */
/*            error estimate, where W(X) = 1/(X-C) */
/* ***LIBRARY   SLATEC (QUADPACK) */
/* ***CATEGORY  H2A2A2, J4 */
/* ***TYPE      DOUBLE PRECISION (QC25C-S, DQC25C-D) */
/* ***KEYWORDS  25-POINT CLENSHAW-CURTIS INTEGRATION, QUADPACK, QUADRATURE */
/* ***AUTHOR  Piessens, Robert */
/*             Applied Mathematics and Programming Division */
/*             K. U. Leuven */
/*           de Doncker, Elise */
/*             Applied Mathematics and Programming Division */
/*             K. U. Leuven */
/* ***DESCRIPTION */

/*        Integration rules for the computation of CAUCHY */
/*        PRINCIPAL VALUE integrals */
/*        Standard fortran subroutine */
/*        Double precision version */

/*        PARAMETERS */
/*           F      - Double precision */
/*                    Function subprogram defining the integrand function */
/*                    F(X). The actual name for F needs to be declared */
/*                    E X T E R N A L  in the driver program. */

/*           A      - Double precision */
/*                    Left end point of the integration interval */

/*           B      - Double precision */
/*                    Right end point of the integration interval, B.GT.A */

/*           C      - Double precision */
/*                    Parameter in the WEIGHT function */

/*           RESULT - Double precision */
/*                    Approximation to the integral */
/*                    result is computed by using a generalized */
/*                    Clenshaw-Curtis method if C lies within ten percent */
/*                    of the integration interval. In the other case the */
/*                    15-point Kronrod rule obtained by optimal addition */
/*                    of abscissae to the 7-point Gauss rule, is applied. */

/*           ABSERR - Double precision */
/*                    Estimate of the modulus of the absolute error, */
/*                    which should equal or exceed ABS(I-RESULT) */

/*           KRUL   - Integer */
/*                    Key which is decreased by 1 if the 15-point */
/*                    Gauss-Kronrod scheme has been used */

/*           NEVAL  - Integer */
/*                    Number of integrand evaluations */

/* ...................................................................... */

/* ***REFERENCES  (NONE) */
/* ***ROUTINES CALLED  DQCHEB, DQK15W, DQWGTC */
/* ***REVISION HISTORY  (YYMMDD) */
/*   810101  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890531  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/* ***END PROLOGUE  DQC25C */




/*           THE VECTOR X CONTAINS THE VALUES COS(K*PI/24), */
/*           K = 1, ..., 11, TO BE USED FOR THE CHEBYSHEV SERIES */
/*           EXPANSION OF F */


/*           LIST OF MAJOR VARIABLES */
/*           ---------------------- */
/*           FVAL   - VALUE OF THE FUNCTION F AT THE POINTS */
/*                    COS(K*PI/24),  K = 0, ..., 24 */
/*           CHEB12 - CHEBYSHEV SERIES EXPANSION COEFFICIENTS, */
/*                    FOR THE FUNCTION F, OF DEGREE 12 */
/*           CHEB24 - CHEBYSHEV SERIES EXPANSION COEFFICIENTS, */
/*                    FOR THE FUNCTION F, OF DEGREE 24 */
/*           RES12  - APPROXIMATION TO THE INTEGRAL CORRESPONDING */
/*                    TO THE USE OF CHEB12 */
/*           RES24  - APPROXIMATION TO THE INTEGRAL CORRESPONDING */
/*                    TO THE USE OF CHEB24 */
/*           DQWGTC - EXTERNAL FUNCTION SUBPROGRAM DEFINING */
/*                    THE WEIGHT FUNCTION */
/*           HLGTH  - HALF-LENGTH OF THE INTERVAL */
/*           CENTR  - MID POINT OF THE INTERVAL */


/*           CHECK THE POSITION OF C. */

/* ***FIRST EXECUTABLE STATEMENT  DQC25C */
    cc = (2. * *c__ - *b - *a) / (*b - *a);
    if (abs(cc) < 1.1) {
	goto L10;
    }

/*           APPLY THE 15-POINT GAUSS-KRONROD SCHEME. */

    --(*krul);
    dqk15w_((D_fp)f, (D_fp)dqwgtc_, c__, &p2, &p3, &p4, &kp, a, b, result, 
	    abserr, &resabs, &resasc);
    *neval = 15;
    if (resasc == *abserr) {
	++(*krul);
    }
    goto L50;

/*           USE THE GENERALIZED CLENSHAW-CURTIS METHOD. */

L10:
    hlgth = (*b - *a) * .5;
    centr = (*b + *a) * .5;
    *neval = 25;
    d__1 = hlgth + centr;
    fval[0] = (*f)(&d__1) * .5;
    fval[12] = (*f)(&centr);
    d__1 = centr - hlgth;
    fval[24] = (*f)(&d__1) * .5;
    for (i__ = 2; i__ <= 12; ++i__) {
	u = hlgth * x[i__ - 2];
	isym = 26 - i__;
	d__1 = u + centr;
	fval[i__ - 1] = (*f)(&d__1);
	d__1 = centr - u;
	fval[isym - 1] = (*f)(&d__1);
/* L20: */
    }

/*           COMPUTE THE CHEBYSHEV SERIES EXPANSION. */

    dqcheb_(x, fval, cheb12, cheb24);

/*           THE MODIFIED CHEBYSHEV MOMENTS ARE COMPUTED BY FORWARD */
/*           RECURSION, USING AMOM0 AND AMOM1 AS STARTING VALUES. */

    amom0 = log((d__1 = (1. - cc) / (cc + 1.), abs(d__1)));
    amom1 = cc * amom0 + 2.;
    res12 = cheb12[0] * amom0 + cheb12[1] * amom1;
    res24 = cheb24[0] * amom0 + cheb24[1] * amom1;
    for (k = 3; k <= 13; ++k) {
	amom2 = cc * 2. * amom1 - amom0;
	ak22 = (doublereal) ((k - 2) * (k - 2));
	if (k / 2 << 1 == k) {
	    amom2 -= 4. / (ak22 - 1.);
	}
	res12 += cheb12[k - 1] * amom2;
	res24 += cheb24[k - 1] * amom2;
	amom0 = amom1;
	amom1 = amom2;
/* L30: */
    }
    for (k = 14; k <= 25; ++k) {
	amom2 = cc * 2. * amom1 - amom0;
	ak22 = (doublereal) ((k - 2) * (k - 2));
	if (k / 2 << 1 == k) {
	    amom2 -= 4. / (ak22 - 1.);
	}
	res24 += cheb24[k - 1] * amom2;
	amom0 = amom1;
	amom1 = amom2;
/* L40: */
    }
    *result = res24;
    *abserr = (d__1 = res24 - res12, abs(d__1));
L50:
    return 0;
} /* dqc25c_ */

