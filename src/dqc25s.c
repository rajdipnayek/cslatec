/* dqc25s.f -- translated by f2c (version 12.02.01).
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

/* DECK DQC25S */
/* Subroutine */ int dqc25s_(D_fp f, doublereal *a, doublereal *b, doublereal 
	*bl, doublereal *br, doublereal *alfa, doublereal *beta, doublereal *
	ri, doublereal *rj, doublereal *rg, doublereal *rh, doublereal *
	result, doublereal *abserr, doublereal *resasc, integer *integr, 
	integer *nev)
{
    /* Initialized data */

    static doublereal x[11] = { .9914448613738104,.9659258262890683,
	    .9238795325112868,.8660254037844386,.7933533402912352,
	    .7071067811865475,.6087614290087206,.5,.3826834323650898,
	    .2588190451025208,.1305261922200516 };

    /* System generated locals */
    doublereal d__1, d__2;

    /* Local variables */
    static integer i__;
    static doublereal u, dc, fix, fval[25], res12, res24;
    static integer isym;
    static doublereal cheb12[13], cheb24[25], hlgth, centr;
    extern /* Subroutine */ int dqk15w_(D_fp, D_fp, doublereal *, doublereal *
	    , doublereal *, doublereal *, integer *, doublereal *, doublereal 
	    *, doublereal *, doublereal *, doublereal *, doublereal *), 
	    dqcheb_(doublereal *, doublereal *, doublereal *, doublereal *);
    static doublereal factor, resabs;
    extern doublereal dqwgts_();

/* ***BEGIN PROLOGUE  DQC25S */
/* ***PURPOSE  To compute I = Integral of F*W over (BL,BR), with error */
/*            estimate, where the weight function W has a singular */
/*            behaviour of ALGEBRAICO-LOGARITHMIC type at the points */
/*            A and/or B. (BL,BR) is a part of (A,B). */
/* ***LIBRARY   SLATEC (QUADPACK) */
/* ***CATEGORY  H2A2A2 */
/* ***TYPE      DOUBLE PRECISION (QC25S-S, DQC25S-D) */
/* ***KEYWORDS  25-POINT CLENSHAW-CURTIS INTEGRATION, QUADPACK, QUADRATURE */
/* ***AUTHOR  Piessens, Robert */
/*             Applied Mathematics and Programming Division */
/*             K. U. Leuven */
/*           de Doncker, Elise */
/*             Applied Mathematics and Programming Division */
/*             K. U. Leuven */
/* ***DESCRIPTION */

/*        Integration rules for integrands having ALGEBRAICO-LOGARITHMIC */
/*        end point singularities */
/*        Standard fortran subroutine */
/*        Double precision version */

/*        PARAMETERS */
/*           F      - Double precision */
/*                    Function subprogram defining the integrand */
/*                    F(X). The actual name for F needs to be declared */
/*                    E X T E R N A L  in the driver program. */

/*           A      - Double precision */
/*                    Left end point of the original interval */

/*           B      - Double precision */
/*                    Right end point of the original interval, B.GT.A */

/*           BL     - Double precision */
/*                    Lower limit of integration, BL.GE.A */

/*           BR     - Double precision */
/*                    Upper limit of integration, BR.LE.B */

/*           ALFA   - Double precision */
/*                    PARAMETER IN THE WEIGHT FUNCTION */

/*           BETA   - Double precision */
/*                    Parameter in the weight function */

/*           RI,RJ,RG,RH - Double precision */
/*                    Modified CHEBYSHEV moments for the application */
/*                    of the generalized CLENSHAW-CURTIS */
/*                    method (computed in subroutine DQMOMO) */

/*           RESULT - Double precision */
/*                    Approximation to the integral */
/*                    RESULT is computed by using a generalized */
/*                    CLENSHAW-CURTIS method if B1 = A or BR = B. */
/*                    in all other cases the 15-POINT KRONROD */
/*                    RULE is applied, obtained by optimal addition of */
/*                    Abscissae to the 7-POINT GAUSS RULE. */

/*           ABSERR - Double precision */
/*                    Estimate of the modulus of the absolute error, */
/*                    which should equal or exceed ABS(I-RESULT) */

/*           RESASC - Double precision */
/*                    Approximation to the integral of ABS(F*W-I/(B-A)) */

/*           INTEGR - Integer */
/*                    Which determines the weight function */
/*                    = 1   W(X) = (X-A)**ALFA*(B-X)**BETA */
/*                    = 2   W(X) = (X-A)**ALFA*(B-X)**BETA*LOG(X-A) */
/*                    = 3   W(X) = (X-A)**ALFA*(B-X)**BETA*LOG(B-X) */
/*                    = 4   W(X) = (X-A)**ALFA*(B-X)**BETA*LOG(X-A)* */
/*                                 LOG(B-X) */

/*           NEV    - Integer */
/*                    Number of integrand evaluations */

/* ***REFERENCES  (NONE) */
/* ***ROUTINES CALLED  DQCHEB, DQK15W, DQWGTS */
/* ***REVISION HISTORY  (YYMMDD) */
/*   810101  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890531  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/* ***END PROLOGUE  DQC25S */




/*           THE VECTOR X CONTAINS THE VALUES COS(K*PI/24) */
/*           K = 1, ..., 11, TO BE USED FOR THE COMPUTATION OF THE */
/*           CHEBYSHEV SERIES EXPANSION OF F. */

    /* Parameter adjustments */
    --rh;
    --rg;
    --rj;
    --ri;

    /* Function Body */

/*           LIST OF MAJOR VARIABLES */
/*           ----------------------- */

/*           FVAL   - VALUE OF THE FUNCTION F AT THE POINTS */
/*                    (BR-BL)*0.5*COS(K*PI/24)+(BR+BL)*0.5 */
/*                    K = 0, ..., 24 */
/*           CHEB12 - COEFFICIENTS OF THE CHEBYSHEV SERIES EXPANSION */
/*                    OF DEGREE 12, FOR THE FUNCTION F, IN THE */
/*                    INTERVAL (BL,BR) */
/*           CHEB24 - COEFFICIENTS OF THE CHEBYSHEV SERIES EXPANSION */
/*                    OF DEGREE 24, FOR THE FUNCTION F, IN THE */
/*                    INTERVAL (BL,BR) */
/*           RES12  - APPROXIMATION TO THE INTEGRAL OBTAINED FROM CHEB12 */
/*           RES24  - APPROXIMATION TO THE INTEGRAL OBTAINED FROM CHEB24 */
/*           DQWGTS - EXTERNAL FUNCTION SUBPROGRAM DEFINING */
/*                    THE FOUR POSSIBLE WEIGHT FUNCTIONS */
/*           HLGTH  - HALF-LENGTH OF THE INTERVAL (BL,BR) */
/*           CENTR  - MID POINT OF THE INTERVAL (BL,BR) */

/* ***FIRST EXECUTABLE STATEMENT  DQC25S */
    *nev = 25;
    if (*bl == *a && (*alfa != 0. || *integr == 2 || *integr == 4)) {
	goto L10;
    }
    if (*br == *b && (*beta != 0. || *integr == 3 || *integr == 4)) {
	goto L140;
    }

/*           IF A.GT.BL AND B.LT.BR, APPLY THE 15-POINT GAUSS-KRONROD */
/*           SCHEME. */


    dqk15w_((D_fp)f, (D_fp)dqwgts_, a, b, alfa, beta, integr, bl, br, result, 
	    abserr, &resabs, resasc);
    *nev = 15;
    goto L270;

/*           THIS PART OF THE PROGRAM IS EXECUTED ONLY IF A = BL. */
/*           ---------------------------------------------------- */

/*           COMPUTE THE CHEBYSHEV SERIES EXPANSION OF THE */
/*           FOLLOWING FUNCTION */
/*           F1 = (0.5*(B+B-BR-A)-0.5*(BR-A)*X)**BETA */
/*                  *F(0.5*(BR-A)*X+0.5*(BR+A)) */

L10:
    hlgth = (*br - *bl) * .5;
    centr = (*br + *bl) * .5;
    fix = *b - centr;
    d__1 = hlgth + centr;
    d__2 = fix - hlgth;
    fval[0] = (*f)(&d__1) * .5 * pow_dd(&d__2, beta);
    fval[12] = (*f)(&centr) * pow_dd(&fix, beta);
    d__1 = centr - hlgth;
    d__2 = fix + hlgth;
    fval[24] = (*f)(&d__1) * .5 * pow_dd(&d__2, beta);
    for (i__ = 2; i__ <= 12; ++i__) {
	u = hlgth * x[i__ - 2];
	isym = 26 - i__;
	d__1 = u + centr;
	d__2 = fix - u;
	fval[i__ - 1] = (*f)(&d__1) * pow_dd(&d__2, beta);
	d__1 = centr - u;
	d__2 = fix + u;
	fval[isym - 1] = (*f)(&d__1) * pow_dd(&d__2, beta);
/* L20: */
    }
    d__1 = *alfa + 1.;
    factor = pow_dd(&hlgth, &d__1);
    *result = 0.;
    *abserr = 0.;
    res12 = 0.;
    res24 = 0.;
    if (*integr > 2) {
	goto L70;
    }
    dqcheb_(x, fval, cheb12, cheb24);

/*           INTEGR = 1  (OR 2) */

    for (i__ = 1; i__ <= 13; ++i__) {
	res12 += cheb12[i__ - 1] * ri[i__];
	res24 += cheb24[i__ - 1] * ri[i__];
/* L30: */
    }
    for (i__ = 14; i__ <= 25; ++i__) {
	res24 += cheb24[i__ - 1] * ri[i__];
/* L40: */
    }
    if (*integr == 1) {
	goto L130;
    }

/*           INTEGR = 2 */

    dc = log(*br - *bl);
    *result = res24 * dc;
    *abserr = (d__1 = (res24 - res12) * dc, abs(d__1));
    res12 = 0.;
    res24 = 0.;
    for (i__ = 1; i__ <= 13; ++i__) {
	res12 += cheb12[i__ - 1] * rg[i__];
	res24 = res12 + cheb24[i__ - 1] * rg[i__];
/* L50: */
    }
    for (i__ = 14; i__ <= 25; ++i__) {
	res24 += cheb24[i__ - 1] * rg[i__];
/* L60: */
    }
    goto L130;

/*           COMPUTE THE CHEBYSHEV SERIES EXPANSION OF THE */
/*           FOLLOWING FUNCTION */
/*           F4 = F1*LOG(0.5*(B+B-BR-A)-0.5*(BR-A)*X) */

L70:
    fval[0] *= log(fix - hlgth);
    fval[12] *= log(fix);
    fval[24] *= log(fix + hlgth);
    for (i__ = 2; i__ <= 12; ++i__) {
	u = hlgth * x[i__ - 2];
	isym = 26 - i__;
	fval[i__ - 1] *= log(fix - u);
	fval[isym - 1] *= log(fix + u);
/* L80: */
    }
    dqcheb_(x, fval, cheb12, cheb24);

/*           INTEGR = 3  (OR 4) */

    for (i__ = 1; i__ <= 13; ++i__) {
	res12 += cheb12[i__ - 1] * ri[i__];
	res24 += cheb24[i__ - 1] * ri[i__];
/* L90: */
    }
    for (i__ = 14; i__ <= 25; ++i__) {
	res24 += cheb24[i__ - 1] * ri[i__];
/* L100: */
    }
    if (*integr == 3) {
	goto L130;
    }

/*           INTEGR = 4 */

    dc = log(*br - *bl);
    *result = res24 * dc;
    *abserr = (d__1 = (res24 - res12) * dc, abs(d__1));
    res12 = 0.;
    res24 = 0.;
    for (i__ = 1; i__ <= 13; ++i__) {
	res12 += cheb12[i__ - 1] * rg[i__];
	res24 += cheb24[i__ - 1] * rg[i__];
/* L110: */
    }
    for (i__ = 14; i__ <= 25; ++i__) {
	res24 += cheb24[i__ - 1] * rg[i__];
/* L120: */
    }
L130:
    *result = (*result + res24) * factor;
    *abserr = (*abserr + (d__1 = res24 - res12, abs(d__1))) * factor;
    goto L270;

/*           THIS PART OF THE PROGRAM IS EXECUTED ONLY IF B = BR. */
/*           ---------------------------------------------------- */

/*           COMPUTE THE CHEBYSHEV SERIES EXPANSION OF THE */
/*           FOLLOWING FUNCTION */
/*           F2 = (0.5*(B+BL-A-A)+0.5*(B-BL)*X)**ALFA */
/*                *F(0.5*(B-BL)*X+0.5*(B+BL)) */

L140:
    hlgth = (*br - *bl) * .5;
    centr = (*br + *bl) * .5;
    fix = centr - *a;
    d__1 = hlgth + centr;
    d__2 = fix + hlgth;
    fval[0] = (*f)(&d__1) * .5 * pow_dd(&d__2, alfa);
    fval[12] = (*f)(&centr) * pow_dd(&fix, alfa);
    d__1 = centr - hlgth;
    d__2 = fix - hlgth;
    fval[24] = (*f)(&d__1) * .5 * pow_dd(&d__2, alfa);
    for (i__ = 2; i__ <= 12; ++i__) {
	u = hlgth * x[i__ - 2];
	isym = 26 - i__;
	d__1 = u + centr;
	d__2 = fix + u;
	fval[i__ - 1] = (*f)(&d__1) * pow_dd(&d__2, alfa);
	d__1 = centr - u;
	d__2 = fix - u;
	fval[isym - 1] = (*f)(&d__1) * pow_dd(&d__2, alfa);
/* L150: */
    }
    d__1 = *beta + 1.;
    factor = pow_dd(&hlgth, &d__1);
    *result = 0.;
    *abserr = 0.;
    res12 = 0.;
    res24 = 0.;
    if (*integr == 2 || *integr == 4) {
	goto L200;
    }

/*           INTEGR = 1  (OR 3) */

    dqcheb_(x, fval, cheb12, cheb24);
    for (i__ = 1; i__ <= 13; ++i__) {
	res12 += cheb12[i__ - 1] * rj[i__];
	res24 += cheb24[i__ - 1] * rj[i__];
/* L160: */
    }
    for (i__ = 14; i__ <= 25; ++i__) {
	res24 += cheb24[i__ - 1] * rj[i__];
/* L170: */
    }
    if (*integr == 1) {
	goto L260;
    }

/*           INTEGR = 3 */

    dc = log(*br - *bl);
    *result = res24 * dc;
    *abserr = (d__1 = (res24 - res12) * dc, abs(d__1));
    res12 = 0.;
    res24 = 0.;
    for (i__ = 1; i__ <= 13; ++i__) {
	res12 += cheb12[i__ - 1] * rh[i__];
	res24 += cheb24[i__ - 1] * rh[i__];
/* L180: */
    }
    for (i__ = 14; i__ <= 25; ++i__) {
	res24 += cheb24[i__ - 1] * rh[i__];
/* L190: */
    }
    goto L260;

/*           COMPUTE THE CHEBYSHEV SERIES EXPANSION OF THE */
/*           FOLLOWING FUNCTION */
/*           F3 = F2*LOG(0.5*(B-BL)*X+0.5*(B+BL-A-A)) */

L200:
    fval[0] *= log(fix + hlgth);
    fval[12] *= log(fix);
    fval[24] *= log(fix - hlgth);
    for (i__ = 2; i__ <= 12; ++i__) {
	u = hlgth * x[i__ - 2];
	isym = 26 - i__;
	fval[i__ - 1] *= log(u + fix);
	fval[isym - 1] *= log(fix - u);
/* L210: */
    }
    dqcheb_(x, fval, cheb12, cheb24);

/*           INTEGR = 2  (OR 4) */

    for (i__ = 1; i__ <= 13; ++i__) {
	res12 += cheb12[i__ - 1] * rj[i__];
	res24 += cheb24[i__ - 1] * rj[i__];
/* L220: */
    }
    for (i__ = 14; i__ <= 25; ++i__) {
	res24 += cheb24[i__ - 1] * rj[i__];
/* L230: */
    }
    if (*integr == 2) {
	goto L260;
    }
    dc = log(*br - *bl);
    *result = res24 * dc;
    *abserr = (d__1 = (res24 - res12) * dc, abs(d__1));
    res12 = 0.;
    res24 = 0.;

/*           INTEGR = 4 */

    for (i__ = 1; i__ <= 13; ++i__) {
	res12 += cheb12[i__ - 1] * rh[i__];
	res24 += cheb24[i__ - 1] * rh[i__];
/* L240: */
    }
    for (i__ = 14; i__ <= 25; ++i__) {
	res24 += cheb24[i__ - 1] * rh[i__];
/* L250: */
    }
L260:
    *result = (*result + res24) * factor;
    *abserr = (*abserr + (d__1 = res24 - res12, abs(d__1))) * factor;
L270:
    return 0;
} /* dqc25s_ */

