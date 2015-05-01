/* qc25s.f -- translated by f2c (version 12.02.01).
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

/* DECK QC25S */
/* Subroutine */ int qc25s_(E_fp f, real *a, real *b, real *bl, real *br, 
	real *alfa, real *beta, real *ri, real *rj, real *rg, real *rh, real *
	result, real *abserr, real *resasc, integer *integr, integer *nev)
{
    /* Initialized data */

    static real x[11] = { .9914448613738104f,.9659258262890683f,
	    .9238795325112868f,.8660254037844386f,.7933533402912352f,
	    .7071067811865475f,.6087614290087206f,.5f,.3826834323650898f,
	    .2588190451025208f,.1305261922200516f };

    /* System generated locals */
    real r__1;
    doublereal d__1, d__2;

    /* Local variables */
    static integer i__;
    static real u, dc, fix, fval[25], res12, res24;
    extern /* Subroutine */ int qk15w_(E_fp, E_fp, real *, real *, real *, 
	    real *, integer *, real *, real *, real *, real *, real *, real *)
	    ;
    static integer isym;
    static real cheb12[13], cheb24[25];
    extern /* Subroutine */ int qcheb_(real *, real *, real *, real *);
    static real hlgth, centr;
    extern doublereal qwgts_();
    static real factor, resabs;

/* ***BEGIN PROLOGUE  QC25S */
/* ***PURPOSE  To compute I = Integral of F*W over (BL,BR), with error */
/*            estimate, where the weight function W has a singular */
/*            behaviour of ALGEBRAICO-LOGARITHMIC type at the points */
/*            A and/or B. (BL,BR) is a part of (A,B). */
/* ***LIBRARY   SLATEC (QUADPACK) */
/* ***CATEGORY  H2A2A2 */
/* ***TYPE      SINGLE PRECISION (QC25S-S, DQC25S-D) */
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
/*        Real version */

/*        PARAMETERS */
/*           F      - Real */
/*                    Function subprogram defining the integrand */
/*                    F(X). The actual name for F needs to be declared */
/*                    E X T E R N A L  in the driver program. */

/*           A      - Real */
/*                    Left end point of the original interval */

/*           B      - Real */
/*                    Right end point of the original interval, B.GT.A */

/*           BL     - Real */
/*                    Lower limit of integration, BL.GE.A */

/*           BR     - Real */
/*                    Upper limit of integration, BR.LE.B */

/*           ALFA   - Real */
/*                    PARAMETER IN THE WEIGHT FUNCTION */

/*           BETA   - Real */
/*                    Parameter in the weight function */

/*           RI,RJ,RG,RH - Real */
/*                    Modified CHEBYSHEV moments for the application */
/*                    of the generalized CLENSHAW-CURTIS */
/*                    method (computed in subroutine DQMOMO) */

/*           RESULT - Real */
/*                    Approximation to the integral */
/*                    RESULT is computed by using a generalized */
/*                    CLENSHAW-CURTIS method if B1 = A or BR = B. */
/*                    in all other cases the 15-POINT KRONROD */
/*                    RULE is applied, obtained by optimal addition of */
/*                    Abscissae to the 7-POINT GAUSS RULE. */

/*           ABSERR - Real */
/*                    Estimate of the modulus of the absolute error, */
/*                    which should equal or exceed ABS(I-RESULT) */

/*           RESASC - Real */
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
/* ***ROUTINES CALLED  QCHEB, QK15W, QWGTS */
/* ***REVISION HISTORY  (YYMMDD) */
/*   810101  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890531  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/* ***END PROLOGUE  QC25S */




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
/*           QWGTS - EXTERNAL FUNCTION SUBPROGRAM DEFINING */
/*                    THE FOUR POSSIBLE WEIGHT FUNCTIONS */
/*           HLGTH  - HALF-LENGTH OF THE INTERVAL (BL,BR) */
/*           CENTR  - MID POINT OF THE INTERVAL (BL,BR) */

/* ***FIRST EXECUTABLE STATEMENT  QC25S */
    *nev = 25;
    if (*bl == *a && (*alfa != 0.f || *integr == 2 || *integr == 4)) {
	goto L10;
    }
    if (*br == *b && (*beta != 0.f || *integr == 3 || *integr == 4)) {
	goto L140;
    }

/*           IF A.GT.BL AND B.LT.BR, APPLY THE 15-POINT GAUSS-KRONROD */
/*           SCHEME. */


    qk15w_((E_fp)f, (E_fp)qwgts_, a, b, alfa, beta, integr, bl, br, result, 
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
    hlgth = (*br - *bl) * .5f;
    centr = (*br + *bl) * .5f;
    fix = *b - centr;
    r__1 = hlgth + centr;
    d__1 = (doublereal) (fix - hlgth);
    d__2 = (doublereal) (*beta);
    fval[0] = (*f)(&r__1) * .5f * pow_dd(&d__1, &d__2);
    d__1 = (doublereal) fix;
    d__2 = (doublereal) (*beta);
    fval[12] = (*f)(&centr) * pow_dd(&d__1, &d__2);
    r__1 = centr - hlgth;
    d__1 = (doublereal) (fix + hlgth);
    d__2 = (doublereal) (*beta);
    fval[24] = (*f)(&r__1) * .5f * pow_dd(&d__1, &d__2);
    for (i__ = 2; i__ <= 12; ++i__) {
	u = hlgth * x[i__ - 2];
	isym = 26 - i__;
	r__1 = u + centr;
	d__1 = (doublereal) (fix - u);
	d__2 = (doublereal) (*beta);
	fval[i__ - 1] = (*f)(&r__1) * pow_dd(&d__1, &d__2);
	r__1 = centr - u;
	d__1 = (doublereal) (fix + u);
	d__2 = (doublereal) (*beta);
	fval[isym - 1] = (*f)(&r__1) * pow_dd(&d__1, &d__2);
/* L20: */
    }
    d__1 = (doublereal) hlgth;
    d__2 = (doublereal) (*alfa + 1.f);
    factor = pow_dd(&d__1, &d__2);
    *result = 0.f;
    *abserr = 0.f;
    res12 = 0.f;
    res24 = 0.f;
    if (*integr > 2) {
	goto L70;
    }
    qcheb_(x, fval, cheb12, cheb24);

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
    *abserr = (r__1 = (res24 - res12) * dc, dabs(r__1));
    res12 = 0.f;
    res24 = 0.f;
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
    qcheb_(x, fval, cheb12, cheb24);

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
    *abserr = (r__1 = (res24 - res12) * dc, dabs(r__1));
    res12 = 0.f;
    res24 = 0.f;
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
    *abserr = (*abserr + (r__1 = res24 - res12, dabs(r__1))) * factor;
    goto L270;

/*           THIS PART OF THE PROGRAM IS EXECUTED ONLY IF B = BR. */
/*           ---------------------------------------------------- */

/*           COMPUTE THE CHEBYSHEV SERIES EXPANSION OF THE */
/*           FOLLOWING FUNCTION */
/*           F2 = (0.5*(B+BL-A-A)+0.5*(B-BL)*X)**ALFA */
/*                *F(0.5*(B-BL)*X+0.5*(B+BL)) */

L140:
    hlgth = (*br - *bl) * .5f;
    centr = (*br + *bl) * .5f;
    fix = centr - *a;
    r__1 = hlgth + centr;
    d__1 = (doublereal) (fix + hlgth);
    d__2 = (doublereal) (*alfa);
    fval[0] = (*f)(&r__1) * .5f * pow_dd(&d__1, &d__2);
    d__1 = (doublereal) fix;
    d__2 = (doublereal) (*alfa);
    fval[12] = (*f)(&centr) * pow_dd(&d__1, &d__2);
    r__1 = centr - hlgth;
    d__1 = (doublereal) (fix - hlgth);
    d__2 = (doublereal) (*alfa);
    fval[24] = (*f)(&r__1) * .5f * pow_dd(&d__1, &d__2);
    for (i__ = 2; i__ <= 12; ++i__) {
	u = hlgth * x[i__ - 2];
	isym = 26 - i__;
	r__1 = u + centr;
	d__1 = (doublereal) (fix + u);
	d__2 = (doublereal) (*alfa);
	fval[i__ - 1] = (*f)(&r__1) * pow_dd(&d__1, &d__2);
	r__1 = centr - u;
	d__1 = (doublereal) (fix - u);
	d__2 = (doublereal) (*alfa);
	fval[isym - 1] = (*f)(&r__1) * pow_dd(&d__1, &d__2);
/* L150: */
    }
    d__1 = (doublereal) hlgth;
    d__2 = (doublereal) (*beta + 1.f);
    factor = pow_dd(&d__1, &d__2);
    *result = 0.f;
    *abserr = 0.f;
    res12 = 0.f;
    res24 = 0.f;
    if (*integr == 2 || *integr == 4) {
	goto L200;
    }

/*           INTEGR = 1  (OR 3) */

    qcheb_(x, fval, cheb12, cheb24);
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
    *abserr = (r__1 = (res24 - res12) * dc, dabs(r__1));
    res12 = 0.f;
    res24 = 0.f;
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
	fval[i__ - 1] *= log(fix + u);
	fval[isym - 1] *= log(fix - u);
/* L210: */
    }
    qcheb_(x, fval, cheb12, cheb24);

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
    *abserr = (r__1 = (res24 - res12) * dc, dabs(r__1));
    res12 = 0.f;
    res24 = 0.f;

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
    *abserr = (*abserr + (r__1 = res24 - res12, dabs(r__1))) * factor;
L270:
    return 0;
} /* qc25s_ */

