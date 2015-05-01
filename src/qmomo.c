/* qmomo.f -- translated by f2c (version 12.02.01).
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

static doublereal c_b2 = 2.;

/* DECK QMOMO */
/* Subroutine */ int qmomo_(real *alfa, real *beta, real *ri, real *rj, real *
	rg, real *rh, integer *integr)
{
    /* System generated locals */
    doublereal d__1;

    /* Local variables */
    static integer i__;
    static real an;
    static integer im1;
    static real anm1, ralf, rbet, alfp1, alfp2, betp1, betp2;

/* ***BEGIN PROLOGUE  QMOMO */
/* ***PURPOSE  This routine computes modified Chebyshev moments.  The K-th */
/*            modified Chebyshev moment is defined as the integral over */
/*            (-1,1) of W(X)*T(K,X), where T(K,X) is the Chebyshev */
/*            polynomial of degree K. */
/* ***LIBRARY   SLATEC (QUADPACK) */
/* ***CATEGORY  H2A2A1, C3A2 */
/* ***TYPE      SINGLE PRECISION (QMOMO-S, DQMOMO-D) */
/* ***KEYWORDS  MODIFIED CHEBYSHEV MOMENTS, QUADPACK, QUADRATURE */
/* ***AUTHOR  Piessens, Robert */
/*             Applied Mathematics and Programming Division */
/*             K. U. Leuven */
/*           de Doncker, Elise */
/*             Applied Mathematics and Programming Division */
/*             K. U. Leuven */
/* ***DESCRIPTION */

/*        MODIFIED CHEBYSHEV MOMENTS */
/*        STANDARD FORTRAN SUBROUTINE */
/*        REAL VERSION */

/*        PARAMETERS */
/*           ALFA   - Real */
/*                    Parameter in the weight function W(X), ALFA.GT.(-1) */

/*           BETA   - Real */
/*                    Parameter in the weight function W(X), BETA.GT.(-1) */

/*           RI     - Real */
/*                    Vector of dimension 25 */
/*                    RI(K) is the integral over (-1,1) of */
/*                    (1+X)**ALFA*T(K-1,X), K = 1, ..., 25. */

/*           RJ     - Real */
/*                    Vector of dimension 25 */
/*                    RJ(K) is the integral over (-1,1) of */
/*                    (1-X)**BETA*T(K-1,X), K = 1, ..., 25. */

/*           RG     - Real */
/*                    Vector of dimension 25 */
/*                    RG(K) is the integral over (-1,1) of */
/*                    (1+X)**ALFA*LOG((1+X)/2)*T(K-1,X), K = 1, ..., 25. */

/*           RH     - Real */
/*                    Vector of dimension 25 */
/*                    RH(K) is the integral over (-1,1) of */
/*                    (1-X)**BETA*LOG((1-X)/2)*T(K-1,X), K = 1, ..., 25. */

/*           INTEGR - Integer */
/*                    Input parameter indicating the modified */
/*                    Moments to be computed */
/*                    INTEGR = 1 compute RI, RJ */
/*                           = 2 compute RI, RJ, RG */
/*                           = 3 compute RI, RJ, RH */
/*                           = 4 compute RI, RJ, RG, RH */

/* ***REFERENCES  (NONE) */
/* ***ROUTINES CALLED  (NONE) */
/* ***REVISION HISTORY  (YYMMDD) */
/*   810101  DATE WRITTEN */
/*   891009  Removed unreferenced statement label.  (WRB) */
/*   891009  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/* ***END PROLOGUE  QMOMO */




/* ***FIRST EXECUTABLE STATEMENT  QMOMO */
    /* Parameter adjustments */
    --rh;
    --rg;
    --rj;
    --ri;

    /* Function Body */
    alfp1 = *alfa + 1.f;
    betp1 = *beta + 1.f;
    alfp2 = *alfa + 2.f;
    betp2 = *beta + 2.f;
    d__1 = (doublereal) alfp1;
    ralf = pow_dd(&c_b2, &d__1);
    d__1 = (doublereal) betp1;
    rbet = pow_dd(&c_b2, &d__1);

/*           COMPUTE RI, RJ USING A FORWARD RECURRENCE RELATION. */

    ri[1] = ralf / alfp1;
    rj[1] = rbet / betp1;
    ri[2] = ri[1] * *alfa / alfp2;
    rj[2] = rj[1] * *beta / betp2;
    an = 2.f;
    anm1 = 1.f;
    for (i__ = 3; i__ <= 25; ++i__) {
	ri[i__] = -(ralf + an * (an - alfp2) * ri[i__ - 1]) / (anm1 * (an + 
		alfp1));
	rj[i__] = -(rbet + an * (an - betp2) * rj[i__ - 1]) / (anm1 * (an + 
		betp1));
	anm1 = an;
	an += 1.f;
/* L20: */
    }
    if (*integr == 1) {
	goto L70;
    }
    if (*integr == 3) {
	goto L40;
    }

/*           COMPUTE RG USING A FORWARD RECURRENCE RELATION. */

    rg[1] = -ri[1] / alfp1;
    rg[2] = -(ralf + ralf) / (alfp2 * alfp2) - rg[1];
    an = 2.f;
    anm1 = 1.f;
    im1 = 2;
    for (i__ = 3; i__ <= 25; ++i__) {
	rg[i__] = -(an * (an - alfp2) * rg[im1] - an * ri[im1] + anm1 * ri[
		i__]) / (anm1 * (an + alfp1));
	anm1 = an;
	an += 1.f;
	im1 = i__;
/* L30: */
    }
    if (*integr == 2) {
	goto L70;
    }

/*           COMPUTE RH USING A FORWARD RECURRENCE RELATION. */

L40:
    rh[1] = -rj[1] / betp1;
    rh[2] = -(rbet + rbet) / (betp2 * betp2) - rh[1];
    an = 2.f;
    anm1 = 1.f;
    im1 = 2;
    for (i__ = 3; i__ <= 25; ++i__) {
	rh[i__] = -(an * (an - betp2) * rh[im1] - an * rj[im1] + anm1 * rj[
		i__]) / (anm1 * (an + betp1));
	anm1 = an;
	an += 1.f;
	im1 = i__;
/* L50: */
    }
    for (i__ = 2; i__ <= 25; i__ += 2) {
	rh[i__] = -rh[i__];
/* L60: */
    }
L70:
    for (i__ = 2; i__ <= 25; i__ += 2) {
	rj[i__] = -rj[i__];
/* L80: */
    }
    return 0;
} /* qmomo_ */

