/* dqcheb.f -- translated by f2c (version 12.02.01).
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

/* DECK DQCHEB */
/* Subroutine */ int dqcheb_(doublereal *x, doublereal *fval, doublereal *
	cheb12, doublereal *cheb24)
{
    static integer i__, j;
    static doublereal v[12], alam, alam1, alam2, part1, part2, part3;

/* ***BEGIN PROLOGUE  DQCHEB */
/* ***SUBSIDIARY */
/* ***PURPOSE  This routine computes the CHEBYSHEV series expansion */
/*            of degrees 12 and 24 of a function using A */
/*            FAST FOURIER TRANSFORM METHOD */
/*            F(X) = SUM(K=1,..,13) (CHEB12(K)*T(K-1,X)), */
/*            F(X) = SUM(K=1,..,25) (CHEB24(K)*T(K-1,X)), */
/*            Where T(K,X) is the CHEBYSHEV POLYNOMIAL OF DEGREE K. */
/* ***LIBRARY   SLATEC */
/* ***TYPE      DOUBLE PRECISION (QCHEB-S, DQCHEB-D) */
/* ***KEYWORDS  CHEBYSHEV SERIES EXPANSION, FAST FOURIER TRANSFORM */
/* ***AUTHOR  Piessens, Robert */
/*             Applied Mathematics and Programming Division */
/*             K. U. Leuven */
/*           de Doncker, Elise */
/*             Applied Mathematics and Programming Division */
/*             K. U. Leuven */
/* ***DESCRIPTION */

/*        Chebyshev Series Expansion */
/*        Standard Fortran Subroutine */
/*        Double precision version */

/*        PARAMETERS */
/*          ON ENTRY */
/*           X      - Double precision */
/*                    Vector of dimension 11 containing the */
/*                    Values COS(K*PI/24), K = 1, ..., 11 */

/*           FVAL   - Double precision */
/*                    Vector of dimension 25 containing the */
/*                    function values at the points */
/*                    (B+A+(B-A)*COS(K*PI/24))/2, K = 0, ...,24, */
/*                    where (A,B) is the approximation interval. */
/*                    FVAL(1) and FVAL(25) are divided by two */
/*                    (these values are destroyed at output). */

/*          ON RETURN */
/*           CHEB12 - Double precision */
/*                    Vector of dimension 13 containing the */
/*                    CHEBYSHEV coefficients for degree 12 */

/*           CHEB24 - Double precision */
/*                    Vector of dimension 25 containing the */
/*                    CHEBYSHEV Coefficients for degree 24 */

/* ***SEE ALSO  DQC25C, DQC25F, DQC25S */
/* ***ROUTINES CALLED  (NONE) */
/* ***REVISION HISTORY  (YYMMDD) */
/*   810101  DATE WRITTEN */
/*   830518  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900328  Added TYPE section.  (WRB) */
/* ***END PROLOGUE  DQCHEB */



/* ***FIRST EXECUTABLE STATEMENT  DQCHEB */
    /* Parameter adjustments */
    --cheb24;
    --cheb12;
    --fval;
    --x;

    /* Function Body */
    for (i__ = 1; i__ <= 12; ++i__) {
	j = 26 - i__;
	v[i__ - 1] = fval[i__] - fval[j];
	fval[i__] += fval[j];
/* L10: */
    }
    alam1 = v[0] - v[8];
    alam2 = x[6] * (v[2] - v[6] - v[10]);
    cheb12[4] = alam1 + alam2;
    cheb12[10] = alam1 - alam2;
    alam1 = v[1] - v[7] - v[9];
    alam2 = v[3] - v[5] - v[11];
    alam = x[3] * alam1 + x[9] * alam2;
    cheb24[4] = cheb12[4] + alam;
    cheb24[22] = cheb12[4] - alam;
    alam = x[9] * alam1 - x[3] * alam2;
    cheb24[10] = cheb12[10] + alam;
    cheb24[16] = cheb12[10] - alam;
    part1 = x[4] * v[4];
    part2 = x[8] * v[8];
    part3 = x[6] * v[6];
    alam1 = v[0] + part1 + part2;
    alam2 = x[2] * v[2] + part3 + x[10] * v[10];
    cheb12[2] = alam1 + alam2;
    cheb12[12] = alam1 - alam2;
    alam = x[1] * v[1] + x[3] * v[3] + x[5] * v[5] + x[7] * v[7] + x[9] * v[9]
	     + x[11] * v[11];
    cheb24[2] = cheb12[2] + alam;
    cheb24[24] = cheb12[2] - alam;
    alam = x[11] * v[1] - x[9] * v[3] + x[7] * v[5] - x[5] * v[7] + x[3] * v[
	    9] - x[1] * v[11];
    cheb24[12] = cheb12[12] + alam;
    cheb24[14] = cheb12[12] - alam;
    alam1 = v[0] - part1 + part2;
    alam2 = x[10] * v[2] - part3 + x[2] * v[10];
    cheb12[6] = alam1 + alam2;
    cheb12[8] = alam1 - alam2;
    alam = x[5] * v[1] - x[9] * v[3] - x[1] * v[5] - x[11] * v[7] + x[3] * v[
	    9] + x[7] * v[11];
    cheb24[6] = cheb12[6] + alam;
    cheb24[20] = cheb12[6] - alam;
    alam = x[7] * v[1] - x[3] * v[3] - x[11] * v[5] + x[1] * v[7] - x[9] * v[
	    9] - x[5] * v[11];
    cheb24[8] = cheb12[8] + alam;
    cheb24[18] = cheb12[8] - alam;
    for (i__ = 1; i__ <= 6; ++i__) {
	j = 14 - i__;
	v[i__ - 1] = fval[i__] - fval[j];
	fval[i__] += fval[j];
/* L20: */
    }
    alam1 = v[0] + x[8] * v[4];
    alam2 = x[4] * v[2];
    cheb12[3] = alam1 + alam2;
    cheb12[11] = alam1 - alam2;
    cheb12[7] = v[0] - v[4];
    alam = x[2] * v[1] + x[6] * v[3] + x[10] * v[5];
    cheb24[3] = cheb12[3] + alam;
    cheb24[23] = cheb12[3] - alam;
    alam = x[6] * (v[1] - v[3] - v[5]);
    cheb24[7] = cheb12[7] + alam;
    cheb24[19] = cheb12[7] - alam;
    alam = x[10] * v[1] - x[6] * v[3] + x[2] * v[5];
    cheb24[11] = cheb12[11] + alam;
    cheb24[15] = cheb12[11] - alam;
    for (i__ = 1; i__ <= 3; ++i__) {
	j = 8 - i__;
	v[i__ - 1] = fval[i__] - fval[j];
	fval[i__] += fval[j];
/* L30: */
    }
    cheb12[5] = v[0] + x[8] * v[2];
    cheb12[9] = fval[1] - x[8] * fval[3];
    alam = x[4] * v[1];
    cheb24[5] = cheb12[5] + alam;
    cheb24[21] = cheb12[5] - alam;
    alam = x[8] * fval[2] - fval[4];
    cheb24[9] = cheb12[9] + alam;
    cheb24[17] = cheb12[9] - alam;
    cheb12[1] = fval[1] + fval[3];
    alam = fval[2] + fval[4];
    cheb24[1] = cheb12[1] + alam;
    cheb24[25] = cheb12[1] - alam;
    cheb12[13] = v[0] - v[2];
    cheb24[13] = cheb12[13];
    alam = .16666666666666666;
    for (i__ = 2; i__ <= 12; ++i__) {
	cheb12[i__] *= alam;
/* L40: */
    }
    alam *= .5;
    cheb12[1] *= alam;
    cheb12[13] *= alam;
    for (i__ = 2; i__ <= 24; ++i__) {
	cheb24[i__] *= alam;
/* L50: */
    }
    cheb24[1] = alam * .5 * cheb24[1];
    cheb24[25] = alam * .5 * cheb24[25];
    return 0;
} /* dqcheb_ */

