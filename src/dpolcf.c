/* dpolcf.f -- translated by f2c (version 12.02.01).
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

/* DECK DPOLCF */
/* Subroutine */ int dpolcf_(doublereal *xx, integer *n, doublereal *x, 
	doublereal *c__, doublereal *d__, doublereal *work)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer i__, k, im1, km1, nm1, km2n;
    static doublereal pone, ptwo;
    static integer km1pi, nmkp1, npkm1, km2npi;

/* ***BEGIN PROLOGUE  DPOLCF */
/* ***PURPOSE  Compute the coefficients of the polynomial fit (including */
/*            Hermite polynomial fits) produced by a previous call to */
/*            POLINT. */
/* ***LIBRARY   SLATEC */
/* ***CATEGORY  E1B */
/* ***TYPE      DOUBLE PRECISION (POLCOF-S, DPOLCF-D) */
/* ***KEYWORDS  COEFFICIENTS, POLYNOMIAL */
/* ***AUTHOR  Huddleston, R. E., (SNLL) */
/* ***DESCRIPTION */

/*     Abstract */
/*        Subroutine DPOLCF computes the coefficients of the polynomial */
/*     fit (including Hermite polynomial fits ) produced by a previous */
/*     call to DPLINT.  The coefficients of the polynomial, expanded */
/*     about XX, are stored in the array D. The expansion is of the form */
/*     P(Z) = D(1) + D(2)*(Z-XX) +D(3)*((Z-XX)**2) + ... + */
/*                                                  D(N)*((Z-XX)**(N-1)). */
/*     Between the call to DPLINT and the call to DPOLCF the variable N */
/*     and the arrays X and C must not be altered. */

/*     *****  INPUT PARAMETERS */
/*      *** All TYPE REAL variables are DOUBLE PRECISION *** */

/*     XX   - The point about which the Taylor expansion is to be made. */

/*     N    - **** */
/*            *     N, X, and C must remain unchanged between the */
/*     X    - *     call to DPLINT and the call to DPOLCF. */
/*     C    - **** */

/*     *****  OUTPUT PARAMETER */
/*      *** All TYPE REAL variables are DOUBLE PRECISION *** */

/*     D    - The array of coefficients for the Taylor expansion as */
/*            explained in the abstract */

/*     *****  STORAGE PARAMETER */

/*     WORK - This is an array to provide internal working storage. It */
/*            must be dimensioned by at least 2*N in the calling program. */


/*     **** Note - There are two methods for evaluating the fit produced */
/*     by DPLINT. You may call DPOLVL to perform the task, or you may */
/*     call DPOLCF to obtain the coefficients of the Taylor expansion and */
/*     then write your own evaluation scheme. Due to the inherent errors */
/*     in the computations of the Taylor expansion from the Newton */
/*     coefficients produced by DPLINT, much more accuracy may be */
/*     expected by calling DPOLVL as opposed to writing your own scheme. */

/* ***REFERENCES  (NONE) */
/* ***ROUTINES CALLED  (NONE) */
/* ***REVISION HISTORY  (YYMMDD) */
/*   890213  DATE WRITTEN */
/*   891006  Cosmetic changes to prologue.  (WRB) */
/*   891024  Corrected KEYWORD section.  (WRB) */
/*   891024  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/* ***END PROLOGUE  DPOLCF */

/* ***FIRST EXECUTABLE STATEMENT  DPOLCF */
    /* Parameter adjustments */
    --work;
    --d__;
    --c__;
    --x;

    /* Function Body */
    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
	d__[k] = c__[k];
/* L10010: */
    }
    if (*n == 1) {
	return 0;
    }
    work[1] = 1.;
    pone = c__[1];
    nm1 = *n - 1;
    i__1 = *n;
    for (k = 2; k <= i__1; ++k) {
	km1 = k - 1;
	npkm1 = *n + k - 1;
	work[npkm1] = *xx - x[km1];
	work[k] = work[npkm1] * work[km1];
	ptwo = pone + work[k] * c__[k];
	pone = ptwo;
/* L10020: */
    }
    d__[1] = ptwo;
    if (*n == 2) {
	return 0;
    }
    i__1 = nm1;
    for (k = 2; k <= i__1; ++k) {
	km1 = k - 1;
	km2n = k - 2 + *n;
	nmkp1 = *n - k + 1;
	i__2 = nmkp1;
	for (i__ = 2; i__ <= i__2; ++i__) {
	    km2npi = km2n + i__;
	    im1 = i__ - 1;
	    km1pi = km1 + i__;
	    work[i__] = work[km2npi] * work[im1] + work[i__];
	    d__[k] += work[i__] * d__[km1pi];
/* L10030: */
	}
    }
    return 0;
} /* dpolcf_ */

