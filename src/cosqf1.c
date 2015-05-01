/* cosqf1.f -- translated by f2c (version 12.02.01).
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

/* DECK COSQF1 */
/* Subroutine */ int cosqf1_(integer *n, real *x, real *w, real *xh)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__, k, kc, np2, ns2;
    static real xim1;
    static integer modn;
    extern /* Subroutine */ int rfftf_(integer *, real *, real *);

/* ***BEGIN PROLOGUE  COSQF1 */
/* ***SUBSIDIARY */
/* ***PURPOSE  Compute the forward cosine transform with odd wave numbers. */
/* ***LIBRARY   SLATEC (FFTPACK) */
/* ***CATEGORY  J1A3 */
/* ***TYPE      SINGLE PRECISION (COSQF1-S) */
/* ***KEYWORDS  FFTPACK, FOURIER TRANSFORM */
/* ***AUTHOR  Swarztrauber, P. N., (NCAR) */
/* ***DESCRIPTION */

/*  Subroutine COSQF1 computes the fast Fourier transform of quarter */
/*  wave data. That is, COSQF1 computes the coefficients in a cosine */
/*  series representation with only odd wave numbers.  The transform */
/*  is defined below at Output Parameter X */

/* ***REFERENCES  P. N. Swarztrauber, Vectorizing the FFTs, in Parallel */
/*                 Computations (G. Rodrigue, ed.), Academic Press, */
/*                 1982, pp. 51-83. */
/* ***ROUTINES CALLED  RFFTF */
/* ***REVISION HISTORY  (YYMMDD) */
/*   790601  DATE WRITTEN */
/*   830401  Modified to use SLATEC library source file format. */
/*   860115  Modified by Ron Boisvert to adhere to Fortran 77 by */
/*           changing dummy array size declarations (1) to (*). */
/*   881128  Modified by Dick Valent to meet prologue standards. */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   920501  Reformatted the REFERENCES section.  (WRB) */
/* ***END PROLOGUE  COSQF1 */
/* ***FIRST EXECUTABLE STATEMENT  COSQF1 */
    /* Parameter adjustments */
    --xh;
    --w;
    --x;

    /* Function Body */
    ns2 = (*n + 1) / 2;
    np2 = *n + 2;
    i__1 = ns2;
    for (k = 2; k <= i__1; ++k) {
	kc = np2 - k;
	xh[k] = x[k] + x[kc];
	xh[kc] = x[k] - x[kc];
/* L101: */
    }
    modn = *n % 2;
    if (modn == 0) {
	xh[ns2 + 1] = x[ns2 + 1] + x[ns2 + 1];
    }
    i__1 = ns2;
    for (k = 2; k <= i__1; ++k) {
	kc = np2 - k;
	x[k] = w[k - 1] * xh[kc] + w[kc - 1] * xh[k];
	x[kc] = w[k - 1] * xh[k] - w[kc - 1] * xh[kc];
/* L102: */
    }
    if (modn == 0) {
	x[ns2 + 1] = w[ns2] * xh[ns2 + 1];
    }
    rfftf_(n, &x[1], &xh[1]);
    i__1 = *n;
    for (i__ = 3; i__ <= i__1; i__ += 2) {
	xim1 = x[i__ - 1] - x[i__];
	x[i__] = x[i__ - 1] + x[i__];
	x[i__ - 1] = xim1;
/* L103: */
    }
    return 0;
} /* cosqf1_ */

