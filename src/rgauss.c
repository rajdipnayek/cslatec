/* rgauss.f -- translated by f2c (version 12.02.01).
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

static real c_b3 = 0.f;

/* DECK RGAUSS */
doublereal rgauss_(real *xmean, real *sd)
{
    /* System generated locals */
    real ret_val;

    /* Local variables */
    static integer i__;
    extern doublereal rand_(real *);

/* ***BEGIN PROLOGUE  RGAUSS */
/* ***PURPOSE  Generate a normally distributed (Gaussian) random number. */
/* ***LIBRARY   SLATEC (FNLIB) */
/* ***CATEGORY  L6A14 */
/* ***TYPE      SINGLE PRECISION (RGAUSS-S) */
/* ***KEYWORDS  FNLIB, GAUSSIAN, NORMAL, RANDOM NUMBER, SPECIAL FUNCTIONS */
/* ***AUTHOR  Fullerton, W., (LANL) */
/* ***DESCRIPTION */

/* Generate a normally distributed random number, i.e., generate random */
/* numbers with a Gaussian distribution.  These random numbers are not */
/* exceptionally good -- especially in the tails of the distribution, */
/* but this implementation is simple and suitable for most applications. */
/* See R. W. Hamming, Numerical Methods for Scientists and Engineers, */
/* McGraw-Hill, 1962, pages 34 and 389. */

/*             Input Arguments -- */
/* XMEAN  the mean of the Guassian distribution. */
/* SD     the standard deviation of the Guassian function */
/*          EXP (-1/2 * (X-XMEAN)**2 / SD**2) */

/* ***REFERENCES  (NONE) */
/* ***ROUTINES CALLED  RAND */
/* ***REVISION HISTORY  (YYMMDD) */
/*   770401  DATE WRITTEN */
/*   861211  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   910819  Added EXTERNAL statement for RAND due to problem on IBM */
/*           RS 6000.  (WRB) */
/* ***END PROLOGUE  RGAUSS */
/* ***FIRST EXECUTABLE STATEMENT  RGAUSS */
    ret_val = -6.f;
    for (i__ = 1; i__ <= 12; ++i__) {
	ret_val += rand_(&c_b3);
/* L10: */
    }

    ret_val = *xmean + *sd * ret_val;

    return ret_val;
} /* rgauss_ */

