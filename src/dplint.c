/* dplint.f -- translated by f2c (version 12.02.01).
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

static integer c__2 = 2;
static integer c__1 = 1;

/* DECK DPLINT */
/* Subroutine */ int dplint_(integer *n, doublereal *x, doublereal *y, 
	doublereal *c__)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer i__, k, km1;
    static doublereal dif;
    extern /* Subroutine */ int xermsg_(char *, char *, char *, integer *, 
	    integer *, ftnlen, ftnlen, ftnlen);

/* ***BEGIN PROLOGUE  DPLINT */
/* ***PURPOSE  Produce the polynomial which interpolates a set of discrete */
/*            data points. */
/* ***LIBRARY   SLATEC */
/* ***CATEGORY  E1B */
/* ***TYPE      DOUBLE PRECISION (POLINT-S, DPLINT-D) */
/* ***KEYWORDS  POLYNOMIAL INTERPOLATION */
/* ***AUTHOR  Huddleston, R. E., (SNLL) */
/* ***DESCRIPTION */

/*     Abstract */
/*        Subroutine DPLINT is designed to produce the polynomial which */
/*     interpolates the data  (X(I),Y(I)), I=1,...,N.  DPLINT sets up */
/*     information in the array C which can be used by subroutine DPOLVL */
/*     to evaluate the polynomial and its derivatives and by subroutine */
/*     DPOLCF to produce the coefficients. */

/*     Formal Parameters */
/*     *** All TYPE REAL variables are DOUBLE PRECISION *** */
/*     N  - the number of data points  (N .GE. 1) */
/*     X  - the array of abscissas (all of which must be distinct) */
/*     Y  - the array of ordinates */
/*     C  - an array of information used by subroutines */
/*     *******  Dimensioning Information  ******* */
/*     Arrays X,Y, and C must be dimensioned at least N in the calling */
/*     program. */

/* ***REFERENCES  L. F. Shampine, S. M. Davenport and R. E. Huddleston, */
/*                 Curve fitting by polynomials in one variable, Report */
/*                 SLA-74-0270, Sandia Laboratories, June 1974. */
/* ***ROUTINES CALLED  XERMSG */
/* ***REVISION HISTORY  (YYMMDD) */
/*   740601  DATE WRITTEN */
/*   891006  Cosmetic changes to prologue.  (WRB) */
/*   891006  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ) */
/*   920501  Reformatted the REFERENCES section.  (WRB) */
/* ***END PROLOGUE  DPLINT */
/* ***FIRST EXECUTABLE STATEMENT  DPLINT */
    /* Parameter adjustments */
    --c__;
    --y;
    --x;

    /* Function Body */
    if (*n <= 0) {
	goto L91;
    }
    c__[1] = y[1];
    if (*n == 1) {
	return 0;
    }
    i__1 = *n;
    for (k = 2; k <= i__1; ++k) {
	c__[k] = y[k];
	km1 = k - 1;
	i__2 = km1;
	for (i__ = 1; i__ <= i__2; ++i__) {
/*     CHECK FOR DISTINCT X VALUES */
	    dif = x[i__] - x[k];
	    if (dif == 0.f) {
		goto L92;
	    }
	    c__[k] = (c__[i__] - c__[k]) / dif;
/* L10010: */
	}
    }
    return 0;
L91:
    xermsg_("SLATEC", "DPLINT", "N IS ZERO OR NEGATIVE.", &c__2, &c__1, (
	    ftnlen)6, (ftnlen)6, (ftnlen)22);
    return 0;
L92:
    xermsg_("SLATEC", "DPLINT", "THE ABSCISSAS ARE NOT DISTINCT.", &c__2, &
	    c__1, (ftnlen)6, (ftnlen)6, (ftnlen)31);
    return 0;
} /* dplint_ */

