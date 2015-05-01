/* avint.f -- translated by f2c (version 12.02.01).
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

/* DECK AVINT */
/* Subroutine */ int avint_(real *x, real *y, integer *n, real *xlo, real *
	xup, real *ans, integer *ierr)
{
    /* System generated locals */
    integer i__1;
    real r__1;
    doublereal d__1, d__2;

    /* Local variables */
    static doublereal a, b, c__;
    static integer i__;
    static doublereal r3, x1, x2, x3, ca, cb, cc;
    static real fl, fr;
    static doublereal x12, x13, x23, rp5, sum, syl, syu, syl2, syl3, syu2, 
	    syu3;
    static integer inrt;
    static doublereal term1, term2, term3;
    static integer inlft;
    static real slope;
    static integer istop;
    extern /* Subroutine */ int xermsg_(char *, char *, char *, integer *, 
	    integer *, ftnlen, ftnlen, ftnlen);
    static integer istart;

/* ***BEGIN PROLOGUE  AVINT */
/* ***PURPOSE  Integrate a function tabulated at arbitrarily spaced */
/*            abscissas using overlapping parabolas. */
/* ***LIBRARY   SLATEC */
/* ***CATEGORY  H2A1B2 */
/* ***TYPE      SINGLE PRECISION (AVINT-S, DAVINT-D) */
/* ***KEYWORDS  INTEGRATION, QUADRATURE, TABULATED DATA */
/* ***AUTHOR  Jones, R. E., (SNLA) */
/* ***DESCRIPTION */

/*     Abstract */
/*         AVINT integrates a function tabulated at arbitrarily spaced */
/*         abscissas.  The limits of integration need not coincide */
/*         with the tabulated abscissas. */

/*         A method of overlapping parabolas fitted to the data is used */
/*         provided that there are at least 3 abscissas between the */
/*         limits of integration.  AVINT also handles two special cases. */
/*         If the limits of integration are equal, AVINT returns a result */
/*         of zero regardless of the number of tabulated values. */
/*         If there are only two function values, AVINT uses the */
/*         trapezoid rule. */

/*     Description of Parameters */
/*         The user must dimension all arrays appearing in the call list */
/*              X(N), Y(N). */

/*         Input-- */
/*         X    - real array of abscissas, which must be in increasing */
/*                order. */
/*         Y    - real array of functional values. i.e., Y(I)=FUNC(X(I)). */
/*         N    - the integer number of function values supplied. */
/*                N .GE. 2 unless XLO = XUP. */
/*         XLO  - real lower limit of integration. */
/*         XUP  - real upper limit of integration. */
/*                Must have XLO .LE. XUP. */

/*         Output-- */
/*         ANS  - computed approximate value of integral */
/*         IERR - a status code */
/*              --normal code */
/*                =1 means the requested integration was performed. */
/*              --abnormal codes */
/*                =2 means XUP was less than XLO. */
/*                =3 means the number of X(I) between XLO and XUP */
/*                   (inclusive) was less than 3 and neither of the two */
/*                   special cases described in the Abstract occurred. */
/*                   No integration was performed. */
/*                =4 means the restriction X(I+1) .GT. X(I) was violated. */
/*                =5 means the number N of function values was .LT. 2. */
/*                ANS is set to zero if IERR=2,3,4,or 5. */

/*     AVINT is documented completely in SC-M-69-335 */
/*     Original program from "Numerical Integration" by Davis & */
/*     Rabinowitz. */
/*     Adaptation and modifications for Sandia Mathematical Program */
/*     Library by Rondall E. Jones. */

/* ***REFERENCES  R. E. Jones, Approximate integrator of functions */
/*                 tabulated at arbitrarily spaced abscissas, */
/*                 Report SC-M-69-335, Sandia Laboratories, 1969. */
/* ***ROUTINES CALLED  XERMSG */
/* ***REVISION HISTORY  (YYMMDD) */
/*   690901  DATE WRITTEN */
/*   890831  Modified array declarations.  (WRB) */
/*   890831  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ) */
/*   900326  Removed duplicate information from DESCRIPTION section. */
/*           (WRB) */
/*   920501  Reformatted the REFERENCES section.  (WRB) */
/* ***END PROLOGUE  AVINT */

/* ***FIRST EXECUTABLE STATEMENT  AVINT */
    /* Parameter adjustments */
    --y;
    --x;

    /* Function Body */
    *ierr = 1;
    *ans = 0.f;
    if ((r__1 = *xlo - *xup) < 0.f) {
	goto L3;
    } else if (r__1 == 0) {
	goto L100;
    } else {
	goto L200;
    }
L3:
    if (*n < 2) {
	goto L215;
    }
    i__1 = *n;
    for (i__ = 2; i__ <= i__1; ++i__) {
	if (x[i__] <= x[i__ - 1]) {
	    goto L210;
	}
	if (x[i__] > *xup) {
	    goto L6;
	}
/* L5: */
    }
L6:
    if (*n >= 3) {
	goto L9;
    }

/*     SPECIAL N=2 CASE */
    slope = (y[2] - y[1]) / (x[2] - x[1]);
    fl = y[1] + slope * (*xlo - x[1]);
    fr = y[2] + slope * (*xup - x[2]);
    *ans = (fl + fr) * .5f * (*xup - *xlo);
    return 0;
L9:
    if (x[*n - 2] < *xlo) {
	goto L205;
    }
    if (x[3] > *xup) {
	goto L205;
    }
    i__ = 1;
L10:
    if (x[i__] >= *xlo) {
	goto L15;
    }
    ++i__;
    goto L10;
L15:
    inlft = i__;
    i__ = *n;
L20:
    if (x[i__] <= *xup) {
	goto L25;
    }
    --i__;
    goto L20;
L25:
    inrt = i__;
    if (inrt - inlft < 2) {
	goto L205;
    }
    istart = inlft;
    if (inlft == 1) {
	istart = 2;
    }
    istop = inrt;
    if (inrt == *n) {
	istop = *n - 1;
    }

    r3 = 3.;
    rp5 = .5;
    sum = 0.f;
    syl = *xlo;
    syl2 = syl * syl;
    syl3 = syl2 * syl;

    i__1 = istop;
    for (i__ = istart; i__ <= i__1; ++i__) {
	x1 = x[i__ - 1];
	x2 = x[i__];
	x3 = x[i__ + 1];
	x12 = x1 - x2;
	x13 = x1 - x3;
	x23 = x2 - x3;
	term1 = (doublereal) y[i__ - 1] / (x12 * x13);
	term2 = -((doublereal) y[i__]) / (x12 * x23);
	term3 = (doublereal) y[i__ + 1] / (x13 * x23);
	a = term1 + term2 + term3;
	b = -(x2 + x3) * term1 - (x1 + x3) * term2 - (x1 + x2) * term3;
	c__ = x2 * x3 * term1 + x1 * x3 * term2 + x1 * x2 * term3;
	if (i__ - istart <= 0) {
	    goto L30;
	} else {
	    goto L35;
	}
L30:
	ca = a;
	cb = b;
	cc = c__;
	goto L40;
L35:
	ca = (a + ca) * .5f;
	cb = (b + cb) * .5f;
	cc = (c__ + cc) * .5f;
L40:
	syu = x2;
	syu2 = syu * syu;
	syu3 = syu2 * syu;
	sum = sum + ca * (syu3 - syl3) / r3 + cb * rp5 * (syu2 - syl2) + cc * 
		(syu - syl);
	ca = a;
	cb = b;
	cc = c__;
	syl = syu;
	syl2 = syu2;
	syl3 = syu3;
/* L50: */
    }
    syu = *xup;
/* Computing 3rd power */
    d__1 = syu;
/* Computing 2nd power */
    d__2 = syu;
    *ans = sum + ca * (d__1 * (d__1 * d__1) - syl3) / r3 + cb * rp5 * (d__2 * 
	    d__2 - syl2) + cc * (syu - syl);
L100:
    return 0;
L200:
    *ierr = 2;
    xermsg_("SLATEC", "AVINT", "THE UPPER LIMIT OF INTEGRATION WAS NOT GREAT"
	    "ER THAN THE LOWER LIMIT.", &c__4, &c__1, (ftnlen)6, (ftnlen)5, (
	    ftnlen)68);
    return 0;
L205:
    *ierr = 3;
    xermsg_("SLATEC", "AVINT", "THERE WERE LESS THAN THREE FUNCTION VALUES B"
	    "ETWEEN THE LIMITS OF INTEGRATION.", &c__4, &c__1, (ftnlen)6, (
	    ftnlen)5, (ftnlen)77);
    return 0;
L210:
    *ierr = 4;
    xermsg_("SLATEC", "AVINT", "THE ABSCISSAS WERE NOT STRICTLY INCREASING. "
	    " MUST HAVE X(I-1) .LT. X(I) FOR ALL I.", &c__4, &c__1, (ftnlen)6, 
	    (ftnlen)5, (ftnlen)82);
    return 0;
L215:
    *ierr = 5;
    xermsg_("SLATEC", "AVINT", "LESS THAN TWO FUNCTION VALUES WERE SUPPLIED.",
	     &c__4, &c__1, (ftnlen)6, (ftnlen)5, (ftnlen)44);
    return 0;
} /* avint_ */

