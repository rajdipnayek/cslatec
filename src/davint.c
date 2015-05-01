/* davint.f -- translated by f2c (version 12.02.01).
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

/* DECK DAVINT */
/* Subroutine */ int davint_(doublereal *x, doublereal *y, integer *n, 
	doublereal *xlo, doublereal *xup, doublereal *ans, integer *ierr)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2;

    /* Local variables */
    static doublereal a, b, c__;
    static integer i__;
    static doublereal r3, x1, x2, x3, ca, cb, cc, fl, fr, x12, x13, x23, rp5, 
	    sum, syl, syu, syl2, syl3, syu2, syu3;
    static integer inrt;
    static doublereal term1, term2, term3;
    static integer inlft;
    static doublereal slope;
    static integer istop;
    extern /* Subroutine */ int xermsg_(char *, char *, char *, integer *, 
	    integer *, ftnlen, ftnlen, ftnlen);
    static integer istart;

/* ***BEGIN PROLOGUE  DAVINT */
/* ***PURPOSE  Integrate a function tabulated at arbitrarily spaced */
/*            abscissas using overlapping parabolas. */
/* ***LIBRARY   SLATEC */
/* ***CATEGORY  H2A1B2 */
/* ***TYPE      DOUBLE PRECISION (AVINT-S, DAVINT-D) */
/* ***KEYWORDS  INTEGRATION, QUADRATURE, TABULATED DATA */
/* ***AUTHOR  Jones, R. E., (SNLA) */
/* ***DESCRIPTION */

/*     Abstract */
/*         DAVINT integrates a function tabulated at arbitrarily spaced */
/*         abscissas.  The limits of integration need not coincide */
/*         with the tabulated abscissas. */

/*         A method of overlapping parabolas fitted to the data is used */
/*         provided that there are at least 3 abscissas between the */
/*         limits of integration.  DAVINT also handles two special cases. */
/*         If the limits of integration are equal, DAVINT returns a */
/*         result of zero regardless of the number of tabulated values. */
/*         If there are only two function values, DAVINT uses the */
/*         trapezoid rule. */

/*     Description of Parameters */
/*         The user must dimension all arrays appearing in the call list */
/*              X(N), Y(N) */

/*         Input-- */
/*      X    - DOUBLE PRECISION array of abscissas, which must be in */
/*             increasing order. */
/*      Y    - DOUBLE PRECISION array of function values. i.e., */
/*                Y(I)=FUNC(X(I)) */
/*      N    - The integer number of function values supplied. */
/*                N .GE. 2 unless XLO = XUP. */
/*      XLO  - DOUBLE PRECISION lower limit of integration */
/*      XUP  - DOUBLE PRECISION upper limit of integration.  Must have */
/*              XLO.LE.XUP */

/*         Output-- */
/*      ANS  - Double Precision computed approximate value of integral */
/*      IERR - A status code */
/*           --Normal Code */
/*                =1 Means the requested integration was performed. */
/*           --Abnormal Codes */
/*                =2 Means XUP was less than XLO. */
/*                =3 Means the number of X(I) between XLO and XUP */
/*                   (inclusive) was less than 3 and neither of the two */
/*                   special cases described in the abstract occurred. */
/*                   No integration was performed. */
/*                =4 Means the restriction X(I+1).GT.X(I) was violated. */
/*                =5 Means the number N of function values was .lt. 2. */
/*                   ANS is set to zero if IERR=2,3,4,or 5. */

/*    DAVINT is documented completely in SC-M-69-335 */
/*    Original program from *Numerical Integration* by Davis & Rabinowitz */
/*    Adaptation and modifications by Rondall E Jones. */

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
/*   920501  Reformatted the REFERENCES section.  (WRB) */
/* ***END PROLOGUE  DAVINT */

/*     BEGIN BLOCK PERMITTING ...EXITS TO 190 */
/*        BEGIN BLOCK PERMITTING ...EXITS TO 180 */
/* ***FIRST EXECUTABLE STATEMENT  DAVINT */
    /* Parameter adjustments */
    --y;
    --x;

    /* Function Body */
    *ierr = 1;
    *ans = 0.;
    if (*xlo > *xup) {
	goto L160;
    }
    if (*xlo == *xup) {
	goto L150;
    }
    if (*n >= 2) {
	goto L10;
    }
    *ierr = 5;
    xermsg_("SLATEC", "DAVINT", "LESS THAN TWO FUNCTION VALUES WERE SUPPLIED."
	    , &c__4, &c__1, (ftnlen)6, (ftnlen)6, (ftnlen)44);
/*     ...............EXIT */
    goto L190;
L10:
    i__1 = *n;
    for (i__ = 2; i__ <= i__1; ++i__) {
/*        ............EXIT */
	if (x[i__] <= x[i__ - 1]) {
	    goto L180;
	}
/*                 ...EXIT */
	if (x[i__] > *xup) {
	    goto L30;
	}
/* L20: */
    }
L30:
    if (*n >= 3) {
	goto L40;
    }

/*                    SPECIAL N=2 CASE */
    slope = (y[2] - y[1]) / (x[2] - x[1]);
    fl = y[1] + slope * (*xlo - x[1]);
    fr = y[2] + slope * (*xup - x[2]);
    *ans = (fl + fr) * .5 * (*xup - *xlo);
/*     ...............EXIT */
    goto L190;
L40:
    if (x[*n - 2] >= *xlo) {
	goto L50;
    }
    *ierr = 3;
    xermsg_("SLATEC", "DAVINT", "THERE WERE LESS THAN THREE FUNCTION VALUES "
	    "BETWEEN THE LIMITS OF INTEGRATION.", &c__4, &c__1, (ftnlen)6, (
	    ftnlen)6, (ftnlen)77);
/*     ...............EXIT */
    goto L190;
L50:
    if (x[3] <= *xup) {
	goto L60;
    }
    *ierr = 3;
    xermsg_("SLATEC", "DAVINT", "THERE WERE LESS THAN THREE FUNCTION VALUES "
	    "BETWEEN THE LIMITS OF INTEGRATION.", &c__4, &c__1, (ftnlen)6, (
	    ftnlen)6, (ftnlen)77);
/*     ...............EXIT */
    goto L190;
L60:
    i__ = 1;
L70:
    if (x[i__] >= *xlo) {
	goto L80;
    }
    ++i__;
    goto L70;
L80:
    inlft = i__;
    i__ = *n;
L90:
    if (x[i__] <= *xup) {
	goto L100;
    }
    --i__;
    goto L90;
L100:
    inrt = i__;
    if (inrt - inlft >= 2) {
	goto L110;
    }
    *ierr = 3;
    xermsg_("SLATEC", "DAVINT", "THERE WERE LESS THAN THREE FUNCTION VALUES "
	    "BETWEEN THE LIMITS OF INTEGRATION.", &c__4, &c__1, (ftnlen)6, (
	    ftnlen)6, (ftnlen)77);
/*     ...............EXIT */
    goto L190;
L110:
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
    sum = 0.;
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
	term1 = y[i__ - 1] / (x12 * x13);
	term2 = -y[i__] / (x12 * x23);
	term3 = y[i__ + 1] / (x13 * x23);
	a = term1 + term2 + term3;
	b = -(x2 + x3) * term1 - (x1 + x3) * term2 - (x1 + x2) * term3;
	c__ = x2 * x3 * term1 + x1 * x3 * term2 + x1 * x2 * term3;
	if (i__ > istart) {
	    goto L120;
	}
	ca = a;
	cb = b;
	cc = c__;
	goto L130;
L120:
	ca = (a + ca) * .5;
	cb = (b + cb) * .5;
	cc = (c__ + cc) * .5;
L130:
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
/* L140: */
    }
    syu = *xup;
/* Computing 3rd power */
    d__1 = syu;
/* Computing 2nd power */
    d__2 = syu;
    *ans = sum + ca * (d__1 * (d__1 * d__1) - syl3) / r3 + cb * rp5 * (d__2 * 
	    d__2 - syl2) + cc * (syu - syl);
L150:
    goto L170;
L160:
    *ierr = 2;
    xermsg_("SLATEC", "DAVINT", "THE UPPER LIMIT OF INTEGRATION WAS NOT GREA"
	    "TER THAN THE LOWER LIMIT.", &c__4, &c__1, (ftnlen)6, (ftnlen)6, (
	    ftnlen)68);
L170:
/*     ......EXIT */
    goto L190;
L180:
    *ierr = 4;
    xermsg_("SLATEC", "DAVINT", "THE ABSCISSAS WERE NOT STRICTLY INCREASING."
	    "  MUST HAVE X(I-1) .LT. X(I) FOR ALL I.", &c__4, &c__1, (ftnlen)6,
	     (ftnlen)6, (ftnlen)82);
L190:
    return 0;
} /* davint_ */

