/* dppval.f -- translated by f2c (version 12.02.01).
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

/* DECK DPPVAL */
doublereal dppval_(integer *ldc, doublereal *c__, doublereal *xi, integer *
	lxi, integer *k, integer *ideriv, doublereal *x, integer *inppv)
{
    /* System generated locals */
    integer c_dim1, c_offset;
    doublereal ret_val;

    /* Local variables */
    static integer i__, j, kk;
    static doublereal dx;
    extern /* Subroutine */ int xermsg_(char *, char *, char *, integer *, 
	    integer *, ftnlen, ftnlen, ftnlen), dintrv_(doublereal *, integer 
	    *, doublereal *, integer *, integer *, integer *);
    static integer ndummy;

/* ***BEGIN PROLOGUE  DPPVAL */
/* ***PURPOSE  Calculate the value of the IDERIV-th derivative of the */
/*            B-spline from the PP-representation. */
/* ***LIBRARY   SLATEC */
/* ***CATEGORY  E3, K6 */
/* ***TYPE      DOUBLE PRECISION (PPVAL-S, DPPVAL-D) */
/* ***KEYWORDS  B-SPLINE, DATA FITTING, INTERPOLATION, SPLINES */
/* ***AUTHOR  Amos, D. E., (SNLA) */
/* ***DESCRIPTION */

/*     Written by Carl de Boor and modified by D. E. Amos */

/*     Abstract    **** a double precision routine **** */
/*         DPPVAL is the PPVALU function of the reference. */

/*         DPPVAL calculates (at X) the value of the IDERIV-th */
/*         derivative of the B-spline from the PP-representation */
/*         (C,XI,LXI,K).  The Taylor expansion about XI(J) for X in */
/*         the interval XI(J) .LE. X .LT. XI(J+1) is evaluated, J=1,LXI. */
/*         Right limiting values at X=XI(J) are obtained.  DPPVAL will */
/*         extrapolate beyond XI(1) and XI(LXI+1). */

/*         To obtain left limiting values (left derivatives) at XI(J) */
/*         replace LXI by J-1 and set X=XI(J),J=2,LXI+1. */

/*     Description of Arguments */

/*         Input      C,XI,X are double precision */
/*          LDC     - leading dimension of C matrix, LDC .GE. K */
/*          C       - matrix of dimension at least (K,LXI) containing */
/*                    right derivatives at break points XI(*). */
/*          XI      - break point vector of length LXI+1 */
/*          LXI     - number of polynomial pieces */
/*          K       - order of B-spline, K .GE. 1 */
/*          IDERIV  - order of the derivative, 0 .LE. IDERIV .LE. K-1 */
/*                    IDERIV=0 gives the B-spline value */
/*          X       - argument, XI(1) .LE. X .LE. XI(LXI+1) */
/*          INPPV   - an initialization parameter which must be set */
/*                    to 1 the first time DPPVAL is called. */

/*         Output     DPPVAL is double precision */
/*          INPPV   - INPPV contains information for efficient process- */
/*                    ing after the initial call and INPPV must not */
/*                    be changed by the user.  Distinct splines require */
/*                    distinct INPPV parameters. */
/*          DPPVAL  - value of the IDERIV-th derivative at X */

/*     Error Conditions */
/*         Improper input is a fatal error */

/* ***REFERENCES  Carl de Boor, Package for calculating with B-splines, */
/*                 SIAM Journal on Numerical Analysis 14, 3 (June 1977), */
/*                 pp. 441-472. */
/* ***ROUTINES CALLED  DINTRV, XERMSG */
/* ***REVISION HISTORY  (YYMMDD) */
/*   800901  DATE WRITTEN */
/*   890831  Modified array declarations.  (WRB) */
/*   890831  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ) */
/*   920501  Reformatted the REFERENCES section.  (WRB) */
/* ***END PROLOGUE  DPPVAL */

/* ***FIRST EXECUTABLE STATEMENT  DPPVAL */
    /* Parameter adjustments */
    c_dim1 = *ldc;
    c_offset = 1 + c_dim1;
    c__ -= c_offset;
    --xi;

    /* Function Body */
    ret_val = 0.;
    if (*k < 1) {
	goto L90;
    }
    if (*ldc < *k) {
	goto L80;
    }
    if (*lxi < 1) {
	goto L85;
    }
    if (*ideriv < 0 || *ideriv >= *k) {
	goto L95;
    }
    i__ = *k - *ideriv;
    kk = i__;
    dintrv_(&xi[1], lxi, x, inppv, &i__, &ndummy);
    dx = *x - xi[i__];
    j = *k;
L10:
    ret_val = ret_val / kk * dx + c__[j + i__ * c_dim1];
    --j;
    --kk;
    if (kk > 0) {
	goto L10;
    }
    return ret_val;


L80:
    xermsg_("SLATEC", "DPPVAL", "LDC DOES NOT SATISFY LDC.GE.K", &c__2, &c__1,
	     (ftnlen)6, (ftnlen)6, (ftnlen)29);
    return ret_val;
L85:
    xermsg_("SLATEC", "DPPVAL", "LXI DOES NOT SATISFY LXI.GE.1", &c__2, &c__1,
	     (ftnlen)6, (ftnlen)6, (ftnlen)29);
    return ret_val;
L90:
    xermsg_("SLATEC", "DPPVAL", "K DOES NOT SATISFY K.GE.1", &c__2, &c__1, (
	    ftnlen)6, (ftnlen)6, (ftnlen)25);
    return ret_val;
L95:
    xermsg_("SLATEC", "DPPVAL", "IDERIV DOES NOT SATISFY 0.LE.IDERIV.LT.K", &
	    c__2, &c__1, (ftnlen)6, (ftnlen)6, (ftnlen)40);
    return ret_val;
} /* dppval_ */

