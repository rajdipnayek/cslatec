/* dbspdr.f -- translated by f2c (version 12.02.01).
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

/* DECK DBSPDR */
/* Subroutine */ int dbspdr_(doublereal *t, doublereal *a, integer *n, 
	integer *k, integer *nderiv, doublereal *ad)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer i__, id, ii, jj, jm;
    static doublereal diff;
    static integer kmid;
    static doublereal fkmid;
    static integer ipkmid;
    extern /* Subroutine */ int xermsg_(char *, char *, char *, integer *, 
	    integer *, ftnlen, ftnlen, ftnlen);

/* ***BEGIN PROLOGUE  DBSPDR */
/* ***PURPOSE  Use the B-representation to construct a divided difference */
/*            table preparatory to a (right) derivative calculation. */
/* ***LIBRARY   SLATEC */
/* ***CATEGORY  E3, K6 */
/* ***TYPE      DOUBLE PRECISION (BSPDR-S, DBSPDR-D) */
/* ***KEYWORDS  B-SPLINE, DATA FITTING, DIFFERENTIATION OF SPLINES, */
/*             INTERPOLATION */
/* ***AUTHOR  Amos, D. E., (SNLA) */
/* ***DESCRIPTION */

/*     Written by Carl de Boor and modified by D. E. Amos */

/*     Abstract     **** a double precision routine **** */
/*         DBSPDR is the BSPLDR routine of the reference. */

/*         DBSPDR uses the B-representation (T,A,N,K) to construct a */
/*         divided difference table ADIF preparatory to a (right) */
/*         derivative calculation in DBSPEV.  The lower triangular matrix */
/*         ADIF is stored in vector AD by columns.  The arrays are */
/*         related by */

/*           ADIF(I,J) = AD(I-J+1 + (2*N-J+2)*(J-1)/2) */

/*         I = J,N   ,   J=1,NDERIV. */

/*     Description of Arguments */

/*         Input      T,A are double precision */
/*          T       - knot vector of length N+K */
/*          A       - B-spline coefficient vector of length N */
/*          N       - number of B-spline coefficients */
/*                    N = sum of knot multiplicities-K */
/*          K       - order of the spline, K .GE. 1 */
/*          NDERIV  - number of derivatives, 1 .LE. NDERIV .LE. K. */
/*                    NDERIV=1 gives the zero-th derivative = */
/*                    function value */

/*         Output     AD is double precision */
/*          AD      - table of differences in a vector of length */
/*                    (2*N-NDERIV+1)*NDERIV/2 for input to DBSPEV */

/*     Error Conditions */
/*         Improper input is a fatal error */

/* ***REFERENCES  Carl de Boor, Package for calculating with B-splines, */
/*                 SIAM Journal on Numerical Analysis 14, 3 (June 1977), */
/*                 pp. 441-472. */
/* ***ROUTINES CALLED  XERMSG */
/* ***REVISION HISTORY  (YYMMDD) */
/*   800901  DATE WRITTEN */
/*   890831  Modified array declarations.  (WRB) */
/*   890911  Removed unnecessary intrinsics.  (WRB) */
/*   890911  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ) */
/*   920501  Reformatted the REFERENCES section.  (WRB) */
/* ***END PROLOGUE  DBSPDR */


/*     DIMENSION T(N+K), AD((2*N-NDERIV+1)*NDERIV/2) */
/* ***FIRST EXECUTABLE STATEMENT  DBSPDR */
    /* Parameter adjustments */
    --ad;
    --a;
    --t;

    /* Function Body */
    if (*k < 1) {
	goto L100;
    }
    if (*n < *k) {
	goto L105;
    }
    if (*nderiv < 1 || *nderiv > *k) {
	goto L110;
    }
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ad[i__] = a[i__];
/* L10: */
    }
    if (*nderiv == 1) {
	return 0;
    }
    kmid = *k;
    jj = *n;
    jm = 0;
    i__1 = *nderiv;
    for (id = 2; id <= i__1; ++id) {
	--kmid;
	fkmid = (doublereal) kmid;
	ii = 1;
	i__2 = *n;
	for (i__ = id; i__ <= i__2; ++i__) {
	    ipkmid = i__ + kmid;
	    diff = t[ipkmid] - t[i__];
	    if (diff != 0.) {
		ad[ii + jj] = (ad[ii + jm + 1] - ad[ii + jm]) / diff * fkmid;
	    }
	    ++ii;
/* L20: */
	}
	jm = jj;
	jj = jj + *n - id + 1;
/* L30: */
    }
    return 0;


L100:
    xermsg_("SLATEC", "DBSPDR", "K DOES NOT SATISFY K.GE.1", &c__2, &c__1, (
	    ftnlen)6, (ftnlen)6, (ftnlen)25);
    return 0;
L105:
    xermsg_("SLATEC", "DBSPDR", "N DOES NOT SATISFY N.GE.K", &c__2, &c__1, (
	    ftnlen)6, (ftnlen)6, (ftnlen)25);
    return 0;
L110:
    xermsg_("SLATEC", "DBSPDR", "NDERIV DOES NOT SATISFY 1.LE.NDERIV.LE.K", &
	    c__2, &c__1, (ftnlen)6, (ftnlen)6, (ftnlen)40);
    return 0;
} /* dbspdr_ */

