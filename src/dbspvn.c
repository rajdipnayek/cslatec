/* dbspvn.f -- translated by f2c (version 12.02.01).
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

/* DECK DBSPVN */
/* Subroutine */ int dbspvn_(doublereal *t, integer *jhigh, integer *k, 
	integer *index, doublereal *x, integer *ileft, doublereal *vnikx, 
	doublereal *work, integer *iwork)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer l;
    static doublereal vm;
    static integer jp1, ipj, imjp1, jp1ml;
    extern /* Subroutine */ int xermsg_(char *, char *, char *, integer *, 
	    integer *, ftnlen, ftnlen, ftnlen);
    static doublereal vmprev;

/* ***BEGIN PROLOGUE  DBSPVN */
/* ***PURPOSE  Calculate the value of all (possibly) nonzero basis */
/*            functions at X. */
/* ***LIBRARY   SLATEC */
/* ***CATEGORY  E3, K6 */
/* ***TYPE      DOUBLE PRECISION (BSPVN-S, DBSPVN-D) */
/* ***KEYWORDS  EVALUATION OF B-SPLINE */
/* ***AUTHOR  Amos, D. E., (SNLA) */
/* ***DESCRIPTION */

/*     Written by Carl de Boor and modified by D. E. Amos */

/*     Abstract    **** a double precision routine **** */
/*         DBSPVN is the BSPLVN routine of the reference. */

/*         DBSPVN calculates the value of all (possibly) nonzero basis */
/*         functions at X of order MAX(JHIGH,(J+1)*(INDEX-1)), where T(K) */
/*         .LE. X .LE. T(N+1) and J=IWORK is set inside the routine on */
/*         the first call when INDEX=1.  ILEFT is such that T(ILEFT) .LE. */
/*         X .LT. T(ILEFT+1).  A call to DINTRV(T,N+1,X,ILO,ILEFT,MFLAG) */
/*         produces the proper ILEFT.  DBSPVN calculates using the basic */
/*         algorithm needed in DBSPVD.  If only basis functions are */
/*         desired, setting JHIGH=K and INDEX=1 can be faster than */
/*         calling DBSPVD, but extra coding is required for derivatives */
/*         (INDEX=2) and DBSPVD is set up for this purpose. */

/*         Left limiting values are set up as described in DBSPVD. */

/*     Description of Arguments */

/*         Input      T,X are double precision */
/*          T       - knot vector of length N+K, where */
/*                    N = number of B-spline basis functions */
/*                    N = sum of knot multiplicities-K */
/*          JHIGH   - order of B-spline, 1 .LE. JHIGH .LE. K */
/*          K       - highest possible order */
/*          INDEX   - INDEX = 1 gives basis functions of order JHIGH */
/*                          = 2 denotes previous entry with work, IWORK */
/*                              values saved for subsequent calls to */
/*                              DBSPVN. */
/*          X       - argument of basis functions, */
/*                    T(K) .LE. X .LE. T(N+1) */
/*          ILEFT   - largest integer such that */
/*                    T(ILEFT) .LE. X .LT.  T(ILEFT+1) */

/*         Output     VNIKX, WORK are double precision */
/*          VNIKX   - vector of length K for spline values. */
/*          WORK    - a work vector of length 2*K */
/*          IWORK   - a work parameter.  Both WORK and IWORK contain */
/*                    information necessary to continue for INDEX = 2. */
/*                    When INDEX = 1 exclusively, these are scratch */
/*                    variables and can be used for other purposes. */

/*     Error Conditions */
/*         Improper input is a fatal error. */

/* ***REFERENCES  Carl de Boor, Package for calculating with B-splines, */
/*                 SIAM Journal on Numerical Analysis 14, 3 (June 1977), */
/*                 pp. 441-472. */
/* ***ROUTINES CALLED  XERMSG */
/* ***REVISION HISTORY  (YYMMDD) */
/*   800901  DATE WRITTEN */
/*   890831  Modified array declarations.  (WRB) */
/*   890831  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ) */
/*   920501  Reformatted the REFERENCES section.  (WRB) */
/* ***END PROLOGUE  DBSPVN */

/*     DIMENSION T(ILEFT+JHIGH) */
/*     CONTENT OF J, DELTAM, DELTAP IS EXPECTED UNCHANGED BETWEEN CALLS. */
/*     WORK(I) = DELTAP(I), WORK(K+I) = DELTAM(I), I = 1,K */
/* ***FIRST EXECUTABLE STATEMENT  DBSPVN */
    /* Parameter adjustments */
    --work;
    --vnikx;
    --t;

    /* Function Body */
    if (*k < 1) {
	goto L90;
    }
    if (*jhigh > *k || *jhigh < 1) {
	goto L100;
    }
    if (*index < 1 || *index > 2) {
	goto L105;
    }
    if (*x < t[*ileft] || *x > t[*ileft + 1]) {
	goto L110;
    }
    switch (*index) {
	case 1:  goto L10;
	case 2:  goto L20;
    }
L10:
    *iwork = 1;
    vnikx[1] = 1.;
    if (*iwork >= *jhigh) {
	goto L40;
    }

L20:
    ipj = *ileft + *iwork;
    work[*iwork] = t[ipj] - *x;
    imjp1 = *ileft - *iwork + 1;
    work[*k + *iwork] = *x - t[imjp1];
    vmprev = 0.;
    jp1 = *iwork + 1;
    i__1 = *iwork;
    for (l = 1; l <= i__1; ++l) {
	jp1ml = jp1 - l;
	vm = vnikx[l] / (work[l] + work[*k + jp1ml]);
	vnikx[l] = vm * work[l] + vmprev;
	vmprev = vm * work[*k + jp1ml];
/* L30: */
    }
    vnikx[jp1] = vmprev;
    *iwork = jp1;
    if (*iwork < *jhigh) {
	goto L20;
    }

L40:
    return 0;


L90:
    xermsg_("SLATEC", "DBSPVN", "K DOES NOT SATISFY K.GE.1", &c__2, &c__1, (
	    ftnlen)6, (ftnlen)6, (ftnlen)25);
    return 0;
L100:
    xermsg_("SLATEC", "DBSPVN", "JHIGH DOES NOT SATISFY 1.LE.JHIGH.LE.K", &
	    c__2, &c__1, (ftnlen)6, (ftnlen)6, (ftnlen)38);
    return 0;
L105:
    xermsg_("SLATEC", "DBSPVN", "INDEX IS NOT 1 OR 2", &c__2, &c__1, (ftnlen)
	    6, (ftnlen)6, (ftnlen)19);
    return 0;
L110:
    xermsg_("SLATEC", "DBSPVN", "X DOES NOT SATISFY T(ILEFT).LE.X.LE.T(ILEFT"
	    "+1)", &c__2, &c__1, (ftnlen)6, (ftnlen)6, (ftnlen)46);
    return 0;
} /* dbspvn_ */

