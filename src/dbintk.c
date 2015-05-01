/* dbintk.f -- translated by f2c (version 12.02.01).
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

static integer c__1 = 1;
static integer c__2 = 2;
static integer c__8 = 8;

/* DECK DBINTK */
/* Subroutine */ int dbintk_(doublereal *x, doublereal *y, doublereal *t, 
	integer *n, integer *k, doublereal *bcoef, doublereal *q, doublereal *
	work)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer i__, j, jj;
    static doublereal xi;
    static integer km1, np1, left, lenq, kpkm2, iflag, iwork, ilp1mx;
    extern /* Subroutine */ int dbnfac_(doublereal *, integer *, integer *, 
	    integer *, integer *, integer *), dbnslv_(doublereal *, integer *,
	     integer *, integer *, integer *, doublereal *), dbspvn_(
	    doublereal *, integer *, integer *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *), xermsg_(char *,
	     char *, char *, integer *, integer *, ftnlen, ftnlen, ftnlen);

/* ***BEGIN PROLOGUE  DBINTK */
/* ***PURPOSE  Compute the B-representation of a spline which interpolates */
/*            given data. */
/* ***LIBRARY   SLATEC */
/* ***CATEGORY  E1A */
/* ***TYPE      DOUBLE PRECISION (BINTK-S, DBINTK-D) */
/* ***KEYWORDS  B-SPLINE, DATA FITTING, INTERPOLATION */
/* ***AUTHOR  Amos, D. E., (SNLA) */
/* ***DESCRIPTION */

/*     Written by Carl de Boor and modified by D. E. Amos */

/*     Abstract    **** a double precision routine **** */

/*         DBINTK is the SPLINT routine of the reference. */

/*         DBINTK produces the B-spline coefficients, BCOEF, of the */
/*         B-spline of order K with knots T(I), I=1,...,N+K, which */
/*         takes on the value Y(I) at X(I), I=1,...,N.  The spline or */
/*         any of its derivatives can be evaluated by calls to DBVALU. */

/*         The I-th equation of the linear system A*BCOEF = B for the */
/*         coefficients of the interpolant enforces interpolation at */
/*         X(I), I=1,...,N.  Hence, B(I) = Y(I), for all I, and A is */
/*         a band matrix with 2K-1 bands if A is invertible.  The matrix */
/*         A is generated row by row and stored, diagonal by diagonal, */
/*         in the rows of Q, with the main diagonal going into row K. */
/*         The banded system is then solved by a call to DBNFAC (which */
/*         constructs the triangular factorization for A and stores it */
/*         again in Q), followed by a call to DBNSLV (which then */
/*         obtains the solution BCOEF by substitution).  DBNFAC does no */
/*         pivoting, since the total positivity of the matrix A makes */
/*         this unnecessary.  The linear system to be solved is */
/*         (theoretically) invertible if and only if */
/*                 T(I) .LT. X(I) .LT. T(I+K),        for all I. */
/*         Equality is permitted on the left for I=1 and on the right */
/*         for I=N when K knots are used at X(1) or X(N).  Otherwise, */
/*         violation of this condition is certain to lead to an error. */

/*     Description of Arguments */

/*         Input       X,Y,T are double precision */
/*           X       - vector of length N containing data point abscissa */
/*                     in strictly increasing order. */
/*           Y       - corresponding vector of length N containing data */
/*                     point ordinates. */
/*           T       - knot vector of length N+K */
/*                     Since T(1),..,T(K) .LE. X(1) and T(N+1),..,T(N+K) */
/*                     .GE. X(N), this leaves only N-K knots (not nec- */
/*                     essarily X(I) values) interior to (X(1),X(N)) */
/*           N       - number of data points, N .GE. K */
/*           K       - order of the spline, K .GE. 1 */

/*         Output      BCOEF,Q,WORK are double precision */
/*           BCOEF   - a vector of length N containing the B-spline */
/*                     coefficients */
/*           Q       - a work vector of length (2*K-1)*N, containing */
/*                     the triangular factorization of the coefficient */
/*                     matrix of the linear system being solved.  The */
/*                     coefficients for the interpolant of an */
/*                     additional data set (X(I),YY(I)), I=1,...,N */
/*                     with the same abscissa can be obtained by loading */
/*                     YY into BCOEF and then executing */
/*                         CALL DBNSLV (Q,2K-1,N,K-1,K-1,BCOEF) */
/*           WORK    - work vector of length 2*K */

/*     Error Conditions */
/*         Improper input is a fatal error */
/*         Singular system of equations is a fatal error */

/* ***REFERENCES  D. E. Amos, Computation with splines and B-splines, */
/*                 Report SAND78-1968, Sandia Laboratories, March 1979. */
/*               Carl de Boor, Package for calculating with B-splines, */
/*                 SIAM Journal on Numerical Analysis 14, 3 (June 1977), */
/*                 pp. 441-472. */
/*               Carl de Boor, A Practical Guide to Splines, Applied */
/*                 Mathematics Series 27, Springer-Verlag, New York, */
/*                 1978. */
/* ***ROUTINES CALLED  DBNFAC, DBNSLV, DBSPVN, XERMSG */
/* ***REVISION HISTORY  (YYMMDD) */
/*   800901  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890831  Modified array declarations.  (WRB) */
/*   890831  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ) */
/*   900326  Removed duplicate information from DESCRIPTION section. */
/*           (WRB) */
/*   920501  Reformatted the REFERENCES section.  (WRB) */
/* ***END PROLOGUE  DBINTK */

/*     DIMENSION Q(2*K-1,N), T(N+K) */
/* ***FIRST EXECUTABLE STATEMENT  DBINTK */
    /* Parameter adjustments */
    --work;
    --q;
    --bcoef;
    --t;
    --y;
    --x;

    /* Function Body */
    if (*k < 1) {
	goto L100;
    }
    if (*n < *k) {
	goto L105;
    }
    jj = *n - 1;
    if (jj == 0) {
	goto L6;
    }
    i__1 = jj;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (x[i__] >= x[i__ + 1]) {
	    goto L110;
	}
/* L5: */
    }
L6:
    np1 = *n + 1;
    km1 = *k - 1;
    kpkm2 = km1 << 1;
    left = *k;
/*                ZERO OUT ALL ENTRIES OF Q */
    lenq = *n * (*k + km1);
    i__1 = lenq;
    for (i__ = 1; i__ <= i__1; ++i__) {
	q[i__] = 0.;
/* L10: */
    }

/*  ***   LOOP OVER I TO CONSTRUCT THE  N  INTERPOLATION EQUATIONS */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	xi = x[i__];
/* Computing MIN */
	i__2 = i__ + *k;
	ilp1mx = min(i__2,np1);
/*        *** FIND  LEFT  IN THE CLOSED INTERVAL (I,I+K-1) SUCH THAT */
/*                T(LEFT) .LE. X(I) .LT. T(LEFT+1) */
/*        MATRIX IS SINGULAR IF THIS IS NOT POSSIBLE */
	left = max(left,i__);
	if (xi < t[left]) {
	    goto L80;
	}
L20:
	if (xi < t[left + 1]) {
	    goto L30;
	}
	++left;
	if (left < ilp1mx) {
	    goto L20;
	}
	--left;
	if (xi > t[left + 1]) {
	    goto L80;
	}
/*        *** THE I-TH EQUATION ENFORCES INTERPOLATION AT XI, HENCE */
/*        A(I,J) = B(J,K,T)(XI), ALL J. ONLY THE  K  ENTRIES WITH  J = */
/*        LEFT-K+1,...,LEFT ACTUALLY MIGHT BE NONZERO. THESE  K  NUMBERS */
/*        ARE RETURNED, IN  BCOEF (USED FOR TEMP. STORAGE HERE), BY THE */
/*        FOLLOWING */
L30:
	dbspvn_(&t[1], k, k, &c__1, &xi, &left, &bcoef[1], &work[1], &iwork);
/*        WE THEREFORE WANT  BCOEF(J) = B(LEFT-K+J)(XI) TO GO INTO */
/*        A(I,LEFT-K+J), I.E., INTO  Q(I-(LEFT+J)+2*K,(LEFT+J)-K) SINCE */
/*        A(I+J,J)  IS TO GO INTO  Q(I+K,J), ALL I,J,  IF WE CONSIDER  Q */
/*        AS A TWO-DIM. ARRAY , WITH  2*K-1  ROWS (SEE COMMENTS IN */
/*        DBNFAC). IN THE PRESENT PROGRAM, WE TREAT  Q  AS AN EQUIVALENT */
/*        ONE-DIMENSIONAL ARRAY (BECAUSE OF FORTRAN RESTRICTIONS ON */
/*        DIMENSION STATEMENTS) . WE THEREFORE WANT  BCOEF(J) TO GO INTO */
/*        ENTRY */
/*            I -(LEFT+J) + 2*K + ((LEFT+J) - K-1)*(2*K-1) */
/*                   =  I-LEFT+1 + (LEFT -K)*(2*K-1) + (2*K-2)*J */
/*        OF  Q . */
	jj = i__ - left + 1 + (left - *k) * (*k + km1);
	i__2 = *k;
	for (j = 1; j <= i__2; ++j) {
	    jj += kpkm2;
	    q[jj] = bcoef[j];
/* L40: */
	}
/* L50: */
    }

/*     ***OBTAIN FACTORIZATION OF  A  , STORED AGAIN IN  Q. */
    i__1 = *k + km1;
    dbnfac_(&q[1], &i__1, n, &km1, &km1, &iflag);
    switch (iflag) {
	case 1:  goto L60;
	case 2:  goto L90;
    }
/*     *** SOLVE  A*BCOEF = Y  BY BACKSUBSTITUTION */
L60:
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	bcoef[i__] = y[i__];
/* L70: */
    }
    i__1 = *k + km1;
    dbnslv_(&q[1], &i__1, n, &km1, &km1, &bcoef[1]);
    return 0;


L80:
    xermsg_("SLATEC", "DBINTK", "SOME ABSCISSA WAS NOT IN THE SUPPORT OF THE"
	    " CORRESPONDING BASIS FUNCTION AND THE SYSTEM IS SINGULAR.", &c__2,
	     &c__1, (ftnlen)6, (ftnlen)6, (ftnlen)100);
    return 0;
L90:
    xermsg_("SLATEC", "DBINTK", "THE SYSTEM OF SOLVER DETECTS A SINGULAR SYS"
	    "TEM ALTHOUGH THE THEORETICAL CONDITIONS FOR A SOLUTION WERE SATI"
	    "SFIED.", &c__8, &c__1, (ftnlen)6, (ftnlen)6, (ftnlen)113);
    return 0;
L100:
    xermsg_("SLATEC", "DBINTK", "K DOES NOT SATISFY K.GE.1", &c__2, &c__1, (
	    ftnlen)6, (ftnlen)6, (ftnlen)25);
    return 0;
L105:
    xermsg_("SLATEC", "DBINTK", "N DOES NOT SATISFY N.GE.K", &c__2, &c__1, (
	    ftnlen)6, (ftnlen)6, (ftnlen)25);
    return 0;
L110:
    xermsg_("SLATEC", "DBINTK", "X(I) DOES NOT SATISFY X(I).LT.X(I+1) FOR SO"
	    "ME I", &c__2, &c__1, (ftnlen)6, (ftnlen)6, (ftnlen)47);
    return 0;
} /* dbintk_ */

