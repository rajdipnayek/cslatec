/* dbint4.f -- translated by f2c (version 12.02.01).
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
static integer c__2 = 2;

/* DECK DBINT4 */
/* Subroutine */ int dbint4_(doublereal *x, doublereal *y, integer *ndata, 
	integer *ibcl, integer *ibcr, doublereal *fbcl, doublereal *fbcr, 
	integer *kntopt, doublereal *t, doublereal *bcoef, integer *n, 
	integer *k, doublereal *w)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer i__, j, it, np, iw, jw;
    static doublereal xl, tx1;
    static integer ilb, ndm, iub;
    static doublereal tol;
    static integer iwp;
    static doublereal txn, work[15];
    static integer iflag, ileft;
    static doublereal wdtol;
    extern doublereal d1mach_(integer *);
    static doublereal vnikx[16]	/* was [4][4] */;
    static integer nwrow;
    extern /* Subroutine */ int dbnfac_(doublereal *, integer *, integer *, 
	    integer *, integer *, integer *), dbspvd_(doublereal *, integer *,
	     integer *, doublereal *, integer *, integer *, doublereal *, 
	    doublereal *), dbnslv_(doublereal *, integer *, integer *, 
	    integer *, integer *, doublereal *), xermsg_(char *, char *, char 
	    *, integer *, integer *, ftnlen, ftnlen, ftnlen);

/* ***BEGIN PROLOGUE  DBINT4 */
/* ***PURPOSE  Compute the B-representation of a cubic spline */
/*            which interpolates given data. */
/* ***LIBRARY   SLATEC */
/* ***CATEGORY  E1A */
/* ***TYPE      DOUBLE PRECISION (BINT4-S, DBINT4-D) */
/* ***KEYWORDS  B-SPLINE, CUBIC SPLINES, DATA FITTING, INTERPOLATION */
/* ***AUTHOR  Amos, D. E., (SNLA) */
/* ***DESCRIPTION */

/*     Abstract    **** a double precision routine **** */

/*         DBINT4 computes the B representation (T,BCOEF,N,K) of a */
/*         cubic spline (K=4) which interpolates data (X(I),Y(I)), */
/*         I=1,NDATA.  Parameters IBCL, IBCR, FBCL, FBCR allow the */
/*         specification of the spline first or second derivative at */
/*         both X(1) and X(NDATA).  When this data is not specified */
/*         by the problem, it is common practice to use a natural */
/*         spline by setting second derivatives at X(1) and X(NDATA) */
/*         to zero (IBCL=IBCR=2,FBCL=FBCR=0.0).  The spline is defined */
/*         on T(4) .LE. X .LE. T(N+1) with (ordered) interior knots at */
/*         X(I) values where N=NDATA+2.  The knots T(1),T(2),T(3) lie to */
/*         the left of T(4)=X(1) and the knots T(N+2), T(N+3), T(N+4) */
/*         lie to the right of T(N+1)=X(NDATA) in increasing order.  If */
/*         no extrapolation outside (X(1),X(NDATA)) is anticipated, the */
/*         knots T(1)=T(2)=T(3)=T(4)=X(1) and T(N+2)=T(N+3)=T(N+4)= */
/*         T(N+1)=X(NDATA) can be specified by KNTOPT=1.  KNTOPT=2 */
/*         selects a knot placement for T(1), T(2), T(3) to make the */
/*         first 7 knots symmetric about T(4)=X(1) and similarly for */
/*         T(N+2), T(N+3), T(N+4) about T(N+1)=X(NDATA).  KNTOPT=3 */
/*         allows the user to make his own selection, in increasing */
/*         order, for T(1), T(2), T(3) to the left of X(1) and T(N+2), */
/*         T(N+3), T(N+4) to the right of X(NDATA) in the work array */
/*         W(1) through W(6).  In any case, the interpolation on */
/*         T(4) .LE. X .LE. T(N+1) by using function DBVALU is unique */
/*         for given boundary conditions. */

/*     Description of Arguments */

/*         Input      X,Y,FBCL,FBCR,W are double precision */
/*           X      - X vector of abscissae of length NDATA, distinct */
/*                    and in increasing order */
/*           Y      - Y vector of ordinates of length NDATA */
/*           NDATA  - number of data points, NDATA .GE. 2 */
/*           IBCL   - selection parameter for left boundary condition */
/*                    IBCL = 1 constrain the first derivative at */
/*                             X(1) to FBCL */
/*                         = 2 constrain the second derivative at */
/*                             X(1) to FBCL */
/*           IBCR   - selection parameter for right boundary condition */
/*                    IBCR = 1 constrain first derivative at */
/*                             X(NDATA) to FBCR */
/*                    IBCR = 2 constrain second derivative at */
/*                             X(NDATA) to FBCR */
/*           FBCL   - left boundary values governed by IBCL */
/*           FBCR   - right boundary values governed by IBCR */
/*           KNTOPT - knot selection parameter */
/*                    KNTOPT = 1 sets knot multiplicity at T(4) and */
/*                               T(N+1) to 4 */
/*                           = 2 sets a symmetric placement of knots */
/*                               about T(4) and T(N+1) */
/*                           = 3 sets T(I)=W(I) and T(N+1+I)=W(3+I),I=1,3 */
/*                               where W(I),I=1,6 is supplied by the user */
/*           W      - work array of dimension at least 5*(NDATA+2) */
/*                    If KNTOPT=3, then W(1),W(2),W(3) are knot values to */
/*                    the left of X(1) and W(4),W(5),W(6) are knot */
/*                    values to the right of X(NDATA) in increasing */
/*                    order to be supplied by the user */

/*         Output     T,BCOEF are double precision */
/*           T      - knot array of length N+4 */
/*           BCOEF  - B spline coefficient array of length N */
/*           N      - number of coefficients, N=NDATA+2 */
/*           K      - order of spline, K=4 */

/*     Error Conditions */
/*         Improper  input is a fatal error */
/*         Singular system of equations is a fatal error */

/* ***REFERENCES  D. E. Amos, Computation with splines and B-splines, */
/*                 Report SAND78-1968, Sandia Laboratories, March 1979. */
/*               Carl de Boor, Package for calculating with B-splines, */
/*                 SIAM Journal on Numerical Analysis 14, 3 (June 1977), */
/*                 pp. 441-472. */
/*               Carl de Boor, A Practical Guide to Splines, Applied */
/*                 Mathematics Series 27, Springer-Verlag, New York, */
/*                 1978. */
/* ***ROUTINES CALLED  D1MACH, DBNFAC, DBNSLV, DBSPVD, XERMSG */
/* ***REVISION HISTORY  (YYMMDD) */
/*   800901  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890531  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ) */
/*   900326  Removed duplicate information from DESCRIPTION section. */
/*           (WRB) */
/*   920501  Reformatted the REFERENCES section.  (WRB) */
/* ***END PROLOGUE  DBINT4 */

/* ***FIRST EXECUTABLE STATEMENT  DBINT4 */
    /* Parameter adjustments */
    w -= 6;
    --bcoef;
    --t;
    --y;
    --x;

    /* Function Body */
    wdtol = d1mach_(&c__4);
    tol = sqrt(wdtol);
    if (*ndata < 2) {
	goto L200;
    }
    ndm = *ndata - 1;
    i__1 = ndm;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (x[i__] >= x[i__ + 1]) {
	    goto L210;
	}
/* L10: */
    }
    if (*ibcl < 1 || *ibcl > 2) {
	goto L220;
    }
    if (*ibcr < 1 || *ibcr > 2) {
	goto L230;
    }
    if (*kntopt < 1 || *kntopt > 3) {
	goto L240;
    }
    *k = 4;
    *n = *ndata + 2;
    np = *n + 1;
    i__1 = *ndata;
    for (i__ = 1; i__ <= i__1; ++i__) {
	t[i__ + 3] = x[i__];
/* L20: */
    }
    switch (*kntopt) {
	case 1:  goto L30;
	case 2:  goto L50;
	case 3:  goto L90;
    }
/*     SET UP KNOT ARRAY WITH MULTIPLICITY 4 AT X(1) AND X(NDATA) */
L30:
    for (i__ = 1; i__ <= 3; ++i__) {
	t[4 - i__] = x[1];
	t[np + i__] = x[*ndata];
/* L40: */
    }
    goto L110;
/*     SET UP KNOT ARRAY WITH SYMMETRIC PLACEMENT ABOUT END POINTS */
L50:
    if (*ndata > 3) {
	goto L70;
    }
    xl = (x[*ndata] - x[1]) / 3.;
    for (i__ = 1; i__ <= 3; ++i__) {
	t[4 - i__] = t[5 - i__] - xl;
	t[np + i__] = t[np + i__ - 1] + xl;
/* L60: */
    }
    goto L110;
L70:
    tx1 = x[1] + x[1];
    txn = x[*ndata] + x[*ndata];
    for (i__ = 1; i__ <= 3; ++i__) {
	t[4 - i__] = tx1 - x[i__ + 1];
	t[np + i__] = txn - x[*ndata - i__];
/* L80: */
    }
    goto L110;
/*     SET UP KNOT ARRAY LESS THAN X(1) AND GREATER THAN X(NDATA) TO BE */
/*     SUPPLIED BY USER IN WORK LOCATIONS W(1) THROUGH W(6) WHEN KNTOPT=3 */
L90:
    for (i__ = 1; i__ <= 3; ++i__) {
	t[4 - i__] = w[4 - i__ + 5];
/* Computing MAX */
	i__1 = 1, i__2 = i__ - 1;
	jw = max(i__1,i__2);
	iw = (i__ + 2) % 5 + 1;
	t[np + i__] = w[iw + jw * 5];
	if (t[4 - i__] > t[5 - i__]) {
	    goto L250;
	}
	if (t[np + i__] < t[np + i__ - 1]) {
	    goto L250;
	}
/* L100: */
    }
L110:

    for (i__ = 1; i__ <= 5; ++i__) {
	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
	    w[i__ + j * 5] = 0.;
/* L120: */
	}
/* L130: */
    }
/*     SET UP LEFT INTERPOLATION POINT AND LEFT BOUNDARY CONDITION FOR */
/*     RIGHT LIMITS */
    it = *ibcl + 1;
    dbspvd_(&t[1], k, &it, &x[1], k, &c__4, vnikx, work);
    iw = 0;
    if (abs(vnikx[2]) < tol) {
	iw = 1;
    }
    for (j = 1; j <= 3; ++j) {
	w[j + 1 + (4 - j) * 5] = vnikx[4 - j + (it << 2) - 5];
	w[j + (4 - j) * 5] = vnikx[4 - j - 1];
/* L140: */
    }
    bcoef[1] = y[1];
    bcoef[2] = *fbcl;
/*     SET UP INTERPOLATION EQUATIONS FOR POINTS I=2 TO I=NDATA-1 */
    ileft = 4;
    if (ndm < 2) {
	goto L170;
    }
    i__1 = ndm;
    for (i__ = 2; i__ <= i__1; ++i__) {
	++ileft;
	dbspvd_(&t[1], k, &c__1, &x[i__], &ileft, &c__4, vnikx, work);
	for (j = 1; j <= 3; ++j) {
	    w[j + 1 + (i__ + 3 - j) * 5] = vnikx[4 - j - 1];
/* L150: */
	}
	bcoef[i__ + 1] = y[i__];
/* L160: */
    }
/*     SET UP RIGHT INTERPOLATION POINT AND RIGHT BOUNDARY CONDITION FOR */
/*     LEFT LIMITS(ILEFT IS ASSOCIATED WITH T(N)=X(NDATA-1)) */
L170:
    it = *ibcr + 1;
    dbspvd_(&t[1], k, &it, &x[*ndata], &ileft, &c__4, vnikx, work);
    jw = 0;
    if (abs(vnikx[1]) < tol) {
	jw = 1;
    }
    for (j = 1; j <= 3; ++j) {
	w[j + 1 + (*ndata + 3 - j) * 5] = vnikx[5 - j + (it << 2) - 5];
	w[j + 2 + (*ndata + 3 - j) * 5] = vnikx[5 - j - 1];
/* L180: */
    }
    bcoef[*n - 1] = *fbcr;
    bcoef[*n] = y[*ndata];
/*     SOLVE SYSTEM OF EQUATIONS */
    ilb = 2 - jw;
    iub = 2 - iw;
    nwrow = 5;
    iwp = iw + 1;
    dbnfac_(&w[iwp + 5], &nwrow, n, &ilb, &iub, &iflag);
    if (iflag == 2) {
	goto L190;
    }
    dbnslv_(&w[iwp + 5], &nwrow, n, &ilb, &iub, &bcoef[1]);
    return 0;


L190:
    xermsg_("SLATEC", "DBINT4", "THE SYSTEM OF EQUATIONS IS SINGULAR", &c__2, 
	    &c__1, (ftnlen)6, (ftnlen)6, (ftnlen)35);
    return 0;
L200:
    xermsg_("SLATEC", "DBINT4", "NDATA IS LESS THAN 2", &c__2, &c__1, (ftnlen)
	    6, (ftnlen)6, (ftnlen)20);
    return 0;
L210:
    xermsg_("SLATEC", "DBINT4", "X VALUES ARE NOT DISTINCT OR NOT ORDERED", &
	    c__2, &c__1, (ftnlen)6, (ftnlen)6, (ftnlen)40);
    return 0;
L220:
    xermsg_("SLATEC", "DBINT4", "IBCL IS NOT 1 OR 2", &c__2, &c__1, (ftnlen)6,
	     (ftnlen)6, (ftnlen)18);
    return 0;
L230:
    xermsg_("SLATEC", "DBINT4", "IBCR IS NOT 1 OR 2", &c__2, &c__1, (ftnlen)6,
	     (ftnlen)6, (ftnlen)18);
    return 0;
L240:
    xermsg_("SLATEC", "DBINT4", "KNTOPT IS NOT 1, 2, OR 3", &c__2, &c__1, (
	    ftnlen)6, (ftnlen)6, (ftnlen)24);
    return 0;
L250:
    xermsg_("SLATEC", "DBINT4", "KNOT INPUT THROUGH W ARRAY IS NOT ORDERED P"
	    "ROPERLY", &c__2, &c__1, (ftnlen)6, (ftnlen)6, (ftnlen)50);
    return 0;
} /* dbint4_ */

