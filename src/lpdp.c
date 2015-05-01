/* lpdp.f -- translated by f2c (version 12.02.01).
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

static integer c__0 = 0;
static integer c__1 = 1;

/* DECK LPDP */
/* Subroutine */ int lpdp_(real *a, integer *mda, integer *m, integer *n1, 
	integer *n2, real *prgopt, real *x, real *wnorm, integer *mode, real *
	ws, integer *is)
{
    /* Initialized data */

    static real zero = 0.f;
    static real one = 1.f;
    static real fac = .1f;

    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;

    /* Local variables */
    static integer i__, j, l, n;
    static real sc;
    static integer iw, ix, np1;
    extern doublereal sdot_(integer *, real *, integer *, real *, integer *), 
	    snrm2_(integer *, real *, integer *);
    extern /* Subroutine */ int sscal_(integer *, real *, real *, integer *);
    static integer modew;
    static real rnorm;
    extern /* Subroutine */ int scopy_(integer *, real *, integer *, real *, 
	    integer *), wnnls_(real *, integer *, integer *, integer *, 
	    integer *, integer *, real *, real *, real *, integer *, integer *
	    , real *);
    static real ynorm;

/* ***BEGIN PROLOGUE  LPDP */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to LSEI */
/* ***LIBRARY   SLATEC */
/* ***TYPE      SINGLE PRECISION (LPDP-S, DLPDP-D) */
/* ***AUTHOR  Hanson, R. J., (SNLA) */
/*           Haskell, K. H., (SNLA) */
/* ***DESCRIPTION */

/*     DIMENSION A(MDA,N+1),PRGOPT(*),X(N),WS((M+2)*(N+7)),IS(M+N+1), */
/*     where N=N1+N2.  This is a slight overestimate for WS(*). */

/*     Determine an N1-vector W, and */
/*               an N2-vector Z */
/*     which minimizes the Euclidean length of W */
/*     subject to G*W+H*Z .GE. Y. */
/*     This is the least projected distance problem, LPDP. */
/*     The matrices G and H are of respective */
/*     dimensions M by N1 and M by N2. */

/*     Called by subprogram LSI( ). */

/*     The matrix */
/*                (G H Y) */

/*     occupies rows 1,...,M and cols 1,...,N1+N2+1 of A(*,*). */

/*     The solution (W) is returned in X(*). */
/*                  (Z) */

/*     The value of MODE indicates the status of */
/*     the computation after returning to the user. */

/*          MODE=1  The solution was successfully obtained. */

/*          MODE=2  The inequalities are inconsistent. */

/* ***SEE ALSO  LSEI */
/* ***ROUTINES CALLED  SCOPY, SDOT, SNRM2, SSCAL, WNNLS */
/* ***REVISION HISTORY  (YYMMDD) */
/*   790701  DATE WRITTEN */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900328  Added TYPE section.  (WRB) */
/*   910408  Updated the AUTHOR section.  (WRB) */
/* ***END PROLOGUE  LPDP */

/*     SUBROUTINES CALLED */

/*     WNNLS         SOLVES A NONNEGATIVELY CONSTRAINED LINEAR LEAST */
/*                   SQUARES PROBLEM WITH LINEAR EQUALITY CONSTRAINTS. */
/*                   PART OF THIS PACKAGE. */

/* ++ */
/*     SDOT,         SUBROUTINES FROM THE BLAS PACKAGE. */
/*     SSCAL,SNRM2,  SEE TRANS. MATH. SOFT., VOL. 5, NO. 3, P. 308. */
/*     SCOPY */

    /* Parameter adjustments */
    a_dim1 = *mda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --prgopt;
    --x;
    --ws;
    --is;

    /* Function Body */
/* ***FIRST EXECUTABLE STATEMENT  LPDP */
    n = *n1 + *n2;
    *mode = 1;
    if (! (*m <= 0)) {
	goto L20;
    }
    if (! (n > 0)) {
	goto L10;
    }
    x[1] = zero;
    scopy_(&n, &x[1], &c__0, &x[1], &c__1);
L10:
    *wnorm = zero;
    return 0;
L20:
    np1 = n + 1;

/*     SCALE NONZERO ROWS OF INEQUALITY MATRIX TO HAVE LENGTH ONE. */
    i__1 = *m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	sc = snrm2_(&n, &a[i__ + a_dim1], mda);
	if (! (sc != zero)) {
	    goto L30;
	}
	sc = one / sc;
	sscal_(&np1, &sc, &a[i__ + a_dim1], mda);
L30:
/* L40: */
	;
    }

/*     SCALE RT.-SIDE VECTOR TO HAVE LENGTH ONE (OR ZERO). */
    ynorm = snrm2_(m, &a[np1 * a_dim1 + 1], &c__1);
    if (! (ynorm != zero)) {
	goto L50;
    }
    sc = one / ynorm;
    sscal_(m, &sc, &a[np1 * a_dim1 + 1], &c__1);

/*     SCALE COLS OF MATRIX H. */
L50:
    j = *n1 + 1;
L60:
    if (! (j <= n)) {
	goto L70;
    }
    sc = snrm2_(m, &a[j * a_dim1 + 1], &c__1);
    if (sc != zero) {
	sc = one / sc;
    }
    sscal_(m, &sc, &a[j * a_dim1 + 1], &c__1);
    x[j] = sc;
    ++j;
    goto L60;
L70:
    if (! (*n1 > 0)) {
	goto L130;
    }

/*     COPY TRANSPOSE OF (H G Y) TO WORK ARRAY WS(*). */
    iw = 0;
    i__1 = *m;
    for (i__ = 1; i__ <= i__1; ++i__) {

/*     MOVE COL OF TRANSPOSE OF H INTO WORK ARRAY. */
	scopy_(n2, &a[i__ + (*n1 + 1) * a_dim1], mda, &ws[iw + 1], &c__1);
	iw += *n2;

/*     MOVE COL OF TRANSPOSE OF G INTO WORK ARRAY. */
	scopy_(n1, &a[i__ + a_dim1], mda, &ws[iw + 1], &c__1);
	iw += *n1;

/*     MOVE COMPONENT OF VECTOR Y INTO WORK ARRAY. */
	ws[iw + 1] = a[i__ + np1 * a_dim1];
	++iw;
/* L80: */
    }
    ws[iw + 1] = zero;
    scopy_(&n, &ws[iw + 1], &c__0, &ws[iw + 1], &c__1);
    iw += n;
    ws[iw + 1] = one;
    ++iw;

/*     SOLVE EU=F SUBJECT TO (TRANSPOSE OF H)U=0, U.GE.0.  THE */
/*     MATRIX E = TRANSPOSE OF (G Y), AND THE (N+1)-VECTOR */
/*     F = TRANSPOSE OF (0,...,0,1). */
    ix = iw + 1;
    iw += *m;

/*     DO NOT CHECK LENGTHS OF WORK ARRAYS IN THIS USAGE OF WNNLS( ). */
    is[1] = 0;
    is[2] = 0;
    i__1 = np1 - *n2;
    wnnls_(&ws[1], &np1, n2, &i__1, m, &c__0, &prgopt[1], &ws[ix], &rnorm, &
	    modew, &is[1], &ws[iw + 1]);

/*     COMPUTE THE COMPONENTS OF THE SOLN DENOTED ABOVE BY W. */
    sc = one - sdot_(m, &a[np1 * a_dim1 + 1], &c__1, &ws[ix], &c__1);
    if (! (one + fac * dabs(sc) != one && rnorm > zero)) {
	goto L110;
    }
    sc = one / sc;
    i__1 = *n1;
    for (j = 1; j <= i__1; ++j) {
	x[j] = sc * sdot_(m, &a[j * a_dim1 + 1], &c__1, &ws[ix], &c__1);
/* L90: */
    }

/*     COMPUTE THE VECTOR Q=Y-GW.  OVERWRITE Y WITH THIS VECTOR. */
    i__1 = *m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	a[i__ + np1 * a_dim1] -= sdot_(n1, &a[i__ + a_dim1], mda, &x[1], &
		c__1);
/* L100: */
    }
    goto L120;
L110:
    *mode = 2;
    return 0;
L120:
L130:
    if (! (*n2 > 0)) {
	goto L180;
    }

/*     COPY TRANSPOSE OF (H Q) TO WORK ARRAY WS(*). */
    iw = 0;
    i__1 = *m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	scopy_(n2, &a[i__ + (*n1 + 1) * a_dim1], mda, &ws[iw + 1], &c__1);
	iw += *n2;
	ws[iw + 1] = a[i__ + np1 * a_dim1];
	++iw;
/* L140: */
    }
    ws[iw + 1] = zero;
    scopy_(n2, &ws[iw + 1], &c__0, &ws[iw + 1], &c__1);
    iw += *n2;
    ws[iw + 1] = one;
    ++iw;
    ix = iw + 1;
    iw += *m;

/*     SOLVE RV=S SUBJECT TO V.GE.0.  THE MATRIX R =(TRANSPOSE */
/*     OF (H Q)), WHERE Q=Y-GW.  THE (N2+1)-VECTOR S =(TRANSPOSE */
/*     OF (0,...,0,1)). */

/*     DO NOT CHECK LENGTHS OF WORK ARRAYS IN THIS USAGE OF WNNLS( ). */
    is[1] = 0;
    is[2] = 0;
    i__1 = *n2 + 1;
    i__2 = *n2 + 1;
    wnnls_(&ws[1], &i__1, &c__0, &i__2, m, &c__0, &prgopt[1], &ws[ix], &rnorm,
	     &modew, &is[1], &ws[iw + 1]);

/*     COMPUTE THE COMPONENTS OF THE SOLN DENOTED ABOVE BY Z. */
    sc = one - sdot_(m, &a[np1 * a_dim1 + 1], &c__1, &ws[ix], &c__1);
    if (! (one + fac * dabs(sc) != one && rnorm > zero)) {
	goto L160;
    }
    sc = one / sc;
    i__1 = *n2;
    for (j = 1; j <= i__1; ++j) {
	l = *n1 + j;
	x[l] = sc * sdot_(m, &a[l * a_dim1 + 1], &c__1, &ws[ix], &c__1) * x[l]
		;
/* L150: */
    }
    goto L170;
L160:
    *mode = 2;
    return 0;
L170:

/*     ACCOUNT FOR SCALING OF RT.-SIDE VECTOR IN SOLUTION. */
L180:
    sscal_(&n, &ynorm, &x[1], &c__1);
    *wnorm = snrm2_(n1, &x[1], &c__1);
    return 0;
} /* lpdp_ */

