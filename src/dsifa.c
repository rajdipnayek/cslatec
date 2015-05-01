/* dsifa.f -- translated by f2c (version 12.02.01).
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

/* DECK DSIFA */
/* Subroutine */ int dsifa_(doublereal *a, integer *lda, integer *n, integer *
	kpvt, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1;
    doublereal d__1, d__2, d__3;

    /* Local variables */
    static integer j, k;
    static doublereal t, ak, bk;
    static integer jj, km1, km2;
    static doublereal akm1, bkm1;
    static integer imax, jmax;
    static doublereal mulk;
    static logical swap;
    static doublereal alpha, denom;
    extern /* Subroutine */ int dswap_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), daxpy_(integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *);
    static integer kstep, imaxp1;
    static doublereal mulkm1, absakk;
    extern integer idamax_(integer *, doublereal *, integer *);
    static doublereal colmax, rowmax;

/* ***BEGIN PROLOGUE  DSIFA */
/* ***PURPOSE  Factor a real symmetric matrix by elimination with */
/*            symmetric pivoting. */
/* ***LIBRARY   SLATEC (LINPACK) */
/* ***CATEGORY  D2B1A */
/* ***TYPE      DOUBLE PRECISION (SSIFA-S, DSIFA-D, CHIFA-C, CSIFA-C) */
/* ***KEYWORDS  LINEAR ALGEBRA, LINPACK, MATRIX FACTORIZATION, SYMMETRIC */
/* ***AUTHOR  Bunch, J., (UCSD) */
/* ***DESCRIPTION */

/*     DSIFA factors a double precision symmetric matrix by elimination */
/*     with symmetric pivoting. */

/*     To solve  A*X = B , follow DSIFA by DSISL. */
/*     To compute  INVERSE(A)*C , follow DSIFA by DSISL. */
/*     To compute  DETERMINANT(A) , follow DSIFA by DSIDI. */
/*     To compute  INERTIA(A) , follow DSIFA by DSIDI. */
/*     To compute  INVERSE(A) , follow DSIFA by DSIDI. */

/*     On Entry */

/*        A       DOUBLE PRECISION(LDA,N) */
/*                the symmetric matrix to be factored. */
/*                Only the diagonal and upper triangle are used. */

/*        LDA     INTEGER */
/*                the leading dimension of the array  A . */

/*        N       INTEGER */
/*                the order of the matrix  A . */

/*     On Return */

/*        A       a block diagonal matrix and the multipliers which */
/*                were used to obtain it. */
/*                The factorization can be written  A = U*D*TRANS(U) */
/*                where  U  is a product of permutation and unit */
/*                upper triangular matrices, TRANS(U) is the */
/*                transpose of  U , and  D  is block diagonal */
/*                with 1 by 1 and 2 by 2 blocks. */

/*        KPVT    INTEGER(N) */
/*                an integer vector of pivot indices. */

/*        INFO    INTEGER */
/*                = 0  normal value. */
/*                = K  if the K-th pivot block is singular.  This is */
/*                     not an error condition for this subroutine, */
/*                     but it does indicate that DSISL or DSIDI may */
/*                     divide by zero if called. */

/* ***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W. */
/*                 Stewart, LINPACK Users' Guide, SIAM, 1979. */
/* ***ROUTINES CALLED  DAXPY, DSWAP, IDAMAX */
/* ***REVISION HISTORY  (YYMMDD) */
/*   780814  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890831  Modified array declarations.  (WRB) */
/*   891107  Modified routine equivalence list.  (WRB) */
/*   891107  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900326  Removed duplicate information from DESCRIPTION section. */
/*           (WRB) */
/*   920501  Reformatted the REFERENCES section.  (WRB) */
/* ***END PROLOGUE  DSIFA */

/* ***FIRST EXECUTABLE STATEMENT  DSIFA */

/*     INITIALIZE */

/*     ALPHA IS USED IN CHOOSING PIVOT BLOCK SIZE. */

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --kpvt;

    /* Function Body */
    alpha = (sqrt(17.) + 1.) / 8.;

    *info = 0;

/*     MAIN LOOP ON K, WHICH GOES FROM N TO 1. */

    k = *n;
L10:

/*        LEAVE THE LOOP IF K=0 OR K=1. */

    if (k == 0) {
	goto L200;
    }
    if (k > 1) {
	goto L20;
    }
    kpvt[1] = 1;
    if (a[a_dim1 + 1] == 0.) {
	*info = 1;
    }
    goto L200;
L20:

/*        THIS SECTION OF CODE DETERMINES THE KIND OF */
/*        ELIMINATION TO BE PERFORMED.  WHEN IT IS COMPLETED, */
/*        KSTEP WILL BE SET TO THE SIZE OF THE PIVOT BLOCK, AND */
/*        SWAP WILL BE SET TO .TRUE. IF AN INTERCHANGE IS */
/*        REQUIRED. */

    km1 = k - 1;
    absakk = (d__1 = a[k + k * a_dim1], abs(d__1));

/*        DETERMINE THE LARGEST OFF-DIAGONAL ELEMENT IN */
/*        COLUMN K. */

    i__1 = k - 1;
    imax = idamax_(&i__1, &a[k * a_dim1 + 1], &c__1);
    colmax = (d__1 = a[imax + k * a_dim1], abs(d__1));
    if (absakk < alpha * colmax) {
	goto L30;
    }
    kstep = 1;
    swap = FALSE_;
    goto L90;
L30:

/*           DETERMINE THE LARGEST OFF-DIAGONAL ELEMENT IN */
/*           ROW IMAX. */

    rowmax = 0.;
    imaxp1 = imax + 1;
    i__1 = k;
    for (j = imaxp1; j <= i__1; ++j) {
/* Computing MAX */
	d__2 = rowmax, d__3 = (d__1 = a[imax + j * a_dim1], abs(d__1));
	rowmax = max(d__2,d__3);
/* L40: */
    }
    if (imax == 1) {
	goto L50;
    }
    i__1 = imax - 1;
    jmax = idamax_(&i__1, &a[imax * a_dim1 + 1], &c__1);
/* Computing MAX */
    d__2 = rowmax, d__3 = (d__1 = a[jmax + imax * a_dim1], abs(d__1));
    rowmax = max(d__2,d__3);
L50:
    if ((d__1 = a[imax + imax * a_dim1], abs(d__1)) < alpha * rowmax) {
	goto L60;
    }
    kstep = 1;
    swap = TRUE_;
    goto L80;
L60:
    if (absakk < alpha * colmax * (colmax / rowmax)) {
	goto L70;
    }
    kstep = 1;
    swap = FALSE_;
    goto L80;
L70:
    kstep = 2;
    swap = imax != km1;
L80:
L90:
    if (max(absakk,colmax) != 0.) {
	goto L100;
    }

/*           COLUMN K IS ZERO.  SET INFO AND ITERATE THE LOOP. */

    kpvt[k] = k;
    *info = k;
    goto L190;
L100:
    if (kstep == 2) {
	goto L140;
    }

/*           1 X 1 PIVOT BLOCK. */

    if (! swap) {
	goto L120;
    }

/*              PERFORM AN INTERCHANGE. */

    dswap_(&imax, &a[imax * a_dim1 + 1], &c__1, &a[k * a_dim1 + 1], &c__1);
    i__1 = k;
    for (jj = imax; jj <= i__1; ++jj) {
	j = k + imax - jj;
	t = a[j + k * a_dim1];
	a[j + k * a_dim1] = a[imax + j * a_dim1];
	a[imax + j * a_dim1] = t;
/* L110: */
    }
L120:

/*           PERFORM THE ELIMINATION. */

    i__1 = km1;
    for (jj = 1; jj <= i__1; ++jj) {
	j = k - jj;
	mulk = -a[j + k * a_dim1] / a[k + k * a_dim1];
	t = mulk;
	daxpy_(&j, &t, &a[k * a_dim1 + 1], &c__1, &a[j * a_dim1 + 1], &c__1);
	a[j + k * a_dim1] = mulk;
/* L130: */
    }

/*           SET THE PIVOT ARRAY. */

    kpvt[k] = k;
    if (swap) {
	kpvt[k] = imax;
    }
    goto L190;
L140:

/*           2 X 2 PIVOT BLOCK. */

    if (! swap) {
	goto L160;
    }

/*              PERFORM AN INTERCHANGE. */

    dswap_(&imax, &a[imax * a_dim1 + 1], &c__1, &a[(k - 1) * a_dim1 + 1], &
	    c__1);
    i__1 = km1;
    for (jj = imax; jj <= i__1; ++jj) {
	j = km1 + imax - jj;
	t = a[j + (k - 1) * a_dim1];
	a[j + (k - 1) * a_dim1] = a[imax + j * a_dim1];
	a[imax + j * a_dim1] = t;
/* L150: */
    }
    t = a[k - 1 + k * a_dim1];
    a[k - 1 + k * a_dim1] = a[imax + k * a_dim1];
    a[imax + k * a_dim1] = t;
L160:

/*           PERFORM THE ELIMINATION. */

    km2 = k - 2;
    if (km2 == 0) {
	goto L180;
    }
    ak = a[k + k * a_dim1] / a[k - 1 + k * a_dim1];
    akm1 = a[k - 1 + (k - 1) * a_dim1] / a[k - 1 + k * a_dim1];
    denom = 1. - ak * akm1;
    i__1 = km2;
    for (jj = 1; jj <= i__1; ++jj) {
	j = km1 - jj;
	bk = a[j + k * a_dim1] / a[k - 1 + k * a_dim1];
	bkm1 = a[j + (k - 1) * a_dim1] / a[k - 1 + k * a_dim1];
	mulk = (akm1 * bk - bkm1) / denom;
	mulkm1 = (ak * bkm1 - bk) / denom;
	t = mulk;
	daxpy_(&j, &t, &a[k * a_dim1 + 1], &c__1, &a[j * a_dim1 + 1], &c__1);
	t = mulkm1;
	daxpy_(&j, &t, &a[(k - 1) * a_dim1 + 1], &c__1, &a[j * a_dim1 + 1], &
		c__1);
	a[j + k * a_dim1] = mulk;
	a[j + (k - 1) * a_dim1] = mulkm1;
/* L170: */
    }
L180:

/*           SET THE PIVOT ARRAY. */

    kpvt[k] = 1 - k;
    if (swap) {
	kpvt[k] = -imax;
    }
    kpvt[k - 1] = kpvt[k];
L190:
    k -= kstep;
    goto L10;
L200:
    return 0;
} /* dsifa_ */

