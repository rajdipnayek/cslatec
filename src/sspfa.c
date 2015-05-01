/* sspfa.f -- translated by f2c (version 12.02.01).
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

/* DECK SSPFA */
/* Subroutine */ int sspfa_(real *ap, integer *n, integer *kpvt, integer *
	info)
{
    /* System generated locals */
    integer i__1;
    real r__1, r__2, r__3;

    /* Local variables */
    static integer j, k;
    static real t, ak, bk;
    static integer ij, ik, jj, im, jk, kk, km1, km2, imj, imk;
    static real akm1, bkm1;
    static integer ikm1, jkm1, km1k, imim, jmim, imax, jmax;
    static real mulk;
    static logical swap;
    static real alpha;
    static integer km1km1;
    static real denom;
    static integer kstep;
    extern /* Subroutine */ int sswap_(integer *, real *, integer *, real *, 
	    integer *), saxpy_(integer *, real *, real *, integer *, real *, 
	    integer *);
    static integer imaxp1;
    static real mulkm1, absakk;
    extern integer isamax_(integer *, real *, integer *);
    static real colmax, rowmax;

/* ***BEGIN PROLOGUE  SSPFA */
/* ***PURPOSE  Factor a real symmetric matrix stored in packed form by */
/*            elimination with symmetric pivoting. */
/* ***LIBRARY   SLATEC (LINPACK) */
/* ***CATEGORY  D2B1A */
/* ***TYPE      SINGLE PRECISION (SSPFA-S, DSPFA-D, CHPFA-C, CSPFA-C) */
/* ***KEYWORDS  LINEAR ALGEBRA, LINPACK, MATRIX FACTORIZATION, PACKED, */
/*             SYMMETRIC */
/* ***AUTHOR  Bunch, J., (UCSD) */
/* ***DESCRIPTION */

/*     SSPFA factors a real symmetric matrix stored in */
/*     packed form by elimination with symmetric pivoting. */

/*     To solve  A*X = B , follow SSPFA by SSPSL. */
/*     To compute  INVERSE(A)*C , follow SSPFA by SSPSL. */
/*     To compute  DETERMINANT(A) , follow SSPFA by SSPDI. */
/*     To compute  INERTIA(A) , follow SSPFA by SSPDI. */
/*     To compute  INVERSE(A) , follow SSPFA by SSPDI. */

/*     On Entry */

/*        AP      REAL (N*(N+1)/2) */
/*                the packed form of a symmetric matrix  A .  The */
/*                columns of the upper triangle are stored sequentially */
/*                in a one-dimensional array of length  N*(N+1)/2 . */
/*                See comments below for details. */

/*        N       INTEGER */
/*                the order of the matrix  A . */

/*     Output */

/*        AP      a block diagonal matrix and the multipliers which */
/*                were used to obtain it stored in packed form. */
/*                The factorization can be written  A = U*D*TRANS(U) */
/*                where  U  is a product of permutation and unit */
/*                upper triangular matrices , TRANS(U) is the */
/*                transpose of  U , and  D  is block diagonal */
/*                with 1 by 1 and 2 by 2 blocks. */

/*        KPVT    INTEGER(N) */
/*                an integer vector of pivot indices. */

/*        INFO    INTEGER */
/*                = 0  normal value. */
/*                = K  if the K-th pivot block is singular.  This is */
/*                     not an error condition for this subroutine, */
/*                     but it does indicate that SSPSL or SSPDI may */
/*                     divide by zero if called. */

/*     Packed Storage */

/*          The following program segment will pack the upper */
/*          triangle of a symmetric matrix. */

/*                K = 0 */
/*                DO 20 J = 1, N */
/*                   DO 10 I = 1, J */
/*                      K = K + 1 */
/*                      AP(K)  = A(I,J) */
/*             10    CONTINUE */
/*             20 CONTINUE */

/* ***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W. */
/*                 Stewart, LINPACK Users' Guide, SIAM, 1979. */
/* ***ROUTINES CALLED  ISAMAX, SAXPY, SSWAP */
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
/* ***END PROLOGUE  SSPFA */

/* ***FIRST EXECUTABLE STATEMENT  SSPFA */

/*     INITIALIZE */

/*     ALPHA IS USED IN CHOOSING PIVOT BLOCK SIZE. */

    /* Parameter adjustments */
    --kpvt;
    --ap;

    /* Function Body */
    alpha = (sqrt(17.f) + 1.f) / 8.f;

    *info = 0;

/*     MAIN LOOP ON K, WHICH GOES FROM N TO 1. */

    k = *n;
    ik = *n * (*n - 1) / 2;
L10:

/*        LEAVE THE LOOP IF K=0 OR K=1. */

    if (k == 0) {
	goto L200;
    }
    if (k > 1) {
	goto L20;
    }
    kpvt[1] = 1;
    if (ap[1] == 0.f) {
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
    kk = ik + k;
    absakk = (r__1 = ap[kk], dabs(r__1));

/*        DETERMINE THE LARGEST OFF-DIAGONAL ELEMENT IN */
/*        COLUMN K. */

    i__1 = k - 1;
    imax = isamax_(&i__1, &ap[ik + 1], &c__1);
    imk = ik + imax;
    colmax = (r__1 = ap[imk], dabs(r__1));
    if (absakk < alpha * colmax) {
	goto L30;
    }
    kstep = 1;
    swap = FALSE_;
    goto L90;
L30:

/*           DETERMINE THE LARGEST OFF-DIAGONAL ELEMENT IN */
/*           ROW IMAX. */

    rowmax = 0.f;
    imaxp1 = imax + 1;
    im = imax * (imax - 1) / 2;
    imj = im + (imax << 1);
    i__1 = k;
    for (j = imaxp1; j <= i__1; ++j) {
/* Computing MAX */
	r__2 = rowmax, r__3 = (r__1 = ap[imj], dabs(r__1));
	rowmax = dmax(r__2,r__3);
	imj += j;
/* L40: */
    }
    if (imax == 1) {
	goto L50;
    }
    i__1 = imax - 1;
    jmax = isamax_(&i__1, &ap[im + 1], &c__1);
    jmim = jmax + im;
/* Computing MAX */
    r__2 = rowmax, r__3 = (r__1 = ap[jmim], dabs(r__1));
    rowmax = dmax(r__2,r__3);
L50:
    imim = imax + im;
    if ((r__1 = ap[imim], dabs(r__1)) < alpha * rowmax) {
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
    if (dmax(absakk,colmax) != 0.f) {
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

    sswap_(&imax, &ap[im + 1], &c__1, &ap[ik + 1], &c__1);
    imj = ik + imax;
    i__1 = k;
    for (jj = imax; jj <= i__1; ++jj) {
	j = k + imax - jj;
	jk = ik + j;
	t = ap[jk];
	ap[jk] = ap[imj];
	ap[imj] = t;
	imj -= j - 1;
/* L110: */
    }
L120:

/*           PERFORM THE ELIMINATION. */

    ij = ik - (k - 1);
    i__1 = km1;
    for (jj = 1; jj <= i__1; ++jj) {
	j = k - jj;
	jk = ik + j;
	mulk = -ap[jk] / ap[kk];
	t = mulk;
	saxpy_(&j, &t, &ap[ik + 1], &c__1, &ap[ij + 1], &c__1);
	ap[jk] = mulk;
	ij -= j - 1;
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

    km1k = ik + k - 1;
    ikm1 = ik - (k - 1);
    if (! swap) {
	goto L160;
    }

/*              PERFORM AN INTERCHANGE. */

    sswap_(&imax, &ap[im + 1], &c__1, &ap[ikm1 + 1], &c__1);
    imj = ikm1 + imax;
    i__1 = km1;
    for (jj = imax; jj <= i__1; ++jj) {
	j = km1 + imax - jj;
	jkm1 = ikm1 + j;
	t = ap[jkm1];
	ap[jkm1] = ap[imj];
	ap[imj] = t;
	imj -= j - 1;
/* L150: */
    }
    t = ap[km1k];
    ap[km1k] = ap[imk];
    ap[imk] = t;
L160:

/*           PERFORM THE ELIMINATION. */

    km2 = k - 2;
    if (km2 == 0) {
	goto L180;
    }
    ak = ap[kk] / ap[km1k];
    km1km1 = ikm1 + k - 1;
    akm1 = ap[km1km1] / ap[km1k];
    denom = 1.f - ak * akm1;
    ij = ik - (k - 1) - (k - 2);
    i__1 = km2;
    for (jj = 1; jj <= i__1; ++jj) {
	j = km1 - jj;
	jk = ik + j;
	bk = ap[jk] / ap[km1k];
	jkm1 = ikm1 + j;
	bkm1 = ap[jkm1] / ap[km1k];
	mulk = (akm1 * bk - bkm1) / denom;
	mulkm1 = (ak * bkm1 - bk) / denom;
	t = mulk;
	saxpy_(&j, &t, &ap[ik + 1], &c__1, &ap[ij + 1], &c__1);
	t = mulkm1;
	saxpy_(&j, &t, &ap[ikm1 + 1], &c__1, &ap[ij + 1], &c__1);
	ap[jk] = mulk;
	ap[jkm1] = mulkm1;
	ij -= j - 1;
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
    ik -= k - 1;
    if (kstep == 2) {
	ik -= k - 2;
    }
    k -= kstep;
    goto L10;
L200:
    return 0;
} /* sspfa_ */

