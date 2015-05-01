/* cspfa.f -- translated by f2c (version 12.02.01).
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

/* DECK CSPFA */
/* Subroutine */ int cspfa_(complex *ap, integer *n, integer *kpvt, integer *
	info)
{
    /* System generated locals */
    integer i__1, i__2, i__3;
    real r__1, r__2, r__3, r__4;
    complex q__1, q__2, q__3;

    /* Local variables */
    static integer j, k;
    static complex t, ak, bk;
    static integer ij, ik, jj, im, jk, kk, km1, km2, imj, imk;
    static complex akm1, bkm1;
    static integer ikm1, jkm1, km1k, imim, jmim, imax, jmax;
    static complex mulk;
    static logical swap;
    static real alpha;
    static integer km1km1;
    static complex denom;
    extern /* Subroutine */ int cswap_(integer *, complex *, integer *, 
	    complex *, integer *), caxpy_(integer *, complex *, complex *, 
	    integer *, complex *, integer *);
    static integer kstep, imaxp1;
    static complex mulkm1;
    static real absakk;
    extern integer icamax_(integer *, complex *, integer *);
    static real colmax, rowmax;

/* ***BEGIN PROLOGUE  CSPFA */
/* ***PURPOSE  Factor a complex symmetric matrix stored in packed form by */
/*            elimination with symmetric pivoting. */
/* ***LIBRARY   SLATEC (LINPACK) */
/* ***CATEGORY  D2C1 */
/* ***TYPE      COMPLEX (SSPFA-S, DSPFA-D, CHPFA-C, CSPFA-C) */
/* ***KEYWORDS  LINEAR ALGEBRA, LINPACK, MATRIX FACTORIZATION, PACKED, */
/*             SYMMETRIC */
/* ***AUTHOR  Bunch, J., (UCSD) */
/* ***DESCRIPTION */

/*     CSPFA factors a complex symmetric matrix stored in */
/*     packed form by elimination with symmetric pivoting. */

/*     To solve  A*X = B , follow CSPFA by CSPSL. */
/*     To compute  INVERSE(A)*C , follow CSPFA by CSPSL. */
/*     To compute  DETERMINANT(A) , follow CSPFA by CSPDI. */
/*     To compute  INVERSE(A) , follow CSPFA by CSPDI. */

/*     On Entry */

/*        AP      COMPLEX (N*(N+1)/2) */
/*                the packed form of a symmetric matrix  A .  The */
/*                columns of the upper triangle are stored sequentially */
/*                in a one-dimensional array of length  N*(N+1)/2 . */
/*                See comments below for details. */

/*        N       INTEGER */
/*                the order of the matrix  A . */

/*     On Return */

/*        AP      a block diagonal matrix and the multipliers which */
/*                were used to obtain it stored in packed form. */
/*                The factorization can be written  A = U*D*TRANS(U) */
/*                where  U  is a product of permutation and unit */
/*                upper triangular matrices , TRANS(U) is the */
/*                transpose of  U , and  D  is block diagonal */
/*                with 1 by 1 and 2 by 2 blocks. */

/*        KVPT    INTEGER(N) */
/*                an integer vector of pivot indices. */

/*        INFO    INTEGER */
/*                = 0  normal value. */
/*                = K  if the K-th pivot block is singular.  This is */
/*                     not an error condition for this subroutine, */
/*                     but it does indicate that CSPSL or CSPDI may */
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
/* ***ROUTINES CALLED  CAXPY, CSWAP, ICAMAX */
/* ***REVISION HISTORY  (YYMMDD) */
/*   780814  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890831  Modified array declarations.  (WRB) */
/*   891107  Corrected category and modified routine equivalence */
/*           list.  (WRB) */
/*   891107  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900326  Removed duplicate information from DESCRIPTION section. */
/*           (WRB) */
/*   920501  Reformatted the REFERENCES section.  (WRB) */
/* ***END PROLOGUE  CSPFA */

/* ***FIRST EXECUTABLE STATEMENT  CSPFA */

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
    if ((r__1 = ap[1].r, dabs(r__1)) + (r__2 = r_imag(&ap[1]), dabs(r__2)) == 
	    0.f) {
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
    i__1 = kk;
    absakk = (r__1 = ap[i__1].r, dabs(r__1)) + (r__2 = r_imag(&ap[kk]), dabs(
	    r__2));

/*        DETERMINE THE LARGEST OFF-DIAGONAL ELEMENT IN */
/*        COLUMN K. */

    i__1 = k - 1;
    imax = icamax_(&i__1, &ap[ik + 1], &c__1);
    imk = ik + imax;
    i__1 = imk;
    colmax = (r__1 = ap[i__1].r, dabs(r__1)) + (r__2 = r_imag(&ap[imk]), dabs(
	    r__2));
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
	i__2 = imj;
	r__3 = rowmax, r__4 = (r__1 = ap[i__2].r, dabs(r__1)) + (r__2 = 
		r_imag(&ap[imj]), dabs(r__2));
	rowmax = dmax(r__3,r__4);
	imj += j;
/* L40: */
    }
    if (imax == 1) {
	goto L50;
    }
    i__1 = imax - 1;
    jmax = icamax_(&i__1, &ap[im + 1], &c__1);
    jmim = jmax + im;
/* Computing MAX */
    i__1 = jmim;
    r__3 = rowmax, r__4 = (r__1 = ap[i__1].r, dabs(r__1)) + (r__2 = r_imag(&
	    ap[jmim]), dabs(r__2));
    rowmax = dmax(r__3,r__4);
L50:
    imim = imax + im;
    i__1 = imim;
    if ((r__1 = ap[i__1].r, dabs(r__1)) + (r__2 = r_imag(&ap[imim]), dabs(
	    r__2)) < alpha * rowmax) {
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

    cswap_(&imax, &ap[im + 1], &c__1, &ap[ik + 1], &c__1);
    imj = ik + imax;
    i__1 = k;
    for (jj = imax; jj <= i__1; ++jj) {
	j = k + imax - jj;
	jk = ik + j;
	i__2 = jk;
	t.r = ap[i__2].r, t.i = ap[i__2].i;
	i__2 = jk;
	i__3 = imj;
	ap[i__2].r = ap[i__3].r, ap[i__2].i = ap[i__3].i;
	i__2 = imj;
	ap[i__2].r = t.r, ap[i__2].i = t.i;
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
	i__2 = jk;
	q__2.r = -ap[i__2].r, q__2.i = -ap[i__2].i;
	c_div(&q__1, &q__2, &ap[kk]);
	mulk.r = q__1.r, mulk.i = q__1.i;
	t.r = mulk.r, t.i = mulk.i;
	caxpy_(&j, &t, &ap[ik + 1], &c__1, &ap[ij + 1], &c__1);
	i__2 = jk;
	ap[i__2].r = mulk.r, ap[i__2].i = mulk.i;
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

    cswap_(&imax, &ap[im + 1], &c__1, &ap[ikm1 + 1], &c__1);
    imj = ikm1 + imax;
    i__1 = km1;
    for (jj = imax; jj <= i__1; ++jj) {
	j = km1 + imax - jj;
	jkm1 = ikm1 + j;
	i__2 = jkm1;
	t.r = ap[i__2].r, t.i = ap[i__2].i;
	i__2 = jkm1;
	i__3 = imj;
	ap[i__2].r = ap[i__3].r, ap[i__2].i = ap[i__3].i;
	i__2 = imj;
	ap[i__2].r = t.r, ap[i__2].i = t.i;
	imj -= j - 1;
/* L150: */
    }
    i__1 = km1k;
    t.r = ap[i__1].r, t.i = ap[i__1].i;
    i__1 = km1k;
    i__2 = imk;
    ap[i__1].r = ap[i__2].r, ap[i__1].i = ap[i__2].i;
    i__1 = imk;
    ap[i__1].r = t.r, ap[i__1].i = t.i;
L160:

/*           PERFORM THE ELIMINATION. */

    km2 = k - 2;
    if (km2 == 0) {
	goto L180;
    }
    c_div(&q__1, &ap[kk], &ap[km1k]);
    ak.r = q__1.r, ak.i = q__1.i;
    km1km1 = ikm1 + k - 1;
    c_div(&q__1, &ap[km1km1], &ap[km1k]);
    akm1.r = q__1.r, akm1.i = q__1.i;
    q__2.r = ak.r * akm1.r - ak.i * akm1.i, q__2.i = ak.r * akm1.i + ak.i * 
	    akm1.r;
    q__1.r = 1.f - q__2.r, q__1.i = -q__2.i;
    denom.r = q__1.r, denom.i = q__1.i;
    ij = ik - (k - 1) - (k - 2);
    i__1 = km2;
    for (jj = 1; jj <= i__1; ++jj) {
	j = km1 - jj;
	jk = ik + j;
	c_div(&q__1, &ap[jk], &ap[km1k]);
	bk.r = q__1.r, bk.i = q__1.i;
	jkm1 = ikm1 + j;
	c_div(&q__1, &ap[jkm1], &ap[km1k]);
	bkm1.r = q__1.r, bkm1.i = q__1.i;
	q__3.r = akm1.r * bk.r - akm1.i * bk.i, q__3.i = akm1.r * bk.i + 
		akm1.i * bk.r;
	q__2.r = q__3.r - bkm1.r, q__2.i = q__3.i - bkm1.i;
	c_div(&q__1, &q__2, &denom);
	mulk.r = q__1.r, mulk.i = q__1.i;
	q__3.r = ak.r * bkm1.r - ak.i * bkm1.i, q__3.i = ak.r * bkm1.i + ak.i 
		* bkm1.r;
	q__2.r = q__3.r - bk.r, q__2.i = q__3.i - bk.i;
	c_div(&q__1, &q__2, &denom);
	mulkm1.r = q__1.r, mulkm1.i = q__1.i;
	t.r = mulk.r, t.i = mulk.i;
	caxpy_(&j, &t, &ap[ik + 1], &c__1, &ap[ij + 1], &c__1);
	t.r = mulkm1.r, t.i = mulkm1.i;
	caxpy_(&j, &t, &ap[ikm1 + 1], &c__1, &ap[ij + 1], &c__1);
	i__2 = jk;
	ap[i__2].r = mulk.r, ap[i__2].i = mulk.i;
	i__2 = jkm1;
	ap[i__2].r = mulkm1.r, ap[i__2].i = mulkm1.i;
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
} /* cspfa_ */

