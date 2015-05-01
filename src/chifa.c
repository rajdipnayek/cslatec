/* chifa.f -- translated by f2c (version 12.02.01).
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

/* DECK CHIFA */
/* Subroutine */ int chifa_(complex *a, integer *lda, integer *n, integer *
	kpvt, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3;
    real r__1, r__2, r__3, r__4;
    complex q__1, q__2, q__3;

    /* Local variables */
    static integer j, k;
    static complex t, ak, bk;
    static integer jj, km1, km2;
    static complex akm1, bkm1;
    static integer imax, jmax;
    static complex mulk;
    static logical swap;
    static real alpha;
    static complex denom;
    extern /* Subroutine */ int cswap_(integer *, complex *, integer *, 
	    complex *, integer *), caxpy_(integer *, complex *, complex *, 
	    integer *, complex *, integer *);
    static integer kstep, imaxp1;
    static complex mulkm1;
    static real absakk;
    extern integer icamax_(integer *, complex *, integer *);
    static real colmax, rowmax;

/* ***BEGIN PROLOGUE  CHIFA */
/* ***PURPOSE  Factor a complex Hermitian matrix by elimination */
/*            (symmetric pivoting). */
/* ***LIBRARY   SLATEC (LINPACK) */
/* ***CATEGORY  D2D1A */
/* ***TYPE      COMPLEX (SSIFA-S, DSIFA-D, CHIFA-C, CSIFA-C) */
/* ***KEYWORDS  HERMITIAN, LINEAR ALGEBRA, LINPACK, MATRIX FACTORIZATION */
/* ***AUTHOR  Bunch, J., (UCSD) */
/* ***DESCRIPTION */

/*     CHIFA factors a complex Hermitian matrix by elimination */
/*     with symmetric pivoting. */

/*     To solve  A*X = B , follow CHIFA by CHISL. */
/*     To compute  INVERSE(A)*C , follow CHIFA by CHISL. */
/*     To compute  DETERMINANT(A) , follow CHIFA by CHIDI. */
/*     To compute  INERTIA(A) , follow CHIFA by CHIDI. */
/*     To compute  INVERSE(A) , follow CHIFA by CHIDI. */

/*     On Entry */

/*        A       COMPLEX(LDA,N) */
/*                the Hermitian matrix to be factored. */
/*                Only the diagonal and upper triangle are used. */

/*        LDA     INTEGER */
/*                the leading dimension of the array  A . */

/*        N       INTEGER */
/*                the order of the matrix  A . */

/*     On Return */

/*        A       a block diagonal matrix and the multipliers which */
/*                were used to obtain it. */
/*                The factorization can be written  A = U*D*CTRANS(U) */
/*                where  U  is a product of permutation and unit */
/*                upper triangular matrices , CTRANS(U) is the */
/*                conjugate transpose of  U , and  D  is block diagonal */
/*                with 1 by 1 and 2 by 2 blocks. */

/*        KVPT    INTEGER(N) */
/*                an integer vector of pivot indices. */

/*        INFO    INTEGER */
/*                = 0  normal value. */
/*                = K  if the K-th pivot block is singular.  This is */
/*                     not an error condition for this subroutine, */
/*                     but it does indicate that CHISL or CHIDI may */
/*                     divide by zero if called. */

/* ***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W. */
/*                 Stewart, LINPACK Users' Guide, SIAM, 1979. */
/* ***ROUTINES CALLED  CAXPY, CSWAP, ICAMAX */
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
/* ***END PROLOGUE  CHIFA */

/* ***FIRST EXECUTABLE STATEMENT  CHIFA */

/*     INITIALIZE */

/*     ALPHA IS USED IN CHOOSING PIVOT BLOCK SIZE. */

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --kpvt;

    /* Function Body */
    alpha = (sqrt(17.f) + 1.f) / 8.f;

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
    i__1 = a_dim1 + 1;
    if ((r__1 = a[i__1].r, dabs(r__1)) + (r__2 = r_imag(&a[a_dim1 + 1]), dabs(
	    r__2)) == 0.f) {
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
    i__1 = k + k * a_dim1;
    absakk = (r__1 = a[i__1].r, dabs(r__1)) + (r__2 = r_imag(&a[k + k * 
	    a_dim1]), dabs(r__2));

/*        DETERMINE THE LARGEST OFF-DIAGONAL ELEMENT IN */
/*        COLUMN K. */

    i__1 = k - 1;
    imax = icamax_(&i__1, &a[k * a_dim1 + 1], &c__1);
    i__1 = imax + k * a_dim1;
    colmax = (r__1 = a[i__1].r, dabs(r__1)) + (r__2 = r_imag(&a[imax + k * 
	    a_dim1]), dabs(r__2));
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
    i__1 = k;
    for (j = imaxp1; j <= i__1; ++j) {
/* Computing MAX */
	i__2 = imax + j * a_dim1;
	r__3 = rowmax, r__4 = (r__1 = a[i__2].r, dabs(r__1)) + (r__2 = r_imag(
		&a[imax + j * a_dim1]), dabs(r__2));
	rowmax = dmax(r__3,r__4);
/* L40: */
    }
    if (imax == 1) {
	goto L50;
    }
    i__1 = imax - 1;
    jmax = icamax_(&i__1, &a[imax * a_dim1 + 1], &c__1);
/* Computing MAX */
    i__1 = jmax + imax * a_dim1;
    r__3 = rowmax, r__4 = (r__1 = a[i__1].r, dabs(r__1)) + (r__2 = r_imag(&a[
	    jmax + imax * a_dim1]), dabs(r__2));
    rowmax = dmax(r__3,r__4);
L50:
    i__1 = imax + imax * a_dim1;
    if ((r__1 = a[i__1].r, dabs(r__1)) + (r__2 = r_imag(&a[imax + imax * 
	    a_dim1]), dabs(r__2)) < alpha * rowmax) {
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

    cswap_(&imax, &a[imax * a_dim1 + 1], &c__1, &a[k * a_dim1 + 1], &c__1);
    i__1 = k;
    for (jj = imax; jj <= i__1; ++jj) {
	j = k + imax - jj;
	r_cnjg(&q__1, &a[j + k * a_dim1]);
	t.r = q__1.r, t.i = q__1.i;
	i__2 = j + k * a_dim1;
	r_cnjg(&q__1, &a[imax + j * a_dim1]);
	a[i__2].r = q__1.r, a[i__2].i = q__1.i;
	i__2 = imax + j * a_dim1;
	a[i__2].r = t.r, a[i__2].i = t.i;
/* L110: */
    }
L120:

/*           PERFORM THE ELIMINATION. */

    i__1 = km1;
    for (jj = 1; jj <= i__1; ++jj) {
	j = k - jj;
	i__2 = j + k * a_dim1;
	q__2.r = -a[i__2].r, q__2.i = -a[i__2].i;
	c_div(&q__1, &q__2, &a[k + k * a_dim1]);
	mulk.r = q__1.r, mulk.i = q__1.i;
	r_cnjg(&q__1, &mulk);
	t.r = q__1.r, t.i = q__1.i;
	caxpy_(&j, &t, &a[k * a_dim1 + 1], &c__1, &a[j * a_dim1 + 1], &c__1);
	i__2 = j + j * a_dim1;
	i__3 = j + j * a_dim1;
	r__1 = a[i__3].r;
	q__1.r = r__1, q__1.i = 0.f;
	a[i__2].r = q__1.r, a[i__2].i = q__1.i;
	i__2 = j + k * a_dim1;
	a[i__2].r = mulk.r, a[i__2].i = mulk.i;
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

    cswap_(&imax, &a[imax * a_dim1 + 1], &c__1, &a[(k - 1) * a_dim1 + 1], &
	    c__1);
    i__1 = km1;
    for (jj = imax; jj <= i__1; ++jj) {
	j = km1 + imax - jj;
	r_cnjg(&q__1, &a[j + (k - 1) * a_dim1]);
	t.r = q__1.r, t.i = q__1.i;
	i__2 = j + (k - 1) * a_dim1;
	r_cnjg(&q__1, &a[imax + j * a_dim1]);
	a[i__2].r = q__1.r, a[i__2].i = q__1.i;
	i__2 = imax + j * a_dim1;
	a[i__2].r = t.r, a[i__2].i = t.i;
/* L150: */
    }
    i__1 = k - 1 + k * a_dim1;
    t.r = a[i__1].r, t.i = a[i__1].i;
    i__1 = k - 1 + k * a_dim1;
    i__2 = imax + k * a_dim1;
    a[i__1].r = a[i__2].r, a[i__1].i = a[i__2].i;
    i__1 = imax + k * a_dim1;
    a[i__1].r = t.r, a[i__1].i = t.i;
L160:

/*           PERFORM THE ELIMINATION. */

    km2 = k - 2;
    if (km2 == 0) {
	goto L180;
    }
    c_div(&q__1, &a[k + k * a_dim1], &a[k - 1 + k * a_dim1]);
    ak.r = q__1.r, ak.i = q__1.i;
    r_cnjg(&q__2, &a[k - 1 + k * a_dim1]);
    c_div(&q__1, &a[k - 1 + (k - 1) * a_dim1], &q__2);
    akm1.r = q__1.r, akm1.i = q__1.i;
    q__2.r = ak.r * akm1.r - ak.i * akm1.i, q__2.i = ak.r * akm1.i + ak.i * 
	    akm1.r;
    q__1.r = 1.f - q__2.r, q__1.i = -q__2.i;
    denom.r = q__1.r, denom.i = q__1.i;
    i__1 = km2;
    for (jj = 1; jj <= i__1; ++jj) {
	j = km1 - jj;
	c_div(&q__1, &a[j + k * a_dim1], &a[k - 1 + k * a_dim1]);
	bk.r = q__1.r, bk.i = q__1.i;
	r_cnjg(&q__2, &a[k - 1 + k * a_dim1]);
	c_div(&q__1, &a[j + (k - 1) * a_dim1], &q__2);
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
	r_cnjg(&q__1, &mulk);
	t.r = q__1.r, t.i = q__1.i;
	caxpy_(&j, &t, &a[k * a_dim1 + 1], &c__1, &a[j * a_dim1 + 1], &c__1);
	r_cnjg(&q__1, &mulkm1);
	t.r = q__1.r, t.i = q__1.i;
	caxpy_(&j, &t, &a[(k - 1) * a_dim1 + 1], &c__1, &a[j * a_dim1 + 1], &
		c__1);
	i__2 = j + k * a_dim1;
	a[i__2].r = mulk.r, a[i__2].i = mulk.i;
	i__2 = j + (k - 1) * a_dim1;
	a[i__2].r = mulkm1.r, a[i__2].i = mulkm1.i;
	i__2 = j + j * a_dim1;
	i__3 = j + j * a_dim1;
	r__1 = a[i__3].r;
	q__1.r = r__1, q__1.i = 0.f;
	a[i__2].r = q__1.r, a[i__2].i = q__1.i;
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
} /* chifa_ */

