/* csidi.f -- translated by f2c (version 12.02.01).
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

static complex c_b3 = {1.f,0.f};
static integer c__1 = 1;

/* DECK CSIDI */
/* Subroutine */ int csidi_(complex *a, integer *lda, integer *n, integer *
	kpvt, complex *det, complex *work, integer *job)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3;
    real r__1, r__2;
    complex q__1, q__2, q__3;

    /* Local variables */
    static complex d__;
    static integer j, k;
    static complex t, ak;
    static integer jb, ks, km1;
    static real ten;
    static complex akp1, temp, akkp1;
    static logical nodet;
    extern /* Subroutine */ int ccopy_(integer *, complex *, integer *, 
	    complex *, integer *);
    extern /* Complex */ void cdotu_(complex *, integer *, complex *, integer 
	    *, complex *, integer *);
    extern /* Subroutine */ int cswap_(integer *, complex *, integer *, 
	    complex *, integer *), caxpy_(integer *, complex *, complex *, 
	    integer *, complex *, integer *);
    static integer kstep;
    static logical noinv;

/* ***BEGIN PROLOGUE  CSIDI */
/* ***PURPOSE  Compute the determinant and inverse of a complex symmetric */
/*            matrix using the factors from CSIFA. */
/* ***LIBRARY   SLATEC (LINPACK) */
/* ***CATEGORY  D2C1, D3C1 */
/* ***TYPE      COMPLEX (SSIDI-S, DSIDI-D, CHIDI-C, CSIDI-C) */
/* ***KEYWORDS  DETERMINANT, INVERSE, LINEAR ALGEBRA, LINPACK, MATRIX, */
/*             SYMMETRIC */
/* ***AUTHOR  Bunch, J., (UCSD) */
/* ***DESCRIPTION */

/*     CSIDI computes the determinant and inverse */
/*     of a complex symmetric matrix using the factors from CSIFA. */

/*     On Entry */

/*        A       COMPLEX(LDA,N) */
/*                the output from CSIFA. */

/*        LDA     INTEGER */
/*                the leading dimension of the array A . */

/*        N       INTEGER */
/*                the order of the matrix A . */

/*        KVPT    INTEGER(N) */
/*                the pivot vector from CSIFA. */

/*        WORK    COMPLEX(N) */
/*                work vector.  Contents destroyed. */

/*        JOB     INTEGER */
/*                JOB has the decimal expansion  AB  where */
/*                   If  B .NE. 0, the inverse is computed, */
/*                   If  A .NE. 0, the determinant is computed, */

/*                For example, JOB = 11  gives both. */

/*     On Return */

/*        Variables not requested by JOB are not used. */

/*        A      contains the upper triangle of the inverse of */
/*               the original matrix.  The strict lower triangle */
/*               is never referenced. */

/*        DET    COMPLEX(2) */
/*               determinant of original matrix. */
/*               Determinant = DET(1) * 10.0**DET(2) */
/*               with 1.0 .LE. ABS(DET(1)) .LT. 10.0 */
/*               or DET(1) = 0.0. */

/*     Error Condition */

/*        A division by zero may occur if the inverse is requested */
/*        and  CSICO  has set RCOND .EQ. 0.0 */
/*        or  CSIFA  has set  INFO .NE. 0 . */

/* ***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W. */
/*                 Stewart, LINPACK Users' Guide, SIAM, 1979. */
/* ***ROUTINES CALLED  CAXPY, CCOPY, CDOTU, CSWAP */
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
/* ***END PROLOGUE  CSIDI */


/* ***FIRST EXECUTABLE STATEMENT  CSIDI */
    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --kpvt;
    --det;
    --work;

    /* Function Body */
    noinv = *job % 10 == 0;
    nodet = *job % 100 / 10 == 0;

    if (nodet) {
	goto L100;
    }
    det[1].r = 1.f, det[1].i = 0.f;
    det[2].r = 0.f, det[2].i = 0.f;
    ten = 10.f;
    t.r = 0.f, t.i = 0.f;
    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
	i__2 = k + k * a_dim1;
	d__.r = a[i__2].r, d__.i = a[i__2].i;

/*           CHECK IF 1 BY 1 */

	if (kpvt[k] > 0) {
	    goto L30;
	}

/*              2 BY 2 BLOCK */
/*              USE DET (D  T)  =  (D/T * C - T) * T */
/*                      (T  C) */
/*              TO AVOID UNDERFLOW/OVERFLOW TROUBLES. */
/*              TAKE TWO PASSES THROUGH SCALING.  USE  T  FOR FLAG. */

	if ((r__1 = t.r, dabs(r__1)) + (r__2 = r_imag(&t), dabs(r__2)) != 0.f)
		 {
	    goto L10;
	}
	i__2 = k + (k + 1) * a_dim1;
	t.r = a[i__2].r, t.i = a[i__2].i;
	c_div(&q__3, &d__, &t);
	i__2 = k + 1 + (k + 1) * a_dim1;
	q__2.r = q__3.r * a[i__2].r - q__3.i * a[i__2].i, q__2.i = q__3.r * a[
		i__2].i + q__3.i * a[i__2].r;
	q__1.r = q__2.r - t.r, q__1.i = q__2.i - t.i;
	d__.r = q__1.r, d__.i = q__1.i;
	goto L20;
L10:
	d__.r = t.r, d__.i = t.i;
	t.r = 0.f, t.i = 0.f;
L20:
L30:

	q__1.r = d__.r * det[1].r - d__.i * det[1].i, q__1.i = d__.r * det[1]
		.i + d__.i * det[1].r;
	det[1].r = q__1.r, det[1].i = q__1.i;
	if ((r__1 = det[1].r, dabs(r__1)) + (r__2 = r_imag(&det[1]), dabs(
		r__2)) == 0.f) {
	    goto L80;
	}
L40:
	if ((r__1 = det[1].r, dabs(r__1)) + (r__2 = r_imag(&det[1]), dabs(
		r__2)) >= 1.f) {
	    goto L50;
	}
	q__2.r = ten, q__2.i = 0.f;
	q__1.r = q__2.r * det[1].r - q__2.i * det[1].i, q__1.i = q__2.r * det[
		1].i + q__2.i * det[1].r;
	det[1].r = q__1.r, det[1].i = q__1.i;
	q__1.r = det[2].r - 1.f, q__1.i = det[2].i - 0.f;
	det[2].r = q__1.r, det[2].i = q__1.i;
	goto L40;
L50:
L60:
	if ((r__1 = det[1].r, dabs(r__1)) + (r__2 = r_imag(&det[1]), dabs(
		r__2)) < ten) {
	    goto L70;
	}
	q__2.r = ten, q__2.i = 0.f;
	c_div(&q__1, &det[1], &q__2);
	det[1].r = q__1.r, det[1].i = q__1.i;
	q__1.r = det[2].r + 1.f, q__1.i = det[2].i + 0.f;
	det[2].r = q__1.r, det[2].i = q__1.i;
	goto L60;
L70:
L80:
/* L90: */
	;
    }
L100:

/*     COMPUTE INVERSE(A) */

    if (noinv) {
	goto L230;
    }
    k = 1;
L110:
    if (k > *n) {
	goto L220;
    }
    km1 = k - 1;
    if (kpvt[k] < 0) {
	goto L140;
    }

/*              1 BY 1 */

    i__1 = k + k * a_dim1;
    c_div(&q__1, &c_b3, &a[k + k * a_dim1]);
    a[i__1].r = q__1.r, a[i__1].i = q__1.i;
    if (km1 < 1) {
	goto L130;
    }
    ccopy_(&km1, &a[k * a_dim1 + 1], &c__1, &work[1], &c__1);
    i__1 = km1;
    for (j = 1; j <= i__1; ++j) {
	i__2 = j + k * a_dim1;
	cdotu_(&q__1, &j, &a[j * a_dim1 + 1], &c__1, &work[1], &c__1);
	a[i__2].r = q__1.r, a[i__2].i = q__1.i;
	i__2 = j - 1;
	caxpy_(&i__2, &work[j], &a[j * a_dim1 + 1], &c__1, &a[k * a_dim1 + 1],
		 &c__1);
/* L120: */
    }
    i__1 = k + k * a_dim1;
    i__2 = k + k * a_dim1;
    cdotu_(&q__2, &km1, &work[1], &c__1, &a[k * a_dim1 + 1], &c__1);
    q__1.r = a[i__2].r + q__2.r, q__1.i = a[i__2].i + q__2.i;
    a[i__1].r = q__1.r, a[i__1].i = q__1.i;
L130:
    kstep = 1;
    goto L180;
L140:

/*              2 BY 2 */

    i__1 = k + (k + 1) * a_dim1;
    t.r = a[i__1].r, t.i = a[i__1].i;
    c_div(&q__1, &a[k + k * a_dim1], &t);
    ak.r = q__1.r, ak.i = q__1.i;
    c_div(&q__1, &a[k + 1 + (k + 1) * a_dim1], &t);
    akp1.r = q__1.r, akp1.i = q__1.i;
    c_div(&q__1, &a[k + (k + 1) * a_dim1], &t);
    akkp1.r = q__1.r, akkp1.i = q__1.i;
    q__3.r = ak.r * akp1.r - ak.i * akp1.i, q__3.i = ak.r * akp1.i + ak.i * 
	    akp1.r;
    q__2.r = q__3.r - 1.f, q__2.i = q__3.i - 0.f;
    q__1.r = t.r * q__2.r - t.i * q__2.i, q__1.i = t.r * q__2.i + t.i * 
	    q__2.r;
    d__.r = q__1.r, d__.i = q__1.i;
    i__1 = k + k * a_dim1;
    c_div(&q__1, &akp1, &d__);
    a[i__1].r = q__1.r, a[i__1].i = q__1.i;
    i__1 = k + 1 + (k + 1) * a_dim1;
    c_div(&q__1, &ak, &d__);
    a[i__1].r = q__1.r, a[i__1].i = q__1.i;
    i__1 = k + (k + 1) * a_dim1;
    q__2.r = -akkp1.r, q__2.i = -akkp1.i;
    c_div(&q__1, &q__2, &d__);
    a[i__1].r = q__1.r, a[i__1].i = q__1.i;
    if (km1 < 1) {
	goto L170;
    }
    ccopy_(&km1, &a[(k + 1) * a_dim1 + 1], &c__1, &work[1], &c__1);
    i__1 = km1;
    for (j = 1; j <= i__1; ++j) {
	i__2 = j + (k + 1) * a_dim1;
	cdotu_(&q__1, &j, &a[j * a_dim1 + 1], &c__1, &work[1], &c__1);
	a[i__2].r = q__1.r, a[i__2].i = q__1.i;
	i__2 = j - 1;
	caxpy_(&i__2, &work[j], &a[j * a_dim1 + 1], &c__1, &a[(k + 1) * 
		a_dim1 + 1], &c__1);
/* L150: */
    }
    i__1 = k + 1 + (k + 1) * a_dim1;
    i__2 = k + 1 + (k + 1) * a_dim1;
    cdotu_(&q__2, &km1, &work[1], &c__1, &a[(k + 1) * a_dim1 + 1], &c__1);
    q__1.r = a[i__2].r + q__2.r, q__1.i = a[i__2].i + q__2.i;
    a[i__1].r = q__1.r, a[i__1].i = q__1.i;
    i__1 = k + (k + 1) * a_dim1;
    i__2 = k + (k + 1) * a_dim1;
    cdotu_(&q__2, &km1, &a[k * a_dim1 + 1], &c__1, &a[(k + 1) * a_dim1 + 1], &
	    c__1);
    q__1.r = a[i__2].r + q__2.r, q__1.i = a[i__2].i + q__2.i;
    a[i__1].r = q__1.r, a[i__1].i = q__1.i;
    ccopy_(&km1, &a[k * a_dim1 + 1], &c__1, &work[1], &c__1);
    i__1 = km1;
    for (j = 1; j <= i__1; ++j) {
	i__2 = j + k * a_dim1;
	cdotu_(&q__1, &j, &a[j * a_dim1 + 1], &c__1, &work[1], &c__1);
	a[i__2].r = q__1.r, a[i__2].i = q__1.i;
	i__2 = j - 1;
	caxpy_(&i__2, &work[j], &a[j * a_dim1 + 1], &c__1, &a[k * a_dim1 + 1],
		 &c__1);
/* L160: */
    }
    i__1 = k + k * a_dim1;
    i__2 = k + k * a_dim1;
    cdotu_(&q__2, &km1, &work[1], &c__1, &a[k * a_dim1 + 1], &c__1);
    q__1.r = a[i__2].r + q__2.r, q__1.i = a[i__2].i + q__2.i;
    a[i__1].r = q__1.r, a[i__1].i = q__1.i;
L170:
    kstep = 2;
L180:

/*           SWAP */

    ks = (i__1 = kpvt[k], abs(i__1));
    if (ks == k) {
	goto L210;
    }
    cswap_(&ks, &a[ks * a_dim1 + 1], &c__1, &a[k * a_dim1 + 1], &c__1);
    i__1 = k;
    for (jb = ks; jb <= i__1; ++jb) {
	j = k + ks - jb;
	i__2 = j + k * a_dim1;
	temp.r = a[i__2].r, temp.i = a[i__2].i;
	i__2 = j + k * a_dim1;
	i__3 = ks + j * a_dim1;
	a[i__2].r = a[i__3].r, a[i__2].i = a[i__3].i;
	i__2 = ks + j * a_dim1;
	a[i__2].r = temp.r, a[i__2].i = temp.i;
/* L190: */
    }
    if (kstep == 1) {
	goto L200;
    }
    i__1 = ks + (k + 1) * a_dim1;
    temp.r = a[i__1].r, temp.i = a[i__1].i;
    i__1 = ks + (k + 1) * a_dim1;
    i__2 = k + (k + 1) * a_dim1;
    a[i__1].r = a[i__2].r, a[i__1].i = a[i__2].i;
    i__1 = k + (k + 1) * a_dim1;
    a[i__1].r = temp.r, a[i__1].i = temp.i;
L200:
L210:
    k += kstep;
    goto L110;
L220:
L230:
    return 0;
} /* csidi_ */

