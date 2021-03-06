/* ssidi.f -- translated by f2c (version 12.02.01).
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

/* DECK SSIDI */
/* Subroutine */ int ssidi_(real *a, integer *lda, integer *n, integer *kpvt, 
	real *det, integer *inert, real *work, integer *job)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;
    real r__1;

    /* Local variables */
    static real d__;
    static integer j, k;
    static real t, ak;
    static integer jb, ks, km1;
    static real ten, akp1, temp;
    extern doublereal sdot_(integer *, real *, integer *, real *, integer *);
    static real akkp1;
    static logical nodet;
    static integer kstep;
    static logical noert, noinv;
    extern /* Subroutine */ int scopy_(integer *, real *, integer *, real *, 
	    integer *), sswap_(integer *, real *, integer *, real *, integer *
	    ), saxpy_(integer *, real *, real *, integer *, real *, integer *)
	    ;

/* ***BEGIN PROLOGUE  SSIDI */
/* ***PURPOSE  Compute the determinant, inertia and inverse of a real */
/*            symmetric matrix using the factors from SSIFA. */
/* ***LIBRARY   SLATEC (LINPACK) */
/* ***CATEGORY  D2B1A, D3B1A */
/* ***TYPE      SINGLE PRECISION (SSIDI-S, DSIDI-D, CHIDI-C, CSIDI-C) */
/* ***KEYWORDS  DETERMINANT, INVERSE, LINEAR ALGEBRA, LINPACK, MATRIX, */
/*             SYMMETRIC */
/* ***AUTHOR  Bunch, J., (UCSD) */
/* ***DESCRIPTION */

/*     SSIDI computes the determinant, inertia and inverse */
/*     of a real symmetric matrix using the factors from SSIFA. */

/*     On Entry */

/*        A       REAL(LDA,N) */
/*                the output from SSIFA. */

/*        LDA     INTEGER */
/*                the leading dimension of the array A. */

/*        N       INTEGER */
/*                the order of the matrix A. */

/*        KPVT    INTEGER(N) */
/*                the pivot vector from SSIFA. */

/*        WORK    REAL(N) */
/*                work vector.  Contents destroyed. */

/*        JOB     INTEGER */
/*                JOB has the decimal expansion  ABC  where */
/*                   If  C .NE. 0, the inverse is computed, */
/*                   If  B .NE. 0, the determinant is computed, */
/*                   If  A .NE. 0, the inertia is computed. */

/*                For example, JOB = 111  gives all three. */

/*     On Return */

/*        Variables not requested by JOB are not used. */

/*        A      contains the upper triangle of the inverse of */
/*               the original matrix.  The strict lower triangle */
/*               is never referenced. */

/*        DET    REAL(2) */
/*               determinant of original matrix. */
/*               Determinant = DET(1) * 10.0**DET(2) */
/*               with 1.0 .LE. ABS(DET(1)) .LT. 10.0 */
/*               or DET(1) = 0.0. */

/*        INERT  INTEGER(3) */
/*               the inertia of the original matrix. */
/*               INERT(1)  =  number of positive eigenvalues. */
/*               INERT(2)  =  number of negative eigenvalues. */
/*               INERT(3)  =  number of zero eigenvalues. */

/*     Error Condition */

/*        A division by zero may occur if the inverse is requested */
/*        and  SSICO  has set RCOND .EQ. 0.0 */
/*        or  SSIFA  has set  INFO .NE. 0 . */

/* ***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W. */
/*                 Stewart, LINPACK Users' Guide, SIAM, 1979. */
/* ***ROUTINES CALLED  SAXPY, SCOPY, SDOT, SSWAP */
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
/* ***END PROLOGUE  SSIDI */

/* ***FIRST EXECUTABLE STATEMENT  SSIDI */
    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --kpvt;
    --det;
    --inert;
    --work;

    /* Function Body */
    noinv = *job % 10 == 0;
    nodet = *job % 100 / 10 == 0;
    noert = *job % 1000 / 100 == 0;

    if (nodet && noert) {
	goto L140;
    }
    if (noert) {
	goto L10;
    }
    inert[1] = 0;
    inert[2] = 0;
    inert[3] = 0;
L10:
    if (nodet) {
	goto L20;
    }
    det[1] = 1.f;
    det[2] = 0.f;
    ten = 10.f;
L20:
    t = 0.f;
    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
	d__ = a[k + k * a_dim1];

/*           CHECK IF 1 BY 1 */

	if (kpvt[k] > 0) {
	    goto L50;
	}

/*              2 BY 2 BLOCK */
/*              USE DET (D  S)  =  (D/T * C - T) * T  ,  T = ABS(S) */
/*                      (S  C) */
/*              TO AVOID UNDERFLOW/OVERFLOW TROUBLES. */
/*              TAKE TWO PASSES THROUGH SCALING.  USE  T  FOR FLAG. */

	if (t != 0.f) {
	    goto L30;
	}
	t = (r__1 = a[k + (k + 1) * a_dim1], dabs(r__1));
	d__ = d__ / t * a[k + 1 + (k + 1) * a_dim1] - t;
	goto L40;
L30:
	d__ = t;
	t = 0.f;
L40:
L50:

	if (noert) {
	    goto L60;
	}
	if (d__ > 0.f) {
	    ++inert[1];
	}
	if (d__ < 0.f) {
	    ++inert[2];
	}
	if (d__ == 0.f) {
	    ++inert[3];
	}
L60:

	if (nodet) {
	    goto L120;
	}
	det[1] = d__ * det[1];
	if (det[1] == 0.f) {
	    goto L110;
	}
L70:
	if (dabs(det[1]) >= 1.f) {
	    goto L80;
	}
	det[1] = ten * det[1];
	det[2] += -1.f;
	goto L70;
L80:
L90:
	if (dabs(det[1]) < ten) {
	    goto L100;
	}
	det[1] /= ten;
	det[2] += 1.f;
	goto L90;
L100:
L110:
L120:
/* L130: */
	;
    }
L140:

/*     COMPUTE INVERSE(A) */

    if (noinv) {
	goto L270;
    }
    k = 1;
L150:
    if (k > *n) {
	goto L260;
    }
    km1 = k - 1;
    if (kpvt[k] < 0) {
	goto L180;
    }

/*              1 BY 1 */

    a[k + k * a_dim1] = 1.f / a[k + k * a_dim1];
    if (km1 < 1) {
	goto L170;
    }
    scopy_(&km1, &a[k * a_dim1 + 1], &c__1, &work[1], &c__1);
    i__1 = km1;
    for (j = 1; j <= i__1; ++j) {
	a[j + k * a_dim1] = sdot_(&j, &a[j * a_dim1 + 1], &c__1, &work[1], &
		c__1);
	i__2 = j - 1;
	saxpy_(&i__2, &work[j], &a[j * a_dim1 + 1], &c__1, &a[k * a_dim1 + 1],
		 &c__1);
/* L160: */
    }
    a[k + k * a_dim1] += sdot_(&km1, &work[1], &c__1, &a[k * a_dim1 + 1], &
	    c__1);
L170:
    kstep = 1;
    goto L220;
L180:

/*              2 BY 2 */

    t = (r__1 = a[k + (k + 1) * a_dim1], dabs(r__1));
    ak = a[k + k * a_dim1] / t;
    akp1 = a[k + 1 + (k + 1) * a_dim1] / t;
    akkp1 = a[k + (k + 1) * a_dim1] / t;
    d__ = t * (ak * akp1 - 1.f);
    a[k + k * a_dim1] = akp1 / d__;
    a[k + 1 + (k + 1) * a_dim1] = ak / d__;
    a[k + (k + 1) * a_dim1] = -akkp1 / d__;
    if (km1 < 1) {
	goto L210;
    }
    scopy_(&km1, &a[(k + 1) * a_dim1 + 1], &c__1, &work[1], &c__1);
    i__1 = km1;
    for (j = 1; j <= i__1; ++j) {
	a[j + (k + 1) * a_dim1] = sdot_(&j, &a[j * a_dim1 + 1], &c__1, &work[
		1], &c__1);
	i__2 = j - 1;
	saxpy_(&i__2, &work[j], &a[j * a_dim1 + 1], &c__1, &a[(k + 1) * 
		a_dim1 + 1], &c__1);
/* L190: */
    }
    a[k + 1 + (k + 1) * a_dim1] += sdot_(&km1, &work[1], &c__1, &a[(k + 1) * 
	    a_dim1 + 1], &c__1);
    a[k + (k + 1) * a_dim1] += sdot_(&km1, &a[k * a_dim1 + 1], &c__1, &a[(k + 
	    1) * a_dim1 + 1], &c__1);
    scopy_(&km1, &a[k * a_dim1 + 1], &c__1, &work[1], &c__1);
    i__1 = km1;
    for (j = 1; j <= i__1; ++j) {
	a[j + k * a_dim1] = sdot_(&j, &a[j * a_dim1 + 1], &c__1, &work[1], &
		c__1);
	i__2 = j - 1;
	saxpy_(&i__2, &work[j], &a[j * a_dim1 + 1], &c__1, &a[k * a_dim1 + 1],
		 &c__1);
/* L200: */
    }
    a[k + k * a_dim1] += sdot_(&km1, &work[1], &c__1, &a[k * a_dim1 + 1], &
	    c__1);
L210:
    kstep = 2;
L220:

/*           SWAP */

    ks = (i__1 = kpvt[k], abs(i__1));
    if (ks == k) {
	goto L250;
    }
    sswap_(&ks, &a[ks * a_dim1 + 1], &c__1, &a[k * a_dim1 + 1], &c__1);
    i__1 = k;
    for (jb = ks; jb <= i__1; ++jb) {
	j = k + ks - jb;
	temp = a[j + k * a_dim1];
	a[j + k * a_dim1] = a[ks + j * a_dim1];
	a[ks + j * a_dim1] = temp;
/* L230: */
    }
    if (kstep == 1) {
	goto L240;
    }
    temp = a[ks + (k + 1) * a_dim1];
    a[ks + (k + 1) * a_dim1] = a[k + (k + 1) * a_dim1];
    a[k + (k + 1) * a_dim1] = temp;
L240:
L250:
    k += kstep;
    goto L150;
L260:
L270:
    return 0;
} /* ssidi_ */

