/* dsidi.f -- translated by f2c (version 12.02.01).
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

/* DECK DSIDI */
/* Subroutine */ int dsidi_(doublereal *a, integer *lda, integer *n, integer *
	kpvt, doublereal *det, integer *inert, doublereal *work, integer *job)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;
    doublereal d__1;

    /* Local variables */
    static doublereal d__;
    static integer j, k;
    static doublereal t, ak;
    static integer jb, ks, km1;
    static doublereal ten, akp1;
    extern doublereal ddot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    static doublereal temp, akkp1;
    static logical nodet;
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), dswap_(integer *, doublereal *, integer 
	    *, doublereal *, integer *), daxpy_(integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *);
    static integer kstep;
    static logical noert, noinv;

/* ***BEGIN PROLOGUE  DSIDI */
/* ***PURPOSE  Compute the determinant, inertia and inverse of a real */
/*            symmetric matrix using the factors from DSIFA. */
/* ***LIBRARY   SLATEC (LINPACK) */
/* ***CATEGORY  D2B1A, D3B1A */
/* ***TYPE      DOUBLE PRECISION (SSIDI-S, DSIDI-D, CHIDI-C, CSIDI-C) */
/* ***KEYWORDS  DETERMINANT, INVERSE, LINEAR ALGEBRA, LINPACK, MATRIX, */
/*             SYMMETRIC */
/* ***AUTHOR  Bunch, J., (UCSD) */
/* ***DESCRIPTION */

/*     DSIDI computes the determinant, inertia and inverse */
/*     of a double precision symmetric matrix using the factors from */
/*     DSIFA. */

/*     On Entry */

/*        A       DOUBLE PRECISION(LDA,N) */
/*                the output from DSIFA. */

/*        LDA     INTEGER */
/*                the leading dimension of the array A. */

/*        N       INTEGER */
/*                the order of the matrix A. */

/*        KPVT    INTEGER(N) */
/*                the pivot vector from DSIFA. */

/*        WORK    DOUBLE PRECISION(N) */
/*                work vector.  Contents destroyed. */

/*        JOB     INTEGER */
/*                JOB has the decimal expansion  ABC  where */
/*                   if  C .NE. 0, the inverse is computed, */
/*                   if  B .NE. 0, the determinant is computed, */
/*                   if  A .NE. 0, the inertia is computed. */

/*                For example, JOB = 111  gives all three. */

/*     On Return */

/*        Variables not requested by JOB are not used. */

/*        A      contains the upper triangle of the inverse of */
/*               the original matrix.  The strict lower triangle */
/*               is never referenced. */

/*        DET    DOUBLE PRECISION(2) */
/*               determinant of original matrix. */
/*               DETERMINANT = DET(1) * 10.0**DET(2) */
/*               with 1.0 .LE. ABS(DET(1)) .LT. 10.0 */
/*               or DET(1) = 0.0. */

/*        INERT  INTEGER(3) */
/*               the inertia of the original matrix. */
/*               INERT(1)  =  number of positive eigenvalues. */
/*               INERT(2)  =  number of negative eigenvalues. */
/*               INERT(3)  =  number of zero eigenvalues. */

/*     Error Condition */

/*        A division by zero may occur if the inverse is requested */
/*        and  DSICO  has set RCOND .EQ. 0.0 */
/*        or  DSIFA  has set  INFO .NE. 0 . */

/* ***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W. */
/*                 Stewart, LINPACK Users' Guide, SIAM, 1979. */
/* ***ROUTINES CALLED  DAXPY, DCOPY, DDOT, DSWAP */
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
/* ***END PROLOGUE  DSIDI */

/* ***FIRST EXECUTABLE STATEMENT  DSIDI */
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
    det[1] = 1.;
    det[2] = 0.;
    ten = 10.;
L20:
    t = 0.;
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

	if (t != 0.) {
	    goto L30;
	}
	t = (d__1 = a[k + (k + 1) * a_dim1], abs(d__1));
	d__ = d__ / t * a[k + 1 + (k + 1) * a_dim1] - t;
	goto L40;
L30:
	d__ = t;
	t = 0.;
L40:
L50:

	if (noert) {
	    goto L60;
	}
	if (d__ > 0.) {
	    ++inert[1];
	}
	if (d__ < 0.) {
	    ++inert[2];
	}
	if (d__ == 0.) {
	    ++inert[3];
	}
L60:

	if (nodet) {
	    goto L120;
	}
	det[1] = d__ * det[1];
	if (det[1] == 0.) {
	    goto L110;
	}
L70:
	if (abs(det[1]) >= 1.) {
	    goto L80;
	}
	det[1] = ten * det[1];
	det[2] += -1.;
	goto L70;
L80:
L90:
	if (abs(det[1]) < ten) {
	    goto L100;
	}
	det[1] /= ten;
	det[2] += 1.;
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

    a[k + k * a_dim1] = 1. / a[k + k * a_dim1];
    if (km1 < 1) {
	goto L170;
    }
    dcopy_(&km1, &a[k * a_dim1 + 1], &c__1, &work[1], &c__1);
    i__1 = km1;
    for (j = 1; j <= i__1; ++j) {
	a[j + k * a_dim1] = ddot_(&j, &a[j * a_dim1 + 1], &c__1, &work[1], &
		c__1);
	i__2 = j - 1;
	daxpy_(&i__2, &work[j], &a[j * a_dim1 + 1], &c__1, &a[k * a_dim1 + 1],
		 &c__1);
/* L160: */
    }
    a[k + k * a_dim1] += ddot_(&km1, &work[1], &c__1, &a[k * a_dim1 + 1], &
	    c__1);
L170:
    kstep = 1;
    goto L220;
L180:

/*              2 BY 2 */

    t = (d__1 = a[k + (k + 1) * a_dim1], abs(d__1));
    ak = a[k + k * a_dim1] / t;
    akp1 = a[k + 1 + (k + 1) * a_dim1] / t;
    akkp1 = a[k + (k + 1) * a_dim1] / t;
    d__ = t * (ak * akp1 - 1.);
    a[k + k * a_dim1] = akp1 / d__;
    a[k + 1 + (k + 1) * a_dim1] = ak / d__;
    a[k + (k + 1) * a_dim1] = -akkp1 / d__;
    if (km1 < 1) {
	goto L210;
    }
    dcopy_(&km1, &a[(k + 1) * a_dim1 + 1], &c__1, &work[1], &c__1);
    i__1 = km1;
    for (j = 1; j <= i__1; ++j) {
	a[j + (k + 1) * a_dim1] = ddot_(&j, &a[j * a_dim1 + 1], &c__1, &work[
		1], &c__1);
	i__2 = j - 1;
	daxpy_(&i__2, &work[j], &a[j * a_dim1 + 1], &c__1, &a[(k + 1) * 
		a_dim1 + 1], &c__1);
/* L190: */
    }
    a[k + 1 + (k + 1) * a_dim1] += ddot_(&km1, &work[1], &c__1, &a[(k + 1) * 
	    a_dim1 + 1], &c__1);
    a[k + (k + 1) * a_dim1] += ddot_(&km1, &a[k * a_dim1 + 1], &c__1, &a[(k + 
	    1) * a_dim1 + 1], &c__1);
    dcopy_(&km1, &a[k * a_dim1 + 1], &c__1, &work[1], &c__1);
    i__1 = km1;
    for (j = 1; j <= i__1; ++j) {
	a[j + k * a_dim1] = ddot_(&j, &a[j * a_dim1 + 1], &c__1, &work[1], &
		c__1);
	i__2 = j - 1;
	daxpy_(&i__2, &work[j], &a[j * a_dim1 + 1], &c__1, &a[k * a_dim1 + 1],
		 &c__1);
/* L200: */
    }
    a[k + k * a_dim1] += ddot_(&km1, &work[1], &c__1, &a[k * a_dim1 + 1], &
	    c__1);
L210:
    kstep = 2;
L220:

/*           SWAP */

    ks = (i__1 = kpvt[k], abs(i__1));
    if (ks == k) {
	goto L250;
    }
    dswap_(&ks, &a[ks * a_dim1 + 1], &c__1, &a[k * a_dim1 + 1], &c__1);
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
} /* dsidi_ */

