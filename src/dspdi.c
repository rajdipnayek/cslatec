/* dspdi.f -- translated by f2c (version 12.02.01).
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

/* DECK DSPDI */
/* Subroutine */ int dspdi_(doublereal *ap, integer *n, integer *kpvt, 
	doublereal *det, integer *inert, doublereal *work, integer *job)
{
    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1;

    /* Local variables */
    static doublereal d__;
    static integer j, k;
    static doublereal t, ak;
    static integer jb, ij, ik, jk, kk, ks, km1;
    static doublereal ten;
    static integer iks, ksj;
    static doublereal akp1;
    static integer ikp1, jkp1, kkp1;
    extern doublereal ddot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    static doublereal temp, akkp1;
    static integer kskp1;
    static logical nodet;
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), dswap_(integer *, doublereal *, integer 
	    *, doublereal *, integer *), daxpy_(integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *);
    static integer kstep;
    static logical noert, noinv;

/* ***BEGIN PROLOGUE  DSPDI */
/* ***PURPOSE  Compute the determinant, inertia, inverse of a real */
/*            symmetric matrix stored in packed form using the factors */
/*            from DSPFA. */
/* ***LIBRARY   SLATEC (LINPACK) */
/* ***CATEGORY  D2B1A, D3B1A */
/* ***TYPE      DOUBLE PRECISION (SSPDI-S, DSPDI-D, CHPDI-C, CSPDI-C) */
/* ***KEYWORDS  DETERMINANT, INVERSE, LINEAR ALGEBRA, LINPACK, MATRIX, */
/*             PACKED, SYMMETRIC */
/* ***AUTHOR  Bunch, J., (UCSD) */
/* ***DESCRIPTION */

/*     DSPDI computes the determinant, inertia and inverse */
/*     of a double precision symmetric matrix using the factors from */
/*     DSPFA, where the matrix is stored in packed form. */

/*     On Entry */

/*        AP      DOUBLE PRECISION (N*(N+1)/2) */
/*                the output from DSPFA. */

/*        N       INTEGER */
/*                the order of the matrix A. */

/*        KPVT    INTEGER(N) */
/*                the pivot vector from DSPFA. */

/*        WORK    DOUBLE PRECISION(N) */
/*                work vector.  Contents ignored. */

/*        JOB     INTEGER */
/*                JOB has the decimal expansion  ABC  where */
/*                   if  C .NE. 0, the inverse is computed, */
/*                   if  B .NE. 0, the determinant is computed, */
/*                   if  A .NE. 0, the inertia is computed. */

/*                For example, JOB = 111  gives all three. */

/*     On Return */

/*        Variables not requested by JOB are not used. */

/*        AP     contains the upper triangle of the inverse of */
/*               the original matrix, stored in packed form. */
/*               The columns of the upper triangle are stored */
/*               sequentially in a one-dimensional array. */

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

/*        A division by zero will occur if the inverse is requested */
/*        and  DSPCO  has set RCOND .EQ. 0.0 */
/*        or  DSPFA  has set  INFO .NE. 0 . */

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
/* ***END PROLOGUE  DSPDI */

/* ***FIRST EXECUTABLE STATEMENT  DSPDI */
    /* Parameter adjustments */
    --work;
    --inert;
    --det;
    --kpvt;
    --ap;

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
    ik = 0;
    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
	kk = ik + k;
	d__ = ap[kk];

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
	ikp1 = ik + k;
	kkp1 = ikp1 + k;
	t = (d__1 = ap[kkp1], abs(d__1));
	d__ = d__ / t * ap[kkp1 + 1] - t;
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
	ik += k;
/* L130: */
    }
L140:

/*     COMPUTE INVERSE(A) */

    if (noinv) {
	goto L270;
    }
    k = 1;
    ik = 0;
L150:
    if (k > *n) {
	goto L260;
    }
    km1 = k - 1;
    kk = ik + k;
    ikp1 = ik + k;
    kkp1 = ikp1 + k;
    if (kpvt[k] < 0) {
	goto L180;
    }

/*              1 BY 1 */

    ap[kk] = 1. / ap[kk];
    if (km1 < 1) {
	goto L170;
    }
    dcopy_(&km1, &ap[ik + 1], &c__1, &work[1], &c__1);
    ij = 0;
    i__1 = km1;
    for (j = 1; j <= i__1; ++j) {
	jk = ik + j;
	ap[jk] = ddot_(&j, &ap[ij + 1], &c__1, &work[1], &c__1);
	i__2 = j - 1;
	daxpy_(&i__2, &work[j], &ap[ij + 1], &c__1, &ap[ik + 1], &c__1);
	ij += j;
/* L160: */
    }
    ap[kk] += ddot_(&km1, &work[1], &c__1, &ap[ik + 1], &c__1);
L170:
    kstep = 1;
    goto L220;
L180:

/*              2 BY 2 */

    t = (d__1 = ap[kkp1], abs(d__1));
    ak = ap[kk] / t;
    akp1 = ap[kkp1 + 1] / t;
    akkp1 = ap[kkp1] / t;
    d__ = t * (ak * akp1 - 1.);
    ap[kk] = akp1 / d__;
    ap[kkp1 + 1] = ak / d__;
    ap[kkp1] = -akkp1 / d__;
    if (km1 < 1) {
	goto L210;
    }
    dcopy_(&km1, &ap[ikp1 + 1], &c__1, &work[1], &c__1);
    ij = 0;
    i__1 = km1;
    for (j = 1; j <= i__1; ++j) {
	jkp1 = ikp1 + j;
	ap[jkp1] = ddot_(&j, &ap[ij + 1], &c__1, &work[1], &c__1);
	i__2 = j - 1;
	daxpy_(&i__2, &work[j], &ap[ij + 1], &c__1, &ap[ikp1 + 1], &c__1);
	ij += j;
/* L190: */
    }
    ap[kkp1 + 1] += ddot_(&km1, &work[1], &c__1, &ap[ikp1 + 1], &c__1);
    ap[kkp1] += ddot_(&km1, &ap[ik + 1], &c__1, &ap[ikp1 + 1], &c__1);
    dcopy_(&km1, &ap[ik + 1], &c__1, &work[1], &c__1);
    ij = 0;
    i__1 = km1;
    for (j = 1; j <= i__1; ++j) {
	jk = ik + j;
	ap[jk] = ddot_(&j, &ap[ij + 1], &c__1, &work[1], &c__1);
	i__2 = j - 1;
	daxpy_(&i__2, &work[j], &ap[ij + 1], &c__1, &ap[ik + 1], &c__1);
	ij += j;
/* L200: */
    }
    ap[kk] += ddot_(&km1, &work[1], &c__1, &ap[ik + 1], &c__1);
L210:
    kstep = 2;
L220:

/*           SWAP */

    ks = (i__1 = kpvt[k], abs(i__1));
    if (ks == k) {
	goto L250;
    }
    iks = ks * (ks - 1) / 2;
    dswap_(&ks, &ap[iks + 1], &c__1, &ap[ik + 1], &c__1);
    ksj = ik + ks;
    i__1 = k;
    for (jb = ks; jb <= i__1; ++jb) {
	j = k + ks - jb;
	jk = ik + j;
	temp = ap[jk];
	ap[jk] = ap[ksj];
	ap[ksj] = temp;
	ksj -= j - 1;
/* L230: */
    }
    if (kstep == 1) {
	goto L240;
    }
    kskp1 = ikp1 + ks;
    temp = ap[kskp1];
    ap[kskp1] = ap[kkp1];
    ap[kkp1] = temp;
L240:
L250:
    ik += k;
    if (kstep == 2) {
	ik = ik + k + 1;
    }
    k += kstep;
    goto L150;
L260:
L270:
    return 0;
} /* dspdi_ */

