/* dqrdc.f -- translated by f2c (version 12.02.01).
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

/* DECK DQRDC */
/* Subroutine */ int dqrdc_(doublereal *x, integer *ldx, integer *n, integer *
	p, doublereal *qraux, integer *jpvt, doublereal *work, integer *job)
{
    /* System generated locals */
    integer x_dim1, x_offset, i__1, i__2, i__3;
    doublereal d__1, d__2;

    /* Local variables */
    static integer j, l;
    static doublereal t;
    static integer jj, jp, pl, pu;
    static doublereal tt;
    static integer lp1, lup;
    static logical negj;
    extern doublereal ddot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    static integer maxj;
    extern doublereal dnrm2_(integer *, doublereal *, integer *);
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *), dswap_(integer *, doublereal *, integer *, doublereal 
	    *, integer *);
    static logical swapj;
    extern /* Subroutine */ int daxpy_(integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *);
    static doublereal nrmxl, maxnrm;

/* ***BEGIN PROLOGUE  DQRDC */
/* ***PURPOSE  Use Householder transformations to compute the QR */
/*            factorization of an N by P matrix.  Column pivoting is a */
/*            users option. */
/* ***LIBRARY   SLATEC (LINPACK) */
/* ***CATEGORY  D5 */
/* ***TYPE      DOUBLE PRECISION (SQRDC-S, DQRDC-D, CQRDC-C) */
/* ***KEYWORDS  LINEAR ALGEBRA, LINPACK, MATRIX, ORTHOGONAL TRIANGULAR, */
/*             QR DECOMPOSITION */
/* ***AUTHOR  Stewart, G. W., (U. of Maryland) */
/* ***DESCRIPTION */

/*     DQRDC uses Householder transformations to compute the QR */
/*     factorization of an N by P matrix X.  Column pivoting */
/*     based on the 2-norms of the reduced columns may be */
/*     performed at the user's option. */

/*     On Entry */

/*        X       DOUBLE PRECISION(LDX,P), where LDX .GE. N. */
/*                X contains the matrix whose decomposition is to be */
/*                computed. */

/*        LDX     INTEGER. */
/*                LDX is the leading dimension of the array X. */

/*        N       INTEGER. */
/*                N is the number of rows of the matrix X. */

/*        P       INTEGER. */
/*                P is the number of columns of the matrix X. */

/*        JPVT    INTEGER(P). */
/*                JPVT contains integers that control the selection */
/*                of the pivot columns.  The K-th column X(K) of X */
/*                is placed in one of three classes according to the */
/*                value of JPVT(K). */

/*                   If JPVT(K) .GT. 0, then X(K) is an initial */
/*                                      column. */

/*                   If JPVT(K) .EQ. 0, then X(K) is a free column. */

/*                   If JPVT(K) .LT. 0, then X(K) is a final column. */

/*                Before the decomposition is computed, initial columns */
/*                are moved to the beginning of the array X and final */
/*                columns to the end.  Both initial and final columns */
/*                are frozen in place during the computation and only */
/*                free columns are moved.  At the K-th stage of the */
/*                reduction, if X(K) is occupied by a free column */
/*                it is interchanged with the free column of largest */
/*                reduced norm.  JPVT is not referenced if */
/*                JOB .EQ. 0. */

/*        WORK    DOUBLE PRECISION(P). */
/*                WORK is a work array.  WORK is not referenced if */
/*                JOB .EQ. 0. */

/*        JOB     INTEGER. */
/*                JOB is an integer that initiates column pivoting. */
/*                If JOB .EQ. 0, no pivoting is done. */
/*                If JOB .NE. 0, pivoting is done. */

/*     On Return */

/*        X       X contains in its upper triangle the upper */
/*                triangular matrix R of the QR factorization. */
/*                Below its diagonal X contains information from */
/*                which the orthogonal part of the decomposition */
/*                can be recovered.  Note that if pivoting has */
/*                been requested, the decomposition is not that */
/*                of the original matrix X but that of X */
/*                with its columns permuted as described by JPVT. */

/*        QRAUX   DOUBLE PRECISION(P). */
/*                QRAUX contains further information required to recover */
/*                the orthogonal part of the decomposition. */

/*        JPVT    JPVT(K) contains the index of the column of the */
/*                original matrix that has been interchanged into */
/*                the K-th column, if pivoting was requested. */

/* ***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W. */
/*                 Stewart, LINPACK Users' Guide, SIAM, 1979. */
/* ***ROUTINES CALLED  DAXPY, DDOT, DNRM2, DSCAL, DSWAP */
/* ***REVISION HISTORY  (YYMMDD) */
/*   780814  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890831  Modified array declarations.  (WRB) */
/*   890831  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900326  Removed duplicate information from DESCRIPTION section. */
/*           (WRB) */
/*   920501  Reformatted the REFERENCES section.  (WRB) */
/* ***END PROLOGUE  DQRDC */


/* ***FIRST EXECUTABLE STATEMENT  DQRDC */
    /* Parameter adjustments */
    x_dim1 = *ldx;
    x_offset = 1 + x_dim1;
    x -= x_offset;
    --qraux;
    --jpvt;
    --work;

    /* Function Body */
    pl = 1;
    pu = 0;
    if (*job == 0) {
	goto L60;
    }

/*        PIVOTING HAS BEEN REQUESTED.  REARRANGE THE COLUMNS */
/*        ACCORDING TO JPVT. */

    i__1 = *p;
    for (j = 1; j <= i__1; ++j) {
	swapj = jpvt[j] > 0;
	negj = jpvt[j] < 0;
	jpvt[j] = j;
	if (negj) {
	    jpvt[j] = -j;
	}
	if (! swapj) {
	    goto L10;
	}
	if (j != pl) {
	    dswap_(n, &x[pl * x_dim1 + 1], &c__1, &x[j * x_dim1 + 1], &c__1);
	}
	jpvt[j] = jpvt[pl];
	jpvt[pl] = j;
	++pl;
L10:
/* L20: */
	;
    }
    pu = *p;
    i__1 = *p;
    for (jj = 1; jj <= i__1; ++jj) {
	j = *p - jj + 1;
	if (jpvt[j] >= 0) {
	    goto L40;
	}
	jpvt[j] = -jpvt[j];
	if (j == pu) {
	    goto L30;
	}
	dswap_(n, &x[pu * x_dim1 + 1], &c__1, &x[j * x_dim1 + 1], &c__1);
	jp = jpvt[pu];
	jpvt[pu] = jpvt[j];
	jpvt[j] = jp;
L30:
	--pu;
L40:
/* L50: */
	;
    }
L60:

/*     COMPUTE THE NORMS OF THE FREE COLUMNS. */

    if (pu < pl) {
	goto L80;
    }
    i__1 = pu;
    for (j = pl; j <= i__1; ++j) {
	qraux[j] = dnrm2_(n, &x[j * x_dim1 + 1], &c__1);
	work[j] = qraux[j];
/* L70: */
    }
L80:

/*     PERFORM THE HOUSEHOLDER REDUCTION OF X. */

    lup = min(*n,*p);
    i__1 = lup;
    for (l = 1; l <= i__1; ++l) {
	if (l < pl || l >= pu) {
	    goto L120;
	}

/*           LOCATE THE COLUMN OF LARGEST NORM AND BRING IT */
/*           INTO THE PIVOT POSITION. */

	maxnrm = 0.;
	maxj = l;
	i__2 = pu;
	for (j = l; j <= i__2; ++j) {
	    if (qraux[j] <= maxnrm) {
		goto L90;
	    }
	    maxnrm = qraux[j];
	    maxj = j;
L90:
/* L100: */
	    ;
	}
	if (maxj == l) {
	    goto L110;
	}
	dswap_(n, &x[l * x_dim1 + 1], &c__1, &x[maxj * x_dim1 + 1], &c__1);
	qraux[maxj] = qraux[l];
	work[maxj] = work[l];
	jp = jpvt[maxj];
	jpvt[maxj] = jpvt[l];
	jpvt[l] = jp;
L110:
L120:
	qraux[l] = 0.;
	if (l == *n) {
	    goto L190;
	}

/*           COMPUTE THE HOUSEHOLDER TRANSFORMATION FOR COLUMN L. */

	i__2 = *n - l + 1;
	nrmxl = dnrm2_(&i__2, &x[l + l * x_dim1], &c__1);
	if (nrmxl == 0.) {
	    goto L180;
	}
	if (x[l + l * x_dim1] != 0.) {
	    nrmxl = d_sign(&nrmxl, &x[l + l * x_dim1]);
	}
	i__2 = *n - l + 1;
	d__1 = 1. / nrmxl;
	dscal_(&i__2, &d__1, &x[l + l * x_dim1], &c__1);
	x[l + l * x_dim1] += 1.;

/*              APPLY THE TRANSFORMATION TO THE REMAINING COLUMNS, */
/*              UPDATING THE NORMS. */

	lp1 = l + 1;
	if (*p < lp1) {
	    goto L170;
	}
	i__2 = *p;
	for (j = lp1; j <= i__2; ++j) {
	    i__3 = *n - l + 1;
	    t = -ddot_(&i__3, &x[l + l * x_dim1], &c__1, &x[l + j * x_dim1], &
		    c__1) / x[l + l * x_dim1];
	    i__3 = *n - l + 1;
	    daxpy_(&i__3, &t, &x[l + l * x_dim1], &c__1, &x[l + j * x_dim1], &
		    c__1);
	    if (j < pl || j > pu) {
		goto L150;
	    }
	    if (qraux[j] == 0.) {
		goto L150;
	    }
/* Computing 2nd power */
	    d__2 = (d__1 = x[l + j * x_dim1], abs(d__1)) / qraux[j];
	    tt = 1. - d__2 * d__2;
	    tt = max(tt,0.);
	    t = tt;
/* Computing 2nd power */
	    d__1 = qraux[j] / work[j];
	    tt = tt * .05 * (d__1 * d__1) + 1.;
	    if (tt == 1.) {
		goto L130;
	    }
	    qraux[j] *= sqrt(t);
	    goto L140;
L130:
	    i__3 = *n - l;
	    qraux[j] = dnrm2_(&i__3, &x[l + 1 + j * x_dim1], &c__1);
	    work[j] = qraux[j];
L140:
L150:
/* L160: */
	    ;
	}
L170:

/*              SAVE THE TRANSFORMATION. */

	qraux[l] = x[l + l * x_dim1];
	x[l + l * x_dim1] = -nrmxl;
L180:
L190:
/* L200: */
	;
    }
    return 0;
} /* dqrdc_ */

