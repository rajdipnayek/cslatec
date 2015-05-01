/* cqrdc.f -- translated by f2c (version 12.02.01).
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
static complex c_b26 = {1.f,0.f};

/* DECK CQRDC */
/* Subroutine */ int cqrdc_(complex *x, integer *ldx, integer *n, integer *p, 
	complex *qraux, integer *jpvt, complex *work, integer *job)
{
    /* System generated locals */
    integer x_dim1, x_offset, i__1, i__2, i__3, i__4;
    real r__1, r__2, r__3, r__4;
    complex q__1, q__2, q__3;

    /* Local variables */
    static integer j, l;
    static complex t;
    static integer jj, jp, pl, pu;
    static real tt;
    static integer lp1, lup;
    static logical negj;
    static integer maxj;
    extern /* Subroutine */ int cscal_(integer *, complex *, complex *, 
	    integer *);
    extern /* Complex */ void cdotc_(complex *, integer *, complex *, integer 
	    *, complex *, integer *);
    extern /* Subroutine */ int cswap_(integer *, complex *, integer *, 
	    complex *, integer *);
    static logical swapj;
    extern /* Subroutine */ int caxpy_(integer *, complex *, complex *, 
	    integer *, complex *, integer *);
    static complex nrmxl;
    extern doublereal scnrm2_(integer *, complex *, integer *);
    static real maxnrm;

/* ***BEGIN PROLOGUE  CQRDC */
/* ***PURPOSE  Use Householder transformations to compute the QR */
/*            factorization of an N by P matrix.  Column pivoting is a */
/*            users option. */
/* ***LIBRARY   SLATEC (LINPACK) */
/* ***CATEGORY  D5 */
/* ***TYPE      COMPLEX (SQRDC-S, DQRDC-D, CQRDC-C) */
/* ***KEYWORDS  LINEAR ALGEBRA, LINPACK, MATRIX, ORTHOGONAL TRIANGULAR, */
/*             QR DECOMPOSITION */
/* ***AUTHOR  Stewart, G. W., (U. of Maryland) */
/* ***DESCRIPTION */

/*     CQRDC uses Householder transformations to compute the QR */
/*     factorization of an N by P matrix X.  Column pivoting */
/*     based on the 2-norms of the reduced columns may be */
/*     performed at the users option. */

/*     On Entry */

/*        X       COMPLEX(LDX,P), where LDX .GE. N. */
/*                X contains the matrix whose decomposition is to be */
/*                computed. */

/*        LDX     INTEGER. */
/*                LDX is the leading dimension of the array X. */

/*        N       INTEGER. */
/*                N is the number of rows of the matrix X. */

/*        P       INTEGER. */
/*                P is the number of columns of the matrix X. */

/*        JVPT    INTEGER(P). */
/*                JVPT contains integers that control the selection */
/*                of the pivot columns.  The K-th column X(K) of X */
/*                is placed in one of three classes according to the */
/*                value of JVPT(K). */

/*                   If JVPT(K) .GT. 0, then X(K) is an initial */
/*                                      column. */

/*                   If JVPT(K) .EQ. 0, then X(K) is a free column. */

/*                   If JVPT(K) .LT. 0, then X(K) is a final column. */

/*                Before the decomposition is computed, initial columns */
/*                are moved to the beginning of the array X and final */
/*                columns to the end.  Both initial and final columns */
/*                are frozen in place during the computation and only */
/*                free columns are moved.  At the K-th stage of the */
/*                reduction, if X(K) is occupied by a free column */
/*                it is interchanged with the free column of largest */
/*                reduced norm.  JVPT is not referenced if */
/*                JOB .EQ. 0. */

/*        WORK    COMPLEX(P). */
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
/*                which the unitary part of the decomposition */
/*                can be recovered.  Note that if pivoting has */
/*                been requested, the decomposition is not that */
/*                of the original matrix X but that of X */
/*                with its columns permuted as described by JVPT. */

/*        QRAUX   COMPLEX(P). */
/*                QRAUX contains further information required to recover */
/*                the unitary part of the decomposition. */

/*        JVPT    JVPT(K) contains the index of the column of the */
/*                original matrix that has been interchanged into */
/*                the K-th column, if pivoting was requested. */

/* ***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W. */
/*                 Stewart, LINPACK Users' Guide, SIAM, 1979. */
/* ***ROUTINES CALLED  CAXPY, CDOTC, CSCAL, CSWAP, SCNRM2 */
/* ***REVISION HISTORY  (YYMMDD) */
/*   780814  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890831  Modified array declarations.  (WRB) */
/*   890831  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900326  Removed duplicate information from DESCRIPTION section. */
/*           (WRB) */
/*   920501  Reformatted the REFERENCES section.  (WRB) */
/* ***END PROLOGUE  CQRDC */


/* ***FIRST EXECUTABLE STATEMENT  CQRDC */
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
	    cswap_(n, &x[pl * x_dim1 + 1], &c__1, &x[j * x_dim1 + 1], &c__1);
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
	cswap_(n, &x[pu * x_dim1 + 1], &c__1, &x[j * x_dim1 + 1], &c__1);
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
	i__2 = j;
	r__1 = scnrm2_(n, &x[j * x_dim1 + 1], &c__1);
	q__1.r = r__1, q__1.i = 0.f;
	qraux[i__2].r = q__1.r, qraux[i__2].i = q__1.i;
	i__2 = j;
	i__3 = j;
	work[i__2].r = qraux[i__3].r, work[i__2].i = qraux[i__3].i;
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

	maxnrm = 0.f;
	maxj = l;
	i__2 = pu;
	for (j = l; j <= i__2; ++j) {
	    i__3 = j;
	    if (qraux[i__3].r <= maxnrm) {
		goto L90;
	    }
	    i__3 = j;
	    maxnrm = qraux[i__3].r;
	    maxj = j;
L90:
/* L100: */
	    ;
	}
	if (maxj == l) {
	    goto L110;
	}
	cswap_(n, &x[l * x_dim1 + 1], &c__1, &x[maxj * x_dim1 + 1], &c__1);
	i__2 = maxj;
	i__3 = l;
	qraux[i__2].r = qraux[i__3].r, qraux[i__2].i = qraux[i__3].i;
	i__2 = maxj;
	i__3 = l;
	work[i__2].r = work[i__3].r, work[i__2].i = work[i__3].i;
	jp = jpvt[maxj];
	jpvt[maxj] = jpvt[l];
	jpvt[l] = jp;
L110:
L120:
	i__2 = l;
	qraux[i__2].r = 0.f, qraux[i__2].i = 0.f;
	if (l == *n) {
	    goto L190;
	}

/*           COMPUTE THE HOUSEHOLDER TRANSFORMATION FOR COLUMN L. */

	i__2 = *n - l + 1;
	r__1 = scnrm2_(&i__2, &x[l + l * x_dim1], &c__1);
	q__1.r = r__1, q__1.i = 0.f;
	nrmxl.r = q__1.r, nrmxl.i = q__1.i;
	if ((r__1 = nrmxl.r, dabs(r__1)) + (r__2 = r_imag(&nrmxl), dabs(r__2))
		 == 0.f) {
	    goto L180;
	}
	i__2 = l + l * x_dim1;
	if ((r__1 = x[i__2].r, dabs(r__1)) + (r__2 = r_imag(&x[l + l * x_dim1]
		), dabs(r__2)) != 0.f) {
	    r__3 = c_abs(&nrmxl);
	    i__3 = l + l * x_dim1;
	    r__4 = c_abs(&x[l + l * x_dim1]);
	    q__2.r = x[i__3].r / r__4, q__2.i = x[i__3].i / r__4;
	    q__1.r = r__3 * q__2.r, q__1.i = r__3 * q__2.i;
	    nrmxl.r = q__1.r, nrmxl.i = q__1.i;
	}
	i__2 = *n - l + 1;
	c_div(&q__1, &c_b26, &nrmxl);
	cscal_(&i__2, &q__1, &x[l + l * x_dim1], &c__1);
	i__2 = l + l * x_dim1;
	i__3 = l + l * x_dim1;
	q__1.r = x[i__3].r + 1.f, q__1.i = x[i__3].i + 0.f;
	x[i__2].r = q__1.r, x[i__2].i = q__1.i;

/*              APPLY THE TRANSFORMATION TO THE REMAINING COLUMNS, */
/*              UPDATING THE NORMS. */

	lp1 = l + 1;
	if (*p < lp1) {
	    goto L170;
	}
	i__2 = *p;
	for (j = lp1; j <= i__2; ++j) {
	    i__3 = *n - l + 1;
	    cdotc_(&q__3, &i__3, &x[l + l * x_dim1], &c__1, &x[l + j * x_dim1]
		    , &c__1);
	    q__2.r = -q__3.r, q__2.i = -q__3.i;
	    c_div(&q__1, &q__2, &x[l + l * x_dim1]);
	    t.r = q__1.r, t.i = q__1.i;
	    i__3 = *n - l + 1;
	    caxpy_(&i__3, &t, &x[l + l * x_dim1], &c__1, &x[l + j * x_dim1], &
		    c__1);
	    if (j < pl || j > pu) {
		goto L150;
	    }
	    i__3 = j;
	    if ((r__1 = qraux[i__3].r, dabs(r__1)) + (r__2 = r_imag(&qraux[j])
		    , dabs(r__2)) == 0.f) {
		goto L150;
	    }
	    i__3 = j;
/* Computing 2nd power */
	    r__1 = c_abs(&x[l + j * x_dim1]) / qraux[i__3].r;
	    tt = 1.f - r__1 * r__1;
	    tt = dmax(tt,0.f);
	    q__1.r = tt, q__1.i = 0.f;
	    t.r = q__1.r, t.i = q__1.i;
	    i__3 = j;
	    i__4 = j;
/* Computing 2nd power */
	    r__1 = qraux[i__3].r / work[i__4].r;
	    tt = tt * .05f * (r__1 * r__1) + 1.f;
	    if (tt == 1.f) {
		goto L130;
	    }
	    i__3 = j;
	    i__4 = j;
	    c_sqrt(&q__2, &t);
	    q__1.r = qraux[i__4].r * q__2.r - qraux[i__4].i * q__2.i, q__1.i =
		     qraux[i__4].r * q__2.i + qraux[i__4].i * q__2.r;
	    qraux[i__3].r = q__1.r, qraux[i__3].i = q__1.i;
	    goto L140;
L130:
	    i__3 = j;
	    i__4 = *n - l;
	    r__1 = scnrm2_(&i__4, &x[l + 1 + j * x_dim1], &c__1);
	    q__1.r = r__1, q__1.i = 0.f;
	    qraux[i__3].r = q__1.r, qraux[i__3].i = q__1.i;
	    i__3 = j;
	    i__4 = j;
	    work[i__3].r = qraux[i__4].r, work[i__3].i = qraux[i__4].i;
L140:
L150:
/* L160: */
	    ;
	}
L170:

/*              SAVE THE TRANSFORMATION. */

	i__2 = l;
	i__3 = l + l * x_dim1;
	qraux[i__2].r = x[i__3].r, qraux[i__2].i = x[i__3].i;
	i__2 = l + l * x_dim1;
	q__1.r = -nrmxl.r, q__1.i = -nrmxl.i;
	x[i__2].r = q__1.r, x[i__2].i = q__1.i;
L180:
L190:
/* L200: */
	;
    }
    return 0;
} /* cqrdc_ */

