/* dchdc.f -- translated by f2c (version 12.02.01).
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

/* DECK DCHDC */
/* Subroutine */ int dchdc_(doublereal *a, integer *lda, integer *p, 
	doublereal *work, integer *jpvt, integer *job, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3;

    /* Local variables */
    static integer j, k, l, kb, jp, pl, jt, pu, km1, kp1, plp1;
    static logical negk;
    static integer maxl;
    static doublereal temp;
    extern /* Subroutine */ int dswap_(integer *, doublereal *, integer *, 
	    doublereal *, integer *);
    static logical swapk;
    extern /* Subroutine */ int daxpy_(integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *);
    static doublereal maxdia;

/* ***BEGIN PROLOGUE  DCHDC */
/* ***PURPOSE  Compute the Cholesky decomposition of a positive definite */
/*            matrix.  A pivoting option allows the user to estimate the */
/*            condition number of a positive definite matrix or determine */
/*            the rank of a positive semidefinite matrix. */
/* ***LIBRARY   SLATEC (LINPACK) */
/* ***CATEGORY  D2B1B */
/* ***TYPE      DOUBLE PRECISION (SCHDC-S, DCHDC-D, CCHDC-C) */
/* ***KEYWORDS  CHOLESKY DECOMPOSITION, LINEAR ALGEBRA, LINPACK, MATRIX, */
/*             POSITIVE DEFINITE */
/* ***AUTHOR  Dongarra, J., (ANL) */
/*           Stewart, G. W., (U. of Maryland) */
/* ***DESCRIPTION */

/*     DCHDC computes the Cholesky decomposition of a positive definite */
/*     matrix.  A pivoting option allows the user to estimate the */
/*     condition of a positive definite matrix or determine the rank */
/*     of a positive semidefinite matrix. */

/*     On Entry */

/*         A      DOUBLE PRECISION(LDA,P). */
/*                A contains the matrix whose decomposition is to */
/*                be computed.  Only the upper half of A need be stored. */
/*                The lower part of the array A is not referenced. */

/*         LDA    INTEGER. */
/*                LDA is the leading dimension of the array A. */

/*         P      INTEGER. */
/*                P is the order of the matrix. */

/*         WORK   DOUBLE PRECISION. */
/*                WORK is a work array. */

/*         JPVT   INTEGER(P). */
/*                JPVT contains integers that control the selection */
/*                of the pivot elements, if pivoting has been requested. */
/*                Each diagonal element A(K,K) */
/*                is placed in one of three classes according to the */
/*                value of JPVT(K). */

/*                   If JPVT(K) .GT. 0, then X(K) is an initial */
/*                                      element. */

/*                   If JPVT(K) .EQ. 0, then X(K) is a free element. */

/*                   If JPVT(K) .LT. 0, then X(K) is a final element. */

/*                Before the decomposition is computed, initial elements */
/*                are moved by symmetric row and column interchanges to */
/*                the beginning of the array A and final */
/*                elements to the end.  Both initial and final elements */
/*                are frozen in place during the computation and only */
/*                free elements are moved.  At the K-th stage of the */
/*                reduction, if A(K,K) is occupied by a free element */
/*                it is interchanged with the largest free element */
/*                A(L,L) with L .GE. K.  JPVT is not referenced if */
/*                JOB .EQ. 0. */

/*        JOB     INTEGER. */
/*                JOB is an integer that initiates column pivoting. */
/*                If JOB .EQ. 0, no pivoting is done. */
/*                If JOB .NE. 0, pivoting is done. */

/*     On Return */

/*         A      A contains in its upper half the Cholesky factor */
/*                of the matrix A as it has been permuted by pivoting. */

/*         JPVT   JPVT(J) contains the index of the diagonal element */
/*                of a that was moved into the J-th position, */
/*                provided pivoting was requested. */

/*         INFO   contains the index of the last positive diagonal */
/*                element of the Cholesky factor. */

/*     For positive definite matrices INFO = P is the normal return. */
/*     For pivoting with positive semidefinite matrices INFO will */
/*     in general be less than P.  However, INFO may be greater than */
/*     the rank of A, since rounding error can cause an otherwise zero */
/*     element to be positive.  Indefinite systems will always cause */
/*     INFO to be less than P. */

/* ***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W. */
/*                 Stewart, LINPACK Users' Guide, SIAM, 1979. */
/* ***ROUTINES CALLED  DAXPY, DSWAP */
/* ***REVISION HISTORY  (YYMMDD) */
/*   790319  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890831  Modified array declarations.  (WRB) */
/*   890831  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900326  Removed duplicate information from DESCRIPTION section. */
/*           (WRB) */
/*   920501  Reformatted the REFERENCES section.  (WRB) */
/* ***END PROLOGUE  DCHDC */

/* ***FIRST EXECUTABLE STATEMENT  DCHDC */
    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --work;
    --jpvt;

    /* Function Body */
    pl = 1;
    pu = 0;
    *info = *p;
    if (*job == 0) {
	goto L160;
    }

/*        PIVOTING HAS BEEN REQUESTED. REARRANGE THE */
/*        THE ELEMENTS ACCORDING TO JPVT. */

    i__1 = *p;
    for (k = 1; k <= i__1; ++k) {
	swapk = jpvt[k] > 0;
	negk = jpvt[k] < 0;
	jpvt[k] = k;
	if (negk) {
	    jpvt[k] = -jpvt[k];
	}
	if (! swapk) {
	    goto L60;
	}
	if (k == pl) {
	    goto L50;
	}
	i__2 = pl - 1;
	dswap_(&i__2, &a[k * a_dim1 + 1], &c__1, &a[pl * a_dim1 + 1], &c__1);
	temp = a[k + k * a_dim1];
	a[k + k * a_dim1] = a[pl + pl * a_dim1];
	a[pl + pl * a_dim1] = temp;
	plp1 = pl + 1;
	if (*p < plp1) {
	    goto L40;
	}
	i__2 = *p;
	for (j = plp1; j <= i__2; ++j) {
	    if (j >= k) {
		goto L10;
	    }
	    temp = a[pl + j * a_dim1];
	    a[pl + j * a_dim1] = a[j + k * a_dim1];
	    a[j + k * a_dim1] = temp;
	    goto L20;
L10:
	    if (j == k) {
		goto L20;
	    }
	    temp = a[k + j * a_dim1];
	    a[k + j * a_dim1] = a[pl + j * a_dim1];
	    a[pl + j * a_dim1] = temp;
L20:
/* L30: */
	    ;
	}
L40:
	jpvt[k] = jpvt[pl];
	jpvt[pl] = k;
L50:
	++pl;
L60:
/* L70: */
	;
    }
    pu = *p;
    if (*p < pl) {
	goto L150;
    }
    i__1 = *p;
    for (kb = pl; kb <= i__1; ++kb) {
	k = *p - kb + pl;
	if (jpvt[k] >= 0) {
	    goto L130;
	}
	jpvt[k] = -jpvt[k];
	if (pu == k) {
	    goto L120;
	}
	i__2 = k - 1;
	dswap_(&i__2, &a[k * a_dim1 + 1], &c__1, &a[pu * a_dim1 + 1], &c__1);
	temp = a[k + k * a_dim1];
	a[k + k * a_dim1] = a[pu + pu * a_dim1];
	a[pu + pu * a_dim1] = temp;
	kp1 = k + 1;
	if (*p < kp1) {
	    goto L110;
	}
	i__2 = *p;
	for (j = kp1; j <= i__2; ++j) {
	    if (j >= pu) {
		goto L80;
	    }
	    temp = a[k + j * a_dim1];
	    a[k + j * a_dim1] = a[j + pu * a_dim1];
	    a[j + pu * a_dim1] = temp;
	    goto L90;
L80:
	    if (j == pu) {
		goto L90;
	    }
	    temp = a[k + j * a_dim1];
	    a[k + j * a_dim1] = a[pu + j * a_dim1];
	    a[pu + j * a_dim1] = temp;
L90:
/* L100: */
	    ;
	}
L110:
	jt = jpvt[k];
	jpvt[k] = jpvt[pu];
	jpvt[pu] = jt;
L120:
	--pu;
L130:
/* L140: */
	;
    }
L150:
L160:
    i__1 = *p;
    for (k = 1; k <= i__1; ++k) {

/*        REDUCTION LOOP. */

	maxdia = a[k + k * a_dim1];
	kp1 = k + 1;
	maxl = k;

/*        DETERMINE THE PIVOT ELEMENT. */

	if (k < pl || k >= pu) {
	    goto L190;
	}
	i__2 = pu;
	for (l = kp1; l <= i__2; ++l) {
	    if (a[l + l * a_dim1] <= maxdia) {
		goto L170;
	    }
	    maxdia = a[l + l * a_dim1];
	    maxl = l;
L170:
/* L180: */
	    ;
	}
L190:

/*        QUIT IF THE PIVOT ELEMENT IS NOT POSITIVE. */

	if (maxdia > 0.) {
	    goto L200;
	}
	*info = k - 1;
	goto L280;
L200:
	if (k == maxl) {
	    goto L210;
	}

/*           START THE PIVOTING AND UPDATE JPVT. */

	km1 = k - 1;
	dswap_(&km1, &a[k * a_dim1 + 1], &c__1, &a[maxl * a_dim1 + 1], &c__1);
	a[maxl + maxl * a_dim1] = a[k + k * a_dim1];
	a[k + k * a_dim1] = maxdia;
	jp = jpvt[maxl];
	jpvt[maxl] = jpvt[k];
	jpvt[k] = jp;
L210:

/*        REDUCTION STEP. PIVOTING IS CONTAINED ACROSS THE ROWS. */

	work[k] = sqrt(a[k + k * a_dim1]);
	a[k + k * a_dim1] = work[k];
	if (*p < kp1) {
	    goto L260;
	}
	i__2 = *p;
	for (j = kp1; j <= i__2; ++j) {
	    if (k == maxl) {
		goto L240;
	    }
	    if (j >= maxl) {
		goto L220;
	    }
	    temp = a[k + j * a_dim1];
	    a[k + j * a_dim1] = a[j + maxl * a_dim1];
	    a[j + maxl * a_dim1] = temp;
	    goto L230;
L220:
	    if (j == maxl) {
		goto L230;
	    }
	    temp = a[k + j * a_dim1];
	    a[k + j * a_dim1] = a[maxl + j * a_dim1];
	    a[maxl + j * a_dim1] = temp;
L230:
L240:
	    a[k + j * a_dim1] /= work[k];
	    work[j] = a[k + j * a_dim1];
	    temp = -a[k + j * a_dim1];
	    i__3 = j - k;
	    daxpy_(&i__3, &temp, &work[kp1], &c__1, &a[kp1 + j * a_dim1], &
		    c__1);
/* L250: */
	}
L260:
/* L270: */
	;
    }
L280:
    return 0;
} /* dchdc_ */

