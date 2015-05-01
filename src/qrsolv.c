/* qrsolv.f -- translated by f2c (version 12.02.01).
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

/* DECK QRSOLV */
/* Subroutine */ int qrsolv_(integer *n, real *r__, integer *ldr, integer *
	ipvt, real *diag, real *qtb, real *x, real *sigma, real *wa)
{
    /* Initialized data */

    static real p5 = .5f;
    static real p25 = .25f;
    static real zero = 0.f;

    /* System generated locals */
    integer r_dim1, r_offset, i__1, i__2, i__3;
    real r__1, r__2;

    /* Local variables */
    static integer i__, j, k, l, jp1, kp1;
    static real tan__, cos__, sin__, sum, temp, cotan;
    static integer nsing;
    static real qtbpj;

/* ***BEGIN PROLOGUE  QRSOLV */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to SNLS1 and SNLS1E */
/* ***LIBRARY   SLATEC */
/* ***TYPE      SINGLE PRECISION (QRSOLV-S, DQRSLV-D) */
/* ***AUTHOR  (UNKNOWN) */
/* ***DESCRIPTION */

/*     Given an M by N matrix A, an N by N diagonal matrix D, */
/*     and an M-vector B, the problem is to determine an X which */
/*     solves the system */

/*           A*X = B ,     D*X = 0 , */

/*     in the least squares sense. */

/*     This subroutine completes the solution of the problem */
/*     if it is provided with the necessary information from the */
/*     QR factorization, with column pivoting, of A. That is, if */
/*     A*P = Q*R, where P is a permutation matrix, Q has orthogonal */
/*     columns, and R is an upper triangular matrix with diagonal */
/*     elements of nonincreasing magnitude, then QRSOLV expects */
/*     the full upper triangle of R, the permutation matrix P, */
/*     and the first N components of (Q TRANSPOSE)*B. The system */
/*     A*X = B, D*X = 0, is then equivalent to */

/*                  T       T */
/*           R*Z = Q *B ,  P *D*P*Z = 0 , */

/*     where X = P*Z. If this system does not have full rank, */
/*     then a least squares solution is obtained. On output QRSOLV */
/*     also provides an upper triangular matrix S such that */

/*            T   T               T */
/*           P *(A *A + D*D)*P = S *S . */

/*     S is computed within QRSOLV and may be of separate interest. */

/*     The subroutine statement is */

/*       SUBROUTINE QRSOLV(N,R,LDR,IPVT,DIAG,QTB,X,SIGMA,WA) */

/*     where */

/*       N is a positive integer input variable set to the order of R. */

/*       R is an N by N array. On input the full upper triangle */
/*         must contain the full upper triangle of the matrix R. */
/*         On output the full upper triangle is unaltered, and the */
/*         strict lower triangle contains the strict upper triangle */
/*         (transposed) of the upper triangular matrix S. */

/*       LDR is a positive integer input variable not less than N */
/*         which specifies the leading dimension of the array R. */

/*       IPVT is an integer input array of length N which defines the */
/*         permutation matrix P such that A*P = Q*R. Column J of P */
/*         is column IPVT(J) of the identity matrix. */

/*       DIAG is an input array of length N which must contain the */
/*         diagonal elements of the matrix D. */

/*       QTB is an input array of length N which must contain the first */
/*         N elements of the vector (Q TRANSPOSE)*B. */

/*       X is an output array of length N which contains the least */
/*         squares solution of the system A*X = B, D*X = 0. */

/*       SIGMA is an output array of length N which contains the */
/*         diagonal elements of the upper triangular matrix S. */

/*       WA is a work array of length N. */

/* ***SEE ALSO  SNLS1, SNLS1E */
/* ***ROUTINES CALLED  (NONE) */
/* ***REVISION HISTORY  (YYMMDD) */
/*   800301  DATE WRITTEN */
/*   890831  Modified array declarations.  (WRB) */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900326  Removed duplicate information from DESCRIPTION section. */
/*           (WRB) */
/*   900328  Added TYPE section.  (WRB) */
/* ***END PROLOGUE  QRSOLV */
    /* Parameter adjustments */
    r_dim1 = *ldr;
    r_offset = 1 + r_dim1;
    r__ -= r_offset;
    --ipvt;
    --diag;
    --qtb;
    --x;
    --sigma;
    --wa;

    /* Function Body */
/* ***FIRST EXECUTABLE STATEMENT  QRSOLV */
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	i__2 = *n;
	for (i__ = j; i__ <= i__2; ++i__) {
	    r__[i__ + j * r_dim1] = r__[j + i__ * r_dim1];
/* L10: */
	}
	x[j] = r__[j + j * r_dim1];
	wa[j] = qtb[j];
/* L20: */
    }

/*     ELIMINATE THE DIAGONAL MATRIX D USING A GIVENS ROTATION. */

    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {

/*        PREPARE THE ROW OF D TO BE ELIMINATED, LOCATING THE */
/*        DIAGONAL ELEMENT USING P FROM THE QR FACTORIZATION. */

	l = ipvt[j];
	if (diag[l] == zero) {
	    goto L90;
	}
	i__2 = *n;
	for (k = j; k <= i__2; ++k) {
	    sigma[k] = zero;
/* L30: */
	}
	sigma[j] = diag[l];

/*        THE TRANSFORMATIONS TO ELIMINATE THE ROW OF D */
/*        MODIFY ONLY A SINGLE ELEMENT OF (Q TRANSPOSE)*B */
/*        BEYOND THE FIRST N, WHICH IS INITIALLY ZERO. */

	qtbpj = zero;
	i__2 = *n;
	for (k = j; k <= i__2; ++k) {

/*           DETERMINE A GIVENS ROTATION WHICH ELIMINATES THE */
/*           APPROPRIATE ELEMENT IN THE CURRENT ROW OF D. */

	    if (sigma[k] == zero) {
		goto L70;
	    }
	    if ((r__1 = r__[k + k * r_dim1], dabs(r__1)) >= (r__2 = sigma[k], 
		    dabs(r__2))) {
		goto L40;
	    }
	    cotan = r__[k + k * r_dim1] / sigma[k];
/* Computing 2nd power */
	    r__1 = cotan;
	    sin__ = p5 / sqrt(p25 + p25 * (r__1 * r__1));
	    cos__ = sin__ * cotan;
	    goto L50;
L40:
	    tan__ = sigma[k] / r__[k + k * r_dim1];
/* Computing 2nd power */
	    r__1 = tan__;
	    cos__ = p5 / sqrt(p25 + p25 * (r__1 * r__1));
	    sin__ = cos__ * tan__;
L50:

/*           COMPUTE THE MODIFIED DIAGONAL ELEMENT OF R AND */
/*           THE MODIFIED ELEMENT OF ((Q TRANSPOSE)*B,0). */

	    r__[k + k * r_dim1] = cos__ * r__[k + k * r_dim1] + sin__ * sigma[
		    k];
	    temp = cos__ * wa[k] + sin__ * qtbpj;
	    qtbpj = -sin__ * wa[k] + cos__ * qtbpj;
	    wa[k] = temp;

/*           ACCUMULATE THE TRANSFORMATION IN THE ROW OF S. */

	    kp1 = k + 1;
	    if (*n < kp1) {
		goto L70;
	    }
	    i__3 = *n;
	    for (i__ = kp1; i__ <= i__3; ++i__) {
		temp = cos__ * r__[i__ + k * r_dim1] + sin__ * sigma[i__];
		sigma[i__] = -sin__ * r__[i__ + k * r_dim1] + cos__ * sigma[
			i__];
		r__[i__ + k * r_dim1] = temp;
/* L60: */
	    }
L70:
/* L80: */
	    ;
	}
L90:

/*        STORE THE DIAGONAL ELEMENT OF S AND RESTORE */
/*        THE CORRESPONDING DIAGONAL ELEMENT OF R. */

	sigma[j] = r__[j + j * r_dim1];
	r__[j + j * r_dim1] = x[j];
/* L100: */
    }

/*     SOLVE THE TRIANGULAR SYSTEM FOR Z. IF THE SYSTEM IS */
/*     SINGULAR, THEN OBTAIN A LEAST SQUARES SOLUTION. */

    nsing = *n;
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	if (sigma[j] == zero && nsing == *n) {
	    nsing = j - 1;
	}
	if (nsing < *n) {
	    wa[j] = zero;
	}
/* L110: */
    }
    if (nsing < 1) {
	goto L150;
    }
    i__1 = nsing;
    for (k = 1; k <= i__1; ++k) {
	j = nsing - k + 1;
	sum = zero;
	jp1 = j + 1;
	if (nsing < jp1) {
	    goto L130;
	}
	i__2 = nsing;
	for (i__ = jp1; i__ <= i__2; ++i__) {
	    sum += r__[i__ + j * r_dim1] * wa[i__];
/* L120: */
	}
L130:
	wa[j] = (wa[j] - sum) / sigma[j];
/* L140: */
    }
L150:

/*     PERMUTE THE COMPONENTS OF Z BACK TO COMPONENTS OF X. */

    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	l = ipvt[j];
	x[l] = wa[j];
/* L160: */
    }
    return 0;

/*     LAST CARD OF SUBROUTINE QRSOLV. */

} /* qrsolv_ */

