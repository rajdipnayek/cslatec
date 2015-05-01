/* rwupdt.f -- translated by f2c (version 12.02.01).
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

/* DECK RWUPDT */
/* Subroutine */ int rwupdt_(integer *n, real *r__, integer *ldr, real *w, 
	real *b, real *alpha, real *cos__, real *sin__)
{
    /* Initialized data */

    static real one = 1.f;
    static real p5 = .5f;
    static real p25 = .25f;
    static real zero = 0.f;

    /* System generated locals */
    integer r_dim1, r_offset, i__1, i__2;
    real r__1;

    /* Local variables */
    static integer i__, j, jm1;
    static real tan__, temp, rowj, cotan;

/* ***BEGIN PROLOGUE  RWUPDT */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to SNLS1 and SNLS1E */
/* ***LIBRARY   SLATEC */
/* ***TYPE      SINGLE PRECISION (RWUPDT-S, DWUPDT-D) */
/* ***AUTHOR  (UNKNOWN) */
/* ***DESCRIPTION */

/*     Given an N by N upper triangular matrix R, this subroutine */
/*     computes the QR decomposition of the matrix formed when a row */
/*     is added to R. If the row is specified by the vector W, then */
/*     RWUPDT determines an orthogonal matrix Q such that when the */
/*     N+1 by N matrix composed of R augmented by W is premultiplied */
/*     by (Q TRANSPOSE), the resulting matrix is upper trapezoidal. */
/*     The orthogonal matrix Q is the product of N transformations */

/*           G(1)*G(2)* ... *G(N) */

/*     where G(I) is a Givens rotation in the (I,N+1) plane which */
/*     eliminates elements in the I-th plane. RWUPDT also */
/*     computes the product (Q TRANSPOSE)*C where C is the */
/*     (N+1)-vector (b,alpha). Q itself is not accumulated, rather */
/*     the information to recover the G rotations is supplied. */

/*     The subroutine statement is */

/*       SUBROUTINE RWUPDT(N,R,LDR,W,B,ALPHA,COS,SIN) */

/*     where */

/*       N is a positive integer input variable set to the order of R. */

/*       R is an N by N array. On input the upper triangular part of */
/*         R must contain the matrix to be updated. On output R */
/*         contains the updated triangular matrix. */

/*       LDR is a positive integer input variable not less than N */
/*         which specifies the leading dimension of the array R. */

/*       W is an input array of length N which must contain the row */
/*         vector to be added to R. */

/*       B is an array of length N. On input B must contain the */
/*         first N elements of the vector C. On output B contains */
/*         the first N elements of the vector (Q TRANSPOSE)*C. */

/*       ALPHA is a variable. On input ALPHA must contain the */
/*         (N+1)-st element of the vector C. On output ALPHA contains */
/*         the (N+1)-st element of the vector (Q TRANSPOSE)*C. */

/*       COS is an output array of length N which contains the */
/*         cosines of the transforming Givens rotations. */

/*       SIN is an output array of length N which contains the */
/*         sines of the transforming Givens rotations. */

/* ***SEE ALSO  SNLS1, SNLS1E */
/* ***ROUTINES CALLED  (NONE) */
/* ***REVISION HISTORY  (YYMMDD) */
/*   800301  DATE WRITTEN */
/*   890831  Modified array declarations.  (WRB) */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900326  Removed duplicate information from DESCRIPTION section. */
/*           (WRB) */
/*   900328  Added TYPE section.  (WRB) */
/* ***END PROLOGUE  RWUPDT */
    /* Parameter adjustments */
    r_dim1 = *ldr;
    r_offset = 1 + r_dim1;
    r__ -= r_offset;
    --w;
    --b;
    --cos__;
    --sin__;

    /* Function Body */
/* ***FIRST EXECUTABLE STATEMENT  RWUPDT */
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	rowj = w[j];
	jm1 = j - 1;

/*        APPLY THE PREVIOUS TRANSFORMATIONS TO */
/*        R(I,J), I=1,2,...,J-1, AND TO W(J). */

	if (jm1 < 1) {
	    goto L20;
	}
	i__2 = jm1;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    temp = cos__[i__] * r__[i__ + j * r_dim1] + sin__[i__] * rowj;
	    rowj = -sin__[i__] * r__[i__ + j * r_dim1] + cos__[i__] * rowj;
	    r__[i__ + j * r_dim1] = temp;
/* L10: */
	}
L20:

/*        DETERMINE A GIVENS ROTATION WHICH ELIMINATES W(J). */

	cos__[j] = one;
	sin__[j] = zero;
	if (rowj == zero) {
	    goto L50;
	}
	if ((r__1 = r__[j + j * r_dim1], dabs(r__1)) >= dabs(rowj)) {
	    goto L30;
	}
	cotan = r__[j + j * r_dim1] / rowj;
/* Computing 2nd power */
	r__1 = cotan;
	sin__[j] = p5 / sqrt(p25 + p25 * (r__1 * r__1));
	cos__[j] = sin__[j] * cotan;
	goto L40;
L30:
	tan__ = rowj / r__[j + j * r_dim1];
/* Computing 2nd power */
	r__1 = tan__;
	cos__[j] = p5 / sqrt(p25 + p25 * (r__1 * r__1));
	sin__[j] = cos__[j] * tan__;
L40:

/*        APPLY THE CURRENT TRANSFORMATION TO R(J,J), B(J), AND ALPHA. */

	r__[j + j * r_dim1] = cos__[j] * r__[j + j * r_dim1] + sin__[j] * 
		rowj;
	temp = cos__[j] * b[j] + sin__[j] * *alpha;
	*alpha = -sin__[j] * b[j] + cos__[j] * *alpha;
	b[j] = temp;
L50:
/* L60: */
	;
    }
    return 0;

/*     LAST CARD OF SUBROUTINE RWUPDT. */

} /* rwupdt_ */

