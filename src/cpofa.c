/* cpofa.f -- translated by f2c (version 12.02.01).
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

/* DECK CPOFA */
/* Subroutine */ int cpofa_(complex *a, integer *lda, integer *n, integer *
	info)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3, i__4;
    real r__1;
    complex q__1, q__2;

    /* Local variables */
    static integer j, k;
    static real s;
    static complex t;
    static integer jm1;
    extern /* Complex */ void cdotc_(complex *, integer *, complex *, integer 
	    *, complex *, integer *);

/* ***BEGIN PROLOGUE  CPOFA */
/* ***PURPOSE  Factor a complex Hermitian positive definite matrix. */
/* ***LIBRARY   SLATEC (LINPACK) */
/* ***CATEGORY  D2D1B */
/* ***TYPE      COMPLEX (SPOFA-S, DPOFA-D, CPOFA-C) */
/* ***KEYWORDS  LINEAR ALGEBRA, LINPACK, MATRIX FACTORIZATION, */
/*             POSITIVE DEFINITE */
/* ***AUTHOR  Moler, C. B., (U. of New Mexico) */
/* ***DESCRIPTION */

/*     CPOFA factors a complex Hermitian positive definite matrix. */

/*     CPOFA is usually called by CPOCO, but it can be called */
/*     directly with a saving in time if  RCOND  is not needed. */
/*     (Time for CPOCO) = (1 + 18/N)*(Time for CPOFA) . */

/*     On Entry */

/*        A       COMPLEX(LDA, N) */
/*                the Hermitian matrix to be factored.  Only the */
/*                diagonal and upper triangle are used. */

/*        LDA     INTEGER */
/*                the leading dimension of the array  A . */

/*        N       INTEGER */
/*                the order of the matrix  A . */

/*     On Return */

/*        A       an upper triangular matrix  R  so that  A = */
/*                CTRANS(R)*R where  CTRANS(R)  is the conjugate */
/*                transpose.  The strict lower triangle is unaltered. */
/*                If  INFO .NE. 0 , the factorization is not complete. */

/*        INFO    INTEGER */
/*                = 0  for normal return. */
/*                = K  signals an error condition.  The leading minor */
/*                     of order  K  is not positive definite. */

/* ***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W. */
/*                 Stewart, LINPACK Users' Guide, SIAM, 1979. */
/* ***ROUTINES CALLED  CDOTC */
/* ***REVISION HISTORY  (YYMMDD) */
/*   780814  DATE WRITTEN */
/*   890831  Modified array declarations.  (WRB) */
/*   890831  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900326  Removed duplicate information from DESCRIPTION section. */
/*           (WRB) */
/*   920501  Reformatted the REFERENCES section.  (WRB) */
/* ***END PROLOGUE  CPOFA */

/* ***FIRST EXECUTABLE STATEMENT  CPOFA */
    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;

    /* Function Body */
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	*info = j;
	s = 0.f;
	jm1 = j - 1;
	if (jm1 < 1) {
	    goto L20;
	}
	i__2 = jm1;
	for (k = 1; k <= i__2; ++k) {
	    i__3 = k + j * a_dim1;
	    i__4 = k - 1;
	    cdotc_(&q__2, &i__4, &a[k * a_dim1 + 1], &c__1, &a[j * a_dim1 + 1]
		    , &c__1);
	    q__1.r = a[i__3].r - q__2.r, q__1.i = a[i__3].i - q__2.i;
	    t.r = q__1.r, t.i = q__1.i;
	    c_div(&q__1, &t, &a[k + k * a_dim1]);
	    t.r = q__1.r, t.i = q__1.i;
	    i__3 = k + j * a_dim1;
	    a[i__3].r = t.r, a[i__3].i = t.i;
	    r_cnjg(&q__2, &t);
	    q__1.r = t.r * q__2.r - t.i * q__2.i, q__1.i = t.r * q__2.i + t.i 
		    * q__2.r;
	    s += q__1.r;
/* L10: */
	}
L20:
	i__2 = j + j * a_dim1;
	s = a[i__2].r - s;
	if (s <= 0.f || r_imag(&a[j + j * a_dim1]) != 0.f) {
	    goto L40;
	}
	i__2 = j + j * a_dim1;
	r__1 = sqrt(s);
	q__1.r = r__1, q__1.i = 0.f;
	a[i__2].r = q__1.r, a[i__2].i = q__1.i;
/* L30: */
    }
    *info = 0;
L40:
    return 0;
} /* cpofa_ */

