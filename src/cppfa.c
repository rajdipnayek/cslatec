/* cppfa.f -- translated by f2c (version 12.02.01).
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

/* DECK CPPFA */
/* Subroutine */ int cppfa_(complex *ap, integer *n, integer *info)
{
    /* System generated locals */
    integer i__1, i__2, i__3, i__4;
    real r__1;
    complex q__1, q__2;

    /* Local variables */
    static integer j, k;
    static real s;
    static complex t;
    static integer jj, kj, kk, jm1;
    extern /* Complex */ void cdotc_(complex *, integer *, complex *, integer 
	    *, complex *, integer *);

/* ***BEGIN PROLOGUE  CPPFA */
/* ***PURPOSE  Factor a complex Hermitian positive definite matrix stored */
/*            in packed form. */
/* ***LIBRARY   SLATEC (LINPACK) */
/* ***CATEGORY  D2D1B */
/* ***TYPE      COMPLEX (SPPFA-S, DPPFA-D, CPPFA-C) */
/* ***KEYWORDS  LINEAR ALGEBRA, LINPACK, MATRIX FACTORIZATION, PACKED, */
/*             POSITIVE DEFINITE */
/* ***AUTHOR  Moler, C. B., (U. of New Mexico) */
/* ***DESCRIPTION */

/*     CPPFA factors a complex Hermitian positive definite matrix */
/*     stored in packed form. */

/*     CPPFA is usually called by CPPCO, but it can be called */
/*     directly with a saving in time if  RCOND  is not needed. */
/*     (Time for CPPCO) = (1 + 18/N)*(Time for CPPFA) . */

/*     On Entry */

/*        AP      COMPLEX (N*(N+1)/2) */
/*                the packed form of a Hermitian matrix  A .  The */
/*                columns of the upper triangle are stored sequentially */
/*                in a one-dimensional array of length  N*(N+1)/2 . */
/*                See comments below for details. */

/*        N       INTEGER */
/*                the order of the matrix  A . */

/*     On Return */

/*        AP      an upper triangular matrix  R , stored in packed */
/*                form, so that  A = CTRANS(R)*R . */

/*        INFO    INTEGER */
/*                = 0  for normal return. */
/*                = K  If the leading minor of order  K  is not */
/*                     positive definite. */


/*     Packed Storage */

/*          The following program segment will pack the upper */
/*          triangle of a Hermitian matrix. */

/*                K = 0 */
/*                DO 20 J = 1, N */
/*                   DO 10 I = 1, J */
/*                      K = K + 1 */
/*                      AP(K) = A(I,J) */
/*             10    CONTINUE */
/*             20 CONTINUE */

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
/* ***END PROLOGUE  CPPFA */

/* ***FIRST EXECUTABLE STATEMENT  CPPFA */
    /* Parameter adjustments */
    --ap;

    /* Function Body */
    jj = 0;
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	*info = j;
	s = 0.f;
	jm1 = j - 1;
	kj = jj;
	kk = 0;
	if (jm1 < 1) {
	    goto L20;
	}
	i__2 = jm1;
	for (k = 1; k <= i__2; ++k) {
	    ++kj;
	    i__3 = kj;
	    i__4 = k - 1;
	    cdotc_(&q__2, &i__4, &ap[kk + 1], &c__1, &ap[jj + 1], &c__1);
	    q__1.r = ap[i__3].r - q__2.r, q__1.i = ap[i__3].i - q__2.i;
	    t.r = q__1.r, t.i = q__1.i;
	    kk += k;
	    c_div(&q__1, &t, &ap[kk]);
	    t.r = q__1.r, t.i = q__1.i;
	    i__3 = kj;
	    ap[i__3].r = t.r, ap[i__3].i = t.i;
	    r_cnjg(&q__2, &t);
	    q__1.r = t.r * q__2.r - t.i * q__2.i, q__1.i = t.r * q__2.i + t.i 
		    * q__2.r;
	    s += q__1.r;
/* L10: */
	}
L20:
	jj += j;
	i__2 = jj;
	s = ap[i__2].r - s;
	if (s <= 0.f || r_imag(&ap[jj]) != 0.f) {
	    goto L40;
	}
	i__2 = jj;
	r__1 = sqrt(s);
	q__1.r = r__1, q__1.i = 0.f;
	ap[i__2].r = q__1.r, ap[i__2].i = q__1.i;
/* L30: */
    }
    *info = 0;
L40:
    return 0;
} /* cppfa_ */

