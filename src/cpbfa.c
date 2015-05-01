/* cpbfa.f -- translated by f2c (version 12.02.01).
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

/* DECK CPBFA */
/* Subroutine */ int cpbfa_(complex *abd, integer *lda, integer *n, integer *
	m, integer *info)
{
    /* System generated locals */
    integer abd_dim1, abd_offset, i__1, i__2, i__3, i__4;
    real r__1;
    complex q__1, q__2;

    /* Local variables */
    static integer j, k;
    static real s;
    static complex t;
    static integer ik, jk, mu;
    extern /* Complex */ void cdotc_(complex *, integer *, complex *, integer 
	    *, complex *, integer *);

/* ***BEGIN PROLOGUE  CPBFA */
/* ***PURPOSE  Factor a complex Hermitian positive definite matrix stored */
/*            in band form. */
/* ***LIBRARY   SLATEC (LINPACK) */
/* ***CATEGORY  D2D2 */
/* ***TYPE      COMPLEX (SPBFA-S, DPBFA-D, CPBFA-C) */
/* ***KEYWORDS  BANDED, LINEAR ALGEBRA, LINPACK, MATRIX FACTORIZATION, */
/*             POSITIVE DEFINITE */
/* ***AUTHOR  Moler, C. B., (U. of New Mexico) */
/* ***DESCRIPTION */

/*     CPBFA factors a complex Hermitian positive definite matrix */
/*     stored in band form. */

/*     CPBFA is usually called by CPBCO, but it can be called */
/*     directly with a saving in time if  RCOND  is not needed. */

/*     On Entry */

/*        ABD     COMPLEX(LDA, N) */
/*                the matrix to be factored.  The columns of the upper */
/*                triangle are stored in the columns of ABD and the */
/*                diagonals of the upper triangle are stored in the */
/*                rows of ABD .  See the comments below for details. */

/*        LDA     INTEGER */
/*                the leading dimension of the array  ABD . */
/*                LDA must be .GE. M + 1 . */

/*        N       INTEGER */
/*                the order of the matrix  A . */

/*        M       INTEGER */
/*                the number of diagonals above the main diagonal. */
/*                0 .LE. M .LT. N . */

/*     On Return */

/*        ABD     an upper triangular matrix  R , stored in band */
/*                form, so that  A = CTRANS(R)*R . */

/*        INFO    INTEGER */
/*                = 0  for normal return. */
/*                = K  if the leading minor of order  K  is not */
/*                     positive definite. */

/*     Band Storage */

/*           If  A  is a Hermitian positive definite band matrix, */
/*           the following program segment will set up the input. */

/*                   M = (band width above diagonal) */
/*                   DO 20 J = 1, N */
/*                      I1 = MAX(1, J-M) */
/*                      DO 10 I = I1, J */
/*                         K = I-J+M+1 */
/*                         ABD(K,J) = A(I,J) */
/*                10    CONTINUE */
/*                20 CONTINUE */

/* ***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W. */
/*                 Stewart, LINPACK Users' Guide, SIAM, 1979. */
/* ***ROUTINES CALLED  CDOTC */
/* ***REVISION HISTORY  (YYMMDD) */
/*   780814  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890831  Modified array declarations.  (WRB) */
/*   890831  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900326  Removed duplicate information from DESCRIPTION section. */
/*           (WRB) */
/*   920501  Reformatted the REFERENCES section.  (WRB) */
/* ***END PROLOGUE  CPBFA */

/* ***FIRST EXECUTABLE STATEMENT  CPBFA */
    /* Parameter adjustments */
    abd_dim1 = *lda;
    abd_offset = 1 + abd_dim1;
    abd -= abd_offset;

    /* Function Body */
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	*info = j;
	s = 0.f;
	ik = *m + 1;
/* Computing MAX */
	i__2 = j - *m;
	jk = max(i__2,1);
/* Computing MAX */
	i__2 = *m + 2 - j;
	mu = max(i__2,1);
	if (*m < mu) {
	    goto L20;
	}
	i__2 = *m;
	for (k = mu; k <= i__2; ++k) {
	    i__3 = k + j * abd_dim1;
	    i__4 = k - mu;
	    cdotc_(&q__2, &i__4, &abd[ik + jk * abd_dim1], &c__1, &abd[mu + j 
		    * abd_dim1], &c__1);
	    q__1.r = abd[i__3].r - q__2.r, q__1.i = abd[i__3].i - q__2.i;
	    t.r = q__1.r, t.i = q__1.i;
	    c_div(&q__1, &t, &abd[*m + 1 + jk * abd_dim1]);
	    t.r = q__1.r, t.i = q__1.i;
	    i__3 = k + j * abd_dim1;
	    abd[i__3].r = t.r, abd[i__3].i = t.i;
	    r_cnjg(&q__2, &t);
	    q__1.r = t.r * q__2.r - t.i * q__2.i, q__1.i = t.r * q__2.i + t.i 
		    * q__2.r;
	    s += q__1.r;
	    --ik;
	    ++jk;
/* L10: */
	}
L20:
	i__2 = *m + 1 + j * abd_dim1;
	s = abd[i__2].r - s;
	if (s <= 0.f || r_imag(&abd[*m + 1 + j * abd_dim1]) != 0.f) {
	    goto L40;
	}
	i__2 = *m + 1 + j * abd_dim1;
	r__1 = sqrt(s);
	q__1.r = r__1, q__1.i = 0.f;
	abd[i__2].r = q__1.r, abd[i__2].i = q__1.i;
/* L30: */
    }
    *info = 0;
L40:
    return 0;
} /* cpbfa_ */

