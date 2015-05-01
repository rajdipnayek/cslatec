/* schex.f -- translated by f2c (version 12.02.01).
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

/* DECK SCHEX */
/* Subroutine */ int schex_(real *r__, integer *ldr, integer *p, integer *k, 
	integer *l, real *z__, integer *ldz, integer *nz, real *c__, real *s, 
	integer *job)
{
    /* System generated locals */
    integer r_dim1, r_offset, z_dim1, z_offset, i__1, i__2, i__3;

    /* Local variables */
    static integer i__, j;
    static real t;
    static integer ii, jj, il, iu, km1, lm1, kp1, lmk;
    extern /* Subroutine */ int srotg_(real *, real *, real *, real *);

/* ***BEGIN PROLOGUE  SCHEX */
/* ***PURPOSE  Update the Cholesky factorization  A=TRANS(R)*R  of A */
/*            positive definite matrix A of order P under diagonal */
/*            permutations of the form TRANS(E)*A*E, where E is a */
/*            permutation matrix. */
/* ***LIBRARY   SLATEC (LINPACK) */
/* ***CATEGORY  D7B */
/* ***TYPE      SINGLE PRECISION (SCHEX-S, DCHEX-D, CCHEX-C) */
/* ***KEYWORDS  CHOLESKY DECOMPOSITION, EXCHANGE, LINEAR ALGEBRA, LINPACK, */
/*             MATRIX, POSITIVE DEFINITE */
/* ***AUTHOR  Stewart, G. W., (U. of Maryland) */
/* ***DESCRIPTION */

/*     SCHEX updates the Cholesky factorization */

/*                   A = TRANS(R)*R */

/*     of a positive definite matrix A of order P under diagonal */
/*     permutations of the form */

/*                   TRANS(E)*A*E */

/*     where E is a permutation matrix.  Specifically, given */
/*     an upper triangular matrix R and a permutation matrix */
/*     E (which is specified by K, L, and JOB), SCHEX determines */
/*     an orthogonal matrix U such that */

/*                           U*R*E = RR, */

/*     where RR is upper triangular.  At the users option, the */
/*     transformation U will be multiplied into the array Z. */
/*     If A = TRANS(X)*X, so that R is the triangular part of the */
/*     QR factorization of X, then RR is the triangular part of the */
/*     QR factorization of X*E, i.e., X with its columns permuted. */
/*     For a less terse description of what SCHEX does and how */
/*     it may be applied, see the LINPACK guide. */

/*     The matrix Q is determined as the product U(L-K)*...*U(1) */
/*     of plane rotations of the form */

/*                           (    C(I)       S(I) ) */
/*                           (                    ) , */
/*                           (    -S(I)      C(I) ) */

/*     where C(I) is real.  The rows these rotations operate on */
/*     are described below. */

/*     There are two types of permutations, which are determined */
/*     by the value of JOB. */

/*     1. Right circular shift (JOB = 1). */

/*         The columns are rearranged in the following order. */

/*                1,...,K-1,L,K,K+1,...,L-1,L+1,...,P. */

/*         U is the product of L-K rotations U(I), where U(I) */
/*         acts in the (L-I,L-I+1)-plane. */

/*     2. Left circular shift (JOB = 2). */
/*         The columns are rearranged in the following order */

/*                1,...,K-1,K+1,K+2,...,L,K,L+1,...,P. */

/*         U is the product of L-K rotations U(I), where U(I) */
/*         acts in the (K+I-1,K+I)-plane. */

/*     On Entry */

/*         R      REAL(LDR,P), where LDR .GE. P. */
/*                R contains the upper triangular factor */
/*                that is to be updated.  Elements of R */
/*                below the diagonal are not referenced. */

/*         LDR    INTEGER. */
/*                LDR is the leading dimension of the array R. */

/*         P      INTEGER. */
/*                P is the order of the matrix R. */

/*         K      INTEGER. */
/*                K is the first column to be permuted. */

/*         L      INTEGER. */
/*                L is the last column to be permuted. */
/*                L must be strictly greater than K. */

/*         Z      REAL(LDZ,NZ), where LDZ.GE.P. */
/*                Z is an array of NZ P-vectors into which the */
/*                transformation U is multiplied.  Z is */
/*                not referenced if NZ = 0. */

/*         LDZ    INTEGER. */
/*                LDZ is the leading dimension of the array Z. */

/*         NZ     INTEGER. */
/*                NZ is the number of columns of the matrix Z. */

/*         JOB    INTEGER. */
/*                JOB determines the type of permutation. */
/*                       JOB = 1  right circular shift. */
/*                       JOB = 2  left circular shift. */

/*     On Return */

/*         R      contains the updated factor. */

/*         Z      contains the updated matrix Z. */

/*         C      REAL(P). */
/*                C contains the cosines of the transforming rotations. */

/*         S      REAL(P). */
/*                S contains the sines of the transforming rotations. */

/* ***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W. */
/*                 Stewart, LINPACK Users' Guide, SIAM, 1979. */
/* ***ROUTINES CALLED  SROTG */
/* ***REVISION HISTORY  (YYMMDD) */
/*   780814  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890531  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900326  Removed duplicate information from DESCRIPTION section. */
/*           (WRB) */
/*   920501  Reformatted the REFERENCES section.  (WRB) */
/* ***END PROLOGUE  SCHEX */


/*     INITIALIZE */

/* ***FIRST EXECUTABLE STATEMENT  SCHEX */
    /* Parameter adjustments */
    r_dim1 = *ldr;
    r_offset = 1 + r_dim1;
    r__ -= r_offset;
    z_dim1 = *ldz;
    z_offset = 1 + z_dim1;
    z__ -= z_offset;
    --c__;
    --s;

    /* Function Body */
    km1 = *k - 1;
    kp1 = *k + 1;
    lmk = *l - *k;
    lm1 = *l - 1;

/*     PERFORM THE APPROPRIATE TASK. */

    switch (*job) {
	case 1:  goto L10;
	case 2:  goto L130;
    }

/*     RIGHT CIRCULAR SHIFT. */

L10:

/*        REORDER THE COLUMNS. */

    i__1 = *l;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ii = *l - i__ + 1;
	s[i__] = r__[ii + *l * r_dim1];
/* L20: */
    }
    i__1 = lm1;
    for (jj = *k; jj <= i__1; ++jj) {
	j = lm1 - jj + *k;
	i__2 = j;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    r__[i__ + (j + 1) * r_dim1] = r__[i__ + j * r_dim1];
/* L30: */
	}
	r__[j + 1 + (j + 1) * r_dim1] = 0.f;
/* L40: */
    }
    if (*k == 1) {
	goto L60;
    }
    i__1 = km1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ii = *l - i__ + 1;
	r__[i__ + *k * r_dim1] = s[ii];
/* L50: */
    }
L60:

/*        CALCULATE THE ROTATIONS. */

    t = s[1];
    i__1 = lmk;
    for (i__ = 1; i__ <= i__1; ++i__) {
	srotg_(&s[i__ + 1], &t, &c__[i__], &s[i__]);
	t = s[i__ + 1];
/* L70: */
    }
    r__[*k + *k * r_dim1] = t;
    i__1 = *p;
    for (j = kp1; j <= i__1; ++j) {
/* Computing MAX */
	i__2 = 1, i__3 = *l - j + 1;
	il = max(i__2,i__3);
	i__2 = lmk;
	for (ii = il; ii <= i__2; ++ii) {
	    i__ = *l - ii;
	    t = c__[ii] * r__[i__ + j * r_dim1] + s[ii] * r__[i__ + 1 + j * 
		    r_dim1];
	    r__[i__ + 1 + j * r_dim1] = c__[ii] * r__[i__ + 1 + j * r_dim1] - 
		    s[ii] * r__[i__ + j * r_dim1];
	    r__[i__ + j * r_dim1] = t;
/* L80: */
	}
/* L90: */
    }

/*        IF REQUIRED, APPLY THE TRANSFORMATIONS TO Z. */

    if (*nz < 1) {
	goto L120;
    }
    i__1 = *nz;
    for (j = 1; j <= i__1; ++j) {
	i__2 = lmk;
	for (ii = 1; ii <= i__2; ++ii) {
	    i__ = *l - ii;
	    t = c__[ii] * z__[i__ + j * z_dim1] + s[ii] * z__[i__ + 1 + j * 
		    z_dim1];
	    z__[i__ + 1 + j * z_dim1] = c__[ii] * z__[i__ + 1 + j * z_dim1] - 
		    s[ii] * z__[i__ + j * z_dim1];
	    z__[i__ + j * z_dim1] = t;
/* L100: */
	}
/* L110: */
    }
L120:
    goto L260;

/*     LEFT CIRCULAR SHIFT */

L130:

/*        REORDER THE COLUMNS */

    i__1 = *k;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ii = lmk + i__;
	s[ii] = r__[i__ + *k * r_dim1];
/* L140: */
    }
    i__1 = lm1;
    for (j = *k; j <= i__1; ++j) {
	i__2 = j;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    r__[i__ + j * r_dim1] = r__[i__ + (j + 1) * r_dim1];
/* L150: */
	}
	jj = j - km1;
	s[jj] = r__[j + 1 + (j + 1) * r_dim1];
/* L160: */
    }
    i__1 = *k;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ii = lmk + i__;
	r__[i__ + *l * r_dim1] = s[ii];
/* L170: */
    }
    i__1 = *l;
    for (i__ = kp1; i__ <= i__1; ++i__) {
	r__[i__ + *l * r_dim1] = 0.f;
/* L180: */
    }

/*        REDUCTION LOOP. */

    i__1 = *p;
    for (j = *k; j <= i__1; ++j) {
	if (j == *k) {
	    goto L200;
	}

/*              APPLY THE ROTATIONS. */

/* Computing MIN */
	i__2 = j - 1, i__3 = *l - 1;
	iu = min(i__2,i__3);
	i__2 = iu;
	for (i__ = *k; i__ <= i__2; ++i__) {
	    ii = i__ - *k + 1;
	    t = c__[ii] * r__[i__ + j * r_dim1] + s[ii] * r__[i__ + 1 + j * 
		    r_dim1];
	    r__[i__ + 1 + j * r_dim1] = c__[ii] * r__[i__ + 1 + j * r_dim1] - 
		    s[ii] * r__[i__ + j * r_dim1];
	    r__[i__ + j * r_dim1] = t;
/* L190: */
	}
L200:
	if (j >= *l) {
	    goto L210;
	}
	jj = j - *k + 1;
	t = s[jj];
	srotg_(&r__[j + j * r_dim1], &t, &c__[jj], &s[jj]);
L210:
/* L220: */
	;
    }

/*        APPLY THE ROTATIONS TO Z. */

    if (*nz < 1) {
	goto L250;
    }
    i__1 = *nz;
    for (j = 1; j <= i__1; ++j) {
	i__2 = lm1;
	for (i__ = *k; i__ <= i__2; ++i__) {
	    ii = i__ - km1;
	    t = c__[ii] * z__[i__ + j * z_dim1] + s[ii] * z__[i__ + 1 + j * 
		    z_dim1];
	    z__[i__ + 1 + j * z_dim1] = c__[ii] * z__[i__ + 1 + j * z_dim1] - 
		    s[ii] * z__[i__ + j * z_dim1];
	    z__[i__ + j * z_dim1] = t;
/* L230: */
	}
/* L240: */
    }
L250:
L260:
    return 0;
} /* schex_ */

