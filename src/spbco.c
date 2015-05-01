/* spbco.f -- translated by f2c (version 12.02.01).
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

/* DECK SPBCO */
/* Subroutine */ int spbco_(real *abd, integer *lda, integer *n, integer *m, 
	real *rcond, real *z__, integer *info)
{
    /* System generated locals */
    integer abd_dim1, abd_offset, i__1, i__2, i__3;
    real r__1, r__2;

    /* Local variables */
    static integer i__, j, k, l;
    static real s, t;
    static integer j2, kb, la, lb;
    static real ek;
    static integer lm;
    static real sm, wk;
    static integer mu, kp1;
    static real wkm;
    extern doublereal sdot_(integer *, real *, integer *, real *, integer *);
    extern /* Subroutine */ int spbfa_(real *, integer *, integer *, integer *
	    , integer *), sscal_(integer *, real *, real *, integer *);
    static real anorm;
    extern doublereal sasum_(integer *, real *, integer *);
    static real ynorm;
    extern /* Subroutine */ int saxpy_(integer *, real *, real *, integer *, 
	    real *, integer *);

/* ***BEGIN PROLOGUE  SPBCO */
/* ***PURPOSE  Factor a real symmetric positive definite matrix stored in */
/*            band form and estimate the condition number of the matrix. */
/* ***LIBRARY   SLATEC (LINPACK) */
/* ***CATEGORY  D2B2 */
/* ***TYPE      SINGLE PRECISION (SPBCO-S, DPBCO-D, CPBCO-C) */
/* ***KEYWORDS  BANDED, CONDITION NUMBER, LINEAR ALGEBRA, LINPACK, */
/*             MATRIX FACTORIZATION, POSITIVE DEFINITE */
/* ***AUTHOR  Moler, C. B., (U. of New Mexico) */
/* ***DESCRIPTION */

/*     SPBCO factors a real symmetric positive definite matrix */
/*     stored in band form and estimates the condition of the matrix. */

/*     If  RCOND  is not needed, SPBFA is slightly faster. */
/*     To solve  A*X = B , follow SPBCO by SPBSL. */
/*     To compute  INVERSE(A)*C , follow SPBCO by SPBSL. */
/*     To compute  DETERMINANT(A) , follow SPBCO by SPBDI. */

/*     On Entry */

/*        ABD     REAL(LDA, N) */
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
/*                form, so that  A = TRANS(R)*R . */
/*                If  INFO .NE. 0 , the factorization is not complete. */

/*        RCOND   REAL */
/*                an estimate of the reciprocal condition of  A . */
/*                For the system  A*X = B , relative perturbations */
/*                in  A  and  B  of size  EPSILON  may cause */
/*                relative perturbations in  X  of size  EPSILON/RCOND . */
/*                If  RCOND  is so small that the logical expression */
/*                           1.0 + RCOND .EQ. 1.0 */
/*                is true, then  A  may be singular to working */
/*                precision.  In particular,  RCOND  is zero  if */
/*                exact singularity is detected or the estimate */
/*                underflows.  If INFO .NE. 0 , RCOND is unchanged. */

/*        Z       REAL(N) */
/*                a work vector whose contents are usually unimportant. */
/*                If  A  is singular to working precision, then  Z  is */
/*                an approximate null vector in the sense that */
/*                NORM(A*Z) = RCOND*NORM(A)*NORM(Z) . */
/*                If  INFO .NE. 0 , Z  is unchanged. */

/*        INFO    INTEGER */
/*                = 0  for normal return. */
/*                = K  signals an error condition.  The leading minor */
/*                     of order  K  is not positive definite. */

/*     Band Storage */

/*           If  A  is a symmetric positive definite band matrix, */
/*           the following program segment will set up the input. */

/*                   M = (band width above diagonal) */
/*                   DO 20 J = 1, N */
/*                      I1 = MAX(1, J-M) */
/*                      DO 10 I = I1, J */
/*                         K = I-J+M+1 */
/*                         ABD(K,J) = A(I,J) */
/*                10    CONTINUE */
/*                20 CONTINUE */

/*           This uses  M + 1  rows of  A , except for the  M by M */
/*           upper left triangle, which is ignored. */

/*     Example:  If the original matrix is */

/*           11 12 13  0  0  0 */
/*           12 22 23 24  0  0 */
/*           13 23 33 34 35  0 */
/*            0 24 34 44 45 46 */
/*            0  0 35 45 55 56 */
/*            0  0  0 46 56 66 */

/*     then  N = 6 , M = 2  and  ABD  should contain */

/*            *  * 13 24 35 46 */
/*            * 12 23 34 45 56 */
/*           11 22 33 44 55 66 */

/* ***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W. */
/*                 Stewart, LINPACK Users' Guide, SIAM, 1979. */
/* ***ROUTINES CALLED  SASUM, SAXPY, SDOT, SPBFA, SSCAL */
/* ***REVISION HISTORY  (YYMMDD) */
/*   780814  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890831  Modified array declarations.  (WRB) */
/*   890831  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900326  Removed duplicate information from DESCRIPTION section. */
/*           (WRB) */
/*   920501  Reformatted the REFERENCES section.  (WRB) */
/* ***END PROLOGUE  SPBCO */


/*     FIND NORM OF A */

/* ***FIRST EXECUTABLE STATEMENT  SPBCO */
    /* Parameter adjustments */
    abd_dim1 = *lda;
    abd_offset = 1 + abd_dim1;
    abd -= abd_offset;
    --z__;

    /* Function Body */
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
/* Computing MIN */
	i__2 = j, i__3 = *m + 1;
	l = min(i__2,i__3);
/* Computing MAX */
	i__2 = *m + 2 - j;
	mu = max(i__2,1);
	z__[j] = sasum_(&l, &abd[mu + j * abd_dim1], &c__1);
	k = j - l;
	if (*m < mu) {
	    goto L20;
	}
	i__2 = *m;
	for (i__ = mu; i__ <= i__2; ++i__) {
	    ++k;
	    z__[k] += (r__1 = abd[i__ + j * abd_dim1], dabs(r__1));
/* L10: */
	}
L20:
/* L30: */
	;
    }
    anorm = 0.f;
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
/* Computing MAX */
	r__1 = anorm, r__2 = z__[j];
	anorm = dmax(r__1,r__2);
/* L40: */
    }

/*     FACTOR */

    spbfa_(&abd[abd_offset], lda, n, m, info);
    if (*info != 0) {
	goto L180;
    }

/*        RCOND = 1/(NORM(A)*(ESTIMATE OF NORM(INVERSE(A)))) . */
/*        ESTIMATE = NORM(Z)/NORM(Y) WHERE  A*Z = Y  AND  A*Y = E . */
/*        THE COMPONENTS OF  E  ARE CHOSEN TO CAUSE MAXIMUM LOCAL */
/*        GROWTH IN THE ELEMENTS OF W  WHERE  TRANS(R)*W = E . */
/*        THE VECTORS ARE FREQUENTLY RESCALED TO AVOID OVERFLOW. */

/*        SOLVE TRANS(R)*W = E */

    ek = 1.f;
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	z__[j] = 0.f;
/* L50: */
    }
    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
	if (z__[k] != 0.f) {
	    r__1 = -z__[k];
	    ek = r_sign(&ek, &r__1);
	}
	if ((r__1 = ek - z__[k], dabs(r__1)) <= abd[*m + 1 + k * abd_dim1]) {
	    goto L60;
	}
	s = abd[*m + 1 + k * abd_dim1] / (r__1 = ek - z__[k], dabs(r__1));
	sscal_(n, &s, &z__[1], &c__1);
	ek = s * ek;
L60:
	wk = ek - z__[k];
	wkm = -ek - z__[k];
	s = dabs(wk);
	sm = dabs(wkm);
	wk /= abd[*m + 1 + k * abd_dim1];
	wkm /= abd[*m + 1 + k * abd_dim1];
	kp1 = k + 1;
/* Computing MIN */
	i__2 = k + *m;
	j2 = min(i__2,*n);
	i__ = *m + 1;
	if (kp1 > j2) {
	    goto L100;
	}
	i__2 = j2;
	for (j = kp1; j <= i__2; ++j) {
	    --i__;
	    sm += (r__1 = z__[j] + wkm * abd[i__ + j * abd_dim1], dabs(r__1));
	    z__[j] += wk * abd[i__ + j * abd_dim1];
	    s += (r__1 = z__[j], dabs(r__1));
/* L70: */
	}
	if (s >= sm) {
	    goto L90;
	}
	t = wkm - wk;
	wk = wkm;
	i__ = *m + 1;
	i__2 = j2;
	for (j = kp1; j <= i__2; ++j) {
	    --i__;
	    z__[j] += t * abd[i__ + j * abd_dim1];
/* L80: */
	}
L90:
L100:
	z__[k] = wk;
/* L110: */
    }
    s = 1.f / sasum_(n, &z__[1], &c__1);
    sscal_(n, &s, &z__[1], &c__1);

/*        SOLVE  R*Y = W */

    i__1 = *n;
    for (kb = 1; kb <= i__1; ++kb) {
	k = *n + 1 - kb;
	if ((r__1 = z__[k], dabs(r__1)) <= abd[*m + 1 + k * abd_dim1]) {
	    goto L120;
	}
	s = abd[*m + 1 + k * abd_dim1] / (r__1 = z__[k], dabs(r__1));
	sscal_(n, &s, &z__[1], &c__1);
L120:
	z__[k] /= abd[*m + 1 + k * abd_dim1];
/* Computing MIN */
	i__2 = k - 1;
	lm = min(i__2,*m);
	la = *m + 1 - lm;
	lb = k - lm;
	t = -z__[k];
	saxpy_(&lm, &t, &abd[la + k * abd_dim1], &c__1, &z__[lb], &c__1);
/* L130: */
    }
    s = 1.f / sasum_(n, &z__[1], &c__1);
    sscal_(n, &s, &z__[1], &c__1);

    ynorm = 1.f;

/*        SOLVE TRANS(R)*V = Y */

    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
/* Computing MIN */
	i__2 = k - 1;
	lm = min(i__2,*m);
	la = *m + 1 - lm;
	lb = k - lm;
	z__[k] -= sdot_(&lm, &abd[la + k * abd_dim1], &c__1, &z__[lb], &c__1);
	if ((r__1 = z__[k], dabs(r__1)) <= abd[*m + 1 + k * abd_dim1]) {
	    goto L140;
	}
	s = abd[*m + 1 + k * abd_dim1] / (r__1 = z__[k], dabs(r__1));
	sscal_(n, &s, &z__[1], &c__1);
	ynorm = s * ynorm;
L140:
	z__[k] /= abd[*m + 1 + k * abd_dim1];
/* L150: */
    }
    s = 1.f / sasum_(n, &z__[1], &c__1);
    sscal_(n, &s, &z__[1], &c__1);
    ynorm = s * ynorm;

/*        SOLVE  R*Z = W */

    i__1 = *n;
    for (kb = 1; kb <= i__1; ++kb) {
	k = *n + 1 - kb;
	if ((r__1 = z__[k], dabs(r__1)) <= abd[*m + 1 + k * abd_dim1]) {
	    goto L160;
	}
	s = abd[*m + 1 + k * abd_dim1] / (r__1 = z__[k], dabs(r__1));
	sscal_(n, &s, &z__[1], &c__1);
	ynorm = s * ynorm;
L160:
	z__[k] /= abd[*m + 1 + k * abd_dim1];
/* Computing MIN */
	i__2 = k - 1;
	lm = min(i__2,*m);
	la = *m + 1 - lm;
	lb = k - lm;
	t = -z__[k];
	saxpy_(&lm, &t, &abd[la + k * abd_dim1], &c__1, &z__[lb], &c__1);
/* L170: */
    }
/*        MAKE ZNORM = 1.0 */
    s = 1.f / sasum_(n, &z__[1], &c__1);
    sscal_(n, &s, &z__[1], &c__1);
    ynorm = s * ynorm;

    if (anorm != 0.f) {
	*rcond = ynorm / anorm;
    }
    if (anorm == 0.f) {
	*rcond = 0.f;
    }
L180:
    return 0;
} /* spbco_ */

