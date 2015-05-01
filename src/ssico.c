/* ssico.f -- translated by f2c (version 12.02.01).
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

/* DECK SSICO */
/* Subroutine */ int ssico_(real *a, integer *lda, integer *n, integer *kpvt, 
	real *rcond, real *z__)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;
    real r__1, r__2;

    /* Local variables */
    static integer i__, j, k;
    static real s, t, ak, bk, ek;
    static integer kp, ks, jm1, kps;
    static real akm1, bkm1;
    static integer info;
    extern doublereal sdot_(integer *, real *, integer *, real *, integer *);
    static real denom;
    extern /* Subroutine */ int ssifa_(real *, integer *, integer *, integer *
	    , integer *), sscal_(integer *, real *, real *, integer *);
    static real anorm;
    extern doublereal sasum_(integer *, real *, integer *);
    static real ynorm;
    extern /* Subroutine */ int saxpy_(integer *, real *, real *, integer *, 
	    real *, integer *);

/* ***BEGIN PROLOGUE  SSICO */
/* ***PURPOSE  Factor a symmetric matrix by elimination with symmetric */
/*            pivoting and estimate the condition number of the matrix. */
/* ***LIBRARY   SLATEC (LINPACK) */
/* ***CATEGORY  D2B1A */
/* ***TYPE      SINGLE PRECISION (SSICO-S, DSICO-D, CHICO-C, CSICO-C) */
/* ***KEYWORDS  CONDITION NUMBER, LINEAR ALGEBRA, LINPACK, */
/*             MATRIX FACTORIZATION, SYMMETRIC */
/* ***AUTHOR  Moler, C. B., (U. of New Mexico) */
/* ***DESCRIPTION */

/*     SSICO factors a real symmetric matrix by elimination with */
/*     symmetric pivoting and estimates the condition of the matrix. */

/*     If  RCOND  is not needed, SSIFA is slightly faster. */
/*     To solve  A*X = B , follow SSICO by SSISL. */
/*     To compute  INVERSE(A)*C , follow SSICO by SSISL. */
/*     To compute  INVERSE(A) , follow SSICO by SSIDI. */
/*     To compute  DETERMINANT(A) , follow SSICO by SSIDI. */
/*     To compute  INERTIA(A), follow SSICO by SSIDI. */

/*     On Entry */

/*        A       REAL(LDA, N) */
/*                the symmetric matrix to be factored. */
/*                Only the diagonal and upper triangle are used. */

/*        LDA     INTEGER */
/*                the leading dimension of the array  A . */

/*        N       INTEGER */
/*                the order of the matrix  A . */

/*     Output */

/*        A       a block diagonal matrix and the multipliers which */
/*                were used to obtain it. */
/*                The factorization can be written  A = U*D*TRANS(U) */
/*                where  U  is a product of permutation and unit */
/*                upper triangular matrices , TRANS(U) is the */
/*                transpose of  U , and  D  is block diagonal */
/*                with 1 by 1 and 2 by 2 blocks. */

/*        KPVT    INTEGER(N) */
/*                an integer vector of pivot indices. */

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
/*                underflows. */

/*        Z       REAL(N) */
/*                a work vector whose contents are usually unimportant. */
/*                If  A  is close to a singular matrix, then  Z  is */
/*                an approximate null vector in the sense that */
/*                NORM(A*Z) = RCOND*NORM(A)*NORM(Z) . */

/* ***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W. */
/*                 Stewart, LINPACK Users' Guide, SIAM, 1979. */
/* ***ROUTINES CALLED  SASUM, SAXPY, SDOT, SSCAL, SSIFA */
/* ***REVISION HISTORY  (YYMMDD) */
/*   780814  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890831  Modified array declarations.  (WRB) */
/*   891107  Modified routine equivalence list.  (WRB) */
/*   891107  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900326  Removed duplicate information from DESCRIPTION section. */
/*           (WRB) */
/*   920501  Reformatted the REFERENCES section.  (WRB) */
/* ***END PROLOGUE  SSICO */


/*     FIND NORM OF A USING ONLY UPPER HALF */

/* ***FIRST EXECUTABLE STATEMENT  SSICO */
    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --kpvt;
    --z__;

    /* Function Body */
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	z__[j] = sasum_(&j, &a[j * a_dim1 + 1], &c__1);
	jm1 = j - 1;
	if (jm1 < 1) {
	    goto L20;
	}
	i__2 = jm1;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    z__[i__] += (r__1 = a[i__ + j * a_dim1], dabs(r__1));
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

    ssifa_(&a[a_offset], lda, n, &kpvt[1], &info);

/*     RCOND = 1/(NORM(A)*(ESTIMATE OF NORM(INVERSE(A)))) . */
/*     ESTIMATE = NORM(Z)/NORM(Y) WHERE  A*Z = Y  AND  A*Y = E . */
/*     THE COMPONENTS OF  E  ARE CHOSEN TO CAUSE MAXIMUM LOCAL */
/*     GROWTH IN THE ELEMENTS OF W  WHERE  U*D*W = E . */
/*     THE VECTORS ARE FREQUENTLY RESCALED TO AVOID OVERFLOW. */

/*     SOLVE U*D*W = E */

    ek = 1.f;
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	z__[j] = 0.f;
/* L50: */
    }
    k = *n;
L60:
    if (k == 0) {
	goto L120;
    }
    ks = 1;
    if (kpvt[k] < 0) {
	ks = 2;
    }
    kp = (i__1 = kpvt[k], abs(i__1));
    kps = k + 1 - ks;
    if (kp == kps) {
	goto L70;
    }
    t = z__[kps];
    z__[kps] = z__[kp];
    z__[kp] = t;
L70:
    if (z__[k] != 0.f) {
	ek = r_sign(&ek, &z__[k]);
    }
    z__[k] += ek;
    i__1 = k - ks;
    saxpy_(&i__1, &z__[k], &a[k * a_dim1 + 1], &c__1, &z__[1], &c__1);
    if (ks == 1) {
	goto L80;
    }
    if (z__[k - 1] != 0.f) {
	ek = r_sign(&ek, &z__[k - 1]);
    }
    z__[k - 1] += ek;
    i__1 = k - ks;
    saxpy_(&i__1, &z__[k - 1], &a[(k - 1) * a_dim1 + 1], &c__1, &z__[1], &
	    c__1);
L80:
    if (ks == 2) {
	goto L100;
    }
    if ((r__1 = z__[k], dabs(r__1)) <= (r__2 = a[k + k * a_dim1], dabs(r__2)))
	     {
	goto L90;
    }
    s = (r__1 = a[k + k * a_dim1], dabs(r__1)) / (r__2 = z__[k], dabs(r__2));
    sscal_(n, &s, &z__[1], &c__1);
    ek = s * ek;
L90:
    if (a[k + k * a_dim1] != 0.f) {
	z__[k] /= a[k + k * a_dim1];
    }
    if (a[k + k * a_dim1] == 0.f) {
	z__[k] = 1.f;
    }
    goto L110;
L100:
    ak = a[k + k * a_dim1] / a[k - 1 + k * a_dim1];
    akm1 = a[k - 1 + (k - 1) * a_dim1] / a[k - 1 + k * a_dim1];
    bk = z__[k] / a[k - 1 + k * a_dim1];
    bkm1 = z__[k - 1] / a[k - 1 + k * a_dim1];
    denom = ak * akm1 - 1.f;
    z__[k] = (akm1 * bk - bkm1) / denom;
    z__[k - 1] = (ak * bkm1 - bk) / denom;
L110:
    k -= ks;
    goto L60;
L120:
    s = 1.f / sasum_(n, &z__[1], &c__1);
    sscal_(n, &s, &z__[1], &c__1);

/*     SOLVE TRANS(U)*Y = W */

    k = 1;
L130:
    if (k > *n) {
	goto L160;
    }
    ks = 1;
    if (kpvt[k] < 0) {
	ks = 2;
    }
    if (k == 1) {
	goto L150;
    }
    i__1 = k - 1;
    z__[k] += sdot_(&i__1, &a[k * a_dim1 + 1], &c__1, &z__[1], &c__1);
    if (ks == 2) {
	i__1 = k - 1;
	z__[k + 1] += sdot_(&i__1, &a[(k + 1) * a_dim1 + 1], &c__1, &z__[1], &
		c__1);
    }
    kp = (i__1 = kpvt[k], abs(i__1));
    if (kp == k) {
	goto L140;
    }
    t = z__[k];
    z__[k] = z__[kp];
    z__[kp] = t;
L140:
L150:
    k += ks;
    goto L130;
L160:
    s = 1.f / sasum_(n, &z__[1], &c__1);
    sscal_(n, &s, &z__[1], &c__1);

    ynorm = 1.f;

/*     SOLVE U*D*V = Y */

    k = *n;
L170:
    if (k == 0) {
	goto L230;
    }
    ks = 1;
    if (kpvt[k] < 0) {
	ks = 2;
    }
    if (k == ks) {
	goto L190;
    }
    kp = (i__1 = kpvt[k], abs(i__1));
    kps = k + 1 - ks;
    if (kp == kps) {
	goto L180;
    }
    t = z__[kps];
    z__[kps] = z__[kp];
    z__[kp] = t;
L180:
    i__1 = k - ks;
    saxpy_(&i__1, &z__[k], &a[k * a_dim1 + 1], &c__1, &z__[1], &c__1);
    if (ks == 2) {
	i__1 = k - ks;
	saxpy_(&i__1, &z__[k - 1], &a[(k - 1) * a_dim1 + 1], &c__1, &z__[1], &
		c__1);
    }
L190:
    if (ks == 2) {
	goto L210;
    }
    if ((r__1 = z__[k], dabs(r__1)) <= (r__2 = a[k + k * a_dim1], dabs(r__2)))
	     {
	goto L200;
    }
    s = (r__1 = a[k + k * a_dim1], dabs(r__1)) / (r__2 = z__[k], dabs(r__2));
    sscal_(n, &s, &z__[1], &c__1);
    ynorm = s * ynorm;
L200:
    if (a[k + k * a_dim1] != 0.f) {
	z__[k] /= a[k + k * a_dim1];
    }
    if (a[k + k * a_dim1] == 0.f) {
	z__[k] = 1.f;
    }
    goto L220;
L210:
    ak = a[k + k * a_dim1] / a[k - 1 + k * a_dim1];
    akm1 = a[k - 1 + (k - 1) * a_dim1] / a[k - 1 + k * a_dim1];
    bk = z__[k] / a[k - 1 + k * a_dim1];
    bkm1 = z__[k - 1] / a[k - 1 + k * a_dim1];
    denom = ak * akm1 - 1.f;
    z__[k] = (akm1 * bk - bkm1) / denom;
    z__[k - 1] = (ak * bkm1 - bk) / denom;
L220:
    k -= ks;
    goto L170;
L230:
    s = 1.f / sasum_(n, &z__[1], &c__1);
    sscal_(n, &s, &z__[1], &c__1);
    ynorm = s * ynorm;

/*     SOLVE TRANS(U)*Z = V */

    k = 1;
L240:
    if (k > *n) {
	goto L270;
    }
    ks = 1;
    if (kpvt[k] < 0) {
	ks = 2;
    }
    if (k == 1) {
	goto L260;
    }
    i__1 = k - 1;
    z__[k] += sdot_(&i__1, &a[k * a_dim1 + 1], &c__1, &z__[1], &c__1);
    if (ks == 2) {
	i__1 = k - 1;
	z__[k + 1] += sdot_(&i__1, &a[(k + 1) * a_dim1 + 1], &c__1, &z__[1], &
		c__1);
    }
    kp = (i__1 = kpvt[k], abs(i__1));
    if (kp == k) {
	goto L250;
    }
    t = z__[k];
    z__[k] = z__[kp];
    z__[kp] = t;
L250:
L260:
    k += ks;
    goto L240;
L270:
/*     MAKE ZNORM = 1.0 */
    s = 1.f / sasum_(n, &z__[1], &c__1);
    sscal_(n, &s, &z__[1], &c__1);
    ynorm = s * ynorm;

    if (anorm != 0.f) {
	*rcond = ynorm / anorm;
    }
    if (anorm == 0.f) {
	*rcond = 0.f;
    }
    return 0;
} /* ssico_ */

