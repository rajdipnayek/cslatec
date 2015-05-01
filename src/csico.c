/* csico.f -- translated by f2c (version 12.02.01).
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

/* DECK CSICO */
/* Subroutine */ int csico_(complex *a, integer *lda, integer *n, integer *
	kpvt, real *rcond, complex *z__)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3, i__4, i__5;
    real r__1, r__2, r__3, r__4, r__5, r__6, r__7, r__8;
    complex q__1, q__2, q__3;

    /* Local variables */
    static integer i__, j, k;
    static real s;
    static complex t, ak, bk, ek;
    static integer kp, ks, jm1, kps;
    static complex akm1, bkm1;
    static integer info;
    extern /* Subroutine */ int csifa_(complex *, integer *, integer *, 
	    integer *, integer *);
    static complex denom;
    static real anorm;
    extern /* Complex */ void cdotu_(complex *, integer *, complex *, integer 
	    *, complex *, integer *);
    extern /* Subroutine */ int caxpy_(integer *, complex *, complex *, 
	    integer *, complex *, integer *);
    static real ynorm;
    extern /* Subroutine */ int csscal_(integer *, real *, complex *, integer 
	    *);
    extern doublereal scasum_(integer *, complex *, integer *);

/* ***BEGIN PROLOGUE  CSICO */
/* ***PURPOSE  Factor a complex symmetric matrix by elimination with */
/*            symmetric pivoting and estimate the condition number of the */
/*            matrix. */
/* ***LIBRARY   SLATEC (LINPACK) */
/* ***CATEGORY  D2C1 */
/* ***TYPE      COMPLEX (SSICO-S, DSICO-D, CHICO-C, CSICO-C) */
/* ***KEYWORDS  CONDITION NUMBER, LINEAR ALGEBRA, LINPACK, */
/*             MATRIX FACTORIZATION, SYMMETRIC */
/* ***AUTHOR  Moler, C. B., (U. of New Mexico) */
/* ***DESCRIPTION */

/*     CSICO factors a complex symmetric matrix by elimination with */
/*     symmetric pivoting and estimates the condition of the matrix. */

/*     If  RCOND  is not needed, CSIFA is slightly faster. */
/*     To solve  A*X = B , follow CSICO by CSISL. */
/*     To compute  INVERSE(A)*C , follow CSICO by CSISL. */
/*     To compute  INVERSE(A) , follow CSICO by CSIDI. */
/*     To compute  DETERMINANT(A) , follow CSICO by CSIDI. */

/*     On Entry */

/*        A       COMPLEX(LDA, N) */
/*                the symmetric matrix to be factored. */
/*                Only the diagonal and upper triangle are used. */

/*        LDA     INTEGER */
/*                the leading dimension of the array  A . */

/*        N       INTEGER */
/*                the order of the matrix  A . */

/*     On Return */

/*        A       a block diagonal matrix and the multipliers which */
/*                were used to obtain it. */
/*                The factorization can be written  A = U*D*TRANS(U) */
/*                where  U  is a product of permutation and unit */
/*                upper triangular matrices , TRANS(U) is the */
/*                transpose of  U , and  D  is block diagonal */
/*                with 1 by 1 and 2 by 2 blocks. */

/*        KVPT    INTEGER(N) */
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

/*        Z       COMPLEX(N) */
/*                a work vector whose contents are usually unimportant. */
/*                If  A  is close to a singular matrix, then  Z  is */
/*                an approximate null vector in the sense that */
/*                NORM(A*Z) = RCOND*NORM(A)*NORM(Z) . */

/* ***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W. */
/*                 Stewart, LINPACK Users' Guide, SIAM, 1979. */
/* ***ROUTINES CALLED  CAXPY, CDOTU, CSIFA, CSSCAL, SCASUM */
/* ***REVISION HISTORY  (YYMMDD) */
/*   780814  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890831  Modified array declarations.  (WRB) */
/*   891107  Corrected category and modified routine equivalence */
/*           list.  (WRB) */
/*   891107  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900326  Removed duplicate information from DESCRIPTION section. */
/*           (WRB) */
/*   920501  Reformatted the REFERENCES section.  (WRB) */
/* ***END PROLOGUE  CSICO */


/*     FIND NORM OF A USING ONLY UPPER HALF */

/* ***FIRST EXECUTABLE STATEMENT  CSICO */
    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --kpvt;
    --z__;

    /* Function Body */
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	i__2 = j;
	r__1 = scasum_(&j, &a[j * a_dim1 + 1], &c__1);
	q__1.r = r__1, q__1.i = 0.f;
	z__[i__2].r = q__1.r, z__[i__2].i = q__1.i;
	jm1 = j - 1;
	if (jm1 < 1) {
	    goto L20;
	}
	i__2 = jm1;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    i__3 = i__;
	    i__4 = i__;
	    i__5 = i__ + j * a_dim1;
	    r__3 = z__[i__4].r + ((r__1 = a[i__5].r, dabs(r__1)) + (r__2 = 
		    r_imag(&a[i__ + j * a_dim1]), dabs(r__2)));
	    q__1.r = r__3, q__1.i = 0.f;
	    z__[i__3].r = q__1.r, z__[i__3].i = q__1.i;
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
	i__2 = j;
	r__1 = anorm, r__2 = z__[i__2].r;
	anorm = dmax(r__1,r__2);
/* L40: */
    }

/*     FACTOR */

    csifa_(&a[a_offset], lda, n, &kpvt[1], &info);

/*     RCOND = 1/(NORM(A)*(ESTIMATE OF NORM(INVERSE(A)))) . */
/*     ESTIMATE = NORM(Z)/NORM(Y) WHERE  A*Z = Y  AND  A*Y = E . */
/*     THE COMPONENTS OF  E  ARE CHOSEN TO CAUSE MAXIMUM LOCAL */
/*     GROWTH IN THE ELEMENTS OF W  WHERE  U*D*W = E . */
/*     THE VECTORS ARE FREQUENTLY RESCALED TO AVOID OVERFLOW. */

/*     SOLVE U*D*W = E */

    ek.r = 1.f, ek.i = 0.f;
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	i__2 = j;
	z__[i__2].r = 0.f, z__[i__2].i = 0.f;
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
    i__1 = kps;
    t.r = z__[i__1].r, t.i = z__[i__1].i;
    i__1 = kps;
    i__2 = kp;
    z__[i__1].r = z__[i__2].r, z__[i__1].i = z__[i__2].i;
    i__1 = kp;
    z__[i__1].r = t.r, z__[i__1].i = t.i;
L70:
    i__1 = k;
    if ((r__1 = z__[i__1].r, dabs(r__1)) + (r__2 = r_imag(&z__[k]), dabs(r__2)
	    ) != 0.f) {
	r__7 = (r__3 = ek.r, dabs(r__3)) + (r__4 = r_imag(&ek), dabs(r__4));
	i__2 = k;
	i__3 = k;
	r__8 = (r__5 = z__[i__3].r, dabs(r__5)) + (r__6 = r_imag(&z__[k]), 
		dabs(r__6));
	q__2.r = z__[i__2].r / r__8, q__2.i = z__[i__2].i / r__8;
	q__1.r = r__7 * q__2.r, q__1.i = r__7 * q__2.i;
	ek.r = q__1.r, ek.i = q__1.i;
    }
    i__1 = k;
    i__2 = k;
    q__1.r = z__[i__2].r + ek.r, q__1.i = z__[i__2].i + ek.i;
    z__[i__1].r = q__1.r, z__[i__1].i = q__1.i;
    i__1 = k - ks;
    caxpy_(&i__1, &z__[k], &a[k * a_dim1 + 1], &c__1, &z__[1], &c__1);
    if (ks == 1) {
	goto L80;
    }
    i__1 = k - 1;
    if ((r__1 = z__[i__1].r, dabs(r__1)) + (r__2 = r_imag(&z__[k - 1]), dabs(
	    r__2)) != 0.f) {
	r__7 = (r__3 = ek.r, dabs(r__3)) + (r__4 = r_imag(&ek), dabs(r__4));
	i__2 = k - 1;
	i__3 = k - 1;
	r__8 = (r__5 = z__[i__3].r, dabs(r__5)) + (r__6 = r_imag(&z__[k - 1]),
		 dabs(r__6));
	q__2.r = z__[i__2].r / r__8, q__2.i = z__[i__2].i / r__8;
	q__1.r = r__7 * q__2.r, q__1.i = r__7 * q__2.i;
	ek.r = q__1.r, ek.i = q__1.i;
    }
    i__1 = k - 1;
    i__2 = k - 1;
    q__1.r = z__[i__2].r + ek.r, q__1.i = z__[i__2].i + ek.i;
    z__[i__1].r = q__1.r, z__[i__1].i = q__1.i;
    i__1 = k - ks;
    caxpy_(&i__1, &z__[k - 1], &a[(k - 1) * a_dim1 + 1], &c__1, &z__[1], &
	    c__1);
L80:
    if (ks == 2) {
	goto L100;
    }
    i__1 = k;
    i__2 = k + k * a_dim1;
    if ((r__1 = z__[i__1].r, dabs(r__1)) + (r__2 = r_imag(&z__[k]), dabs(r__2)
	    ) <= (r__3 = a[i__2].r, dabs(r__3)) + (r__4 = r_imag(&a[k + k * 
	    a_dim1]), dabs(r__4))) {
	goto L90;
    }
    i__1 = k + k * a_dim1;
    i__2 = k;
    s = ((r__1 = a[i__1].r, dabs(r__1)) + (r__2 = r_imag(&a[k + k * a_dim1]), 
	    dabs(r__2))) / ((r__3 = z__[i__2].r, dabs(r__3)) + (r__4 = r_imag(
	    &z__[k]), dabs(r__4)));
    csscal_(n, &s, &z__[1], &c__1);
    q__2.r = s, q__2.i = 0.f;
    q__1.r = q__2.r * ek.r - q__2.i * ek.i, q__1.i = q__2.r * ek.i + q__2.i * 
	    ek.r;
    ek.r = q__1.r, ek.i = q__1.i;
L90:
    i__1 = k + k * a_dim1;
    if ((r__1 = a[i__1].r, dabs(r__1)) + (r__2 = r_imag(&a[k + k * a_dim1]), 
	    dabs(r__2)) != 0.f) {
	i__2 = k;
	c_div(&q__1, &z__[k], &a[k + k * a_dim1]);
	z__[i__2].r = q__1.r, z__[i__2].i = q__1.i;
    }
    i__1 = k + k * a_dim1;
    if ((r__1 = a[i__1].r, dabs(r__1)) + (r__2 = r_imag(&a[k + k * a_dim1]), 
	    dabs(r__2)) == 0.f) {
	i__2 = k;
	z__[i__2].r = 1.f, z__[i__2].i = 0.f;
    }
    goto L110;
L100:
    c_div(&q__1, &a[k + k * a_dim1], &a[k - 1 + k * a_dim1]);
    ak.r = q__1.r, ak.i = q__1.i;
    c_div(&q__1, &a[k - 1 + (k - 1) * a_dim1], &a[k - 1 + k * a_dim1]);
    akm1.r = q__1.r, akm1.i = q__1.i;
    c_div(&q__1, &z__[k], &a[k - 1 + k * a_dim1]);
    bk.r = q__1.r, bk.i = q__1.i;
    c_div(&q__1, &z__[k - 1], &a[k - 1 + k * a_dim1]);
    bkm1.r = q__1.r, bkm1.i = q__1.i;
    q__2.r = ak.r * akm1.r - ak.i * akm1.i, q__2.i = ak.r * akm1.i + ak.i * 
	    akm1.r;
    q__1.r = q__2.r - 1.f, q__1.i = q__2.i;
    denom.r = q__1.r, denom.i = q__1.i;
    i__1 = k;
    q__3.r = akm1.r * bk.r - akm1.i * bk.i, q__3.i = akm1.r * bk.i + akm1.i * 
	    bk.r;
    q__2.r = q__3.r - bkm1.r, q__2.i = q__3.i - bkm1.i;
    c_div(&q__1, &q__2, &denom);
    z__[i__1].r = q__1.r, z__[i__1].i = q__1.i;
    i__1 = k - 1;
    q__3.r = ak.r * bkm1.r - ak.i * bkm1.i, q__3.i = ak.r * bkm1.i + ak.i * 
	    bkm1.r;
    q__2.r = q__3.r - bk.r, q__2.i = q__3.i - bk.i;
    c_div(&q__1, &q__2, &denom);
    z__[i__1].r = q__1.r, z__[i__1].i = q__1.i;
L110:
    k -= ks;
    goto L60;
L120:
    s = 1.f / scasum_(n, &z__[1], &c__1);
    csscal_(n, &s, &z__[1], &c__1);

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
    i__1 = k;
    i__2 = k;
    i__3 = k - 1;
    cdotu_(&q__2, &i__3, &a[k * a_dim1 + 1], &c__1, &z__[1], &c__1);
    q__1.r = z__[i__2].r + q__2.r, q__1.i = z__[i__2].i + q__2.i;
    z__[i__1].r = q__1.r, z__[i__1].i = q__1.i;
    if (ks == 2) {
	i__1 = k + 1;
	i__2 = k + 1;
	i__3 = k - 1;
	cdotu_(&q__2, &i__3, &a[(k + 1) * a_dim1 + 1], &c__1, &z__[1], &c__1);
	q__1.r = z__[i__2].r + q__2.r, q__1.i = z__[i__2].i + q__2.i;
	z__[i__1].r = q__1.r, z__[i__1].i = q__1.i;
    }
    kp = (i__1 = kpvt[k], abs(i__1));
    if (kp == k) {
	goto L140;
    }
    i__1 = k;
    t.r = z__[i__1].r, t.i = z__[i__1].i;
    i__1 = k;
    i__2 = kp;
    z__[i__1].r = z__[i__2].r, z__[i__1].i = z__[i__2].i;
    i__1 = kp;
    z__[i__1].r = t.r, z__[i__1].i = t.i;
L140:
L150:
    k += ks;
    goto L130;
L160:
    s = 1.f / scasum_(n, &z__[1], &c__1);
    csscal_(n, &s, &z__[1], &c__1);

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
    i__1 = kps;
    t.r = z__[i__1].r, t.i = z__[i__1].i;
    i__1 = kps;
    i__2 = kp;
    z__[i__1].r = z__[i__2].r, z__[i__1].i = z__[i__2].i;
    i__1 = kp;
    z__[i__1].r = t.r, z__[i__1].i = t.i;
L180:
    i__1 = k - ks;
    caxpy_(&i__1, &z__[k], &a[k * a_dim1 + 1], &c__1, &z__[1], &c__1);
    if (ks == 2) {
	i__1 = k - ks;
	caxpy_(&i__1, &z__[k - 1], &a[(k - 1) * a_dim1 + 1], &c__1, &z__[1], &
		c__1);
    }
L190:
    if (ks == 2) {
	goto L210;
    }
    i__1 = k;
    i__2 = k + k * a_dim1;
    if ((r__1 = z__[i__1].r, dabs(r__1)) + (r__2 = r_imag(&z__[k]), dabs(r__2)
	    ) <= (r__3 = a[i__2].r, dabs(r__3)) + (r__4 = r_imag(&a[k + k * 
	    a_dim1]), dabs(r__4))) {
	goto L200;
    }
    i__1 = k + k * a_dim1;
    i__2 = k;
    s = ((r__1 = a[i__1].r, dabs(r__1)) + (r__2 = r_imag(&a[k + k * a_dim1]), 
	    dabs(r__2))) / ((r__3 = z__[i__2].r, dabs(r__3)) + (r__4 = r_imag(
	    &z__[k]), dabs(r__4)));
    csscal_(n, &s, &z__[1], &c__1);
    ynorm = s * ynorm;
L200:
    i__1 = k + k * a_dim1;
    if ((r__1 = a[i__1].r, dabs(r__1)) + (r__2 = r_imag(&a[k + k * a_dim1]), 
	    dabs(r__2)) != 0.f) {
	i__2 = k;
	c_div(&q__1, &z__[k], &a[k + k * a_dim1]);
	z__[i__2].r = q__1.r, z__[i__2].i = q__1.i;
    }
    i__1 = k + k * a_dim1;
    if ((r__1 = a[i__1].r, dabs(r__1)) + (r__2 = r_imag(&a[k + k * a_dim1]), 
	    dabs(r__2)) == 0.f) {
	i__2 = k;
	z__[i__2].r = 1.f, z__[i__2].i = 0.f;
    }
    goto L220;
L210:
    c_div(&q__1, &a[k + k * a_dim1], &a[k - 1 + k * a_dim1]);
    ak.r = q__1.r, ak.i = q__1.i;
    c_div(&q__1, &a[k - 1 + (k - 1) * a_dim1], &a[k - 1 + k * a_dim1]);
    akm1.r = q__1.r, akm1.i = q__1.i;
    c_div(&q__1, &z__[k], &a[k - 1 + k * a_dim1]);
    bk.r = q__1.r, bk.i = q__1.i;
    c_div(&q__1, &z__[k - 1], &a[k - 1 + k * a_dim1]);
    bkm1.r = q__1.r, bkm1.i = q__1.i;
    q__2.r = ak.r * akm1.r - ak.i * akm1.i, q__2.i = ak.r * akm1.i + ak.i * 
	    akm1.r;
    q__1.r = q__2.r - 1.f, q__1.i = q__2.i;
    denom.r = q__1.r, denom.i = q__1.i;
    i__1 = k;
    q__3.r = akm1.r * bk.r - akm1.i * bk.i, q__3.i = akm1.r * bk.i + akm1.i * 
	    bk.r;
    q__2.r = q__3.r - bkm1.r, q__2.i = q__3.i - bkm1.i;
    c_div(&q__1, &q__2, &denom);
    z__[i__1].r = q__1.r, z__[i__1].i = q__1.i;
    i__1 = k - 1;
    q__3.r = ak.r * bkm1.r - ak.i * bkm1.i, q__3.i = ak.r * bkm1.i + ak.i * 
	    bkm1.r;
    q__2.r = q__3.r - bk.r, q__2.i = q__3.i - bk.i;
    c_div(&q__1, &q__2, &denom);
    z__[i__1].r = q__1.r, z__[i__1].i = q__1.i;
L220:
    k -= ks;
    goto L170;
L230:
    s = 1.f / scasum_(n, &z__[1], &c__1);
    csscal_(n, &s, &z__[1], &c__1);
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
    i__1 = k;
    i__2 = k;
    i__3 = k - 1;
    cdotu_(&q__2, &i__3, &a[k * a_dim1 + 1], &c__1, &z__[1], &c__1);
    q__1.r = z__[i__2].r + q__2.r, q__1.i = z__[i__2].i + q__2.i;
    z__[i__1].r = q__1.r, z__[i__1].i = q__1.i;
    if (ks == 2) {
	i__1 = k + 1;
	i__2 = k + 1;
	i__3 = k - 1;
	cdotu_(&q__2, &i__3, &a[(k + 1) * a_dim1 + 1], &c__1, &z__[1], &c__1);
	q__1.r = z__[i__2].r + q__2.r, q__1.i = z__[i__2].i + q__2.i;
	z__[i__1].r = q__1.r, z__[i__1].i = q__1.i;
    }
    kp = (i__1 = kpvt[k], abs(i__1));
    if (kp == k) {
	goto L250;
    }
    i__1 = k;
    t.r = z__[i__1].r, t.i = z__[i__1].i;
    i__1 = k;
    i__2 = kp;
    z__[i__1].r = z__[i__2].r, z__[i__1].i = z__[i__2].i;
    i__1 = kp;
    z__[i__1].r = t.r, z__[i__1].i = t.i;
L250:
L260:
    k += ks;
    goto L240;
L270:
/*     MAKE ZNORM = 1.0 */
    s = 1.f / scasum_(n, &z__[1], &c__1);
    csscal_(n, &s, &z__[1], &c__1);
    ynorm = s * ynorm;

    if (anorm != 0.f) {
	*rcond = ynorm / anorm;
    }
    if (anorm == 0.f) {
	*rcond = 0.f;
    }
    return 0;
} /* csico_ */

