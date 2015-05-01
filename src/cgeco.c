/* cgeco.f -- translated by f2c (version 12.02.01).
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

/* DECK CGECO */
/* Subroutine */ int cgeco_(complex *a, integer *lda, integer *n, integer *
	ipvt, real *rcond, complex *z__)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3, i__4;
    real r__1, r__2, r__3, r__4, r__5, r__6, r__7, r__8;
    complex q__1, q__2, q__3, q__4;

    /* Local variables */
    static integer j, k, l;
    static real s;
    static complex t;
    static integer kb;
    static complex ek;
    static real sm;
    static complex wk;
    static integer kp1;
    static complex wkm;
    static integer info;
    extern /* Subroutine */ int cgefa_(complex *, integer *, integer *, 
	    integer *, integer *);
    extern /* Complex */ void cdotc_(complex *, integer *, complex *, integer 
	    *, complex *, integer *);
    static real anorm;
    extern /* Subroutine */ int caxpy_(integer *, complex *, complex *, 
	    integer *, complex *, integer *);
    static real ynorm;
    extern /* Subroutine */ int csscal_(integer *, real *, complex *, integer 
	    *);
    extern doublereal scasum_(integer *, complex *, integer *);

/* ***BEGIN PROLOGUE  CGECO */
/* ***PURPOSE  Factor a matrix using Gaussian elimination and estimate */
/*            the condition number of the matrix. */
/* ***LIBRARY   SLATEC (LINPACK) */
/* ***CATEGORY  D2C1 */
/* ***TYPE      COMPLEX (SGECO-S, DGECO-D, CGECO-C) */
/* ***KEYWORDS  CONDITION NUMBER, GENERAL MATRIX, LINEAR ALGEBRA, LINPACK, */
/*             MATRIX FACTORIZATION */
/* ***AUTHOR  Moler, C. B., (U. of New Mexico) */
/* ***DESCRIPTION */

/*     CGECO factors a complex matrix by Gaussian elimination */
/*     and estimates the condition of the matrix. */

/*     If  RCOND  is not needed, CGEFA is slightly faster. */
/*     To solve  A*X = B , follow CGECO By CGESL. */
/*     To Compute  INVERSE(A)*C , follow CGECO by CGESL. */
/*     To compute  DETERMINANT(A) , follow CGECO by CGEDI. */
/*     To compute  INVERSE(A) , follow CGECO by CGEDI. */

/*     On Entry */

/*        A       COMPLEX(LDA, N) */
/*                the matrix to be factored. */

/*        LDA     INTEGER */
/*                the leading dimension of the array  A . */

/*        N       INTEGER */
/*                the order of the matrix  A . */

/*     On Return */

/*        A       an upper triangular matrix and the multipliers */
/*                which were used to obtain it. */
/*                The factorization can be written  A = L*U  where */
/*                L  is a product of permutation and unit lower */
/*                triangular matrices and  U  is upper triangular. */

/*        IPVT    INTEGER(N) */
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
/* ***ROUTINES CALLED  CAXPY, CDOTC, CGEFA, CSSCAL, SCASUM */
/* ***REVISION HISTORY  (YYMMDD) */
/*   780814  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890831  Modified array declarations.  (WRB) */
/*   890831  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900326  Removed duplicate information from DESCRIPTION section. */
/*           (WRB) */
/*   920501  Reformatted the REFERENCES section.  (WRB) */
/* ***END PROLOGUE  CGECO */


/*     COMPUTE 1-NORM OF A */

/* ***FIRST EXECUTABLE STATEMENT  CGECO */
    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --ipvt;
    --z__;

    /* Function Body */
    anorm = 0.f;
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
/* Computing MAX */
	r__1 = anorm, r__2 = scasum_(n, &a[j * a_dim1 + 1], &c__1);
	anorm = dmax(r__1,r__2);
/* L10: */
    }

/*     FACTOR */

    cgefa_(&a[a_offset], lda, n, &ipvt[1], &info);

/*     RCOND = 1/(NORM(A)*(ESTIMATE OF NORM(INVERSE(A)))) . */
/*     ESTIMATE = NORM(Z)/NORM(Y) WHERE  A*Z = Y  AND  CTRANS(A)*Y = E . */
/*     CTRANS(A)  IS THE CONJUGATE TRANSPOSE OF A . */
/*     THE COMPONENTS OF  E  ARE CHOSEN TO CAUSE MAXIMUM LOCAL */
/*     GROWTH IN THE ELEMENTS OF W  WHERE  CTRANS(U)*W = E . */
/*     THE VECTORS ARE FREQUENTLY RESCALED TO AVOID OVERFLOW. */

/*     SOLVE CTRANS(U)*W = E */

    ek.r = 1.f, ek.i = 0.f;
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	i__2 = j;
	z__[i__2].r = 0.f, z__[i__2].i = 0.f;
/* L20: */
    }
    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
	i__2 = k;
	if ((r__1 = z__[i__2].r, dabs(r__1)) + (r__2 = r_imag(&z__[k]), dabs(
		r__2)) != 0.f) {
	    i__3 = k;
	    q__2.r = -z__[i__3].r, q__2.i = -z__[i__3].i;
	    q__1.r = q__2.r, q__1.i = q__2.i;
	    r__7 = (r__3 = ek.r, dabs(r__3)) + (r__4 = r_imag(&ek), dabs(r__4)
		    );
	    r__8 = (r__5 = q__1.r, dabs(r__5)) + (r__6 = r_imag(&q__1), dabs(
		    r__6));
	    q__4.r = q__1.r / r__8, q__4.i = q__1.i / r__8;
	    q__3.r = r__7 * q__4.r, q__3.i = r__7 * q__4.i;
	    ek.r = q__3.r, ek.i = q__3.i;
	}
	i__2 = k;
	q__2.r = ek.r - z__[i__2].r, q__2.i = ek.i - z__[i__2].i;
	q__1.r = q__2.r, q__1.i = q__2.i;
	i__3 = k + k * a_dim1;
	if ((r__1 = q__1.r, dabs(r__1)) + (r__2 = r_imag(&q__1), dabs(r__2)) 
		<= (r__3 = a[i__3].r, dabs(r__3)) + (r__4 = r_imag(&a[k + k * 
		a_dim1]), dabs(r__4))) {
	    goto L30;
	}
	i__2 = k;
	q__2.r = ek.r - z__[i__2].r, q__2.i = ek.i - z__[i__2].i;
	q__1.r = q__2.r, q__1.i = q__2.i;
	i__3 = k + k * a_dim1;
	s = ((r__1 = a[i__3].r, dabs(r__1)) + (r__2 = r_imag(&a[k + k * 
		a_dim1]), dabs(r__2))) / ((r__3 = q__1.r, dabs(r__3)) + (r__4 
		= r_imag(&q__1), dabs(r__4)));
	csscal_(n, &s, &z__[1], &c__1);
	q__2.r = s, q__2.i = 0.f;
	q__1.r = q__2.r * ek.r - q__2.i * ek.i, q__1.i = q__2.r * ek.i + 
		q__2.i * ek.r;
	ek.r = q__1.r, ek.i = q__1.i;
L30:
	i__2 = k;
	q__1.r = ek.r - z__[i__2].r, q__1.i = ek.i - z__[i__2].i;
	wk.r = q__1.r, wk.i = q__1.i;
	q__2.r = -ek.r, q__2.i = -ek.i;
	i__2 = k;
	q__1.r = q__2.r - z__[i__2].r, q__1.i = q__2.i - z__[i__2].i;
	wkm.r = q__1.r, wkm.i = q__1.i;
	s = (r__1 = wk.r, dabs(r__1)) + (r__2 = r_imag(&wk), dabs(r__2));
	sm = (r__1 = wkm.r, dabs(r__1)) + (r__2 = r_imag(&wkm), dabs(r__2));
	i__2 = k + k * a_dim1;
	if ((r__1 = a[i__2].r, dabs(r__1)) + (r__2 = r_imag(&a[k + k * a_dim1]
		), dabs(r__2)) == 0.f) {
	    goto L40;
	}
	r_cnjg(&q__2, &a[k + k * a_dim1]);
	c_div(&q__1, &wk, &q__2);
	wk.r = q__1.r, wk.i = q__1.i;
	r_cnjg(&q__2, &a[k + k * a_dim1]);
	c_div(&q__1, &wkm, &q__2);
	wkm.r = q__1.r, wkm.i = q__1.i;
	goto L50;
L40:
	wk.r = 1.f, wk.i = 0.f;
	wkm.r = 1.f, wkm.i = 0.f;
L50:
	kp1 = k + 1;
	if (kp1 > *n) {
	    goto L90;
	}
	i__2 = *n;
	for (j = kp1; j <= i__2; ++j) {
	    i__3 = j;
	    r_cnjg(&q__4, &a[k + j * a_dim1]);
	    q__3.r = wkm.r * q__4.r - wkm.i * q__4.i, q__3.i = wkm.r * q__4.i 
		    + wkm.i * q__4.r;
	    q__2.r = z__[i__3].r + q__3.r, q__2.i = z__[i__3].i + q__3.i;
	    q__1.r = q__2.r, q__1.i = q__2.i;
	    sm += (r__1 = q__1.r, dabs(r__1)) + (r__2 = r_imag(&q__1), dabs(
		    r__2));
	    i__3 = j;
	    i__4 = j;
	    r_cnjg(&q__3, &a[k + j * a_dim1]);
	    q__2.r = wk.r * q__3.r - wk.i * q__3.i, q__2.i = wk.r * q__3.i + 
		    wk.i * q__3.r;
	    q__1.r = z__[i__4].r + q__2.r, q__1.i = z__[i__4].i + q__2.i;
	    z__[i__3].r = q__1.r, z__[i__3].i = q__1.i;
	    i__3 = j;
	    s += (r__1 = z__[i__3].r, dabs(r__1)) + (r__2 = r_imag(&z__[j]), 
		    dabs(r__2));
/* L60: */
	}
	if (s >= sm) {
	    goto L80;
	}
	q__1.r = wkm.r - wk.r, q__1.i = wkm.i - wk.i;
	t.r = q__1.r, t.i = q__1.i;
	wk.r = wkm.r, wk.i = wkm.i;
	i__2 = *n;
	for (j = kp1; j <= i__2; ++j) {
	    i__3 = j;
	    i__4 = j;
	    r_cnjg(&q__3, &a[k + j * a_dim1]);
	    q__2.r = t.r * q__3.r - t.i * q__3.i, q__2.i = t.r * q__3.i + t.i 
		    * q__3.r;
	    q__1.r = z__[i__4].r + q__2.r, q__1.i = z__[i__4].i + q__2.i;
	    z__[i__3].r = q__1.r, z__[i__3].i = q__1.i;
/* L70: */
	}
L80:
L90:
	i__2 = k;
	z__[i__2].r = wk.r, z__[i__2].i = wk.i;
/* L100: */
    }
    s = 1.f / scasum_(n, &z__[1], &c__1);
    csscal_(n, &s, &z__[1], &c__1);

/*     SOLVE CTRANS(L)*Y = W */

    i__1 = *n;
    for (kb = 1; kb <= i__1; ++kb) {
	k = *n + 1 - kb;
	if (k < *n) {
	    i__2 = k;
	    i__3 = k;
	    i__4 = *n - k;
	    cdotc_(&q__2, &i__4, &a[k + 1 + k * a_dim1], &c__1, &z__[k + 1], &
		    c__1);
	    q__1.r = z__[i__3].r + q__2.r, q__1.i = z__[i__3].i + q__2.i;
	    z__[i__2].r = q__1.r, z__[i__2].i = q__1.i;
	}
	i__2 = k;
	if ((r__1 = z__[i__2].r, dabs(r__1)) + (r__2 = r_imag(&z__[k]), dabs(
		r__2)) <= 1.f) {
	    goto L110;
	}
	i__2 = k;
	s = 1.f / ((r__1 = z__[i__2].r, dabs(r__1)) + (r__2 = r_imag(&z__[k]),
		 dabs(r__2)));
	csscal_(n, &s, &z__[1], &c__1);
L110:
	l = ipvt[k];
	i__2 = l;
	t.r = z__[i__2].r, t.i = z__[i__2].i;
	i__2 = l;
	i__3 = k;
	z__[i__2].r = z__[i__3].r, z__[i__2].i = z__[i__3].i;
	i__2 = k;
	z__[i__2].r = t.r, z__[i__2].i = t.i;
/* L120: */
    }
    s = 1.f / scasum_(n, &z__[1], &c__1);
    csscal_(n, &s, &z__[1], &c__1);

    ynorm = 1.f;

/*     SOLVE L*V = Y */

    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
	l = ipvt[k];
	i__2 = l;
	t.r = z__[i__2].r, t.i = z__[i__2].i;
	i__2 = l;
	i__3 = k;
	z__[i__2].r = z__[i__3].r, z__[i__2].i = z__[i__3].i;
	i__2 = k;
	z__[i__2].r = t.r, z__[i__2].i = t.i;
	if (k < *n) {
	    i__2 = *n - k;
	    caxpy_(&i__2, &t, &a[k + 1 + k * a_dim1], &c__1, &z__[k + 1], &
		    c__1);
	}
	i__2 = k;
	if ((r__1 = z__[i__2].r, dabs(r__1)) + (r__2 = r_imag(&z__[k]), dabs(
		r__2)) <= 1.f) {
	    goto L130;
	}
	i__2 = k;
	s = 1.f / ((r__1 = z__[i__2].r, dabs(r__1)) + (r__2 = r_imag(&z__[k]),
		 dabs(r__2)));
	csscal_(n, &s, &z__[1], &c__1);
	ynorm = s * ynorm;
L130:
/* L140: */
	;
    }
    s = 1.f / scasum_(n, &z__[1], &c__1);
    csscal_(n, &s, &z__[1], &c__1);
    ynorm = s * ynorm;

/*     SOLVE  U*Z = V */

    i__1 = *n;
    for (kb = 1; kb <= i__1; ++kb) {
	k = *n + 1 - kb;
	i__2 = k;
	i__3 = k + k * a_dim1;
	if ((r__1 = z__[i__2].r, dabs(r__1)) + (r__2 = r_imag(&z__[k]), dabs(
		r__2)) <= (r__3 = a[i__3].r, dabs(r__3)) + (r__4 = r_imag(&a[
		k + k * a_dim1]), dabs(r__4))) {
	    goto L150;
	}
	i__2 = k + k * a_dim1;
	i__3 = k;
	s = ((r__1 = a[i__2].r, dabs(r__1)) + (r__2 = r_imag(&a[k + k * 
		a_dim1]), dabs(r__2))) / ((r__3 = z__[i__3].r, dabs(r__3)) + (
		r__4 = r_imag(&z__[k]), dabs(r__4)));
	csscal_(n, &s, &z__[1], &c__1);
	ynorm = s * ynorm;
L150:
	i__2 = k + k * a_dim1;
	if ((r__1 = a[i__2].r, dabs(r__1)) + (r__2 = r_imag(&a[k + k * a_dim1]
		), dabs(r__2)) != 0.f) {
	    i__3 = k;
	    c_div(&q__1, &z__[k], &a[k + k * a_dim1]);
	    z__[i__3].r = q__1.r, z__[i__3].i = q__1.i;
	}
	i__2 = k + k * a_dim1;
	if ((r__1 = a[i__2].r, dabs(r__1)) + (r__2 = r_imag(&a[k + k * a_dim1]
		), dabs(r__2)) == 0.f) {
	    i__3 = k;
	    z__[i__3].r = 1.f, z__[i__3].i = 0.f;
	}
	i__2 = k;
	q__1.r = -z__[i__2].r, q__1.i = -z__[i__2].i;
	t.r = q__1.r, t.i = q__1.i;
	i__2 = k - 1;
	caxpy_(&i__2, &t, &a[k * a_dim1 + 1], &c__1, &z__[1], &c__1);
/* L160: */
    }
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
} /* cgeco_ */

