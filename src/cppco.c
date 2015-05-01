/* cppco.f -- translated by f2c (version 12.02.01).
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

/* DECK CPPCO */
/* Subroutine */ int cppco_(complex *ap, integer *n, real *rcond, complex *
	z__, integer *info)
{
    /* System generated locals */
    integer i__1, i__2, i__3, i__4, i__5;
    real r__1, r__2, r__3, r__4, r__5, r__6, r__7, r__8;
    complex q__1, q__2, q__3, q__4;

    /* Local variables */
    static integer i__, j, k;
    static real s;
    static complex t;
    static integer j1, kb;
    static complex ek;
    static integer ij, kj, kk;
    static real sm;
    static complex wk;
    static integer jm1, kp1;
    static complex wkm;
    extern /* Subroutine */ int cppfa_(complex *, integer *, integer *);
    extern /* Complex */ void cdotc_(complex *, integer *, complex *, integer 
	    *, complex *, integer *);
    static real anorm;
    extern /* Subroutine */ int caxpy_(integer *, complex *, complex *, 
	    integer *, complex *, integer *);
    static real ynorm;
    extern /* Subroutine */ int csscal_(integer *, real *, complex *, integer 
	    *);
    extern doublereal scasum_(integer *, complex *, integer *);

/* ***BEGIN PROLOGUE  CPPCO */
/* ***PURPOSE  Factor a complex Hermitian positive definite matrix stored */
/*            in packed form and estimate the condition number of the */
/*            matrix. */
/* ***LIBRARY   SLATEC (LINPACK) */
/* ***CATEGORY  D2D1B */
/* ***TYPE      COMPLEX (SPPCO-S, DPPCO-D, CPPCO-C) */
/* ***KEYWORDS  CONDITION NUMBER, LINEAR ALGEBRA, LINPACK, */
/*             MATRIX FACTORIZATION, PACKED, POSITIVE DEFINITE */
/* ***AUTHOR  Moler, C. B., (U. of New Mexico) */
/* ***DESCRIPTION */

/*     CPPCO factors a complex Hermitian positive definite matrix */
/*     stored in packed form and estimates the condition of the matrix. */

/*     If  RCOND  is not needed, CPPFA is slightly faster. */
/*     To solve  A*X = B , follow CPPCO by CPPSL. */
/*     To compute  INVERSE(A)*C , follow CPPCO by CPPSL. */
/*     To compute  DETERMINANT(A) , follow CPPCO by CPPDI. */
/*     To compute  INVERSE(A) , follow CPPCO by CPPDI. */

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

/*        Z       COMPLEX(N) */
/*                a work vector whose contents are usually unimportant. */
/*                If  A  is singular to working precision, then  Z  is */
/*                an approximate null vector in the sense that */
/*                NORM(A*Z) = RCOND*NORM(A)*NORM(Z) . */
/*                If  INFO .NE. 0 , Z  is unchanged. */

/*        INFO    INTEGER */
/*                = 0  for normal return. */
/*                = K  signals an error condition.  The leading minor */
/*                     of order  K  is not positive definite. */

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
/* ***ROUTINES CALLED  CAXPY, CDOTC, CPPFA, CSSCAL, SCASUM */
/* ***REVISION HISTORY  (YYMMDD) */
/*   780814  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890831  Modified array declarations.  (WRB) */
/*   890831  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900326  Removed duplicate information from DESCRIPTION section. */
/*           (WRB) */
/*   920501  Reformatted the REFERENCES section.  (WRB) */
/* ***END PROLOGUE  CPPCO */


/*     FIND NORM OF A */

/* ***FIRST EXECUTABLE STATEMENT  CPPCO */
    /* Parameter adjustments */
    --z__;
    --ap;

    /* Function Body */
    j1 = 1;
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	i__2 = j;
	r__1 = scasum_(&j, &ap[j1], &c__1);
	q__1.r = r__1, q__1.i = 0.f;
	z__[i__2].r = q__1.r, z__[i__2].i = q__1.i;
	ij = j1;
	j1 += j;
	jm1 = j - 1;
	if (jm1 < 1) {
	    goto L20;
	}
	i__2 = jm1;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    i__3 = i__;
	    i__4 = i__;
	    i__5 = ij;
	    r__3 = z__[i__4].r + ((r__1 = ap[i__5].r, dabs(r__1)) + (r__2 = 
		    r_imag(&ap[ij]), dabs(r__2)));
	    q__1.r = r__3, q__1.i = 0.f;
	    z__[i__3].r = q__1.r, z__[i__3].i = q__1.i;
	    ++ij;
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

    cppfa_(&ap[1], n, info);
    if (*info != 0) {
	goto L180;
    }

/*        RCOND = 1/(NORM(A)*(ESTIMATE OF NORM(INVERSE(A)))) . */
/*        ESTIMATE = NORM(Z)/NORM(Y) WHERE  A*Z = Y  AND  A*Y = E . */
/*        THE COMPONENTS OF  E  ARE CHOSEN TO CAUSE MAXIMUM LOCAL */
/*        GROWTH IN THE ELEMENTS OF W  WHERE  CTRANS(R)*W = E . */
/*        THE VECTORS ARE FREQUENTLY RESCALED TO AVOID OVERFLOW. */

/*        SOLVE CTRANS(R)*W = E */

    ek.r = 1.f, ek.i = 0.f;
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	i__2 = j;
	z__[i__2].r = 0.f, z__[i__2].i = 0.f;
/* L50: */
    }
    kk = 0;
    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
	kk += k;
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
	i__3 = kk;
	if ((r__1 = q__1.r, dabs(r__1)) + (r__2 = r_imag(&q__1), dabs(r__2)) 
		<= ap[i__3].r) {
	    goto L60;
	}
	i__2 = k;
	q__2.r = ek.r - z__[i__2].r, q__2.i = ek.i - z__[i__2].i;
	q__1.r = q__2.r, q__1.i = q__2.i;
	i__3 = kk;
	s = ap[i__3].r / ((r__1 = q__1.r, dabs(r__1)) + (r__2 = r_imag(&q__1),
		 dabs(r__2)));
	csscal_(n, &s, &z__[1], &c__1);
	q__2.r = s, q__2.i = 0.f;
	q__1.r = q__2.r * ek.r - q__2.i * ek.i, q__1.i = q__2.r * ek.i + 
		q__2.i * ek.r;
	ek.r = q__1.r, ek.i = q__1.i;
L60:
	i__2 = k;
	q__1.r = ek.r - z__[i__2].r, q__1.i = ek.i - z__[i__2].i;
	wk.r = q__1.r, wk.i = q__1.i;
	q__2.r = -ek.r, q__2.i = -ek.i;
	i__2 = k;
	q__1.r = q__2.r - z__[i__2].r, q__1.i = q__2.i - z__[i__2].i;
	wkm.r = q__1.r, wkm.i = q__1.i;
	s = (r__1 = wk.r, dabs(r__1)) + (r__2 = r_imag(&wk), dabs(r__2));
	sm = (r__1 = wkm.r, dabs(r__1)) + (r__2 = r_imag(&wkm), dabs(r__2));
	c_div(&q__1, &wk, &ap[kk]);
	wk.r = q__1.r, wk.i = q__1.i;
	c_div(&q__1, &wkm, &ap[kk]);
	wkm.r = q__1.r, wkm.i = q__1.i;
	kp1 = k + 1;
	kj = kk + k;
	if (kp1 > *n) {
	    goto L100;
	}
	i__2 = *n;
	for (j = kp1; j <= i__2; ++j) {
	    i__3 = j;
	    r_cnjg(&q__4, &ap[kj]);
	    q__3.r = wkm.r * q__4.r - wkm.i * q__4.i, q__3.i = wkm.r * q__4.i 
		    + wkm.i * q__4.r;
	    q__2.r = z__[i__3].r + q__3.r, q__2.i = z__[i__3].i + q__3.i;
	    q__1.r = q__2.r, q__1.i = q__2.i;
	    sm += (r__1 = q__1.r, dabs(r__1)) + (r__2 = r_imag(&q__1), dabs(
		    r__2));
	    i__3 = j;
	    i__4 = j;
	    r_cnjg(&q__3, &ap[kj]);
	    q__2.r = wk.r * q__3.r - wk.i * q__3.i, q__2.i = wk.r * q__3.i + 
		    wk.i * q__3.r;
	    q__1.r = z__[i__4].r + q__2.r, q__1.i = z__[i__4].i + q__2.i;
	    z__[i__3].r = q__1.r, z__[i__3].i = q__1.i;
	    i__3 = j;
	    s += (r__1 = z__[i__3].r, dabs(r__1)) + (r__2 = r_imag(&z__[j]), 
		    dabs(r__2));
	    kj += j;
/* L70: */
	}
	if (s >= sm) {
	    goto L90;
	}
	q__1.r = wkm.r - wk.r, q__1.i = wkm.i - wk.i;
	t.r = q__1.r, t.i = q__1.i;
	wk.r = wkm.r, wk.i = wkm.i;
	kj = kk + k;
	i__2 = *n;
	for (j = kp1; j <= i__2; ++j) {
	    i__3 = j;
	    i__4 = j;
	    r_cnjg(&q__3, &ap[kj]);
	    q__2.r = t.r * q__3.r - t.i * q__3.i, q__2.i = t.r * q__3.i + t.i 
		    * q__3.r;
	    q__1.r = z__[i__4].r + q__2.r, q__1.i = z__[i__4].i + q__2.i;
	    z__[i__3].r = q__1.r, z__[i__3].i = q__1.i;
	    kj += j;
/* L80: */
	}
L90:
L100:
	i__2 = k;
	z__[i__2].r = wk.r, z__[i__2].i = wk.i;
/* L110: */
    }
    s = 1.f / scasum_(n, &z__[1], &c__1);
    csscal_(n, &s, &z__[1], &c__1);

/*        SOLVE R*Y = W */

    i__1 = *n;
    for (kb = 1; kb <= i__1; ++kb) {
	k = *n + 1 - kb;
	i__2 = k;
	i__3 = kk;
	if ((r__1 = z__[i__2].r, dabs(r__1)) + (r__2 = r_imag(&z__[k]), dabs(
		r__2)) <= ap[i__3].r) {
	    goto L120;
	}
	i__2 = kk;
	i__3 = k;
	s = ap[i__2].r / ((r__1 = z__[i__3].r, dabs(r__1)) + (r__2 = r_imag(&
		z__[k]), dabs(r__2)));
	csscal_(n, &s, &z__[1], &c__1);
L120:
	i__2 = k;
	c_div(&q__1, &z__[k], &ap[kk]);
	z__[i__2].r = q__1.r, z__[i__2].i = q__1.i;
	kk -= k;
	i__2 = k;
	q__1.r = -z__[i__2].r, q__1.i = -z__[i__2].i;
	t.r = q__1.r, t.i = q__1.i;
	i__2 = k - 1;
	caxpy_(&i__2, &t, &ap[kk + 1], &c__1, &z__[1], &c__1);
/* L130: */
    }
    s = 1.f / scasum_(n, &z__[1], &c__1);
    csscal_(n, &s, &z__[1], &c__1);

    ynorm = 1.f;

/*        SOLVE CTRANS(R)*V = Y */

    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
	i__2 = k;
	i__3 = k;
	i__4 = k - 1;
	cdotc_(&q__2, &i__4, &ap[kk + 1], &c__1, &z__[1], &c__1);
	q__1.r = z__[i__3].r - q__2.r, q__1.i = z__[i__3].i - q__2.i;
	z__[i__2].r = q__1.r, z__[i__2].i = q__1.i;
	kk += k;
	i__2 = k;
	i__3 = kk;
	if ((r__1 = z__[i__2].r, dabs(r__1)) + (r__2 = r_imag(&z__[k]), dabs(
		r__2)) <= ap[i__3].r) {
	    goto L140;
	}
	i__2 = kk;
	i__3 = k;
	s = ap[i__2].r / ((r__1 = z__[i__3].r, dabs(r__1)) + (r__2 = r_imag(&
		z__[k]), dabs(r__2)));
	csscal_(n, &s, &z__[1], &c__1);
	ynorm = s * ynorm;
L140:
	i__2 = k;
	c_div(&q__1, &z__[k], &ap[kk]);
	z__[i__2].r = q__1.r, z__[i__2].i = q__1.i;
/* L150: */
    }
    s = 1.f / scasum_(n, &z__[1], &c__1);
    csscal_(n, &s, &z__[1], &c__1);
    ynorm = s * ynorm;

/*        SOLVE R*Z = V */

    i__1 = *n;
    for (kb = 1; kb <= i__1; ++kb) {
	k = *n + 1 - kb;
	i__2 = k;
	i__3 = kk;
	if ((r__1 = z__[i__2].r, dabs(r__1)) + (r__2 = r_imag(&z__[k]), dabs(
		r__2)) <= ap[i__3].r) {
	    goto L160;
	}
	i__2 = kk;
	i__3 = k;
	s = ap[i__2].r / ((r__1 = z__[i__3].r, dabs(r__1)) + (r__2 = r_imag(&
		z__[k]), dabs(r__2)));
	csscal_(n, &s, &z__[1], &c__1);
	ynorm = s * ynorm;
L160:
	i__2 = k;
	c_div(&q__1, &z__[k], &ap[kk]);
	z__[i__2].r = q__1.r, z__[i__2].i = q__1.i;
	kk -= k;
	i__2 = k;
	q__1.r = -z__[i__2].r, q__1.i = -z__[i__2].i;
	t.r = q__1.r, t.i = q__1.i;
	i__2 = k - 1;
	caxpy_(&i__2, &t, &ap[kk + 1], &c__1, &z__[1], &c__1);
/* L170: */
    }
/*        MAKE ZNORM = 1.0 */
    s = 1.f / scasum_(n, &z__[1], &c__1);
    csscal_(n, &s, &z__[1], &c__1);
    ynorm = s * ynorm;

    if (anorm != 0.f) {
	*rcond = ynorm / anorm;
    }
    if (anorm == 0.f) {
	*rcond = 0.f;
    }
L180:
    return 0;
} /* cppco_ */

