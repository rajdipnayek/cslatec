/* cpodi.f -- translated by f2c (version 12.02.01).
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

static complex c_b11 = {1.f,0.f};
static integer c__1 = 1;

/* DECK CPODI */
/* Subroutine */ int cpodi_(complex *a, integer *lda, integer *n, real *det, 
	integer *job)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3;
    real r__1;
    complex q__1;

    /* Local variables */
    static integer i__, j, k;
    static real s;
    static complex t;
    static integer jm1, kp1;
    extern /* Subroutine */ int cscal_(integer *, complex *, complex *, 
	    integer *), caxpy_(integer *, complex *, complex *, integer *, 
	    complex *, integer *);

/* ***BEGIN PROLOGUE  CPODI */
/* ***PURPOSE  Compute the determinant and inverse of a certain complex */
/*            Hermitian positive definite matrix using the factors */
/*            computed by CPOCO, CPOFA, or CQRDC. */
/* ***LIBRARY   SLATEC (LINPACK) */
/* ***CATEGORY  D2D1B, D3D1B */
/* ***TYPE      COMPLEX (SPODI-S, DPODI-D, CPODI-C) */
/* ***KEYWORDS  DETERMINANT, INVERSE, LINEAR ALGEBRA, LINPACK, MATRIX, */
/*             POSITIVE DEFINITE */
/* ***AUTHOR  Moler, C. B., (U. of New Mexico) */
/* ***DESCRIPTION */

/*     CPODI computes the determinant and inverse of a certain */
/*     complex Hermitian positive definite matrix (see below) */
/*     using the factors computed by CPOCO, CPOFA or CQRDC. */

/*     On Entry */

/*        A       COMPLEX(LDA, N) */
/*                the output  A  from CPOCO or CPOFA */
/*                or the output  X  from CQRDC. */

/*        LDA     INTEGER */
/*                the leading dimension of the array  A . */

/*        N       INTEGER */
/*                the order of the matrix  A . */

/*        JOB     INTEGER */
/*                = 11   both determinant and inverse. */
/*                = 01   inverse only. */
/*                = 10   determinant only. */

/*     On Return */

/*        A       If CPOCO or CPOFA was used to factor  A  then */
/*                CPODI produces the upper half of INVERSE(A) . */
/*                If CQRDC was used to decompose  X  then */
/*                CPODI produces the upper half of INVERSE(CTRANS(X)*X) */
/*                where CTRANS(X) is the conjugate transpose. */
/*                Elements of  A  below the diagonal are unchanged. */
/*                If the units digit of JOB is zero,  A  is unchanged. */

/*        DET     REAL(2) */
/*                determinant of  A  or of  CTRANS(X)*X  if requested. */
/*                Otherwise not referenced. */
/*                Determinant = DET(1) * 10.0**DET(2) */
/*                with  1.0 .LE. DET(1) .LT. 10.0 */
/*                or  DET(1) .EQ. 0.0 . */

/*     Error Condition */

/*        a division by zero will occur if the input factor contains */
/*        a zero on the diagonal and the inverse is requested. */
/*        It will not occur if the subroutines are called correctly */
/*        and if CPOCO or CPOFA has set INFO .EQ. 0 . */

/* ***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W. */
/*                 Stewart, LINPACK Users' Guide, SIAM, 1979. */
/* ***ROUTINES CALLED  CAXPY, CSCAL */
/* ***REVISION HISTORY  (YYMMDD) */
/*   780814  DATE WRITTEN */
/*   890831  Modified array declarations.  (WRB) */
/*   890831  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900326  Removed duplicate information from DESCRIPTION section. */
/*           (WRB) */
/*   920501  Reformatted the REFERENCES section.  (WRB) */
/* ***END PROLOGUE  CPODI */

/* ***FIRST EXECUTABLE STATEMENT  CPODI */

/*     COMPUTE DETERMINANT */

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --det;

    /* Function Body */
    if (*job / 10 == 0) {
	goto L70;
    }
    det[1] = 1.f;
    det[2] = 0.f;
    s = 10.f;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = i__ + i__ * a_dim1;
/* Computing 2nd power */
	r__1 = a[i__2].r;
	det[1] = r__1 * r__1 * det[1];
	if (det[1] == 0.f) {
	    goto L60;
	}
L10:
	if (det[1] >= 1.f) {
	    goto L20;
	}
	det[1] = s * det[1];
	det[2] += -1.f;
	goto L10;
L20:
L30:
	if (det[1] < s) {
	    goto L40;
	}
	det[1] /= s;
	det[2] += 1.f;
	goto L30;
L40:
/* L50: */
	;
    }
L60:
L70:

/*     COMPUTE INVERSE(R) */

    if (*job % 10 == 0) {
	goto L140;
    }
    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
	i__2 = k + k * a_dim1;
	c_div(&q__1, &c_b11, &a[k + k * a_dim1]);
	a[i__2].r = q__1.r, a[i__2].i = q__1.i;
	i__2 = k + k * a_dim1;
	q__1.r = -a[i__2].r, q__1.i = -a[i__2].i;
	t.r = q__1.r, t.i = q__1.i;
	i__2 = k - 1;
	cscal_(&i__2, &t, &a[k * a_dim1 + 1], &c__1);
	kp1 = k + 1;
	if (*n < kp1) {
	    goto L90;
	}
	i__2 = *n;
	for (j = kp1; j <= i__2; ++j) {
	    i__3 = k + j * a_dim1;
	    t.r = a[i__3].r, t.i = a[i__3].i;
	    i__3 = k + j * a_dim1;
	    a[i__3].r = 0.f, a[i__3].i = 0.f;
	    caxpy_(&k, &t, &a[k * a_dim1 + 1], &c__1, &a[j * a_dim1 + 1], &
		    c__1);
/* L80: */
	}
L90:
/* L100: */
	;
    }

/*        FORM  INVERSE(R) * CTRANS(INVERSE(R)) */

    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	jm1 = j - 1;
	if (jm1 < 1) {
	    goto L120;
	}
	i__2 = jm1;
	for (k = 1; k <= i__2; ++k) {
	    r_cnjg(&q__1, &a[k + j * a_dim1]);
	    t.r = q__1.r, t.i = q__1.i;
	    caxpy_(&k, &t, &a[j * a_dim1 + 1], &c__1, &a[k * a_dim1 + 1], &
		    c__1);
/* L110: */
	}
L120:
	r_cnjg(&q__1, &a[j + j * a_dim1]);
	t.r = q__1.r, t.i = q__1.i;
	cscal_(&j, &t, &a[j * a_dim1 + 1], &c__1);
/* L130: */
    }
L140:
    return 0;
} /* cpodi_ */

