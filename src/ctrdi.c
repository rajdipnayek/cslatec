/* ctrdi.f -- translated by f2c (version 12.02.01).
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

static complex c_b3 = {1.f,0.f};
static integer c__1 = 1;

/* DECK CTRDI */
/* Subroutine */ int ctrdi_(complex *t, integer *ldt, integer *n, complex *
	det, integer *job, integer *info)
{
    /* System generated locals */
    integer t_dim1, t_offset, i__1, i__2, i__3;
    real r__1, r__2;
    complex q__1, q__2;

    /* Local variables */
    static integer i__, j, k, kb, km1, kp1;
    static real ten;
    static complex temp;
    extern /* Subroutine */ int cscal_(integer *, complex *, complex *, 
	    integer *), caxpy_(integer *, complex *, complex *, integer *, 
	    complex *, integer *);

/* ***BEGIN PROLOGUE  CTRDI */
/* ***PURPOSE  Compute the determinant and inverse of a triangular matrix. */
/* ***LIBRARY   SLATEC (LINPACK) */
/* ***CATEGORY  D2C3, D3C3 */
/* ***TYPE      COMPLEX (STRDI-S, DTRDI-D, CTRDI-C) */
/* ***KEYWORDS  DETERMINANT, INVERSE, LINEAR ALGEBRA, LINPACK, */
/*             TRIANGULAR MATRIX */
/* ***AUTHOR  Moler, C. B., (U. of New Mexico) */
/* ***DESCRIPTION */

/*     CTRDI computes the determinant and inverse of a complex */
/*     triangular matrix. */

/*     On Entry */

/*        T       COMPLEX(LDT,N) */
/*                T contains the triangular matrix.  The zero */
/*                elements of the matrix are not referenced, and */
/*                the corresponding elements of the array can be */
/*                used to store other information. */

/*        LDT     INTEGER */
/*                LDT is the leading dimension of the array T. */

/*        N       INTEGER */
/*                N is the order of the system. */

/*        JOB     INTEGER */
/*                = 010       no det, inverse of lower triangular. */
/*                = 011       no det, inverse of upper triangular. */
/*                = 100       det, no inverse. */
/*                = 110       det, inverse of lower triangular. */
/*                = 111       det, inverse of upper triangular. */

/*     On Return */

/*        T       inverse of original matrix if requested. */
/*                Otherwise unchanged. */

/*        DET     COMPLEX(2) */
/*                determinant of original matrix if requested. */
/*                Otherwise not referenced. */
/*                Determinant = DET(1) * 10.0**DET(2) */
/*                with  1.0 .LE. CABS1(DET(1)) .LT. 10.0 */
/*                or  DET(1) .EQ. 0.0 . */

/*        INFO    INTEGER */
/*                INFO contains zero if the system is nonsingular */
/*                and the inverse is requested. */
/*                Otherwise INFO contains the index of */
/*                a zero diagonal element of T. */


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
/* ***END PROLOGUE  CTRDI */

/* ***FIRST EXECUTABLE STATEMENT  CTRDI */

/*        COMPUTE DETERMINANT */

    /* Parameter adjustments */
    t_dim1 = *ldt;
    t_offset = 1 + t_dim1;
    t -= t_offset;
    --det;

    /* Function Body */
    if (*job / 100 == 0) {
	goto L70;
    }
    det[1].r = 1.f, det[1].i = 0.f;
    det[2].r = 0.f, det[2].i = 0.f;
    ten = 10.f;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = i__ + i__ * t_dim1;
	q__1.r = t[i__2].r * det[1].r - t[i__2].i * det[1].i, q__1.i = t[i__2]
		.r * det[1].i + t[i__2].i * det[1].r;
	det[1].r = q__1.r, det[1].i = q__1.i;
	if ((r__1 = det[1].r, dabs(r__1)) + (r__2 = r_imag(&det[1]), dabs(
		r__2)) == 0.f) {
	    goto L60;
	}
L10:
	if ((r__1 = det[1].r, dabs(r__1)) + (r__2 = r_imag(&det[1]), dabs(
		r__2)) >= 1.f) {
	    goto L20;
	}
	q__2.r = ten, q__2.i = 0.f;
	q__1.r = q__2.r * det[1].r - q__2.i * det[1].i, q__1.i = q__2.r * det[
		1].i + q__2.i * det[1].r;
	det[1].r = q__1.r, det[1].i = q__1.i;
	q__1.r = det[2].r - 1.f, q__1.i = det[2].i - 0.f;
	det[2].r = q__1.r, det[2].i = q__1.i;
	goto L10;
L20:
L30:
	if ((r__1 = det[1].r, dabs(r__1)) + (r__2 = r_imag(&det[1]), dabs(
		r__2)) < ten) {
	    goto L40;
	}
	q__2.r = ten, q__2.i = 0.f;
	c_div(&q__1, &det[1], &q__2);
	det[1].r = q__1.r, det[1].i = q__1.i;
	q__1.r = det[2].r + 1.f, q__1.i = det[2].i + 0.f;
	det[2].r = q__1.r, det[2].i = q__1.i;
	goto L30;
L40:
/* L50: */
	;
    }
L60:
L70:

/*        COMPUTE INVERSE OF UPPER TRIANGULAR */

    if (*job / 10 % 10 == 0) {
	goto L170;
    }
    if (*job % 10 == 0) {
	goto L120;
    }
    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
	*info = k;
	i__2 = k + k * t_dim1;
	if ((r__1 = t[i__2].r, dabs(r__1)) + (r__2 = r_imag(&t[k + k * t_dim1]
		), dabs(r__2)) == 0.f) {
	    goto L110;
	}
	i__2 = k + k * t_dim1;
	c_div(&q__1, &c_b3, &t[k + k * t_dim1]);
	t[i__2].r = q__1.r, t[i__2].i = q__1.i;
	i__2 = k + k * t_dim1;
	q__1.r = -t[i__2].r, q__1.i = -t[i__2].i;
	temp.r = q__1.r, temp.i = q__1.i;
	i__2 = k - 1;
	cscal_(&i__2, &temp, &t[k * t_dim1 + 1], &c__1);
	kp1 = k + 1;
	if (*n < kp1) {
	    goto L90;
	}
	i__2 = *n;
	for (j = kp1; j <= i__2; ++j) {
	    i__3 = k + j * t_dim1;
	    temp.r = t[i__3].r, temp.i = t[i__3].i;
	    i__3 = k + j * t_dim1;
	    t[i__3].r = 0.f, t[i__3].i = 0.f;
	    caxpy_(&k, &temp, &t[k * t_dim1 + 1], &c__1, &t[j * t_dim1 + 1], &
		    c__1);
/* L80: */
	}
L90:
/* L100: */
	;
    }
    *info = 0;
L110:
    goto L160;
L120:

/*              COMPUTE INVERSE OF LOWER TRIANGULAR */

    i__1 = *n;
    for (kb = 1; kb <= i__1; ++kb) {
	k = *n + 1 - kb;
	*info = k;
	i__2 = k + k * t_dim1;
	if ((r__1 = t[i__2].r, dabs(r__1)) + (r__2 = r_imag(&t[k + k * t_dim1]
		), dabs(r__2)) == 0.f) {
	    goto L180;
	}
	i__2 = k + k * t_dim1;
	c_div(&q__1, &c_b3, &t[k + k * t_dim1]);
	t[i__2].r = q__1.r, t[i__2].i = q__1.i;
	i__2 = k + k * t_dim1;
	q__1.r = -t[i__2].r, q__1.i = -t[i__2].i;
	temp.r = q__1.r, temp.i = q__1.i;
	if (k != *n) {
	    i__2 = *n - k;
	    cscal_(&i__2, &temp, &t[k + 1 + k * t_dim1], &c__1);
	}
	km1 = k - 1;
	if (km1 < 1) {
	    goto L140;
	}
	i__2 = km1;
	for (j = 1; j <= i__2; ++j) {
	    i__3 = k + j * t_dim1;
	    temp.r = t[i__3].r, temp.i = t[i__3].i;
	    i__3 = k + j * t_dim1;
	    t[i__3].r = 0.f, t[i__3].i = 0.f;
	    i__3 = *n - k + 1;
	    caxpy_(&i__3, &temp, &t[k + k * t_dim1], &c__1, &t[k + j * t_dim1]
		    , &c__1);
/* L130: */
	}
L140:
/* L150: */
	;
    }
    *info = 0;
L160:
L170:
L180:
    return 0;
} /* ctrdi_ */

