/* sppdi.f -- translated by f2c (version 12.02.01).
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

/* DECK SPPDI */
/* Subroutine */ int sppdi_(real *ap, integer *n, real *det, integer *job)
{
    /* System generated locals */
    integer i__1, i__2;
    real r__1;

    /* Local variables */
    static integer i__, j, k;
    static real s, t;
    static integer j1, k1, ii, jj, kj, kk, jm1, kp1;
    extern /* Subroutine */ int sscal_(integer *, real *, real *, integer *), 
	    saxpy_(integer *, real *, real *, integer *, real *, integer *);

/* ***BEGIN PROLOGUE  SPPDI */
/* ***PURPOSE  Compute the determinant and inverse of a real symmetric */
/*            positive definite matrix using factors from SPPCO or SPPFA. */
/* ***LIBRARY   SLATEC (LINPACK) */
/* ***CATEGORY  D2B1B, D3B1B */
/* ***TYPE      SINGLE PRECISION (SPPDI-S, DPPDI-D, CPPDI-C) */
/* ***KEYWORDS  DETERMINANT, INVERSE, LINEAR ALGEBRA, LINPACK, MATRIX, */
/*             PACKED, POSITIVE DEFINITE */
/* ***AUTHOR  Moler, C. B., (U. of New Mexico) */
/* ***DESCRIPTION */

/*     SPPDI computes the determinant and inverse */
/*     of a real symmetric positive definite matrix */
/*     using the factors computed by SPPCO or SPPFA . */

/*     On Entry */

/*        AP      REAL (N*(N+1)/2) */
/*                the output from SPPCO or SPPFA. */

/*        N       INTEGER */
/*                the order of the matrix  A . */

/*        JOB     INTEGER */
/*                = 11   both determinant and inverse. */
/*                = 01   inverse only. */
/*                = 10   determinant only. */

/*     On Return */

/*        AP      the upper triangular half of the inverse . */
/*                The strict lower triangle is unaltered. */

/*        DET     REAL(2) */
/*                determinant of original matrix if requested. */
/*                Otherwise not referenced. */
/*                Determinant = DET(1) * 10.0**DET(2) */
/*                with  1.0 .LE. DET(1) .LT. 10.0 */
/*                or  DET(1) .EQ. 0.0 . */

/*     Error Condition */

/*        A division by zero will occur if the input factor contains */
/*        a zero on the diagonal and the inverse is requested. */
/*        It will not occur if the subroutines are called correctly */
/*        and if SPOCO or SPOFA has set INFO .EQ. 0 . */

/* ***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W. */
/*                 Stewart, LINPACK Users' Guide, SIAM, 1979. */
/* ***ROUTINES CALLED  SAXPY, SSCAL */
/* ***REVISION HISTORY  (YYMMDD) */
/*   780814  DATE WRITTEN */
/*   890831  Modified array declarations.  (WRB) */
/*   890831  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900326  Removed duplicate information from DESCRIPTION section. */
/*           (WRB) */
/*   920501  Reformatted the REFERENCES section.  (WRB) */
/* ***END PROLOGUE  SPPDI */

/* ***FIRST EXECUTABLE STATEMENT  SPPDI */

/*     COMPUTE DETERMINANT */

    /* Parameter adjustments */
    --det;
    --ap;

    /* Function Body */
    if (*job / 10 == 0) {
	goto L70;
    }
    det[1] = 1.f;
    det[2] = 0.f;
    s = 10.f;
    ii = 0;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ii += i__;
/* Computing 2nd power */
	r__1 = ap[ii];
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
    kk = 0;
    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
	k1 = kk + 1;
	kk += k;
	ap[kk] = 1.f / ap[kk];
	t = -ap[kk];
	i__2 = k - 1;
	sscal_(&i__2, &t, &ap[k1], &c__1);
	kp1 = k + 1;
	j1 = kk + 1;
	kj = kk + k;
	if (*n < kp1) {
	    goto L90;
	}
	i__2 = *n;
	for (j = kp1; j <= i__2; ++j) {
	    t = ap[kj];
	    ap[kj] = 0.f;
	    saxpy_(&k, &t, &ap[k1], &c__1, &ap[j1], &c__1);
	    j1 += j;
	    kj += j;
/* L80: */
	}
L90:
/* L100: */
	;
    }

/*        FORM  INVERSE(R) * TRANS(INVERSE(R)) */

    jj = 0;
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	j1 = jj + 1;
	jj += j;
	jm1 = j - 1;
	k1 = 1;
	kj = j1;
	if (jm1 < 1) {
	    goto L120;
	}
	i__2 = jm1;
	for (k = 1; k <= i__2; ++k) {
	    t = ap[kj];
	    saxpy_(&k, &t, &ap[j1], &c__1, &ap[k1], &c__1);
	    k1 += k;
	    ++kj;
/* L110: */
	}
L120:
	t = ap[jj];
	sscal_(&j, &t, &ap[j1], &c__1);
/* L130: */
    }
L140:
    return 0;
} /* sppdi_ */

