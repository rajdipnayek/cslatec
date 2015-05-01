/* cgefa.f -- translated by f2c (version 12.02.01).
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
static complex c_b7 = {-1.f,-0.f};

/* DECK CGEFA */
/* Subroutine */ int cgefa_(complex *a, integer *lda, integer *n, integer *
	ipvt, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3, i__4;
    real r__1, r__2;
    complex q__1;

    /* Local variables */
    static integer j, k, l;
    static complex t;
    static integer kp1, nm1;
    extern /* Subroutine */ int cscal_(integer *, complex *, complex *, 
	    integer *), caxpy_(integer *, complex *, complex *, integer *, 
	    complex *, integer *);
    extern integer icamax_(integer *, complex *, integer *);

/* ***BEGIN PROLOGUE  CGEFA */
/* ***PURPOSE  Factor a matrix using Gaussian elimination. */
/* ***LIBRARY   SLATEC (LINPACK) */
/* ***CATEGORY  D2C1 */
/* ***TYPE      COMPLEX (SGEFA-S, DGEFA-D, CGEFA-C) */
/* ***KEYWORDS  GENERAL MATRIX, LINEAR ALGEBRA, LINPACK, */
/*             MATRIX FACTORIZATION */
/* ***AUTHOR  Moler, C. B., (U. of New Mexico) */
/* ***DESCRIPTION */

/*     CGEFA factors a complex matrix by Gaussian elimination. */

/*     CGEFA is usually called by CGECO, but it can be called */
/*     directly with a saving in time if  RCOND  is not needed. */
/*     (Time for CGECO) = (1 + 9/N)*(Time for CGEFA) . */

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

/*        INFO    INTEGER */
/*                = 0  normal value. */
/*                = K  if  U(K,K) .EQ. 0.0 .  This is not an error */
/*                     condition for this subroutine, but it does */
/*                     indicate that CGESL or CGEDI will divide by zero */
/*                     if called.  Use  RCOND  in CGECO for a reliable */
/*                     indication of singularity. */

/* ***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W. */
/*                 Stewart, LINPACK Users' Guide, SIAM, 1979. */
/* ***ROUTINES CALLED  CAXPY, CSCAL, ICAMAX */
/* ***REVISION HISTORY  (YYMMDD) */
/*   780814  DATE WRITTEN */
/*   890831  Modified array declarations.  (WRB) */
/*   890831  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900326  Removed duplicate information from DESCRIPTION section. */
/*           (WRB) */
/*   920501  Reformatted the REFERENCES section.  (WRB) */
/* ***END PROLOGUE  CGEFA */


/*     GAUSSIAN ELIMINATION WITH PARTIAL PIVOTING */

/* ***FIRST EXECUTABLE STATEMENT  CGEFA */
    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --ipvt;

    /* Function Body */
    *info = 0;
    nm1 = *n - 1;
    if (nm1 < 1) {
	goto L70;
    }
    i__1 = nm1;
    for (k = 1; k <= i__1; ++k) {
	kp1 = k + 1;

/*        FIND L = PIVOT INDEX */

	i__2 = *n - k + 1;
	l = icamax_(&i__2, &a[k + k * a_dim1], &c__1) + k - 1;
	ipvt[k] = l;

/*        ZERO PIVOT IMPLIES THIS COLUMN ALREADY TRIANGULARIZED */

	i__2 = l + k * a_dim1;
	if ((r__1 = a[i__2].r, dabs(r__1)) + (r__2 = r_imag(&a[l + k * a_dim1]
		), dabs(r__2)) == 0.f) {
	    goto L40;
	}

/*           INTERCHANGE IF NECESSARY */

	if (l == k) {
	    goto L10;
	}
	i__2 = l + k * a_dim1;
	t.r = a[i__2].r, t.i = a[i__2].i;
	i__2 = l + k * a_dim1;
	i__3 = k + k * a_dim1;
	a[i__2].r = a[i__3].r, a[i__2].i = a[i__3].i;
	i__2 = k + k * a_dim1;
	a[i__2].r = t.r, a[i__2].i = t.i;
L10:

/*           COMPUTE MULTIPLIERS */

	c_div(&q__1, &c_b7, &a[k + k * a_dim1]);
	t.r = q__1.r, t.i = q__1.i;
	i__2 = *n - k;
	cscal_(&i__2, &t, &a[k + 1 + k * a_dim1], &c__1);

/*           ROW ELIMINATION WITH COLUMN INDEXING */

	i__2 = *n;
	for (j = kp1; j <= i__2; ++j) {
	    i__3 = l + j * a_dim1;
	    t.r = a[i__3].r, t.i = a[i__3].i;
	    if (l == k) {
		goto L20;
	    }
	    i__3 = l + j * a_dim1;
	    i__4 = k + j * a_dim1;
	    a[i__3].r = a[i__4].r, a[i__3].i = a[i__4].i;
	    i__3 = k + j * a_dim1;
	    a[i__3].r = t.r, a[i__3].i = t.i;
L20:
	    i__3 = *n - k;
	    caxpy_(&i__3, &t, &a[k + 1 + k * a_dim1], &c__1, &a[k + 1 + j * 
		    a_dim1], &c__1);
/* L30: */
	}
	goto L50;
L40:
	*info = k;
L50:
/* L60: */
	;
    }
L70:
    ipvt[*n] = *n;
    i__1 = *n + *n * a_dim1;
    if ((r__1 = a[i__1].r, dabs(r__1)) + (r__2 = r_imag(&a[*n + *n * a_dim1]),
	     dabs(r__2)) == 0.f) {
	*info = *n;
    }
    return 0;
} /* cgefa_ */

