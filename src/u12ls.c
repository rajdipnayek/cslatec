/* u12ls.f -- translated by f2c (version 12.02.01).
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

/* DECK U12LS */
/* Subroutine */ int u12ls_(real *a, integer *mda, integer *m, integer *n, 
	real *b, integer *mdb, integer *nb, integer *mode, integer *krank, 
	real *rnorm, real *h__, real *w, integer *ic, integer *ir)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, i__1, i__2, i__3;
    real r__1;

    /* Local variables */
    static integer i__, j, k;
    static real bb;
    static integer jb, ij;
    static real tt;
    static integer im1, kp1, nmk;
    extern doublereal sdot_(integer *, real *, integer *, real *, integer *), 
	    snrm2_(integer *, real *, integer *);
    extern /* Subroutine */ int sswap_(integer *, real *, integer *, real *, 
	    integer *), saxpy_(integer *, real *, real *, integer *, real *, 
	    integer *);

/* ***BEGIN PROLOGUE  U12LS */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to LLSIA */
/* ***LIBRARY   SLATEC */
/* ***TYPE      SINGLE PRECISION (U12LS-S, DU12LS-D) */
/* ***AUTHOR  (UNKNOWN) */
/* ***DESCRIPTION */

/*        Given the Householder QR factorization of A, this */
/*        subroutine solves the system AX=B. If the system */
/*        is of reduced rank, this routine returns a solution */
/*        according to the selected mode. */

/*       Note - If MODE.NE.2, W is never accessed. */

/* ***SEE ALSO  LLSIA */
/* ***ROUTINES CALLED  SAXPY, SDOT, SNRM2, SSWAP */
/* ***REVISION HISTORY  (YYMMDD) */
/*   810801  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890831  Modified array declarations.  (WRB) */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900328  Added TYPE section.  (WRB) */
/* ***END PROLOGUE  U12LS */
/* ***FIRST EXECUTABLE STATEMENT  U12LS */
    /* Parameter adjustments */
    a_dim1 = *mda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    b_dim1 = *mdb;
    b_offset = 1 + b_dim1;
    b -= b_offset;
    --rnorm;
    --h__;
    --w;
    --ic;
    --ir;

    /* Function Body */
    k = *krank;
    kp1 = k + 1;

/*        RANK=0 */

    if (k > 0) {
	goto L410;
    }
    i__1 = *nb;
    for (jb = 1; jb <= i__1; ++jb) {
	rnorm[jb] = snrm2_(m, &b[jb * b_dim1 + 1], &c__1);
/* L404: */
    }
    i__1 = *nb;
    for (jb = 1; jb <= i__1; ++jb) {
	i__2 = *n;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    b[i__ + jb * b_dim1] = 0.f;
/* L406: */
	}
    }
    return 0;

/*     REORDER B TO REFLECT ROW INTERCHANGES */

L410:
    i__ = 0;
L412:
    ++i__;
    if (i__ == *m) {
	goto L418;
    }
    j = ir[i__];
    if (j == i__) {
	goto L412;
    }
    if (j < 0) {
	goto L412;
    }
    ir[i__] = -ir[i__];
    i__2 = *nb;
    for (jb = 1; jb <= i__2; ++jb) {
	rnorm[jb] = b[i__ + jb * b_dim1];
/* L413: */
    }
    ij = i__;
L414:
    i__2 = *nb;
    for (jb = 1; jb <= i__2; ++jb) {
	b[ij + jb * b_dim1] = b[j + jb * b_dim1];
/* L415: */
    }
    ij = j;
    j = ir[ij];
    ir[ij] = -ir[ij];
    if (j != i__) {
	goto L414;
    }
    i__2 = *nb;
    for (jb = 1; jb <= i__2; ++jb) {
	b[ij + jb * b_dim1] = rnorm[jb];
/* L416: */
    }
    goto L412;
L418:
    i__2 = *m;
    for (i__ = 1; i__ <= i__2; ++i__) {
	ir[i__] = (i__1 = ir[i__], abs(i__1));
/* L420: */
    }

/*     APPLY HOUSEHOLDER TRANSFORMATIONS TO B */

    i__2 = k;
    for (j = 1; j <= i__2; ++j) {
	tt = a[j + j * a_dim1];
	a[j + j * a_dim1] = h__[j];
	i__1 = *nb;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    i__3 = *m - j + 1;
	    bb = -sdot_(&i__3, &a[j + j * a_dim1], &c__1, &b[j + i__ * b_dim1]
		    , &c__1) / h__[j];
	    i__3 = *m - j + 1;
	    saxpy_(&i__3, &bb, &a[j + j * a_dim1], &c__1, &b[j + i__ * b_dim1]
		    , &c__1);
/* L425: */
	}
	a[j + j * a_dim1] = tt;
/* L430: */
    }

/*        FIND NORMS OF RESIDUAL VECTOR(S)..(BEFORE OVERWRITE B) */

    i__2 = *nb;
    for (jb = 1; jb <= i__2; ++jb) {
	i__1 = *m - k;
	rnorm[jb] = snrm2_(&i__1, &b[kp1 + jb * b_dim1], &c__1);
/* L440: */
    }

/*     BACK SOLVE UPPER TRIANGULAR R */

    i__ = k;
L442:
    i__2 = *nb;
    for (jb = 1; jb <= i__2; ++jb) {
	b[i__ + jb * b_dim1] /= a[i__ + i__ * a_dim1];
/* L444: */
    }
    if (i__ == 1) {
	goto L450;
    }
    im1 = i__ - 1;
    i__2 = *nb;
    for (jb = 1; jb <= i__2; ++jb) {
	r__1 = -b[i__ + jb * b_dim1];
	saxpy_(&im1, &r__1, &a[i__ * a_dim1 + 1], &c__1, &b[jb * b_dim1 + 1], 
		&c__1);
/* L448: */
    }
    i__ = im1;
    goto L442;
L450:

/*     RANK LT N */

/*      TRUNCATED SOLUTION */

    if (k == *n) {
	goto L480;
    }
    i__2 = *nb;
    for (jb = 1; jb <= i__2; ++jb) {
	i__1 = *n;
	for (i__ = kp1; i__ <= i__1; ++i__) {
	    b[i__ + jb * b_dim1] = 0.f;
/* L460: */
	}
    }
    if (*mode == 1) {
	goto L480;
    }

/*      MINIMAL LENGTH SOLUTION */

    nmk = *n - k;
    i__1 = *nb;
    for (jb = 1; jb <= i__1; ++jb) {
	i__2 = k;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    tt = -sdot_(&nmk, &a[i__ + kp1 * a_dim1], mda, &b[kp1 + jb * 
		    b_dim1], &c__1) / w[i__];
	    tt -= b[i__ + jb * b_dim1];
	    saxpy_(&nmk, &tt, &a[i__ + kp1 * a_dim1], mda, &b[kp1 + jb * 
		    b_dim1], &c__1);
	    b[i__ + jb * b_dim1] += tt * w[i__];
/* L465: */
	}
/* L470: */
    }


/*     REORDER B TO REFLECT COLUMN INTERCHANGES */

L480:
    i__ = 0;
L482:
    ++i__;
    if (i__ == *n) {
	goto L488;
    }
    j = ic[i__];
    if (j == i__) {
	goto L482;
    }
    if (j < 0) {
	goto L482;
    }
    ic[i__] = -ic[i__];
L484:
    sswap_(nb, &b[j + b_dim1], mdb, &b[i__ + b_dim1], mdb);
    ij = ic[j];
    ic[j] = -ic[j];
    j = ij;
    if (j == i__) {
	goto L482;
    }
    goto L484;
L488:
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ic[i__] = (i__2 = ic[i__], abs(i__2));
/* L490: */
    }

/*        SOLUTION VECTORS ARE IN FIRST N ROWS OF B(,) */

    return 0;
} /* u12ls_ */

