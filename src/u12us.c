/* u12us.f -- translated by f2c (version 12.02.01).
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

/* DECK U12US */
/* Subroutine */ int u12us_(real *a, integer *mda, integer *m, integer *n, 
	real *b, integer *mdb, integer *nb, integer *mode, integer *krank, 
	real *rnorm, real *h__, real *w, integer *ir, integer *ic)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, i__1, i__2, i__3;
    real r__1;

    /* Local variables */
    static integer i__, j, k;
    static real bb;
    static integer jb, ij;
    static real tt;
    static integer ip1, kp1, mmk;
    extern doublereal sdot_(integer *, real *, integer *, real *, integer *), 
	    snrm2_(integer *, real *, integer *);
    extern /* Subroutine */ int sswap_(integer *, real *, integer *, real *, 
	    integer *), saxpy_(integer *, real *, real *, integer *, real *, 
	    integer *);

/* ***BEGIN PROLOGUE  U12US */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to ULSIA */
/* ***LIBRARY   SLATEC */
/* ***TYPE      SINGLE PRECISION (U12US-S, DU12US-D) */
/* ***AUTHOR  (UNKNOWN) */
/* ***DESCRIPTION */

/*        Given the Householder LQ factorization of A, this */
/*        subroutine solves the system AX=B. If the system */
/*        is of reduced rank, this routine returns a solution */
/*        according to the selected mode. */

/*       Note - If MODE.NE.2, W is never accessed. */

/* ***SEE ALSO  ULSIA */
/* ***ROUTINES CALLED  SAXPY, SDOT, SNRM2, SSWAP */
/* ***REVISION HISTORY  (YYMMDD) */
/*   810801  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890831  Modified array declarations.  (WRB) */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900328  Added TYPE section.  (WRB) */
/* ***END PROLOGUE  U12US */
/* ***FIRST EXECUTABLE STATEMENT  U12US */
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
    --ir;
    --ic;

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

/*     IF A IS OF REDUCED RANK AND MODE=2, */
/*     APPLY HOUSEHOLDER TRANSFORMATIONS TO B */

    if (*mode < 2 || k == *m) {
	goto L440;
    }
    mmk = *m - k;
    i__2 = *nb;
    for (jb = 1; jb <= i__2; ++jb) {
	i__1 = k;
	for (j = 1; j <= i__1; ++j) {
	    i__ = kp1 - j;
	    tt = -sdot_(&mmk, &a[kp1 + i__ * a_dim1], &c__1, &b[kp1 + jb * 
		    b_dim1], &c__1) / w[i__];
	    tt -= b[i__ + jb * b_dim1];
	    saxpy_(&mmk, &tt, &a[kp1 + i__ * a_dim1], &c__1, &b[kp1 + jb * 
		    b_dim1], &c__1);
	    b[i__ + jb * b_dim1] += tt * w[i__];
/* L425: */
	}
/* L430: */
    }

/*     FIND NORMS OF RESIDUAL VECTOR(S)..(BEFORE OVERWRITE B) */

L440:
    i__2 = *nb;
    for (jb = 1; jb <= i__2; ++jb) {
	i__1 = *m - k;
	rnorm[jb] = snrm2_(&i__1, &b[kp1 + jb * b_dim1], &c__1);
/* L442: */
    }

/*     BACK SOLVE LOWER TRIANGULAR L */

    i__2 = *nb;
    for (jb = 1; jb <= i__2; ++jb) {
	i__1 = k;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    b[i__ + jb * b_dim1] /= a[i__ + i__ * a_dim1];
	    if (i__ == k) {
		goto L450;
	    }
	    ip1 = i__ + 1;
	    i__3 = k - i__;
	    r__1 = -b[i__ + jb * b_dim1];
	    saxpy_(&i__3, &r__1, &a[ip1 + i__ * a_dim1], &c__1, &b[ip1 + jb * 
		    b_dim1], &c__1);
/* L448: */
	}
L450:
	;
    }


/*      TRUNCATED SOLUTION */

    if (k == *n) {
	goto L462;
    }
    i__2 = *nb;
    for (jb = 1; jb <= i__2; ++jb) {
	i__1 = *n;
	for (i__ = kp1; i__ <= i__1; ++i__) {
	    b[i__ + jb * b_dim1] = 0.f;
/* L460: */
	}
    }

/*     APPLY HOUSEHOLDER TRANSFORMATIONS TO B */

L462:
    i__1 = k;
    for (i__ = 1; i__ <= i__1; ++i__) {
	j = kp1 - i__;
	tt = a[j + j * a_dim1];
	a[j + j * a_dim1] = h__[j];
	i__2 = *nb;
	for (jb = 1; jb <= i__2; ++jb) {
	    i__3 = *n - j + 1;
	    bb = -sdot_(&i__3, &a[j + j * a_dim1], mda, &b[j + jb * b_dim1], &
		    c__1) / h__[j];
	    i__3 = *n - j + 1;
	    saxpy_(&i__3, &bb, &a[j + j * a_dim1], mda, &b[j + jb * b_dim1], &
		    c__1);
/* L465: */
	}
	a[j + j * a_dim1] = tt;
/* L470: */
    }


/*     REORDER B TO REFLECT COLUMN INTERCHANGES */

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
} /* u12us_ */

