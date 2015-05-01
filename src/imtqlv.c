/* imtqlv.f -- translated by f2c (version 12.02.01).
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

static real c_b11 = 1.f;

/* DECK IMTQLV */
/* Subroutine */ int imtqlv_(integer *n, real *d__, real *e, real *e2, real *
	w, integer *ind, integer *ierr, real *rv1)
{
    /* System generated locals */
    integer i__1, i__2;
    real r__1, r__2;

    /* Local variables */
    static real b, c__, f, g;
    static integer i__, j, k, l, m;
    static real p, r__, s, s1, s2;
    static integer ii, tag, mml;
    extern doublereal pythag_(real *, real *);

/* ***BEGIN PROLOGUE  IMTQLV */
/* ***PURPOSE  Compute the eigenvalues of a symmetric tridiagonal matrix */
/*            using the implicit QL method.  Eigenvectors may be computed */
/*            later. */
/* ***LIBRARY   SLATEC (EISPACK) */
/* ***CATEGORY  D4A5, D4C2A */
/* ***TYPE      SINGLE PRECISION (IMTQLV-S) */
/* ***KEYWORDS  EIGENVALUES, EIGENVECTORS, EISPACK */
/* ***AUTHOR  Smith, B. T., et al. */
/* ***DESCRIPTION */

/*     This subroutine is a variant of  IMTQL1  which is a translation of */
/*     ALGOL procedure IMTQL1, NUM. MATH. 12, 377-383(1968) by Martin and */
/*     Wilkinson, as modified in NUM. MATH. 15, 450(1970) by Dubrulle. */
/*     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 241-248(1971). */

/*     This subroutine finds the eigenvalues of a SYMMETRIC TRIDIAGONAL */
/*     matrix by the implicit QL method and associates with them */
/*     their corresponding submatrix indices. */

/*     On INPUT */

/*        N is the order of the matrix.  N is an INTEGER variable. */

/*        D contains the diagonal elements of the symmetric tridiagonal */
/*          matrix.  D is a one-dimensional REAL array, dimensioned D(N). */

/*        E contains the subdiagonal elements of the symmetric */
/*          tridiagonal matrix in its last N-1 positions.  E(1) is */
/*          arbitrary.  E is a one-dimensional REAL array, dimensioned */
/*          E(N). */

/*        E2 contains the squares of the corresponding elements of E in */
/*          its last N-1 positions.  E2(1) is arbitrary.  E2 is a one- */
/*          dimensional REAL array, dimensioned E2(N). */

/*     On OUTPUT */

/*        D and E are unaltered. */

/*        Elements of E2, corresponding to elements of E regarded as */
/*          negligible, have been replaced by zero causing the matrix to */
/*          split into a direct sum of submatrices.  E2(1) is also set */
/*          to zero. */

/*        W contains the eigenvalues in ascending order.  If an error */
/*          exit is made, the eigenvalues are correct and ordered for */
/*          indices 1, 2, ..., IERR-1, but may not be the smallest */
/*          eigenvalues.  W is a one-dimensional REAL array, dimensioned */
/*          W(N). */

/*        IND contains the submatrix indices associated with the */
/*          corresponding eigenvalues in W -- 1 for eigenvalues belonging */
/*          to the first submatrix from the top, 2 for those belonging to */
/*          the second submatrix, etc.  IND is a one-dimensional REAL */
/*          array, dimensioned IND(N). */

/*        IERR is an INTEGER flag set to */
/*          Zero       for normal return, */
/*          J          if the J-th eigenvalue has not been */
/*                     determined after 30 iterations. */
/*                     The eigenvalues should be correct for indices */
/*                     1, 2, ..., IERR-1.  These eigenvalues are */
/*                     ordered, but are not necessarily the smallest. */

/*        RV1 is a one-dimensional REAL array used for temporary storage, */
/*          dimensioned RV1(N). */

/*     Calls PYTHAG(A,B) for sqrt(A**2 + B**2). */

/*     Questions and comments should be directed to B. S. Garbow, */
/*     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY */
/*     ------------------------------------------------------------------ */

/* ***REFERENCES  B. T. Smith, J. M. Boyle, J. J. Dongarra, B. S. Garbow, */
/*                 Y. Ikebe, V. C. Klema and C. B. Moler, Matrix Eigen- */
/*                 system Routines - EISPACK Guide, Springer-Verlag, */
/*                 1976. */
/* ***ROUTINES CALLED  PYTHAG */
/* ***REVISION HISTORY  (YYMMDD) */
/*   760101  DATE WRITTEN */
/*   890831  Modified array declarations.  (WRB) */
/*   890831  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   920501  Reformatted the REFERENCES section.  (WRB) */
/* ***END PROLOGUE  IMTQLV */


/* ***FIRST EXECUTABLE STATEMENT  IMTQLV */
    /* Parameter adjustments */
    --rv1;
    --ind;
    --w;
    --e2;
    --e;
    --d__;

    /* Function Body */
    *ierr = 0;
    k = 0;
    tag = 0;

    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	w[i__] = d__[i__];
	if (i__ != 1) {
	    rv1[i__ - 1] = e[i__];
	}
/* L100: */
    }

    e2[1] = 0.f;
    rv1[*n] = 0.f;

    i__1 = *n;
    for (l = 1; l <= i__1; ++l) {
	j = 0;
/*     .......... LOOK FOR SMALL SUB-DIAGONAL ELEMENT .......... */
L105:
	i__2 = *n;
	for (m = l; m <= i__2; ++m) {
	    if (m == *n) {
		goto L120;
	    }
	    s1 = (r__1 = w[m], dabs(r__1)) + (r__2 = w[m + 1], dabs(r__2));
	    s2 = s1 + (r__1 = rv1[m], dabs(r__1));
	    if (s2 == s1) {
		goto L120;
	    }
/*     .......... GUARD AGAINST UNDERFLOWED ELEMENT OF E2 .......... */
	    if (e2[m + 1] == 0.f) {
		goto L125;
	    }
/* L110: */
	}

L120:
	if (m <= k) {
	    goto L130;
	}
	if (m != *n) {
	    e2[m + 1] = 0.f;
	}
L125:
	k = m;
	++tag;
L130:
	p = w[l];
	if (m == l) {
	    goto L215;
	}
	if (j == 30) {
	    goto L1000;
	}
	++j;
/*     .......... FORM SHIFT .......... */
	g = (w[l + 1] - p) / (rv1[l] * 2.f);
	r__ = pythag_(&g, &c_b11);
	g = w[m] - p + rv1[l] / (g + r_sign(&r__, &g));
	s = 1.f;
	c__ = 1.f;
	p = 0.f;
	mml = m - l;
/*     .......... FOR I=M-1 STEP -1 UNTIL L DO -- .......... */
	i__2 = mml;
	for (ii = 1; ii <= i__2; ++ii) {
	    i__ = m - ii;
	    f = s * rv1[i__];
	    b = c__ * rv1[i__];
	    if (dabs(f) < dabs(g)) {
		goto L150;
	    }
	    c__ = g / f;
	    r__ = sqrt(c__ * c__ + 1.f);
	    rv1[i__ + 1] = f * r__;
	    s = 1.f / r__;
	    c__ *= s;
	    goto L160;
L150:
	    s = f / g;
	    r__ = sqrt(s * s + 1.f);
	    rv1[i__ + 1] = g * r__;
	    c__ = 1.f / r__;
	    s *= c__;
L160:
	    g = w[i__ + 1] - p;
	    r__ = (w[i__] - g) * s + c__ * 2.f * b;
	    p = s * r__;
	    w[i__ + 1] = g + p;
	    g = c__ * r__ - b;
/* L200: */
	}

	w[l] -= p;
	rv1[l] = g;
	rv1[m] = 0.f;
	goto L105;
/*     .......... ORDER EIGENVALUES .......... */
L215:
	if (l == 1) {
	    goto L250;
	}
/*     .......... FOR I=L STEP -1 UNTIL 2 DO -- .......... */
	i__2 = l;
	for (ii = 2; ii <= i__2; ++ii) {
	    i__ = l + 2 - ii;
	    if (p >= w[i__ - 1]) {
		goto L270;
	    }
	    w[i__] = w[i__ - 1];
	    ind[i__] = ind[i__ - 1];
/* L230: */
	}

L250:
	i__ = 1;
L270:
	w[i__] = p;
	ind[i__] = tag;
/* L290: */
    }

    goto L1001;
/*     .......... SET ERROR -- NO CONVERGENCE TO AN */
/*                EIGENVALUE AFTER 30 ITERATIONS .......... */
L1000:
    *ierr = l;
L1001:
    return 0;
} /* imtqlv_ */

