/* tred3.f -- translated by f2c (version 12.02.01).
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

/* DECK TRED3 */
/* Subroutine */ int tred3_(integer *n, integer *nv, real *a, real *d__, real 
	*e, real *e2)
{
    /* System generated locals */
    integer i__1, i__2, i__3;
    real r__1;

    /* Local variables */
    static real f, g, h__;
    static integer i__, j, k, l;
    static real hh;
    static integer ii, jk, iz;
    static real scale;

/* ***BEGIN PROLOGUE  TRED3 */
/* ***PURPOSE  Reduce a real symmetric matrix stored in packed form to */
/*            symmetric tridiagonal matrix using orthogonal */
/*            transformations. */
/* ***LIBRARY   SLATEC (EISPACK) */
/* ***CATEGORY  D4C1B1 */
/* ***TYPE      SINGLE PRECISION (TRED3-S) */
/* ***KEYWORDS  EIGENVALUES, EIGENVECTORS, EISPACK */
/* ***AUTHOR  Smith, B. T., et al. */
/* ***DESCRIPTION */

/*     This subroutine is a translation of the ALGOL procedure TRED3, */
/*     NUM. MATH. 11, 181-195(1968) by Martin, Reinsch, and Wilkinson. */
/*     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 212-226(1971). */

/*     This subroutine reduces a REAL SYMMETRIC matrix, stored as */
/*     a one-dimensional array, to a symmetric tridiagonal matrix */
/*     using orthogonal similarity transformations. */

/*     On Input */

/*        N is the order of the matrix A.  N is an INTEGER variable. */

/*        NV is an INTEGER variable set equal to the dimension of the */
/*          array A as specified in the calling program.  NV must not */
/*          be less than  N*(N+1)/2. */

/*        A contains the lower triangle, stored row-wise, of the real */
/*          symmetric packed matrix.  A is a one-dimensional REAL */
/*          array, dimensioned A(NV). */

/*     On Output */

/*        A contains information about the orthogonal transformations */
/*          used in the reduction in its first N*(N+1)/2 positions. */

/*        D contains the diagonal elements of the symmetric tridiagonal */
/*          matrix.  D is a one-dimensional REAL array, dimensioned D(N). */

/*        E contains the subdiagonal elements of the symmetric */
/*          tridiagonal matrix in its last N-1 positions.  E(1) is set */
/*          to zero.  E is a one-dimensional REAL array, dimensioned */
/*          E(N). */

/*        E2 contains the squares of the corresponding elements of E. */
/*          E2 may coincide with E if the squares are not needed. */
/*          E2 is a one-dimensional REAL array, dimensioned E2(N). */

/*     Questions and comments should be directed to B. S. Garbow, */
/*     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY */
/*     ------------------------------------------------------------------ */

/* ***REFERENCES  B. T. Smith, J. M. Boyle, J. J. Dongarra, B. S. Garbow, */
/*                 Y. Ikebe, V. C. Klema and C. B. Moler, Matrix Eigen- */
/*                 system Routines - EISPACK Guide, Springer-Verlag, */
/*                 1976. */
/* ***ROUTINES CALLED  (NONE) */
/* ***REVISION HISTORY  (YYMMDD) */
/*   760101  DATE WRITTEN */
/*   890831  Modified array declarations.  (WRB) */
/*   890831  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   920501  Reformatted the REFERENCES section.  (WRB) */
/* ***END PROLOGUE  TRED3 */


/*     .......... FOR I=N STEP -1 UNTIL 1 DO -- .......... */
/* ***FIRST EXECUTABLE STATEMENT  TRED3 */
    /* Parameter adjustments */
    --e2;
    --e;
    --d__;
    --a;

    /* Function Body */
    i__1 = *n;
    for (ii = 1; ii <= i__1; ++ii) {
	i__ = *n + 1 - ii;
	l = i__ - 1;
	iz = i__ * l / 2;
	h__ = 0.f;
	scale = 0.f;
	if (l < 1) {
	    goto L130;
	}
/*     .......... SCALE ROW (ALGOL TOL THEN NOT NEEDED) .......... */
	i__2 = l;
	for (k = 1; k <= i__2; ++k) {
	    ++iz;
	    d__[k] = a[iz];
	    scale += (r__1 = d__[k], dabs(r__1));
/* L120: */
	}

	if (scale != 0.f) {
	    goto L140;
	}
L130:
	e[i__] = 0.f;
	e2[i__] = 0.f;
	goto L290;

L140:
	i__2 = l;
	for (k = 1; k <= i__2; ++k) {
	    d__[k] /= scale;
	    h__ += d__[k] * d__[k];
/* L150: */
	}

	e2[i__] = scale * scale * h__;
	f = d__[l];
	r__1 = sqrt(h__);
	g = -r_sign(&r__1, &f);
	e[i__] = scale * g;
	h__ -= f * g;
	d__[l] = f - g;
	a[iz] = scale * d__[l];
	if (l == 1) {
	    goto L290;
	}
	f = 0.f;

	i__2 = l;
	for (j = 1; j <= i__2; ++j) {
	    g = 0.f;
	    jk = j * (j - 1) / 2;
/*     .......... FORM ELEMENT OF A*U .......... */
	    i__3 = l;
	    for (k = 1; k <= i__3; ++k) {
		++jk;
		if (k > j) {
		    jk = jk + k - 2;
		}
		g += a[jk] * d__[k];
/* L180: */
	    }
/*     .......... FORM ELEMENT OF P .......... */
	    e[j] = g / h__;
	    f += e[j] * d__[j];
/* L240: */
	}

	hh = f / (h__ + h__);
	jk = 0;
/*     .......... FORM REDUCED A .......... */
	i__2 = l;
	for (j = 1; j <= i__2; ++j) {
	    f = d__[j];
	    g = e[j] - hh * f;
	    e[j] = g;

	    i__3 = j;
	    for (k = 1; k <= i__3; ++k) {
		++jk;
		a[jk] = a[jk] - f * e[k] - g * d__[k];
/* L260: */
	    }
	}

L290:
	d__[i__] = a[iz + 1];
	a[iz + 1] = scale * sqrt(h__);
/* L300: */
    }

    return 0;
} /* tred3_ */

