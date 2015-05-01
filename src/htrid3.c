/* htrid3.f -- translated by f2c (version 12.02.01).
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

/* DECK HTRID3 */
/* Subroutine */ int htrid3_(integer *nm, integer *n, real *a, real *d__, 
	real *e, real *e2, real *tau)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3;
    real r__1, r__2;

    /* Local variables */
    static real f, g, h__;
    static integer i__, j, k, l;
    static real fi, gi, hh;
    static integer ii;
    static real si;
    static integer jm1, jp1;
    static real scale;
    extern doublereal pythag_(real *, real *);

/* ***BEGIN PROLOGUE  HTRID3 */
/* ***PURPOSE  Reduce a complex Hermitian (packed) matrix to a real */
/*            symmetric tridiagonal matrix by unitary similarity */
/*            transformations. */
/* ***LIBRARY   SLATEC (EISPACK) */
/* ***CATEGORY  D4C1B1 */
/* ***TYPE      SINGLE PRECISION (HTRID3-S) */
/* ***KEYWORDS  EIGENVALUES, EIGENVECTORS, EISPACK */
/* ***AUTHOR  Smith, B. T., et al. */
/* ***DESCRIPTION */

/*     This subroutine is a translation of a complex analogue of */
/*     the ALGOL procedure TRED3, NUM. MATH. 11, 181-195(1968) */
/*     by Martin, Reinsch, and Wilkinson. */
/*     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 212-226(1971). */

/*     This subroutine reduces a COMPLEX HERMITIAN matrix, stored as */
/*     a single square array, to a real symmetric tridiagonal matrix */
/*     using unitary similarity transformations. */

/*     On INPUT */

/*        NM must be set to the row dimension of the two-dimensional */
/*          array parameter, A, as declared in the calling program */
/*          dimension statement.  NM is an INTEGER variable. */

/*        N is the order of the matrix.  N is an INTEGER variable. */
/*          N must be less than or equal to NM. */

/*        A contains the lower triangle of the complex Hermitian input */
/*          matrix.  The real parts of the matrix elements are stored */
/*          in the full lower triangle of A, and the imaginary parts */
/*          are stored in the transposed positions of the strict upper */
/*          triangle of A.  No storage is required for the zero */
/*          imaginary parts of the diagonal elements.  A is a two- */
/*          dimensional REAL array, dimensioned A(NM,N). */

/*     On OUTPUT */

/*        A contains some information about the unitary transformations */
/*          used in the reduction. */

/*        D contains the diagonal elements of the real symmetric */
/*          tridiagonal matrix.  D is a one-dimensional REAL array, */
/*          dimensioned D(N). */

/*        E contains the subdiagonal elements of the real tridiagonal */
/*          matrix in its last N-1 positions.  E(1) is set to zero. */
/*          E is a one-dimensional REAL array, dimensioned E(N). */

/*        E2 contains the squares of the corresponding elements of E. */
/*          E2(1) is set to zero.  E2 may coincide with E if the squares */
/*          are not needed.  E2 is a one-dimensional REAL array, */
/*          dimensioned E2(N). */

/*        TAU contains further information about the transformations. */
/*          TAU is a one-dimensional REAL array, dimensioned TAU(2,N). */

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
/* ***END PROLOGUE  HTRID3 */


/* ***FIRST EXECUTABLE STATEMENT  HTRID3 */
    /* Parameter adjustments */
    a_dim1 = *nm;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --d__;
    --e;
    --e2;
    tau -= 3;

    /* Function Body */
    tau[(*n << 1) + 1] = 1.f;
    tau[(*n << 1) + 2] = 0.f;
/*     .......... FOR I=N STEP -1 UNTIL 1 DO -- .......... */
    i__1 = *n;
    for (ii = 1; ii <= i__1; ++ii) {
	i__ = *n + 1 - ii;
	l = i__ - 1;
	h__ = 0.f;
	scale = 0.f;
	if (l < 1) {
	    goto L130;
	}
/*     .......... SCALE ROW (ALGOL TOL THEN NOT NEEDED) .......... */
	i__2 = l;
	for (k = 1; k <= i__2; ++k) {
/* L120: */
	    scale = scale + (r__1 = a[i__ + k * a_dim1], dabs(r__1)) + (r__2 =
		     a[k + i__ * a_dim1], dabs(r__2));
	}

	if (scale != 0.f) {
	    goto L140;
	}
	tau[(l << 1) + 1] = 1.f;
	tau[(l << 1) + 2] = 0.f;
L130:
	e[i__] = 0.f;
	e2[i__] = 0.f;
	goto L290;

L140:
	i__2 = l;
	for (k = 1; k <= i__2; ++k) {
	    a[i__ + k * a_dim1] /= scale;
	    a[k + i__ * a_dim1] /= scale;
	    h__ = h__ + a[i__ + k * a_dim1] * a[i__ + k * a_dim1] + a[k + i__ 
		    * a_dim1] * a[k + i__ * a_dim1];
/* L150: */
	}

	e2[i__] = scale * scale * h__;
	g = sqrt(h__);
	e[i__] = scale * g;
	f = pythag_(&a[i__ + l * a_dim1], &a[l + i__ * a_dim1]);
/*     .......... FORM NEXT DIAGONAL ELEMENT OF MATRIX T .......... */
	if (f == 0.f) {
	    goto L160;
	}
	tau[(l << 1) + 1] = (a[l + i__ * a_dim1] * tau[(i__ << 1) + 2] - a[
		i__ + l * a_dim1] * tau[(i__ << 1) + 1]) / f;
	si = (a[i__ + l * a_dim1] * tau[(i__ << 1) + 2] + a[l + i__ * a_dim1] 
		* tau[(i__ << 1) + 1]) / f;
	h__ += f * g;
	g = g / f + 1.f;
	a[i__ + l * a_dim1] = g * a[i__ + l * a_dim1];
	a[l + i__ * a_dim1] = g * a[l + i__ * a_dim1];
	if (l == 1) {
	    goto L270;
	}
	goto L170;
L160:
	tau[(l << 1) + 1] = -tau[(i__ << 1) + 1];
	si = tau[(i__ << 1) + 2];
	a[i__ + l * a_dim1] = g;
L170:
	f = 0.f;

	i__2 = l;
	for (j = 1; j <= i__2; ++j) {
	    g = 0.f;
	    gi = 0.f;
	    if (j == 1) {
		goto L190;
	    }
	    jm1 = j - 1;
/*     .......... FORM ELEMENT OF A*U .......... */
	    i__3 = jm1;
	    for (k = 1; k <= i__3; ++k) {
		g = g + a[j + k * a_dim1] * a[i__ + k * a_dim1] + a[k + j * 
			a_dim1] * a[k + i__ * a_dim1];
		gi = gi - a[j + k * a_dim1] * a[k + i__ * a_dim1] + a[k + j * 
			a_dim1] * a[i__ + k * a_dim1];
/* L180: */
	    }

L190:
	    g += a[j + j * a_dim1] * a[i__ + j * a_dim1];
	    gi -= a[j + j * a_dim1] * a[j + i__ * a_dim1];
	    jp1 = j + 1;
	    if (l < jp1) {
		goto L220;
	    }

	    i__3 = l;
	    for (k = jp1; k <= i__3; ++k) {
		g = g + a[k + j * a_dim1] * a[i__ + k * a_dim1] - a[j + k * 
			a_dim1] * a[k + i__ * a_dim1];
		gi = gi - a[k + j * a_dim1] * a[k + i__ * a_dim1] - a[j + k * 
			a_dim1] * a[i__ + k * a_dim1];
/* L200: */
	    }
/*     .......... FORM ELEMENT OF P .......... */
L220:
	    e[j] = g / h__;
	    tau[(j << 1) + 2] = gi / h__;
	    f = f + e[j] * a[i__ + j * a_dim1] - tau[(j << 1) + 2] * a[j + 
		    i__ * a_dim1];
/* L240: */
	}

	hh = f / (h__ + h__);
/*     .......... FORM REDUCED A .......... */
	i__2 = l;
	for (j = 1; j <= i__2; ++j) {
	    f = a[i__ + j * a_dim1];
	    g = e[j] - hh * f;
	    e[j] = g;
	    fi = -a[j + i__ * a_dim1];
	    gi = tau[(j << 1) + 2] - hh * fi;
	    tau[(j << 1) + 2] = -gi;
	    a[j + j * a_dim1] -= (f * g + fi * gi) * 2.f;
	    if (j == 1) {
		goto L260;
	    }
	    jm1 = j - 1;

	    i__3 = jm1;
	    for (k = 1; k <= i__3; ++k) {
		a[j + k * a_dim1] = a[j + k * a_dim1] - f * e[k] - g * a[i__ 
			+ k * a_dim1] + fi * tau[(k << 1) + 2] + gi * a[k + 
			i__ * a_dim1];
		a[k + j * a_dim1] = a[k + j * a_dim1] - f * tau[(k << 1) + 2] 
			- g * a[k + i__ * a_dim1] - fi * e[k] - gi * a[i__ + 
			k * a_dim1];
/* L250: */
	    }

L260:
	    ;
	}

L270:
	i__2 = l;
	for (k = 1; k <= i__2; ++k) {
	    a[i__ + k * a_dim1] = scale * a[i__ + k * a_dim1];
	    a[k + i__ * a_dim1] = scale * a[k + i__ * a_dim1];
/* L280: */
	}

	tau[(l << 1) + 2] = -si;
L290:
	d__[i__] = a[i__ + i__ * a_dim1];
	a[i__ + i__ * a_dim1] = scale * sqrt(h__);
/* L300: */
    }

    return 0;
} /* htrid3_ */

