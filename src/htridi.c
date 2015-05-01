/* htridi.f -- translated by f2c (version 12.02.01).
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

/* DECK HTRIDI */
/* Subroutine */ int htridi_(integer *nm, integer *n, real *ar, real *ai, 
	real *d__, real *e, real *e2, real *tau)
{
    /* System generated locals */
    integer ar_dim1, ar_offset, ai_dim1, ai_offset, i__1, i__2, i__3;
    real r__1, r__2;

    /* Local variables */
    static real f, g, h__;
    static integer i__, j, k, l;
    static real fi, gi, hh;
    static integer ii;
    static real si;
    static integer jp1;
    static real scale;
    extern doublereal pythag_(real *, real *);

/* ***BEGIN PROLOGUE  HTRIDI */
/* ***PURPOSE  Reduce a complex Hermitian matrix to a real symmetric */
/*            tridiagonal matrix using unitary similarity */
/*            transformations. */
/* ***LIBRARY   SLATEC (EISPACK) */
/* ***CATEGORY  D4C1B1 */
/* ***TYPE      SINGLE PRECISION (HTRIDI-S) */
/* ***KEYWORDS  EIGENVALUES, EIGENVECTORS, EISPACK */
/* ***AUTHOR  Smith, B. T., et al. */
/* ***DESCRIPTION */

/*     This subroutine is a translation of a complex analogue of */
/*     the ALGOL procedure TRED1, NUM. MATH. 11, 181-195(1968) */
/*     by Martin, Reinsch, and Wilkinson. */
/*     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 212-226(1971). */

/*     This subroutine reduces a COMPLEX HERMITIAN matrix */
/*     to a real symmetric tridiagonal matrix using */
/*     unitary similarity transformations. */

/*     On INPUT */

/*        NM must be set to the row dimension of the two-dimensional */
/*          array parameters, AR and AI, as declared in the calling */
/*          program dimension statement.  NM is an INTEGER variable. */

/*        N is the order of the matrix A=(AR,AI).  N is an INTEGER */
/*          variable. N must be less than or equal to NM. */

/*        AR and AI contain the real and imaginary parts, respectively, */
/*          of the complex Hermitian input matrix.  Only the lower */
/*          triangle of the matrix need be supplied.  AR and AI are two- */
/*          dimensional REAL arrays, dimensioned AR(NM,N) and AI(NM,N). */

/*     On OUTPUT */

/*        AR and AI contain some information about the unitary trans- */
/*          formations used in the reduction in the strict lower triangle */
/*          of AR and the full lower triangle of AI.  The rest of the */
/*          matrices are unaltered. */

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
/* ***END PROLOGUE  HTRIDI */


/* ***FIRST EXECUTABLE STATEMENT  HTRIDI */
    /* Parameter adjustments */
    ai_dim1 = *nm;
    ai_offset = 1 + ai_dim1;
    ai -= ai_offset;
    ar_dim1 = *nm;
    ar_offset = 1 + ar_dim1;
    ar -= ar_offset;
    --d__;
    --e;
    --e2;
    tau -= 3;

    /* Function Body */
    tau[(*n << 1) + 1] = 1.f;
    tau[(*n << 1) + 2] = 0.f;

    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L100: */
	d__[i__] = ar[i__ + i__ * ar_dim1];
    }
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
	    scale = scale + (r__1 = ar[i__ + k * ar_dim1], dabs(r__1)) + (
		    r__2 = ai[i__ + k * ai_dim1], dabs(r__2));
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
	    ar[i__ + k * ar_dim1] /= scale;
	    ai[i__ + k * ai_dim1] /= scale;
	    h__ = h__ + ar[i__ + k * ar_dim1] * ar[i__ + k * ar_dim1] + ai[
		    i__ + k * ai_dim1] * ai[i__ + k * ai_dim1];
/* L150: */
	}

	e2[i__] = scale * scale * h__;
	g = sqrt(h__);
	e[i__] = scale * g;
	f = pythag_(&ar[i__ + l * ar_dim1], &ai[i__ + l * ai_dim1]);
/*     .......... FORM NEXT DIAGONAL ELEMENT OF MATRIX T .......... */
	if (f == 0.f) {
	    goto L160;
	}
	tau[(l << 1) + 1] = (ai[i__ + l * ai_dim1] * tau[(i__ << 1) + 2] - ar[
		i__ + l * ar_dim1] * tau[(i__ << 1) + 1]) / f;
	si = (ar[i__ + l * ar_dim1] * tau[(i__ << 1) + 2] + ai[i__ + l * 
		ai_dim1] * tau[(i__ << 1) + 1]) / f;
	h__ += f * g;
	g = g / f + 1.f;
	ar[i__ + l * ar_dim1] = g * ar[i__ + l * ar_dim1];
	ai[i__ + l * ai_dim1] = g * ai[i__ + l * ai_dim1];
	if (l == 1) {
	    goto L270;
	}
	goto L170;
L160:
	tau[(l << 1) + 1] = -tau[(i__ << 1) + 1];
	si = tau[(i__ << 1) + 2];
	ar[i__ + l * ar_dim1] = g;
L170:
	f = 0.f;

	i__2 = l;
	for (j = 1; j <= i__2; ++j) {
	    g = 0.f;
	    gi = 0.f;
/*     .......... FORM ELEMENT OF A*U .......... */
	    i__3 = j;
	    for (k = 1; k <= i__3; ++k) {
		g = g + ar[j + k * ar_dim1] * ar[i__ + k * ar_dim1] + ai[j + 
			k * ai_dim1] * ai[i__ + k * ai_dim1];
		gi = gi - ar[j + k * ar_dim1] * ai[i__ + k * ai_dim1] + ai[j 
			+ k * ai_dim1] * ar[i__ + k * ar_dim1];
/* L180: */
	    }

	    jp1 = j + 1;
	    if (l < jp1) {
		goto L220;
	    }

	    i__3 = l;
	    for (k = jp1; k <= i__3; ++k) {
		g = g + ar[k + j * ar_dim1] * ar[i__ + k * ar_dim1] - ai[k + 
			j * ai_dim1] * ai[i__ + k * ai_dim1];
		gi = gi - ar[k + j * ar_dim1] * ai[i__ + k * ai_dim1] - ai[k 
			+ j * ai_dim1] * ar[i__ + k * ar_dim1];
/* L200: */
	    }
/*     .......... FORM ELEMENT OF P .......... */
L220:
	    e[j] = g / h__;
	    tau[(j << 1) + 2] = gi / h__;
	    f = f + e[j] * ar[i__ + j * ar_dim1] - tau[(j << 1) + 2] * ai[i__ 
		    + j * ai_dim1];
/* L240: */
	}

	hh = f / (h__ + h__);
/*     .......... FORM REDUCED A .......... */
	i__2 = l;
	for (j = 1; j <= i__2; ++j) {
	    f = ar[i__ + j * ar_dim1];
	    g = e[j] - hh * f;
	    e[j] = g;
	    fi = -ai[i__ + j * ai_dim1];
	    gi = tau[(j << 1) + 2] - hh * fi;
	    tau[(j << 1) + 2] = -gi;

	    i__3 = j;
	    for (k = 1; k <= i__3; ++k) {
		ar[j + k * ar_dim1] = ar[j + k * ar_dim1] - f * e[k] - g * ar[
			i__ + k * ar_dim1] + fi * tau[(k << 1) + 2] + gi * ai[
			i__ + k * ai_dim1];
		ai[j + k * ai_dim1] = ai[j + k * ai_dim1] - f * tau[(k << 1) 
			+ 2] - g * ai[i__ + k * ai_dim1] - fi * e[k] - gi * 
			ar[i__ + k * ar_dim1];
/* L260: */
	    }
	}

L270:
	i__3 = l;
	for (k = 1; k <= i__3; ++k) {
	    ar[i__ + k * ar_dim1] = scale * ar[i__ + k * ar_dim1];
	    ai[i__ + k * ai_dim1] = scale * ai[i__ + k * ai_dim1];
/* L280: */
	}

	tau[(l << 1) + 2] = -si;
L290:
	hh = d__[i__];
	d__[i__] = ar[i__ + i__ * ar_dim1];
	ar[i__ + i__ * ar_dim1] = hh;
	ai[i__ + i__ * ai_dim1] = scale * sqrt(h__);
/* L300: */
    }

    return 0;
} /* htridi_ */

