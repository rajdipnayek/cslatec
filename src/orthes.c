/* orthes.f -- translated by f2c (version 12.02.01).
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

/* DECK ORTHES */
/* Subroutine */ int orthes_(integer *nm, integer *n, integer *low, integer *
	igh, real *a, real *ort)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3;
    real r__1;

    /* Local variables */
    static real f, g, h__;
    static integer i__, j, m, la, ii, jj, mp, kp1;
    static real scale;

/* ***BEGIN PROLOGUE  ORTHES */
/* ***PURPOSE  Reduce a real general matrix to upper Hessenberg form */
/*            using orthogonal similarity transformations. */
/* ***LIBRARY   SLATEC (EISPACK) */
/* ***CATEGORY  D4C1B2 */
/* ***TYPE      SINGLE PRECISION (ORTHES-S, CORTH-C) */
/* ***KEYWORDS  EIGENVALUES, EIGENVECTORS, EISPACK */
/* ***AUTHOR  Smith, B. T., et al. */
/* ***DESCRIPTION */

/*     This subroutine is a translation of the ALGOL procedure ORTHES, */
/*     NUM. MATH. 12, 349-368(1968) by Martin and Wilkinson. */
/*     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 339-358(1971). */

/*     Given a REAL GENERAL matrix, this subroutine */
/*     reduces a submatrix situated in rows and columns */
/*     LOW through IGH to upper Hessenberg form by */
/*     orthogonal similarity transformations. */

/*     On INPUT */

/*        NM must be set to the row dimension of the two-dimensional */
/*          array parameter, A, as declared in the calling program */
/*          dimension statement.  NM is an INTEGER variable. */

/*        N is the order of the matrix A.  N is an INTEGER variable. */
/*          N must be less than or equal to NM. */

/*        LOW and IGH are two INTEGER variables determined by the */
/*          balancing subroutine  BALANC.  If  BALANC  has not been */
/*          used, set LOW=1 and IGH equal to the order of the matrix, N. */

/*        A contains the general matrix to be reduced to upper */
/*          Hessenberg form.  A is a two-dimensional REAL array, */
/*          dimensioned A(NM,N). */

/*     On OUTPUT */

/*        A contains the upper Hessenberg matrix.  Some information about */
/*          the orthogonal transformations used in the reduction */
/*          is stored in the remaining triangle under the Hessenberg */
/*          matrix. */

/*        ORT contains further information about the orthogonal trans- */
/*          formations used in the reduction.  Only elements LOW+1 */
/*          through IGH are used.  ORT is a one-dimensional REAL array, */
/*          dimensioned ORT(IGH). */

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
/* ***END PROLOGUE  ORTHES */


/* ***FIRST EXECUTABLE STATEMENT  ORTHES */
    /* Parameter adjustments */
    a_dim1 = *nm;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --ort;

    /* Function Body */
    la = *igh - 1;
    kp1 = *low + 1;
    if (la < kp1) {
	goto L200;
    }

    i__1 = la;
    for (m = kp1; m <= i__1; ++m) {
	h__ = 0.f;
	ort[m] = 0.f;
	scale = 0.f;
/*     .......... SCALE COLUMN (ALGOL TOL THEN NOT NEEDED) .......... */
	i__2 = *igh;
	for (i__ = m; i__ <= i__2; ++i__) {
/* L90: */
	    scale += (r__1 = a[i__ + (m - 1) * a_dim1], dabs(r__1));
	}

	if (scale == 0.f) {
	    goto L180;
	}
	mp = m + *igh;
/*     .......... FOR I=IGH STEP -1 UNTIL M DO -- .......... */
	i__2 = *igh;
	for (ii = m; ii <= i__2; ++ii) {
	    i__ = mp - ii;
	    ort[i__] = a[i__ + (m - 1) * a_dim1] / scale;
	    h__ += ort[i__] * ort[i__];
/* L100: */
	}

	r__1 = sqrt(h__);
	g = -r_sign(&r__1, &ort[m]);
	h__ -= ort[m] * g;
	ort[m] -= g;
/*     .......... FORM (I-(U*UT)/H) * A .......... */
	i__2 = *n;
	for (j = m; j <= i__2; ++j) {
	    f = 0.f;
/*     .......... FOR I=IGH STEP -1 UNTIL M DO -- .......... */
	    i__3 = *igh;
	    for (ii = m; ii <= i__3; ++ii) {
		i__ = mp - ii;
		f += ort[i__] * a[i__ + j * a_dim1];
/* L110: */
	    }

	    f /= h__;

	    i__3 = *igh;
	    for (i__ = m; i__ <= i__3; ++i__) {
/* L120: */
		a[i__ + j * a_dim1] -= f * ort[i__];
	    }

/* L130: */
	}
/*     .......... FORM (I-(U*UT)/H)*A*(I-(U*UT)/H) .......... */
	i__2 = *igh;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    f = 0.f;
/*     .......... FOR J=IGH STEP -1 UNTIL M DO -- .......... */
	    i__3 = *igh;
	    for (jj = m; jj <= i__3; ++jj) {
		j = mp - jj;
		f += ort[j] * a[i__ + j * a_dim1];
/* L140: */
	    }

	    f /= h__;

	    i__3 = *igh;
	    for (j = m; j <= i__3; ++j) {
/* L150: */
		a[i__ + j * a_dim1] -= f * ort[j];
	    }

/* L160: */
	}

	ort[m] = scale * ort[m];
	a[m + (m - 1) * a_dim1] = scale * g;
L180:
	;
    }

L200:
    return 0;
} /* orthes_ */

