/* corth.f -- translated by f2c (version 12.02.01).
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

/* DECK CORTH */
/* Subroutine */ int corth_(integer *nm, integer *n, integer *low, integer *
	igh, real *ar, real *ai, real *ortr, real *orti)
{
    /* System generated locals */
    integer ar_dim1, ar_offset, ai_dim1, ai_offset, i__1, i__2, i__3;
    real r__1, r__2;

    /* Local variables */
    static real f, g, h__;
    static integer i__, j, m, la;
    static real fi;
    static integer ii, jj;
    static real fr;
    static integer mp, kp1;
    static real scale;
    extern doublereal pythag_(real *, real *);

/* ***BEGIN PROLOGUE  CORTH */
/* ***PURPOSE  Reduce a complex general matrix to complex upper Hessenberg */
/*            form using unitary similarity transformations. */
/* ***LIBRARY   SLATEC (EISPACK) */
/* ***CATEGORY  D4C1B2 */
/* ***TYPE      COMPLEX (ORTHES-S, CORTH-C) */
/* ***KEYWORDS  EIGENVALUES, EIGENVECTORS, EISPACK */
/* ***AUTHOR  Smith, B. T., et al. */
/* ***DESCRIPTION */

/*     This subroutine is a translation of a complex analogue of */
/*     the ALGOL procedure ORTHES, NUM. MATH. 12, 349-368(1968) */
/*     by Martin and Wilkinson. */
/*     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 339-358(1971). */

/*     Given a COMPLEX GENERAL matrix, this subroutine */
/*     reduces a submatrix situated in rows and columns */
/*     LOW through IGH to upper Hessenberg form by */
/*     unitary similarity transformations. */

/*     On INPUT */

/*        NM must be set to the row dimension of the two-dimensional */
/*          array parameters, AR and AI, as declared in the calling */
/*          program dimension statement.  NM is an INTEGER variable. */

/*        N is the order of the matrix A=(AR,AI).  N is an INTEGER */
/*          variable.  N must be less than or equal to NM. */

/*        LOW and IGH are two INTEGER variables determined by the */
/*          balancing subroutine  CBAL.  If  CBAL  has not been used, */
/*          set LOW=1 and IGH equal to the order of the matrix, N. */

/*        AR and AI contain the real and imaginary parts, respectively, */
/*          of the complex input matrix.  AR and AI are two-dimensional */
/*          REAL arrays, dimensioned AR(NM,N) and AI(NM,N). */

/*     On OUTPUT */

/*        AR and AI contain the real and imaginary parts, respectively, */
/*          of the Hessenberg matrix.  Information about the unitary */
/*          transformations used in the reduction is stored in the */
/*          remaining triangles under the Hessenberg matrix. */

/*        ORTR and ORTI contain further information about the unitary */
/*          transformations.  Only elements LOW through IGH are used. */
/*          ORTR and ORTI are one-dimensional REAL arrays, dimensioned */
/*          ORTR(IGH) and ORTI(IGH). */

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
/* ***END PROLOGUE  CORTH */


/* ***FIRST EXECUTABLE STATEMENT  CORTH */
    /* Parameter adjustments */
    ai_dim1 = *nm;
    ai_offset = 1 + ai_dim1;
    ai -= ai_offset;
    ar_dim1 = *nm;
    ar_offset = 1 + ar_dim1;
    ar -= ar_offset;
    --ortr;
    --orti;

    /* Function Body */
    la = *igh - 1;
    kp1 = *low + 1;
    if (la < kp1) {
	goto L200;
    }

    i__1 = la;
    for (m = kp1; m <= i__1; ++m) {
	h__ = 0.f;
	ortr[m] = 0.f;
	orti[m] = 0.f;
	scale = 0.f;
/*     .......... SCALE COLUMN (ALGOL TOL THEN NOT NEEDED) .......... */
	i__2 = *igh;
	for (i__ = m; i__ <= i__2; ++i__) {
/* L90: */
	    scale = scale + (r__1 = ar[i__ + (m - 1) * ar_dim1], dabs(r__1)) 
		    + (r__2 = ai[i__ + (m - 1) * ai_dim1], dabs(r__2));
	}

	if (scale == 0.f) {
	    goto L180;
	}
	mp = m + *igh;
/*     .......... FOR I=IGH STEP -1 UNTIL M DO -- .......... */
	i__2 = *igh;
	for (ii = m; ii <= i__2; ++ii) {
	    i__ = mp - ii;
	    ortr[i__] = ar[i__ + (m - 1) * ar_dim1] / scale;
	    orti[i__] = ai[i__ + (m - 1) * ai_dim1] / scale;
	    h__ = h__ + ortr[i__] * ortr[i__] + orti[i__] * orti[i__];
/* L100: */
	}

	g = sqrt(h__);
	f = pythag_(&ortr[m], &orti[m]);
	if (f == 0.f) {
	    goto L103;
	}
	h__ += f * g;
	g /= f;
	ortr[m] = (g + 1.f) * ortr[m];
	orti[m] = (g + 1.f) * orti[m];
	goto L105;

L103:
	ortr[m] = g;
	ar[m + (m - 1) * ar_dim1] = scale;
/*     .......... FORM (I-(U*UT)/H) * A .......... */
L105:
	i__2 = *n;
	for (j = m; j <= i__2; ++j) {
	    fr = 0.f;
	    fi = 0.f;
/*     .......... FOR I=IGH STEP -1 UNTIL M DO -- .......... */
	    i__3 = *igh;
	    for (ii = m; ii <= i__3; ++ii) {
		i__ = mp - ii;
		fr = fr + ortr[i__] * ar[i__ + j * ar_dim1] + orti[i__] * ai[
			i__ + j * ai_dim1];
		fi = fi + ortr[i__] * ai[i__ + j * ai_dim1] - orti[i__] * ar[
			i__ + j * ar_dim1];
/* L110: */
	    }

	    fr /= h__;
	    fi /= h__;

	    i__3 = *igh;
	    for (i__ = m; i__ <= i__3; ++i__) {
		ar[i__ + j * ar_dim1] = ar[i__ + j * ar_dim1] - fr * ortr[i__]
			 + fi * orti[i__];
		ai[i__ + j * ai_dim1] = ai[i__ + j * ai_dim1] - fr * orti[i__]
			 - fi * ortr[i__];
/* L120: */
	    }

/* L130: */
	}
/*     .......... FORM (I-(U*UT)/H)*A*(I-(U*UT)/H) .......... */
	i__2 = *igh;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    fr = 0.f;
	    fi = 0.f;
/*     .......... FOR J=IGH STEP -1 UNTIL M DO -- .......... */
	    i__3 = *igh;
	    for (jj = m; jj <= i__3; ++jj) {
		j = mp - jj;
		fr = fr + ortr[j] * ar[i__ + j * ar_dim1] - orti[j] * ai[i__ 
			+ j * ai_dim1];
		fi = fi + ortr[j] * ai[i__ + j * ai_dim1] + orti[j] * ar[i__ 
			+ j * ar_dim1];
/* L140: */
	    }

	    fr /= h__;
	    fi /= h__;

	    i__3 = *igh;
	    for (j = m; j <= i__3; ++j) {
		ar[i__ + j * ar_dim1] = ar[i__ + j * ar_dim1] - fr * ortr[j] 
			- fi * orti[j];
		ai[i__ + j * ai_dim1] = ai[i__ + j * ai_dim1] + fr * orti[j] 
			- fi * ortr[j];
/* L150: */
	    }

/* L160: */
	}

	ortr[m] = scale * ortr[m];
	orti[m] = scale * orti[m];
	ar[m + (m - 1) * ar_dim1] = -g * ar[m + (m - 1) * ar_dim1];
	ai[m + (m - 1) * ai_dim1] = -g * ai[m + (m - 1) * ai_dim1];
L180:
	;
    }

L200:
    return 0;
} /* corth_ */

