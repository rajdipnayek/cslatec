/* imtql2.f -- translated by f2c (version 12.02.01).
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

static real c_b9 = 1.f;

/* DECK IMTQL2 */
/* Subroutine */ int imtql2_(integer *nm, integer *n, real *d__, real *e, 
	real *z__, integer *ierr)
{
    /* System generated locals */
    integer z_dim1, z_offset, i__1, i__2, i__3;
    real r__1, r__2;

    /* Local variables */
    static real b, c__, f, g;
    static integer i__, j, k, l, m;
    static real p, r__, s, s1, s2;
    static integer ii, mml;
    extern doublereal pythag_(real *, real *);

/* ***BEGIN PROLOGUE  IMTQL2 */
/* ***PURPOSE  Compute the eigenvalues and eigenvectors of a symmetric */
/*            tridiagonal matrix using the implicit QL method. */
/* ***LIBRARY   SLATEC (EISPACK) */
/* ***CATEGORY  D4A5, D4C2A */
/* ***TYPE      SINGLE PRECISION (IMTQL2-S) */
/* ***KEYWORDS  EIGENVALUES, EIGENVECTORS, EISPACK */
/* ***AUTHOR  Smith, B. T., et al. */
/* ***DESCRIPTION */

/*     This subroutine is a translation of the ALGOL procedure IMTQL2, */
/*     NUM. MATH. 12, 377-383(1968) by Martin and Wilkinson, */
/*     as modified in NUM. MATH. 15, 450(1970) by Dubrulle. */
/*     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 241-248(1971). */

/*     This subroutine finds the eigenvalues and eigenvectors */
/*     of a SYMMETRIC TRIDIAGONAL matrix by the implicit QL method. */
/*     The eigenvectors of a FULL SYMMETRIC matrix can also */
/*     be found if  TRED2  has been used to reduce this */
/*     full matrix to tridiagonal form. */

/*     On INPUT */

/*        NM must be set to the row dimension of the two-dimensional */
/*          array parameter, Z, as declared in the calling program */
/*          dimension statement.  NM is an INTEGER variable. */

/*        N is the order of the matrix.  N is an INTEGER variable. */
/*          N must be less than or equal to NM. */

/*        D contains the diagonal elements of the symmetric tridiagonal */
/*          matrix.  D is a one-dimensional REAL array, dimensioned D(N). */

/*        E contains the subdiagonal elements of the symmetric */
/*          tridiagonal matrix in its last N-1 positions.  E(1) is */
/*          arbitrary.  E is a one-dimensional REAL array, dimensioned */
/*          E(N). */

/*        Z contains the transformation matrix produced in the reduction */
/*          by  TRED2,  if performed.  This transformation matrix is */
/*          necessary if you want to obtain the eigenvectors of the full */
/*          symmetric matrix.  If the eigenvectors of the symmetric */
/*          tridiagonal matrix are desired, Z must contain the identity */
/*          matrix.  Z is a two-dimensional REAL array, dimensioned */
/*          Z(NM,N). */

/*      On OUTPUT */

/*        D contains the eigenvalues in ascending order.  If an */
/*          error exit is made, the eigenvalues are correct but */
/*          unordered for indices 1, 2, ..., IERR-1. */

/*        E has been destroyed. */

/*        Z contains orthonormal eigenvectors of the full symmetric */
/*          or symmetric tridiagonal matrix, depending on what it */
/*          contained on input.  If an error exit is made,  Z contains */
/*          the eigenvectors associated with the stored eigenvalues. */

/*        IERR is an INTEGER flag set to */
/*          Zero       for normal return, */
/*          J          if the J-th eigenvalue has not been */
/*                     determined after 30 iterations. */
/*                     The eigenvalues and eigenvectors should be correct */
/*                     for indices 1, 2, ..., IERR-1, but the eigenvalues */
/*                     are not ordered. */

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
/* ***END PROLOGUE  IMTQL2 */


/* ***FIRST EXECUTABLE STATEMENT  IMTQL2 */
    /* Parameter adjustments */
    z_dim1 = *nm;
    z_offset = 1 + z_dim1;
    z__ -= z_offset;
    --d__;
    --e;

    /* Function Body */
    *ierr = 0;
    if (*n == 1) {
	goto L1001;
    }

    i__1 = *n;
    for (i__ = 2; i__ <= i__1; ++i__) {
/* L100: */
	e[i__ - 1] = e[i__];
    }

    e[*n] = 0.f;

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
	    s1 = (r__1 = d__[m], dabs(r__1)) + (r__2 = d__[m + 1], dabs(r__2))
		    ;
	    s2 = s1 + (r__1 = e[m], dabs(r__1));
	    if (s2 == s1) {
		goto L120;
	    }
/* L110: */
	}

L120:
	p = d__[l];
	if (m == l) {
	    goto L240;
	}
	if (j == 30) {
	    goto L1000;
	}
	++j;
/*     .......... FORM SHIFT .......... */
	g = (d__[l + 1] - p) / (e[l] * 2.f);
	r__ = pythag_(&g, &c_b9);
	g = d__[m] - p + e[l] / (g + r_sign(&r__, &g));
	s = 1.f;
	c__ = 1.f;
	p = 0.f;
	mml = m - l;
/*     .......... FOR I=M-1 STEP -1 UNTIL L DO -- .......... */
	i__2 = mml;
	for (ii = 1; ii <= i__2; ++ii) {
	    i__ = m - ii;
	    f = s * e[i__];
	    b = c__ * e[i__];
	    if (dabs(f) < dabs(g)) {
		goto L150;
	    }
	    c__ = g / f;
	    r__ = sqrt(c__ * c__ + 1.f);
	    e[i__ + 1] = f * r__;
	    s = 1.f / r__;
	    c__ *= s;
	    goto L160;
L150:
	    s = f / g;
	    r__ = sqrt(s * s + 1.f);
	    e[i__ + 1] = g * r__;
	    c__ = 1.f / r__;
	    s *= c__;
L160:
	    g = d__[i__ + 1] - p;
	    r__ = (d__[i__] - g) * s + c__ * 2.f * b;
	    p = s * r__;
	    d__[i__ + 1] = g + p;
	    g = c__ * r__ - b;
/*     .......... FORM VECTOR .......... */
	    i__3 = *n;
	    for (k = 1; k <= i__3; ++k) {
		f = z__[k + (i__ + 1) * z_dim1];
		z__[k + (i__ + 1) * z_dim1] = s * z__[k + i__ * z_dim1] + c__ 
			* f;
		z__[k + i__ * z_dim1] = c__ * z__[k + i__ * z_dim1] - s * f;
/* L180: */
	    }

/* L200: */
	}

	d__[l] -= p;
	e[l] = g;
	e[m] = 0.f;
	goto L105;
L240:
	;
    }
/*     .......... ORDER EIGENVALUES AND EIGENVECTORS .......... */
    i__1 = *n;
    for (ii = 2; ii <= i__1; ++ii) {
	i__ = ii - 1;
	k = i__;
	p = d__[i__];

	i__2 = *n;
	for (j = ii; j <= i__2; ++j) {
	    if (d__[j] >= p) {
		goto L260;
	    }
	    k = j;
	    p = d__[j];
L260:
	    ;
	}

	if (k == i__) {
	    goto L300;
	}
	d__[k] = d__[i__];
	d__[i__] = p;

	i__2 = *n;
	for (j = 1; j <= i__2; ++j) {
	    p = z__[j + i__ * z_dim1];
	    z__[j + i__ * z_dim1] = z__[j + k * z_dim1];
	    z__[j + k * z_dim1] = p;
/* L280: */
	}

L300:
	;
    }

    goto L1001;
/*     .......... SET ERROR -- NO CONVERGENCE TO AN */
/*                EIGENVALUE AFTER 30 ITERATIONS .......... */
L1000:
    *ierr = l;
L1001:
    return 0;
} /* imtql2_ */

