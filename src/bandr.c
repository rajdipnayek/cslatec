/* bandr.f -- translated by f2c (version 12.02.01).
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

/* DECK BANDR */
/* Subroutine */ int bandr_(integer *nm, integer *n, integer *mb, real *a, 
	real *d__, real *e, real *e2, logical *matz, real *z__)
{
    /* System generated locals */
    integer a_dim1, a_offset, z_dim1, z_offset, i__1, i__2, i__3, i__4, i__5, 
	    i__6;
    real r__1;

    /* Local variables */
    static real g;
    static integer j, k, l, r__;
    static real u, b1, b2, c2, f1, f2;
    static integer i1, i2, j1, j2, m1, n2, r1;
    static real s2;
    static integer kr, mr, ugl;
    static real dmin__;
    static integer maxl, maxr;
    static real dminrt;

/* ***BEGIN PROLOGUE  BANDR */
/* ***PURPOSE  Reduce a real symmetric band matrix to symmetric */
/*            tridiagonal matrix and, optionally, accumulate */
/*            orthogonal similarity transformations. */
/* ***LIBRARY   SLATEC (EISPACK) */
/* ***CATEGORY  D4C1B1 */
/* ***TYPE      SINGLE PRECISION (BANDR-S) */
/* ***KEYWORDS  EIGENVALUES, EIGENVECTORS, EISPACK */
/* ***AUTHOR  Smith, B. T., et al. */
/* ***DESCRIPTION */

/*     This subroutine is a translation of the ALGOL procedure BANDRD, */
/*     NUM. MATH. 12, 231-241(1968) by Schwarz. */
/*     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 273-283(1971). */

/*     This subroutine reduces a REAL SYMMETRIC BAND matrix */
/*     to a symmetric tridiagonal matrix using and optionally */
/*     accumulating orthogonal similarity transformations. */

/*     On INPUT */

/*        NM must be set to the row dimension of the two-dimensional */
/*          array parameters, A and Z, as declared in the calling */
/*          program dimension statement.  NM is an INTEGER variable. */

/*        N is the order of the matrix A.  N is an INTEGER variable. */
/*          N must be less than or equal to NM. */

/*        MB is the (half) band width of the matrix, defined as the */
/*          number of adjacent diagonals, including the principal */
/*          diagonal, required to specify the non-zero portion of the */
/*          lower triangle of the matrix.  MB is less than or equal */
/*          to N.  MB is an INTEGER variable. */

/*        A contains the lower triangle of the real symmetric band */
/*          matrix.  Its lowest subdiagonal is stored in the last */
/*          N+1-MB  positions of the first column, its next subdiagonal */
/*          in the last  N+2-MB  positions of the second column, further */
/*          subdiagonals similarly, and finally its principal diagonal */
/*          in the  N  positions of the last column.  Contents of storage */
/*          locations not part of the matrix are arbitrary.  A is a */
/*          two-dimensional REAL array, dimensioned A(NM,MB). */

/*        MATZ should be set to .TRUE. if the transformation matrix is */
/*          to be accumulated, and to .FALSE. otherwise.  MATZ is a */
/*          LOGICAL variable. */

/*     On OUTPUT */

/*        A has been destroyed, except for its last two columns which */
/*          contain a copy of the tridiagonal matrix. */

/*        D contains the diagonal elements of the tridiagonal matrix. */
/*          D is a one-dimensional REAL array, dimensioned D(N). */

/*        E contains the subdiagonal elements of the tridiagonal */
/*          matrix in its last N-1 positions.  E(1) is set to zero. */
/*          E is a one-dimensional REAL array, dimensioned E(N). */

/*        E2 contains the squares of the corresponding elements of E. */
/*          E2 may coincide with E if the squares are not needed. */
/*          E2 is a one-dimensional REAL array, dimensioned E2(N). */

/*        Z contains the orthogonal transformation matrix produced in */
/*          the reduction if MATZ has been set to .TRUE.  Otherwise, Z */
/*          is not referenced.  Z is a two-dimensional REAL array, */
/*          dimensioned Z(NM,N). */

/*     Questions and comments should be directed to B. S. Garbow, */
/*     Applied Mathematics Division, ARGONNE NATIONAL LABORATORY */
/*     ------------------------------------------------------------------ */

/* ***REFERENCES  B. T. Smith, J. M. Boyle, J. J. Dongarra, B. S. Garbow, */
/*                 Y. Ikebe, V. C. Klema and C. B. Moler, Matrix Eigen- */
/*                 system Routines - EISPACK Guide, Springer-Verlag, */
/*                 1976. */
/* ***ROUTINES CALLED  (NONE) */
/* ***REVISION HISTORY  (YYMMDD) */
/*   760101  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890831  Modified array declarations.  (WRB) */
/*   890831  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   920501  Reformatted the REFERENCES section.  (WRB) */
/* ***END PROLOGUE  BANDR */


/* ***FIRST EXECUTABLE STATEMENT  BANDR */
    /* Parameter adjustments */
    z_dim1 = *nm;
    z_offset = 1 + z_dim1;
    z__ -= z_offset;
    a_dim1 = *nm;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --d__;
    --e;
    --e2;

    /* Function Body */
    dmin__ = 5.4210108624275222e-20f;
    dminrt = 2.3283064365386963e-10f;
/*     .......... INITIALIZE DIAGONAL SCALING MATRIX .......... */
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
/* L30: */
	d__[j] = 1.f;
    }

    if (! (*matz)) {
	goto L60;
    }

    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {

	i__2 = *n;
	for (k = 1; k <= i__2; ++k) {
/* L40: */
	    z__[j + k * z_dim1] = 0.f;
	}

	z__[j + j * z_dim1] = 1.f;
/* L50: */
    }

L60:
    m1 = *mb - 1;
    if ((i__1 = m1 - 1) < 0) {
	goto L900;
    } else if (i__1 == 0) {
	goto L800;
    } else {
	goto L70;
    }
L70:
    n2 = *n - 2;

    i__1 = n2;
    for (k = 1; k <= i__1; ++k) {
/* Computing MIN */
	i__2 = m1, i__3 = *n - k;
	maxr = min(i__2,i__3);
/*     .......... FOR R=MAXR STEP -1 UNTIL 2 DO -- .......... */
	i__2 = maxr;
	for (r1 = 2; r1 <= i__2; ++r1) {
	    r__ = maxr + 2 - r1;
	    kr = k + r__;
	    mr = *mb - r__;
	    g = a[kr + mr * a_dim1];
	    a[kr - 1 + a_dim1] = a[kr - 1 + (mr + 1) * a_dim1];
	    ugl = k;

	    i__3 = *n;
	    i__4 = m1;
	    for (j = kr; i__4 < 0 ? j >= i__3 : j <= i__3; j += i__4) {
		j1 = j - 1;
		j2 = j1 - 1;
		if (g == 0.f) {
		    goto L600;
		}
		b1 = a[j1 + a_dim1] / g;
		b2 = b1 * d__[j1] / d__[j];
		s2 = 1.f / (b1 * b2 + 1.f);
		if (s2 >= .5f) {
		    goto L450;
		}
		b1 = g / a[j1 + a_dim1];
		b2 = b1 * d__[j] / d__[j1];
		c2 = 1.f - s2;
		d__[j1] = c2 * d__[j1];
		d__[j] = c2 * d__[j];
		f1 = a[j + m1 * a_dim1] * 2.f;
		f2 = b1 * a[j1 + *mb * a_dim1];
		a[j + m1 * a_dim1] = -b2 * (b1 * a[j + m1 * a_dim1] - a[j + *
			mb * a_dim1]) - f2 + a[j + m1 * a_dim1];
		a[j1 + *mb * a_dim1] = b2 * (b2 * a[j + *mb * a_dim1] + f1) + 
			a[j1 + *mb * a_dim1];
		a[j + *mb * a_dim1] = b1 * (f2 - f1) + a[j + *mb * a_dim1];

		i__5 = j2;
		for (l = ugl; l <= i__5; ++l) {
		    i2 = *mb - j + l;
		    u = a[j1 + (i2 + 1) * a_dim1] + b2 * a[j + i2 * a_dim1];
		    a[j + i2 * a_dim1] = -b1 * a[j1 + (i2 + 1) * a_dim1] + a[
			    j + i2 * a_dim1];
		    a[j1 + (i2 + 1) * a_dim1] = u;
/* L200: */
		}

		ugl = j;
		a[j1 + a_dim1] += b2 * g;
		if (j == *n) {
		    goto L350;
		}
/* Computing MIN */
		i__5 = m1, i__6 = *n - j1;
		maxl = min(i__5,i__6);

		i__5 = maxl;
		for (l = 2; l <= i__5; ++l) {
		    i1 = j1 + l;
		    i2 = *mb - l;
		    u = a[i1 + i2 * a_dim1] + b2 * a[i1 + (i2 + 1) * a_dim1];
		    a[i1 + (i2 + 1) * a_dim1] = -b1 * a[i1 + i2 * a_dim1] + a[
			    i1 + (i2 + 1) * a_dim1];
		    a[i1 + i2 * a_dim1] = u;
/* L300: */
		}

		i1 = j + m1;
		if (i1 > *n) {
		    goto L350;
		}
		g = b2 * a[i1 + a_dim1];
L350:
		if (! (*matz)) {
		    goto L500;
		}

		i__5 = *n;
		for (l = 1; l <= i__5; ++l) {
		    u = z__[l + j1 * z_dim1] + b2 * z__[l + j * z_dim1];
		    z__[l + j * z_dim1] = -b1 * z__[l + j1 * z_dim1] + z__[l 
			    + j * z_dim1];
		    z__[l + j1 * z_dim1] = u;
/* L400: */
		}

		goto L500;

L450:
		u = d__[j1];
		d__[j1] = s2 * d__[j];
		d__[j] = s2 * u;
		f1 = a[j + m1 * a_dim1] * 2.f;
		f2 = b1 * a[j + *mb * a_dim1];
		u = b1 * (f2 - f1) + a[j1 + *mb * a_dim1];
		a[j + m1 * a_dim1] = b2 * (b1 * a[j + m1 * a_dim1] - a[j1 + *
			mb * a_dim1]) + f2 - a[j + m1 * a_dim1];
		a[j1 + *mb * a_dim1] = b2 * (b2 * a[j1 + *mb * a_dim1] + f1) 
			+ a[j + *mb * a_dim1];
		a[j + *mb * a_dim1] = u;

		i__5 = j2;
		for (l = ugl; l <= i__5; ++l) {
		    i2 = *mb - j + l;
		    u = b2 * a[j1 + (i2 + 1) * a_dim1] + a[j + i2 * a_dim1];
		    a[j + i2 * a_dim1] = -a[j1 + (i2 + 1) * a_dim1] + b1 * a[
			    j + i2 * a_dim1];
		    a[j1 + (i2 + 1) * a_dim1] = u;
/* L460: */
		}

		ugl = j;
		a[j1 + a_dim1] = b2 * a[j1 + a_dim1] + g;
		if (j == *n) {
		    goto L480;
		}
/* Computing MIN */
		i__5 = m1, i__6 = *n - j1;
		maxl = min(i__5,i__6);

		i__5 = maxl;
		for (l = 2; l <= i__5; ++l) {
		    i1 = j1 + l;
		    i2 = *mb - l;
		    u = b2 * a[i1 + i2 * a_dim1] + a[i1 + (i2 + 1) * a_dim1];
		    a[i1 + (i2 + 1) * a_dim1] = -a[i1 + i2 * a_dim1] + b1 * a[
			    i1 + (i2 + 1) * a_dim1];
		    a[i1 + i2 * a_dim1] = u;
/* L470: */
		}

		i1 = j + m1;
		if (i1 > *n) {
		    goto L480;
		}
		g = a[i1 + a_dim1];
		a[i1 + a_dim1] = b1 * a[i1 + a_dim1];
L480:
		if (! (*matz)) {
		    goto L500;
		}

		i__5 = *n;
		for (l = 1; l <= i__5; ++l) {
		    u = b2 * z__[l + j1 * z_dim1] + z__[l + j * z_dim1];
		    z__[l + j * z_dim1] = -z__[l + j1 * z_dim1] + b1 * z__[l 
			    + j * z_dim1];
		    z__[l + j1 * z_dim1] = u;
/* L490: */
		}

L500:
		;
	    }

L600:
	    ;
	}

	if (k % 64 != 0) {
	    goto L700;
	}
/*     .......... RESCALE TO AVOID UNDERFLOW OR OVERFLOW .......... */
	i__2 = *n;
	for (j = k; j <= i__2; ++j) {
	    if (d__[j] >= dmin__) {
		goto L650;
	    }
/* Computing MAX */
	    i__4 = 1, i__3 = *mb + 1 - j;
	    maxl = max(i__4,i__3);

	    i__4 = m1;
	    for (l = maxl; l <= i__4; ++l) {
/* L610: */
		a[j + l * a_dim1] = dminrt * a[j + l * a_dim1];
	    }

	    if (j == *n) {
		goto L630;
	    }
/* Computing MIN */
	    i__4 = m1, i__3 = *n - j;
	    maxl = min(i__4,i__3);

	    i__4 = maxl;
	    for (l = 1; l <= i__4; ++l) {
		i1 = j + l;
		i2 = *mb - l;
		a[i1 + i2 * a_dim1] = dminrt * a[i1 + i2 * a_dim1];
/* L620: */
	    }

L630:
	    if (! (*matz)) {
		goto L645;
	    }

	    i__4 = *n;
	    for (l = 1; l <= i__4; ++l) {
/* L640: */
		z__[l + j * z_dim1] = dminrt * z__[l + j * z_dim1];
	    }

L645:
	    a[j + *mb * a_dim1] = dmin__ * a[j + *mb * a_dim1];
	    d__[j] /= dmin__;
L650:
	    ;
	}

L700:
	;
    }
/*     .......... FORM SQUARE ROOT OF SCALING MATRIX .......... */
L800:
    i__1 = *n;
    for (j = 2; j <= i__1; ++j) {
/* L810: */
	e[j] = sqrt(d__[j]);
    }

    if (! (*matz)) {
	goto L840;
    }

    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {

	i__2 = *n;
	for (k = 2; k <= i__2; ++k) {
/* L820: */
	    z__[j + k * z_dim1] = e[k] * z__[j + k * z_dim1];
	}

/* L830: */
    }

L840:
    u = 1.f;

    i__1 = *n;
    for (j = 2; j <= i__1; ++j) {
	a[j + m1 * a_dim1] = u * e[j] * a[j + m1 * a_dim1];
	u = e[j];
/* Computing 2nd power */
	r__1 = a[j + m1 * a_dim1];
	e2[j] = r__1 * r__1;
	a[j + *mb * a_dim1] = d__[j] * a[j + *mb * a_dim1];
	d__[j] = a[j + *mb * a_dim1];
	e[j] = a[j + m1 * a_dim1];
/* L850: */
    }

    d__[1] = a[*mb * a_dim1 + 1];
    e[1] = 0.f;
    e2[1] = 0.f;
    goto L1001;

L900:
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	d__[j] = a[j + *mb * a_dim1];
	e[j] = 0.f;
	e2[j] = 0.f;
/* L950: */
    }

L1001:
    return 0;
} /* bandr_ */

