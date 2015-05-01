/* bqr.f -- translated by f2c (version 12.02.01).
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

static real c_b8 = 1.f;

/* DECK BQR */
/* Subroutine */ int bqr_(integer *nm, integer *n, integer *mb, real *a, real 
	*t, real *r__, integer *ierr, integer *nv, real *rv)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3;
    real r__1;

    /* Local variables */
    static real f, g;
    static integer i__, j, k, l, m;
    static real q, s;
    static integer m1, m2, m3, m4, m21, m31, ii, ik, jk, kj, jm, kk, km, ll, 
	    mk, mn, ni, mz, kj1, its;
    static real scale;
    static integer imult;
    extern doublereal pythag_(real *, real *);

/* ***BEGIN PROLOGUE  BQR */
/* ***PURPOSE  Compute some of the eigenvalues of a real symmetric */
/*            matrix using the QR method with shifts of origin. */
/* ***LIBRARY   SLATEC (EISPACK) */
/* ***CATEGORY  D4A6 */
/* ***TYPE      SINGLE PRECISION (BQR-S) */
/* ***KEYWORDS  EIGENVALUES, EISPACK */
/* ***AUTHOR  Smith, B. T., et al. */
/* ***DESCRIPTION */

/*     This subroutine is a translation of the ALGOL procedure BQR, */
/*     NUM. MATH. 16, 85-92(1970) by Martin, Reinsch, and Wilkinson. */
/*     HANDBOOK FOR AUTO. COMP., VOL II-LINEAR ALGEBRA, 266-272(1971). */

/*     This subroutine finds the eigenvalue of smallest (usually) */
/*     magnitude of a REAL SYMMETRIC BAND matrix using the */
/*     QR algorithm with shifts of origin.  Consecutive calls */
/*     can be made to find further eigenvalues. */

/*     On INPUT */

/*        NM must be set to the row dimension of the two-dimensional */
/*          array parameter, A, as declared in the calling program */
/*          dimension statement.  NM is an INTEGER variable. */

/*        N is the order of the matrix A.  N is an INTEGER variable. */
/*          N must be less than or equal to NM. */

/*        MB is the (half) band width of the matrix, defined as the */
/*          number of adjacent diagonals, including the principal */
/*          diagonal, required to specify the non-zero portion of the */
/*          lower triangle of the matrix.  MB is an INTEGER variable. */
/*          MB must be less than or equal to N on first call. */

/*        A contains the lower triangle of the symmetric band input */
/*          matrix stored as an N by MB array.  Its lowest subdiagonal */
/*          is stored in the last N+1-MB positions of the first column, */
/*          its next subdiagonal in the last N+2-MB positions of the */
/*          second column, further subdiagonals similarly, and finally */
/*          its principal diagonal in the N positions of the last column. */
/*          Contents of storages not part of the matrix are arbitrary. */
/*          On a subsequent call, its output contents from the previous */
/*          call should be passed.  A is a two-dimensional REAL array, */
/*          dimensioned A(NM,MB). */

/*        T specifies the shift (of eigenvalues) applied to the diagonal */
/*          of A in forming the input matrix. What is actually determined */
/*          is the eigenvalue of A+TI (I is the identity matrix) nearest */
/*          to T.  On a subsequent call, the output value of T from the */
/*          previous call should be passed if the next nearest eigenvalue */
/*          is sought.  T is a REAL variable. */

/*        R should be specified as zero on the first call, and as its */
/*          output value from the previous call on a subsequent call. */
/*          It is used to determine when the last row and column of */
/*          the transformed band matrix can be regarded as negligible. */
/*          R is a REAL variable. */

/*        NV must be set to the dimension of the array parameter RV */
/*          as declared in the calling program dimension statement. */
/*          NV is an INTEGER variable. */

/*     On OUTPUT */

/*        A contains the transformed band matrix.  The matrix A+TI */
/*          derived from the output parameters is similar to the */
/*          input A+TI to within rounding errors.  Its last row and */
/*          column are null (if IERR is zero). */

/*        T contains the computed eigenvalue of A+TI (if IERR is zero), */
/*          where I is the identity matrix. */

/*        R contains the maximum of its input value and the norm of the */
/*          last column of the input matrix A. */

/*        IERR is an INTEGER flag set to */
/*          Zero       for normal return, */
/*          J          if the J-th eigenvalue has not been */
/*                     determined after a total of 30 iterations. */

/*        RV is a one-dimensional REAL array of dimension NV which is */
/*          at least (2*MB**2+4*MB-3), used for temporary storage.  The */
/*          first (3*MB-2) locations correspond to the ALGOL array B, */
/*          the next (2*MB-1) locations correspond to the ALGOL array H, */
/*          and the final (2*MB**2-MB) locations correspond to the MB */
/*          by (2*MB-1) ALGOL array U. */

/*     NOTE. For a subsequent call, N should be replaced by N-1, but */
/*     MB should not be altered even when it exceeds the current N. */

/*     Calls PYTHAG(A,B) for SQRT(A**2 + B**2). */

/*     Questions and comments should be directed to B. S. Garbow, */
/*     Applied Mathematics Division, ARGONNE NATIONAL LABORATORY */
/*     ------------------------------------------------------------------ */

/* ***REFERENCES  B. T. Smith, J. M. Boyle, J. J. Dongarra, B. S. Garbow, */
/*                 Y. Ikebe, V. C. Klema and C. B. Moler, Matrix Eigen- */
/*                 system Routines - EISPACK Guide, Springer-Verlag, */
/*                 1976. */
/* ***ROUTINES CALLED  PYTHAG */
/* ***REVISION HISTORY  (YYMMDD) */
/*   760101  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890831  Modified array declarations.  (WRB) */
/*   890831  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   920501  Reformatted the REFERENCES section.  (WRB) */
/* ***END PROLOGUE  BQR */


/* ***FIRST EXECUTABLE STATEMENT  BQR */
    /* Parameter adjustments */
    a_dim1 = *nm;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --rv;

    /* Function Body */
    *ierr = 0;
    m1 = min(*mb,*n);
    m = m1 - 1;
    m2 = m + m;
    m21 = m2 + 1;
    m3 = m21 + m;
    m31 = m3 + 1;
    m4 = m31 + m2;
    mn = m + *n;
    mz = *mb - m1;
    its = 0;
/*     .......... TEST FOR CONVERGENCE .......... */
L40:
    g = a[*n + *mb * a_dim1];
    if (m == 0) {
	goto L360;
    }
    f = 0.f;

    i__1 = m;
    for (k = 1; k <= i__1; ++k) {
	mk = k + mz;
	f += (r__1 = a[*n + mk * a_dim1], dabs(r__1));
/* L50: */
    }

    if (its == 0 && f > *r__) {
	*r__ = f;
    }
    if (*r__ + f <= *r__) {
	goto L360;
    }
    if (its == 30) {
	goto L1000;
    }
    ++its;
/*     .......... FORM SHIFT FROM BOTTOM 2 BY 2 MINOR .......... */
    if (f > *r__ * .25f && its < 5) {
	goto L90;
    }
    f = a[*n + (*mb - 1) * a_dim1];
    if (f == 0.f) {
	goto L70;
    }
    q = (a[*n - 1 + *mb * a_dim1] - g) / (f * 2.f);
    s = pythag_(&q, &c_b8);
    g -= f / (q + r_sign(&s, &q));
L70:
    *t += g;

    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L80: */
	a[i__ + *mb * a_dim1] -= g;
    }

L90:
    i__1 = m4;
    for (k = m31; k <= i__1; ++k) {
/* L100: */
	rv[k] = 0.f;
    }

    i__1 = mn;
    for (ii = 1; ii <= i__1; ++ii) {
	i__ = ii - m;
	ni = *n - ii;
	if (ni < 0) {
	    goto L230;
	}
/*     .......... FORM COLUMN OF SHIFTED MATRIX A-G*I .......... */
/* Computing MAX */
	i__2 = 1, i__3 = 2 - i__;
	l = max(i__2,i__3);

	i__2 = m3;
	for (k = 1; k <= i__2; ++k) {
/* L110: */
	    rv[k] = 0.f;
	}

	i__2 = m1;
	for (k = l; k <= i__2; ++k) {
	    km = k + m;
	    mk = k + mz;
	    rv[km] = a[ii + mk * a_dim1];
/* L120: */
	}

	ll = min(m,ni);
	if (ll == 0) {
	    goto L135;
	}

	i__2 = ll;
	for (k = 1; k <= i__2; ++k) {
	    km = k + m21;
	    ik = ii + k;
	    mk = *mb - k;
	    rv[km] = a[ik + mk * a_dim1];
/* L130: */
	}
/*     .......... PRE-MULTIPLY WITH HOUSEHOLDER REFLECTIONS .......... */
L135:
	ll = m2;
	imult = 0;
/*     .......... MULTIPLICATION PROCEDURE .......... */
L140:
	kj = m4 - m1;

	i__2 = ll;
	for (j = 1; j <= i__2; ++j) {
	    kj += m1;
	    jm = j + m3;
	    if (rv[jm] == 0.f) {
		goto L170;
	    }
	    f = 0.f;

	    i__3 = m1;
	    for (k = 1; k <= i__3; ++k) {
		++kj;
		jk = j + k - 1;
		f += rv[kj] * rv[jk];
/* L150: */
	    }

	    f /= rv[jm];
	    kj -= m1;

	    i__3 = m1;
	    for (k = 1; k <= i__3; ++k) {
		++kj;
		jk = j + k - 1;
		rv[jk] -= rv[kj] * f;
/* L160: */
	    }

	    kj -= m1;
L170:
	    ;
	}

	if (imult != 0) {
	    goto L280;
	}
/*     .......... HOUSEHOLDER REFLECTION .......... */
	f = rv[m21];
	s = 0.f;
	rv[m4] = 0.f;
	scale = 0.f;

	i__2 = m3;
	for (k = m21; k <= i__2; ++k) {
/* L180: */
	    scale += (r__1 = rv[k], dabs(r__1));
	}

	if (scale == 0.f) {
	    goto L210;
	}

	i__2 = m3;
	for (k = m21; k <= i__2; ++k) {
/* L190: */
/* Computing 2nd power */
	    r__1 = rv[k] / scale;
	    s += r__1 * r__1;
	}

	s = scale * scale * s;
	r__1 = sqrt(s);
	g = -r_sign(&r__1, &f);
	rv[m21] = g;
	rv[m4] = s - f * g;
	kj = m4 + m2 * m1 + 1;
	rv[kj] = f - g;

	i__2 = m1;
	for (k = 2; k <= i__2; ++k) {
	    ++kj;
	    km = k + m2;
	    rv[kj] = rv[km];
/* L200: */
	}
/*     .......... SAVE COLUMN OF TRIANGULAR FACTOR R .......... */
L210:
	i__2 = m1;
	for (k = l; k <= i__2; ++k) {
	    km = k + m;
	    mk = k + mz;
	    a[ii + mk * a_dim1] = rv[km];
/* L220: */
	}

L230:
/* Computing MAX */
	i__2 = 1, i__3 = m1 + 1 - i__;
	l = max(i__2,i__3);
	if (i__ <= 0) {
	    goto L300;
	}
/*     .......... PERFORM ADDITIONAL STEPS .......... */
	i__2 = m21;
	for (k = 1; k <= i__2; ++k) {
/* L240: */
	    rv[k] = 0.f;
	}

/* Computing MIN */
	i__2 = m1, i__3 = ni + m1;
	ll = min(i__2,i__3);
/*     .......... GET ROW OF TRIANGULAR FACTOR R .......... */
	i__2 = ll;
	for (kk = 1; kk <= i__2; ++kk) {
	    k = kk - 1;
	    km = k + m1;
	    ik = i__ + k;
	    mk = *mb - k;
	    rv[km] = a[ik + mk * a_dim1];
/* L250: */
	}
/*     .......... POST-MULTIPLY WITH HOUSEHOLDER REFLECTIONS .......... */
	ll = m1;
	imult = 1;
	goto L140;
/*     .......... STORE COLUMN OF NEW A MATRIX .......... */
L280:
	i__2 = m1;
	for (k = l; k <= i__2; ++k) {
	    mk = k + mz;
	    a[i__ + mk * a_dim1] = rv[k];
/* L290: */
	}
/*     .......... UPDATE HOUSEHOLDER REFLECTIONS .......... */
L300:
	if (l > 1) {
	    --l;
	}
	kj1 = m4 + l * m1;

	i__2 = m2;
	for (j = l; j <= i__2; ++j) {
	    jm = j + m3;
	    rv[jm] = rv[jm + 1];

	    i__3 = m1;
	    for (k = 1; k <= i__3; ++k) {
		++kj1;
		kj = kj1 - m1;
		rv[kj] = rv[kj1];
/* L320: */
	    }
	}

/* L350: */
    }

    goto L40;
/*     .......... CONVERGENCE .......... */
L360:
    *t += g;

    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L380: */
	a[i__ + *mb * a_dim1] -= g;
    }

    i__1 = m1;
    for (k = 1; k <= i__1; ++k) {
	mk = k + mz;
	a[*n + mk * a_dim1] = 0.f;
/* L400: */
    }

    goto L1001;
/*     .......... SET ERROR -- NO CONVERGENCE TO */
/*                EIGENVALUE AFTER 30 ITERATIONS .......... */
L1000:
    *ierr = *n;
L1001:
    return 0;
} /* bqr_ */

