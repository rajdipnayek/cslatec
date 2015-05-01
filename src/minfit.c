/* minfit.f -- translated by f2c (version 12.02.01).
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

static real c_b39 = 1.f;

/* DECK MINFIT */
/* Subroutine */ int minfit_(integer *nm, integer *m, integer *n, real *a, 
	real *w, integer *ip, real *b, integer *ierr, real *rv1)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, i__1, i__2, i__3;
    real r__1, r__2, r__3, r__4;

    /* Local variables */
    static real c__, f, g, h__;
    static integer i__, j, k, l;
    static real s, x, y, z__;
    static integer i1, k1, l1, m1;
    static real s1;
    static integer ii, kk, ll, its;
    static real scale;
    extern doublereal pythag_(real *, real *);

/* ***BEGIN PROLOGUE  MINFIT */
/* ***PURPOSE  Compute the singular value decomposition of a rectangular */
/*            matrix and solve the related linear least squares problem. */
/* ***LIBRARY   SLATEC (EISPACK) */
/* ***CATEGORY  D9 */
/* ***TYPE      SINGLE PRECISION (MINFIT-S) */
/* ***KEYWORDS  EIGENVALUES, EIGENVECTORS, EISPACK */
/* ***AUTHOR  Smith, B. T., et al. */
/* ***DESCRIPTION */

/*     This subroutine is a translation of the ALGOL procedure MINFIT, */
/*     NUM. MATH. 14, 403-420(1970) by Golub and Reinsch. */
/*     HANDBOOK FOR AUTO. COMP., VOL II-LINEAR ALGEBRA, 134-151(1971). */

/*     This subroutine determines, towards the solution of the linear */
/*                                                        T */
/*     system AX=B, the singular value decomposition A=USV  of a real */
/*                                         T */
/*     M by N rectangular matrix, forming U B rather than U.  Householder */
/*     bidiagonalization and a variant of the QR algorithm are used. */

/*     On INPUT */

/*        NM must be set to the row dimension of the two-dimensional */
/*          array parameters, A and B, as declared in the calling */
/*          program dimension statement.  Note that NM must be at least */
/*          as large as the maximum of M and N.  NM is an INTEGER */
/*          variable. */

/*        M is the number of rows of A and B.  M is an INTEGER variable. */

/*        N is the number of columns of A and the order of V.  N is an */
/*          INTEGER variable. */

/*        A contains the rectangular coefficient matrix of the system. */
/*          A is a two-dimensional REAL array, dimensioned A(NM,N). */

/*        IP is the number of columns of B.  IP can be zero. */

/*        B contains the constant column matrix of the system if IP is */
/*          not zero.  Otherwise, B is not referenced.  B is a two- */
/*          dimensional REAL array, dimensioned B(NM,IP). */

/*     On OUTPUT */

/*        A has been overwritten by the matrix V (orthogonal) of the */
/*          decomposition in its first N rows and columns.  If an */
/*          error exit is made, the columns of V corresponding to */
/*          indices of correct singular values should be correct. */

/*        W contains the N (non-negative) singular values of A (the */
/*          diagonal elements of S).  They are unordered.  If an */
/*          error exit is made, the singular values should be correct */
/*          for indices IERR+1, IERR+2, ..., N.  W is a one-dimensional */
/*          REAL array, dimensioned W(N). */

/*                                   T */
/*        B has been overwritten by U B.  If an error exit is made, */
/*                       T */
/*          the rows of U B corresponding to indices of correct singular */
/*          values should be correct. */

/*        IERR is an INTEGER flag set to */
/*          Zero       for normal return, */
/*          K          if the K-th singular value has not been */
/*                     determined after 30 iterations. */
/*                     The singular values should be correct for */
/*                     indices IERR+1, IERR+2, ..., N. */

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
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890831  Modified array declarations.  (WRB) */
/*   890831  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   920501  Reformatted the REFERENCES section.  (WRB) */
/* ***END PROLOGUE  MINFIT */


/* ***FIRST EXECUTABLE STATEMENT  MINFIT */
    /* Parameter adjustments */
    a_dim1 = *nm;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --w;
    b_dim1 = *nm;
    b_offset = 1 + b_dim1;
    b -= b_offset;
    --rv1;

    /* Function Body */
    *ierr = 0;
/*     .......... HOUSEHOLDER REDUCTION TO BIDIAGONAL FORM .......... */
    g = 0.f;
    scale = 0.f;
    s1 = 0.f;

    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	l = i__ + 1;
	rv1[i__] = scale * g;
	g = 0.f;
	s = 0.f;
	scale = 0.f;
	if (i__ > *m) {
	    goto L210;
	}

	i__2 = *m;
	for (k = i__; k <= i__2; ++k) {
/* L120: */
	    scale += (r__1 = a[k + i__ * a_dim1], dabs(r__1));
	}

	if (scale == 0.f) {
	    goto L210;
	}

	i__2 = *m;
	for (k = i__; k <= i__2; ++k) {
	    a[k + i__ * a_dim1] /= scale;
/* Computing 2nd power */
	    r__1 = a[k + i__ * a_dim1];
	    s += r__1 * r__1;
/* L130: */
	}

	f = a[i__ + i__ * a_dim1];
	r__1 = sqrt(s);
	g = -r_sign(&r__1, &f);
	h__ = f * g - s;
	a[i__ + i__ * a_dim1] = f - g;
	if (i__ == *n) {
	    goto L160;
	}

	i__2 = *n;
	for (j = l; j <= i__2; ++j) {
	    s = 0.f;

	    i__3 = *m;
	    for (k = i__; k <= i__3; ++k) {
/* L140: */
		s += a[k + i__ * a_dim1] * a[k + j * a_dim1];
	    }

	    f = s / h__;

	    i__3 = *m;
	    for (k = i__; k <= i__3; ++k) {
		a[k + j * a_dim1] += f * a[k + i__ * a_dim1];
/* L150: */
	    }
	}

L160:
	if (*ip == 0) {
	    goto L190;
	}

	i__3 = *ip;
	for (j = 1; j <= i__3; ++j) {
	    s = 0.f;

	    i__2 = *m;
	    for (k = i__; k <= i__2; ++k) {
/* L170: */
		s += a[k + i__ * a_dim1] * b[k + j * b_dim1];
	    }

	    f = s / h__;

	    i__2 = *m;
	    for (k = i__; k <= i__2; ++k) {
		b[k + j * b_dim1] += f * a[k + i__ * a_dim1];
/* L180: */
	    }
	}

L190:
	i__2 = *m;
	for (k = i__; k <= i__2; ++k) {
/* L200: */
	    a[k + i__ * a_dim1] = scale * a[k + i__ * a_dim1];
	}

L210:
	w[i__] = scale * g;
	g = 0.f;
	s = 0.f;
	scale = 0.f;
	if (i__ > *m || i__ == *n) {
	    goto L290;
	}

	i__2 = *n;
	for (k = l; k <= i__2; ++k) {
/* L220: */
	    scale += (r__1 = a[i__ + k * a_dim1], dabs(r__1));
	}

	if (scale == 0.f) {
	    goto L290;
	}

	i__2 = *n;
	for (k = l; k <= i__2; ++k) {
	    a[i__ + k * a_dim1] /= scale;
/* Computing 2nd power */
	    r__1 = a[i__ + k * a_dim1];
	    s += r__1 * r__1;
/* L230: */
	}

	f = a[i__ + l * a_dim1];
	r__1 = sqrt(s);
	g = -r_sign(&r__1, &f);
	h__ = f * g - s;
	a[i__ + l * a_dim1] = f - g;

	i__2 = *n;
	for (k = l; k <= i__2; ++k) {
/* L240: */
	    rv1[k] = a[i__ + k * a_dim1] / h__;
	}

	if (i__ == *m) {
	    goto L270;
	}

	i__2 = *m;
	for (j = l; j <= i__2; ++j) {
	    s = 0.f;

	    i__3 = *n;
	    for (k = l; k <= i__3; ++k) {
/* L250: */
		s += a[j + k * a_dim1] * a[i__ + k * a_dim1];
	    }

	    i__3 = *n;
	    for (k = l; k <= i__3; ++k) {
		a[j + k * a_dim1] += s * rv1[k];
/* L260: */
	    }
	}

L270:
	i__3 = *n;
	for (k = l; k <= i__3; ++k) {
/* L280: */
	    a[i__ + k * a_dim1] = scale * a[i__ + k * a_dim1];
	}

L290:
/* Computing MAX */
	r__3 = s1, r__4 = (r__1 = w[i__], dabs(r__1)) + (r__2 = rv1[i__], 
		dabs(r__2));
	s1 = dmax(r__3,r__4);
/* L300: */
    }
/*     .......... ACCUMULATION OF RIGHT-HAND TRANSFORMATIONS. */
/*                FOR I=N STEP -1 UNTIL 1 DO -- .......... */
    i__1 = *n;
    for (ii = 1; ii <= i__1; ++ii) {
	i__ = *n + 1 - ii;
	if (i__ == *n) {
	    goto L390;
	}
	if (g == 0.f) {
	    goto L360;
	}

	i__3 = *n;
	for (j = l; j <= i__3; ++j) {
/*     .......... DOUBLE DIVISION AVOIDS POSSIBLE UNDERFLOW .......... */
/* L320: */
	    a[j + i__ * a_dim1] = a[i__ + j * a_dim1] / a[i__ + l * a_dim1] / 
		    g;
	}

	i__3 = *n;
	for (j = l; j <= i__3; ++j) {
	    s = 0.f;

	    i__2 = *n;
	    for (k = l; k <= i__2; ++k) {
/* L340: */
		s += a[i__ + k * a_dim1] * a[k + j * a_dim1];
	    }

	    i__2 = *n;
	    for (k = l; k <= i__2; ++k) {
		a[k + j * a_dim1] += s * a[k + i__ * a_dim1];
/* L350: */
	    }
	}

L360:
	i__2 = *n;
	for (j = l; j <= i__2; ++j) {
	    a[i__ + j * a_dim1] = 0.f;
	    a[j + i__ * a_dim1] = 0.f;
/* L380: */
	}

L390:
	a[i__ + i__ * a_dim1] = 1.f;
	g = rv1[i__];
	l = i__;
/* L400: */
    }

    if (*m >= *n || *ip == 0) {
	goto L510;
    }
    m1 = *m + 1;

    i__1 = *n;
    for (i__ = m1; i__ <= i__1; ++i__) {

	i__2 = *ip;
	for (j = 1; j <= i__2; ++j) {
	    b[i__ + j * b_dim1] = 0.f;
/* L500: */
	}
    }
/*     .......... DIAGONALIZATION OF THE BIDIAGONAL FORM .......... */
L510:
/*     .......... FOR K=N STEP -1 UNTIL 1 DO -- .......... */
    i__2 = *n;
    for (kk = 1; kk <= i__2; ++kk) {
	k1 = *n - kk;
	k = k1 + 1;
	its = 0;
/*     .......... TEST FOR SPLITTING. */
/*                FOR L=K STEP -1 UNTIL 1 DO -- .......... */
L520:
	i__1 = k;
	for (ll = 1; ll <= i__1; ++ll) {
	    l1 = k - ll;
	    l = l1 + 1;
	    if (s1 + (r__1 = rv1[l], dabs(r__1)) == s1) {
		goto L565;
	    }
/*     .......... RV1(1) IS ALWAYS ZERO, SO THERE IS NO EXIT */
/*                THROUGH THE BOTTOM OF THE LOOP .......... */
	    if (s1 + (r__1 = w[l1], dabs(r__1)) == s1) {
		goto L540;
	    }
/* L530: */
	}
/*     .......... CANCELLATION OF RV1(L) IF L GREATER THAN 1 .......... */
L540:
	c__ = 0.f;
	s = 1.f;

	i__1 = k;
	for (i__ = l; i__ <= i__1; ++i__) {
	    f = s * rv1[i__];
	    rv1[i__] = c__ * rv1[i__];
	    if (s1 + dabs(f) == s1) {
		goto L565;
	    }
	    g = w[i__];
	    h__ = pythag_(&f, &g);
	    w[i__] = h__;
	    c__ = g / h__;
	    s = -f / h__;
	    if (*ip == 0) {
		goto L560;
	    }

	    i__3 = *ip;
	    for (j = 1; j <= i__3; ++j) {
		y = b[l1 + j * b_dim1];
		z__ = b[i__ + j * b_dim1];
		b[l1 + j * b_dim1] = y * c__ + z__ * s;
		b[i__ + j * b_dim1] = -y * s + z__ * c__;
/* L550: */
	    }

L560:
	    ;
	}
/*     .......... TEST FOR CONVERGENCE .......... */
L565:
	z__ = w[k];
	if (l == k) {
	    goto L650;
	}
/*     .......... SHIFT FROM BOTTOM 2 BY 2 MINOR .......... */
	if (its == 30) {
	    goto L1000;
	}
	++its;
	x = w[l];
	y = w[k1];
	g = rv1[k1];
	h__ = rv1[k];
	f = ((g + z__) / h__ * ((g - z__) / y) + y / h__ - h__ / y) * .5f;
	g = pythag_(&f, &c_b39);
	f = x - z__ / x * z__ + h__ / x * (y / (f + r_sign(&g, &f)) - h__);
/*     .......... NEXT QR TRANSFORMATION .......... */
	c__ = 1.f;
	s = 1.f;

	i__1 = k1;
	for (i1 = l; i1 <= i__1; ++i1) {
	    i__ = i1 + 1;
	    g = rv1[i__];
	    y = w[i__];
	    h__ = s * g;
	    g = c__ * g;
	    z__ = pythag_(&f, &h__);
	    rv1[i1] = z__;
	    c__ = f / z__;
	    s = h__ / z__;
	    f = x * c__ + g * s;
	    g = -x * s + g * c__;
	    h__ = y * s;
	    y *= c__;

	    i__3 = *n;
	    for (j = 1; j <= i__3; ++j) {
		x = a[j + i1 * a_dim1];
		z__ = a[j + i__ * a_dim1];
		a[j + i1 * a_dim1] = x * c__ + z__ * s;
		a[j + i__ * a_dim1] = -x * s + z__ * c__;
/* L570: */
	    }

	    z__ = pythag_(&f, &h__);
	    w[i1] = z__;
/*     .......... ROTATION CAN BE ARBITRARY IF Z IS ZERO .......... */
	    if (z__ == 0.f) {
		goto L580;
	    }
	    c__ = f / z__;
	    s = h__ / z__;
L580:
	    f = c__ * g + s * y;
	    x = -s * g + c__ * y;
	    if (*ip == 0) {
		goto L600;
	    }

	    i__3 = *ip;
	    for (j = 1; j <= i__3; ++j) {
		y = b[i1 + j * b_dim1];
		z__ = b[i__ + j * b_dim1];
		b[i1 + j * b_dim1] = y * c__ + z__ * s;
		b[i__ + j * b_dim1] = -y * s + z__ * c__;
/* L590: */
	    }

L600:
	    ;
	}

	rv1[l] = 0.f;
	rv1[k] = f;
	w[k] = x;
	goto L520;
/*     .......... CONVERGENCE .......... */
L650:
	if (z__ >= 0.f) {
	    goto L700;
	}
/*     .......... W(K) IS MADE NON-NEGATIVE .......... */
	w[k] = -z__;

	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
/* L690: */
	    a[j + k * a_dim1] = -a[j + k * a_dim1];
	}

L700:
	;
    }

    goto L1001;
/*     .......... SET ERROR -- NO CONVERGENCE TO A */
/*                SINGULAR VALUE AFTER 30 ITERATIONS .......... */
L1000:
    *ierr = k;
L1001:
    return 0;
} /* minfit_ */

