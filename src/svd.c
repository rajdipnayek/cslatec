/* svd.f -- translated by f2c (version 12.02.01).
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

static real c_b47 = 1.f;

/* DECK SVD */
/* Subroutine */ int svd_(integer *nm, integer *m, integer *n, real *a, real *
	w, logical *matu, real *u, logical *matv, real *v, integer *ierr, 
	real *rv1)
{
    /* System generated locals */
    integer a_dim1, a_offset, u_dim1, u_offset, v_dim1, v_offset, i__1, i__2, 
	    i__3;
    real r__1, r__2, r__3, r__4;

    /* Local variables */
    static real c__, f, g, h__;
    static integer i__, j, k, l;
    static real s, x, y, z__;
    static integer i1, k1, l1;
    static real s1;
    static integer ii, kk, ll, mn, its;
    static real scale;
    extern doublereal pythag_(real *, real *);

/* ***BEGIN PROLOGUE  SVD */
/* ***SUBSIDIARY */
/* ***PURPOSE  Perform the singular value decomposition of a rectangular */
/*            matrix. */
/* ***LIBRARY   SLATEC */
/* ***TYPE      SINGLE PRECISION (SVD-S) */
/* ***AUTHOR  (UNKNOWN) */
/* ***DESCRIPTION */

/*     This subroutine is a translation of the ALGOL procedure SVD, */
/*     NUM. MATH. 14, 403-420(1970) by Golub and Reinsch. */
/*     HANDBOOK FOR AUTO. COMP., VOL II-LINEAR ALGEBRA, 134-151(1971). */

/*     This subroutine determines the singular value decomposition */
/*          T */
/*     A=USV  of a REAL M by N rectangular matrix.  Householder */
/*     bidiagonalization and a variant of the QR algorithm are used. */

/*     On Input */

/*        NM must be set to the row dimension of the two-dimensional */
/*          array parameters, A, U and V, as declared in the calling */
/*          program dimension statement.  NM is an INTEGER variable. */
/*          Note that NM must be at least as large as the maximum */
/*          of M and N. */

/*        M is the number of rows of A and U. */

/*        N is the number of columns of A and U and the order of V. */

/*        A contains the rectangular input matrix to be decomposed.  A is */
/*          a two-dimensional REAL array, dimensioned A(NM,N). */

/*        MATU should be set to .TRUE. if the U matrix in the */
/*          decomposition is desired, and to .FALSE. otherwise. */
/*          MATU is a LOGICAL variable. */

/*        MATV should be set to .TRUE. if the V matrix in the */
/*          decomposition is desired, and to .FALSE. otherwise. */
/*          MATV is a LOGICAL variable. */

/*     On Output */

/*        A is unaltered (unless overwritten by U or V). */

/*        W contains the N (non-negative) singular values of A (the */
/*          diagonal elements of S).  They are unordered.  If an */
/*          error exit is made, the singular values should be correct */
/*          for indices IERR+1, IERR+2, ..., N.  W is a one-dimensional */
/*          REAL array, dimensioned W(N). */

/*        U contains the matrix U (orthogonal column vectors) of the */
/*          decomposition if MATU has been set to .TRUE.  Otherwise, */
/*          U is used as a temporary array.  U may coincide with A. */
/*          If an error exit is made, the columns of U corresponding */
/*          to indices of correct singular values should be correct. */
/*          U is a two-dimensional REAL array, dimensioned U(NM,N). */

/*        V contains the matrix V (orthogonal) of the decomposition if */
/*          MATV has been set to .TRUE.  Otherwise, V is not referenced. */
/*          V may also coincide with A if U does not.  If an error */
/*          exit is made, the columns of V corresponding to indices of */
/*          correct singular values should be correct.  V is a two- */
/*          dimensional REAL array, dimensioned V(NM,N). */

/*        IERR is an INTEGER flag set to */
/*          Zero       for normal return, */
/*          K          if the K-th singular value has not been */
/*                     determined after 30 iterations. */

/*        RV1 is a one-dimensional REAL array used for temporary storage, */
/*          dimensioned RV1(N). */

/*     CALLS PYTHAG(A,B) for sqrt(A**2 + B**2). */

/*     Questions and comments should be directed to B. S. Garbow, */
/*     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY */
/*     ------------------------------------------------------------------ */

/* ***SEE ALSO  EISDOC */
/* ***ROUTINES CALLED  PYTHAG */
/* ***REVISION HISTORY  (YYMMDD) */
/*   811101  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890831  Modified array declarations.  (WRB) */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900402  Added TYPE section.  (WRB) */
/* ***END PROLOGUE  SVD */


/* ***FIRST EXECUTABLE STATEMENT  SVD */
    /* Parameter adjustments */
    v_dim1 = *nm;
    v_offset = 1 + v_dim1;
    v -= v_offset;
    u_dim1 = *nm;
    u_offset = 1 + u_dim1;
    u -= u_offset;
    a_dim1 = *nm;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --w;
    --rv1;

    /* Function Body */
    *ierr = 0;

    i__1 = *m;
    for (i__ = 1; i__ <= i__1; ++i__) {

	i__2 = *n;
	for (j = 1; j <= i__2; ++j) {
	    u[i__ + j * u_dim1] = a[i__ + j * a_dim1];
/* L100: */
	}
    }
/*     .......... HOUSEHOLDER REDUCTION TO BIDIAGONAL FORM .......... */
    g = 0.f;
    scale = 0.f;
    s1 = 0.f;

    i__2 = *n;
    for (i__ = 1; i__ <= i__2; ++i__) {
	l = i__ + 1;
	rv1[i__] = scale * g;
	g = 0.f;
	s = 0.f;
	scale = 0.f;
	if (i__ > *m) {
	    goto L210;
	}

	i__1 = *m;
	for (k = i__; k <= i__1; ++k) {
/* L120: */
	    scale += (r__1 = u[k + i__ * u_dim1], dabs(r__1));
	}

	if (scale == 0.f) {
	    goto L210;
	}

	i__1 = *m;
	for (k = i__; k <= i__1; ++k) {
	    u[k + i__ * u_dim1] /= scale;
/* Computing 2nd power */
	    r__1 = u[k + i__ * u_dim1];
	    s += r__1 * r__1;
/* L130: */
	}

	f = u[i__ + i__ * u_dim1];
	r__1 = sqrt(s);
	g = -r_sign(&r__1, &f);
	h__ = f * g - s;
	u[i__ + i__ * u_dim1] = f - g;
	if (i__ == *n) {
	    goto L190;
	}

	i__1 = *n;
	for (j = l; j <= i__1; ++j) {
	    s = 0.f;

	    i__3 = *m;
	    for (k = i__; k <= i__3; ++k) {
/* L140: */
		s += u[k + i__ * u_dim1] * u[k + j * u_dim1];
	    }

	    f = s / h__;

	    i__3 = *m;
	    for (k = i__; k <= i__3; ++k) {
		u[k + j * u_dim1] += f * u[k + i__ * u_dim1];
/* L150: */
	    }
	}

L190:
	i__3 = *m;
	for (k = i__; k <= i__3; ++k) {
/* L200: */
	    u[k + i__ * u_dim1] = scale * u[k + i__ * u_dim1];
	}

L210:
	w[i__] = scale * g;
	g = 0.f;
	s = 0.f;
	scale = 0.f;
	if (i__ > *m || i__ == *n) {
	    goto L290;
	}

	i__3 = *n;
	for (k = l; k <= i__3; ++k) {
/* L220: */
	    scale += (r__1 = u[i__ + k * u_dim1], dabs(r__1));
	}

	if (scale == 0.f) {
	    goto L290;
	}

	i__3 = *n;
	for (k = l; k <= i__3; ++k) {
	    u[i__ + k * u_dim1] /= scale;
/* Computing 2nd power */
	    r__1 = u[i__ + k * u_dim1];
	    s += r__1 * r__1;
/* L230: */
	}

	f = u[i__ + l * u_dim1];
	r__1 = sqrt(s);
	g = -r_sign(&r__1, &f);
	h__ = f * g - s;
	u[i__ + l * u_dim1] = f - g;

	i__3 = *n;
	for (k = l; k <= i__3; ++k) {
/* L240: */
	    rv1[k] = u[i__ + k * u_dim1] / h__;
	}

	if (i__ == *m) {
	    goto L270;
	}

	i__3 = *m;
	for (j = l; j <= i__3; ++j) {
	    s = 0.f;

	    i__1 = *n;
	    for (k = l; k <= i__1; ++k) {
/* L250: */
		s += u[j + k * u_dim1] * u[i__ + k * u_dim1];
	    }

	    i__1 = *n;
	    for (k = l; k <= i__1; ++k) {
		u[j + k * u_dim1] += s * rv1[k];
/* L260: */
	    }
	}

L270:
	i__1 = *n;
	for (k = l; k <= i__1; ++k) {
/* L280: */
	    u[i__ + k * u_dim1] = scale * u[i__ + k * u_dim1];
	}

L290:
/* Computing MAX */
	r__3 = s1, r__4 = (r__1 = w[i__], dabs(r__1)) + (r__2 = rv1[i__], 
		dabs(r__2));
	s1 = dmax(r__3,r__4);
/* L300: */
    }
/*     .......... ACCUMULATION OF RIGHT-HAND TRANSFORMATIONS .......... */
    if (! (*matv)) {
	goto L410;
    }
/*     .......... FOR I=N STEP -1 UNTIL 1 DO -- .......... */
    i__2 = *n;
    for (ii = 1; ii <= i__2; ++ii) {
	i__ = *n + 1 - ii;
	if (i__ == *n) {
	    goto L390;
	}
	if (g == 0.f) {
	    goto L360;
	}

	i__1 = *n;
	for (j = l; j <= i__1; ++j) {
/*     .......... DOUBLE DIVISION AVOIDS POSSIBLE UNDERFLOW .......... */
/* L320: */
	    v[j + i__ * v_dim1] = u[i__ + j * u_dim1] / u[i__ + l * u_dim1] / 
		    g;
	}

	i__1 = *n;
	for (j = l; j <= i__1; ++j) {
	    s = 0.f;

	    i__3 = *n;
	    for (k = l; k <= i__3; ++k) {
/* L340: */
		s += u[i__ + k * u_dim1] * v[k + j * v_dim1];
	    }

	    i__3 = *n;
	    for (k = l; k <= i__3; ++k) {
		v[k + j * v_dim1] += s * v[k + i__ * v_dim1];
/* L350: */
	    }
	}

L360:
	i__3 = *n;
	for (j = l; j <= i__3; ++j) {
	    v[i__ + j * v_dim1] = 0.f;
	    v[j + i__ * v_dim1] = 0.f;
/* L380: */
	}

L390:
	v[i__ + i__ * v_dim1] = 1.f;
	g = rv1[i__];
	l = i__;
/* L400: */
    }
/*     .......... ACCUMULATION OF LEFT-HAND TRANSFORMATIONS .......... */
L410:
    if (! (*matu)) {
	goto L510;
    }
/*     ..........FOR I=MIN(M,N) STEP -1 UNTIL 1 DO -- .......... */
    mn = *n;
    if (*m < *n) {
	mn = *m;
    }

    i__2 = mn;
    for (ii = 1; ii <= i__2; ++ii) {
	i__ = mn + 1 - ii;
	l = i__ + 1;
	g = w[i__];
	if (i__ == *n) {
	    goto L430;
	}

	i__3 = *n;
	for (j = l; j <= i__3; ++j) {
/* L420: */
	    u[i__ + j * u_dim1] = 0.f;
	}

L430:
	if (g == 0.f) {
	    goto L475;
	}
	if (i__ == mn) {
	    goto L460;
	}

	i__3 = *n;
	for (j = l; j <= i__3; ++j) {
	    s = 0.f;

	    i__1 = *m;
	    for (k = l; k <= i__1; ++k) {
/* L440: */
		s += u[k + i__ * u_dim1] * u[k + j * u_dim1];
	    }
/*     .......... DOUBLE DIVISION AVOIDS POSSIBLE UNDERFLOW .......... */
	    f = s / u[i__ + i__ * u_dim1] / g;

	    i__1 = *m;
	    for (k = i__; k <= i__1; ++k) {
		u[k + j * u_dim1] += f * u[k + i__ * u_dim1];
/* L450: */
	    }
	}

L460:
	i__1 = *m;
	for (j = i__; j <= i__1; ++j) {
/* L470: */
	    u[j + i__ * u_dim1] /= g;
	}

	goto L490;

L475:
	i__1 = *m;
	for (j = i__; j <= i__1; ++j) {
/* L480: */
	    u[j + i__ * u_dim1] = 0.f;
	}

L490:
	u[i__ + i__ * u_dim1] += 1.f;
/* L500: */
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
	    if (! (*matu)) {
		goto L560;
	    }

	    i__3 = *m;
	    for (j = 1; j <= i__3; ++j) {
		y = u[j + l1 * u_dim1];
		z__ = u[j + i__ * u_dim1];
		u[j + l1 * u_dim1] = y * c__ + z__ * s;
		u[j + i__ * u_dim1] = -y * s + z__ * c__;
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
	g = pythag_(&f, &c_b47);
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
	    if (! (*matv)) {
		goto L575;
	    }

	    i__3 = *n;
	    for (j = 1; j <= i__3; ++j) {
		x = v[j + i1 * v_dim1];
		z__ = v[j + i__ * v_dim1];
		v[j + i1 * v_dim1] = x * c__ + z__ * s;
		v[j + i__ * v_dim1] = -x * s + z__ * c__;
/* L570: */
	    }

L575:
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
	    if (! (*matu)) {
		goto L600;
	    }

	    i__3 = *m;
	    for (j = 1; j <= i__3; ++j) {
		y = u[j + i1 * u_dim1];
		z__ = u[j + i__ * u_dim1];
		u[j + i1 * u_dim1] = y * c__ + z__ * s;
		u[j + i__ * u_dim1] = -y * s + z__ * c__;
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
	if (! (*matv)) {
	    goto L700;
	}

	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
/* L690: */
	    v[j + k * v_dim1] = -v[j + k * v_dim1];
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
} /* svd_ */

