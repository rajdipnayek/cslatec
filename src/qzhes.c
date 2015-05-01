/* qzhes.f -- translated by f2c (version 12.02.01).
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

/* DECK QZHES */
/* Subroutine */ int qzhes_(integer *nm, integer *n, real *a, real *b, 
	logical *matz, real *z__)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, z_dim1, z_offset, i__1, i__2, 
	    i__3;
    real r__1, r__2;

    /* Local variables */
    static integer i__, j, k, l;
    static real r__, s, t;
    static integer l1;
    static real u1, u2, v1, v2;
    static integer lb, nk1, nm1, nm2;
    static real rho;

/* ***BEGIN PROLOGUE  QZHES */
/* ***PURPOSE  The first step of the QZ algorithm for solving generalized */
/*            matrix eigenproblems.  Accepts a pair of real general */
/*            matrices and reduces one of them to upper Hessenberg */
/*            and the other to upper triangular form using orthogonal */
/*            transformations. Usually followed by QZIT, QZVAL, QZVEC. */
/* ***LIBRARY   SLATEC (EISPACK) */
/* ***CATEGORY  D4C1B3 */
/* ***TYPE      SINGLE PRECISION (QZHES-S) */
/* ***KEYWORDS  EIGENVALUES, EIGENVECTORS, EISPACK */
/* ***AUTHOR  Smith, B. T., et al. */
/* ***DESCRIPTION */

/*     This subroutine is the first step of the QZ algorithm */
/*     for solving generalized matrix eigenvalue problems, */
/*     SIAM J. NUMER. ANAL. 10, 241-256(1973) by MOLER and STEWART. */

/*     This subroutine accepts a pair of REAL GENERAL matrices and */
/*     reduces one of them to upper Hessenberg form and the other */
/*     to upper triangular form using orthogonal transformations. */
/*     It is usually followed by  QZIT,  QZVAL  and, possibly,  QZVEC. */

/*     On Input */

/*        NM must be set to the row dimension of the two-dimensional */
/*          array parameters, A, B, and Z, as declared in the calling */
/*          program dimension statement.  NM is an INTEGER variable. */

/*        N is the order of the matrices A and B.  N is an INTEGER */
/*          variable.  N must be less than or equal to NM. */

/*        A contains a real general matrix.  A is a two-dimensional */
/*          REAL array, dimensioned A(NM,N). */

/*        B contains a real general matrix.  B is a two-dimensional */
/*          REAL array, dimensioned B(NM,N). */

/*        MATZ should be set to .TRUE. if the right hand transformations */
/*          are to be accumulated for later use in computing */
/*          eigenvectors, and to .FALSE. otherwise.  MATZ is a LOGICAL */
/*          variable. */

/*     On Output */

/*        A has been reduced to upper Hessenberg form.  The elements */
/*          below the first subdiagonal have been set to zero. */

/*        B has been reduced to upper triangular form.  The elements */
/*          below the main diagonal have been set to zero. */

/*        Z contains the product of the right hand transformations if */
/*          MATZ has been set to .TRUE.  Otherwise, Z is not referenced. */
/*          Z is a two-dimensional REAL array, dimensioned Z(NM,N). */

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
/* ***END PROLOGUE  QZHES */


/*     .......... INITIALIZE Z .......... */
/* ***FIRST EXECUTABLE STATEMENT  QZHES */
    /* Parameter adjustments */
    z_dim1 = *nm;
    z_offset = 1 + z_dim1;
    z__ -= z_offset;
    b_dim1 = *nm;
    b_offset = 1 + b_dim1;
    b -= b_offset;
    a_dim1 = *nm;
    a_offset = 1 + a_dim1;
    a -= a_offset;

    /* Function Body */
    if (! (*matz)) {
	goto L10;
    }

    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {

	i__2 = *n;
	for (j = 1; j <= i__2; ++j) {
	    z__[i__ + j * z_dim1] = 0.f;
/* L2: */
	}

	z__[i__ + i__ * z_dim1] = 1.f;
/* L3: */
    }
/*     .......... REDUCE B TO UPPER TRIANGULAR FORM .......... */
L10:
    if (*n <= 1) {
	goto L170;
    }
    nm1 = *n - 1;

    i__1 = nm1;
    for (l = 1; l <= i__1; ++l) {
	l1 = l + 1;
	s = 0.f;

	i__2 = *n;
	for (i__ = l1; i__ <= i__2; ++i__) {
	    s += (r__1 = b[i__ + l * b_dim1], dabs(r__1));
/* L20: */
	}

	if (s == 0.f) {
	    goto L100;
	}
	s += (r__1 = b[l + l * b_dim1], dabs(r__1));
	r__ = 0.f;

	i__2 = *n;
	for (i__ = l; i__ <= i__2; ++i__) {
	    b[i__ + l * b_dim1] /= s;
/* Computing 2nd power */
	    r__1 = b[i__ + l * b_dim1];
	    r__ += r__1 * r__1;
/* L25: */
	}

	r__1 = sqrt(r__);
	r__ = r_sign(&r__1, &b[l + l * b_dim1]);
	b[l + l * b_dim1] += r__;
	rho = r__ * b[l + l * b_dim1];

	i__2 = *n;
	for (j = l1; j <= i__2; ++j) {
	    t = 0.f;

	    i__3 = *n;
	    for (i__ = l; i__ <= i__3; ++i__) {
		t += b[i__ + l * b_dim1] * b[i__ + j * b_dim1];
/* L30: */
	    }

	    t = -t / rho;

	    i__3 = *n;
	    for (i__ = l; i__ <= i__3; ++i__) {
		b[i__ + j * b_dim1] += t * b[i__ + l * b_dim1];
/* L40: */
	    }

/* L50: */
	}

	i__2 = *n;
	for (j = 1; j <= i__2; ++j) {
	    t = 0.f;

	    i__3 = *n;
	    for (i__ = l; i__ <= i__3; ++i__) {
		t += b[i__ + l * b_dim1] * a[i__ + j * a_dim1];
/* L60: */
	    }

	    t = -t / rho;

	    i__3 = *n;
	    for (i__ = l; i__ <= i__3; ++i__) {
		a[i__ + j * a_dim1] += t * b[i__ + l * b_dim1];
/* L70: */
	    }

/* L80: */
	}

	b[l + l * b_dim1] = -s * r__;

	i__2 = *n;
	for (i__ = l1; i__ <= i__2; ++i__) {
	    b[i__ + l * b_dim1] = 0.f;
/* L90: */
	}

L100:
	;
    }
/*     .......... REDUCE A TO UPPER HESSENBERG FORM, WHILE */
/*                KEEPING B TRIANGULAR .......... */
    if (*n == 2) {
	goto L170;
    }
    nm2 = *n - 2;

    i__1 = nm2;
    for (k = 1; k <= i__1; ++k) {
	nk1 = nm1 - k;
/*     .......... FOR L=N-1 STEP -1 UNTIL K+1 DO -- .......... */
	i__2 = nk1;
	for (lb = 1; lb <= i__2; ++lb) {
	    l = *n - lb;
	    l1 = l + 1;
/*     .......... ZERO A(L+1,K) .......... */
	    s = (r__1 = a[l + k * a_dim1], dabs(r__1)) + (r__2 = a[l1 + k * 
		    a_dim1], dabs(r__2));
	    if (s == 0.f) {
		goto L150;
	    }
	    u1 = a[l + k * a_dim1] / s;
	    u2 = a[l1 + k * a_dim1] / s;
	    r__1 = sqrt(u1 * u1 + u2 * u2);
	    r__ = r_sign(&r__1, &u1);
	    v1 = -(u1 + r__) / r__;
	    v2 = -u2 / r__;
	    u2 = v2 / v1;

	    i__3 = *n;
	    for (j = k; j <= i__3; ++j) {
		t = a[l + j * a_dim1] + u2 * a[l1 + j * a_dim1];
		a[l + j * a_dim1] += t * v1;
		a[l1 + j * a_dim1] += t * v2;
/* L110: */
	    }

	    a[l1 + k * a_dim1] = 0.f;

	    i__3 = *n;
	    for (j = l; j <= i__3; ++j) {
		t = b[l + j * b_dim1] + u2 * b[l1 + j * b_dim1];
		b[l + j * b_dim1] += t * v1;
		b[l1 + j * b_dim1] += t * v2;
/* L120: */
	    }
/*     .......... ZERO B(L+1,L) .......... */
	    s = (r__1 = b[l1 + l1 * b_dim1], dabs(r__1)) + (r__2 = b[l1 + l * 
		    b_dim1], dabs(r__2));
	    if (s == 0.f) {
		goto L150;
	    }
	    u1 = b[l1 + l1 * b_dim1] / s;
	    u2 = b[l1 + l * b_dim1] / s;
	    r__1 = sqrt(u1 * u1 + u2 * u2);
	    r__ = r_sign(&r__1, &u1);
	    v1 = -(u1 + r__) / r__;
	    v2 = -u2 / r__;
	    u2 = v2 / v1;

	    i__3 = l1;
	    for (i__ = 1; i__ <= i__3; ++i__) {
		t = b[i__ + l1 * b_dim1] + u2 * b[i__ + l * b_dim1];
		b[i__ + l1 * b_dim1] += t * v1;
		b[i__ + l * b_dim1] += t * v2;
/* L130: */
	    }

	    b[l1 + l * b_dim1] = 0.f;

	    i__3 = *n;
	    for (i__ = 1; i__ <= i__3; ++i__) {
		t = a[i__ + l1 * a_dim1] + u2 * a[i__ + l * a_dim1];
		a[i__ + l1 * a_dim1] += t * v1;
		a[i__ + l * a_dim1] += t * v2;
/* L140: */
	    }

	    if (! (*matz)) {
		goto L150;
	    }

	    i__3 = *n;
	    for (i__ = 1; i__ <= i__3; ++i__) {
		t = z__[i__ + l1 * z_dim1] + u2 * z__[i__ + l * z_dim1];
		z__[i__ + l1 * z_dim1] += t * v1;
		z__[i__ + l * z_dim1] += t * v2;
/* L145: */
	    }

L150:
	    ;
	}

/* L160: */
    }

L170:
    return 0;
} /* qzhes_ */

