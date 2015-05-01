/* qzit.f -- translated by f2c (version 12.02.01).
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

/* DECK QZIT */
/* Subroutine */ int qzit_(integer *nm, integer *n, real *a, real *b, real *
	eps1, logical *matz, real *z__, integer *ierr)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, z_dim1, z_offset, i__1, i__2, 
	    i__3;
    real r__1, r__2, r__3;

    /* Local variables */
    static integer i__, j, k, l;
    static real r__, s, t, a1, a2, a3;
    static integer k1, k2, l1;
    static real u1, u2, u3, v1, v2, v3, a11, a12, a21, a22, a33, a34, a43, 
	    a44, b11, b12, b22, b33;
    static integer na, ld;
    static real b34, b44;
    static integer en;
    static real ep;
    static integer ll;
    static real sh;
    static integer km1, lm1;
    static real ani, bni;
    static integer ish, itn, its, enm2, lor1;
    static real epsa, epsb, anorm, bnorm;
    static integer enorn;
    static logical notlas;

/* ***BEGIN PROLOGUE  QZIT */
/* ***PURPOSE  The second step of the QZ algorithm for generalized */
/*            eigenproblems.  Accepts an upper Hessenberg and an upper */
/*            triangular matrix and reduces the former to */
/*            quasi-triangular form while preserving the form of the */
/*            latter.  Usually preceded by QZHES and followed by QZVAL */
/*            and QZVEC. */
/* ***LIBRARY   SLATEC (EISPACK) */
/* ***CATEGORY  D4C1B3 */
/* ***TYPE      SINGLE PRECISION (QZIT-S) */
/* ***KEYWORDS  EIGENVALUES, EIGENVECTORS, EISPACK */
/* ***AUTHOR  Smith, B. T., et al. */
/* ***DESCRIPTION */

/*     This subroutine is the second step of the QZ algorithm */
/*     for solving generalized matrix eigenvalue problems, */
/*     SIAM J. NUMER. ANAL. 10, 241-256(1973) by MOLER and STEWART, */
/*     as modified in technical note NASA TN D-7305(1973) by WARD. */

/*     This subroutine accepts a pair of REAL matrices, one of them */
/*     in upper Hessenberg form and the other in upper triangular form. */
/*     It reduces the Hessenberg matrix to quasi-triangular form using */
/*     orthogonal transformations while maintaining the triangular form */
/*     of the other matrix.  It is usually preceded by  QZHES  and */
/*     followed by  QZVAL  and, possibly,  QZVEC. */

/*     On Input */

/*        NM must be set to the row dimension of the two-dimensional */
/*          array parameters, A, B, and Z, as declared in the calling */
/*          program dimension statement.  NM is an INTEGER variable. */

/*        N is the order of the matrices A and B.  N is an INTEGER */
/*          variable.  N must be less than or equal to NM. */

/*        A contains a real upper Hessenberg matrix.  A is a two- */
/*          dimensional REAL array, dimensioned A(NM,N). */

/*        B contains a real upper triangular matrix.  B is a two- */
/*          dimensional REAL array, dimensioned B(NM,N). */

/*        EPS1 is a tolerance used to determine negligible elements. */
/*          EPS1 = 0.0 (or negative) may be input, in which case an */
/*          element will be neglected only if it is less than roundoff */
/*          error times the norm of its matrix.  If the input EPS1 is */
/*          positive, then an element will be considered negligible */
/*          if it is less than EPS1 times the norm of its matrix.  A */
/*          positive value of EPS1 may result in faster execution, */
/*          but less accurate results.  EPS1 is a REAL variable. */

/*        MATZ should be set to .TRUE. if the right hand transformations */
/*          are to be accumulated for later use in computing */
/*          eigenvectors, and to .FALSE. otherwise.  MATZ is a LOGICAL */
/*          variable. */

/*        Z contains, if MATZ has been set to .TRUE., the transformation */
/*          matrix produced in the reduction by  QZHES, if performed, or */
/*          else the identity matrix.  If MATZ has been set to .FALSE., */
/*          Z is not referenced.  Z is a two-dimensional REAL array, */
/*          dimensioned Z(NM,N). */

/*     On Output */

/*        A has been reduced to quasi-triangular form.  The elements */
/*          below the first subdiagonal are still zero, and no two */
/*          consecutive subdiagonal elements are nonzero. */

/*        B is still in upper triangular form, although its elements */
/*          have been altered.  The location B(N,1) is used to store */
/*          EPS1 times the norm of B for later use by  QZVAL  and  QZVEC. */

/*        Z contains the product of the right hand transformations */
/*          (for both steps) if MATZ has been set to .TRUE. */

/*        IERR is an INTEGER flag set to */
/*          Zero       for normal return, */
/*          J          if neither A(J,J-1) nor A(J-1,J-2) has become */
/*                     zero after a total of 30*N iterations. */

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
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890831  Modified array declarations.  (WRB) */
/*   890831  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   920501  Reformatted the REFERENCES section.  (WRB) */
/* ***END PROLOGUE  QZIT */


/* ***FIRST EXECUTABLE STATEMENT  QZIT */
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
    *ierr = 0;
/*     .......... COMPUTE EPSA,EPSB .......... */
    anorm = 0.f;
    bnorm = 0.f;

    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ani = 0.f;
	if (i__ != 1) {
	    ani = (r__1 = a[i__ + (i__ - 1) * a_dim1], dabs(r__1));
	}
	bni = 0.f;

	i__2 = *n;
	for (j = i__; j <= i__2; ++j) {
	    ani += (r__1 = a[i__ + j * a_dim1], dabs(r__1));
	    bni += (r__1 = b[i__ + j * b_dim1], dabs(r__1));
/* L20: */
	}

	if (ani > anorm) {
	    anorm = ani;
	}
	if (bni > bnorm) {
	    bnorm = bni;
	}
/* L30: */
    }

    if (anorm == 0.f) {
	anorm = 1.f;
    }
    if (bnorm == 0.f) {
	bnorm = 1.f;
    }
    ep = *eps1;
    if (ep > 0.f) {
	goto L50;
    }
/*     .......... COMPUTE ROUNDOFF LEVEL IF EPS1 IS ZERO .......... */
    ep = 1.f;
L40:
    ep /= 2.f;
    if (ep + 1.f > 1.f) {
	goto L40;
    }
L50:
    epsa = ep * anorm;
    epsb = ep * bnorm;
/*     .......... REDUCE A TO QUASI-TRIANGULAR FORM, WHILE */
/*                KEEPING B TRIANGULAR .......... */
    lor1 = 1;
    enorn = *n;
    en = *n;
    itn = *n * 30;
/*     .......... BEGIN QZ STEP .......... */
L60:
    if (en <= 2) {
	goto L1001;
    }
    if (! (*matz)) {
	enorn = en;
    }
    its = 0;
    na = en - 1;
    enm2 = na - 1;
L70:
    ish = 2;
/*     .......... CHECK FOR CONVERGENCE OR REDUCIBILITY. */
/*                FOR L=EN STEP -1 UNTIL 1 DO -- .......... */
    i__1 = en;
    for (ll = 1; ll <= i__1; ++ll) {
	lm1 = en - ll;
	l = lm1 + 1;
	if (l == 1) {
	    goto L95;
	}
	if ((r__1 = a[l + lm1 * a_dim1], dabs(r__1)) <= epsa) {
	    goto L90;
	}
/* L80: */
    }

L90:
    a[l + lm1 * a_dim1] = 0.f;
    if (l < na) {
	goto L95;
    }
/*     .......... 1-BY-1 OR 2-BY-2 BLOCK ISOLATED .......... */
    en = lm1;
    goto L60;
/*     .......... CHECK FOR SMALL TOP OF B .......... */
L95:
    ld = l;
L100:
    l1 = l + 1;
    b11 = b[l + l * b_dim1];
    if (dabs(b11) > epsb) {
	goto L120;
    }
    b[l + l * b_dim1] = 0.f;
    s = (r__1 = a[l + l * a_dim1], dabs(r__1)) + (r__2 = a[l1 + l * a_dim1], 
	    dabs(r__2));
    u1 = a[l + l * a_dim1] / s;
    u2 = a[l1 + l * a_dim1] / s;
    r__1 = sqrt(u1 * u1 + u2 * u2);
    r__ = r_sign(&r__1, &u1);
    v1 = -(u1 + r__) / r__;
    v2 = -u2 / r__;
    u2 = v2 / v1;

    i__1 = enorn;
    for (j = l; j <= i__1; ++j) {
	t = a[l + j * a_dim1] + u2 * a[l1 + j * a_dim1];
	a[l + j * a_dim1] += t * v1;
	a[l1 + j * a_dim1] += t * v2;
	t = b[l + j * b_dim1] + u2 * b[l1 + j * b_dim1];
	b[l + j * b_dim1] += t * v1;
	b[l1 + j * b_dim1] += t * v2;
/* L110: */
    }

    if (l != 1) {
	a[l + lm1 * a_dim1] = -a[l + lm1 * a_dim1];
    }
    lm1 = l;
    l = l1;
    goto L90;
L120:
    a11 = a[l + l * a_dim1] / b11;
    a21 = a[l1 + l * a_dim1] / b11;
    if (ish == 1) {
	goto L140;
    }
/*     .......... ITERATION STRATEGY .......... */
    if (itn == 0) {
	goto L1000;
    }
    if (its == 10) {
	goto L155;
    }
/*     .......... DETERMINE TYPE OF SHIFT .......... */
    b22 = b[l1 + l1 * b_dim1];
    if (dabs(b22) < epsb) {
	b22 = epsb;
    }
    b33 = b[na + na * b_dim1];
    if (dabs(b33) < epsb) {
	b33 = epsb;
    }
    b44 = b[en + en * b_dim1];
    if (dabs(b44) < epsb) {
	b44 = epsb;
    }
    a33 = a[na + na * a_dim1] / b33;
    a34 = a[na + en * a_dim1] / b44;
    a43 = a[en + na * a_dim1] / b33;
    a44 = a[en + en * a_dim1] / b44;
    b34 = b[na + en * b_dim1] / b44;
    t = (a43 * b34 - a33 - a44) * .5f;
    r__ = t * t + a34 * a43 - a33 * a44;
    if (r__ < 0.f) {
	goto L150;
    }
/*     .......... DETERMINE SINGLE SHIFT ZEROTH COLUMN OF A .......... */
    ish = 1;
    r__ = sqrt(r__);
    sh = -t + r__;
    s = -t - r__;
    if ((r__1 = s - a44, dabs(r__1)) < (r__2 = sh - a44, dabs(r__2))) {
	sh = s;
    }
/*     .......... LOOK FOR TWO CONSECUTIVE SMALL */
/*                SUB-DIAGONAL ELEMENTS OF A. */
/*                FOR L=EN-2 STEP -1 UNTIL LD DO -- .......... */
    i__1 = enm2;
    for (ll = ld; ll <= i__1; ++ll) {
	l = enm2 + ld - ll;
	if (l == ld) {
	    goto L140;
	}
	lm1 = l - 1;
	l1 = l + 1;
	t = a[l + l * a_dim1];
	if ((r__1 = b[l + l * b_dim1], dabs(r__1)) > epsb) {
	    t -= sh * b[l + l * b_dim1];
	}
	if ((r__2 = a[l + lm1 * a_dim1], dabs(r__2)) <= (r__1 = t / a[l1 + l *
		 a_dim1], dabs(r__1)) * epsa) {
	    goto L100;
	}
/* L130: */
    }

L140:
    a1 = a11 - sh;
    a2 = a21;
    if (l != ld) {
	a[l + lm1 * a_dim1] = -a[l + lm1 * a_dim1];
    }
    goto L160;
/*     .......... DETERMINE DOUBLE SHIFT ZEROTH COLUMN OF A .......... */
L150:
    a12 = a[l + l1 * a_dim1] / b22;
    a22 = a[l1 + l1 * a_dim1] / b22;
    b12 = b[l + l1 * b_dim1] / b22;
    a1 = ((a33 - a11) * (a44 - a11) - a34 * a43 + a43 * b34 * a11) / a21 + 
	    a12 - a11 * b12;
    a2 = a22 - a11 - a21 * b12 - (a33 - a11) - (a44 - a11) + a43 * b34;
    a3 = a[l1 + 1 + l1 * a_dim1] / b22;
    goto L160;
/*     .......... AD HOC SHIFT .......... */
L155:
    a1 = 0.f;
    a2 = 1.f;
    a3 = 1.1605f;
L160:
    ++its;
    --itn;
    if (! (*matz)) {
	lor1 = ld;
    }
/*     .......... MAIN LOOP .......... */
    i__1 = na;
    for (k = l; k <= i__1; ++k) {
	notlas = k != na && ish == 2;
	k1 = k + 1;
	k2 = k + 2;
/* Computing MAX */
	i__2 = k - 1;
	km1 = max(i__2,l);
/* Computing MIN */
	i__2 = en, i__3 = k1 + ish;
	ll = min(i__2,i__3);
	if (notlas) {
	    goto L190;
	}
/*     .......... ZERO A(K+1,K-1) .......... */
	if (k == l) {
	    goto L170;
	}
	a1 = a[k + km1 * a_dim1];
	a2 = a[k1 + km1 * a_dim1];
L170:
	s = dabs(a1) + dabs(a2);
	if (s == 0.f) {
	    goto L70;
	}
	u1 = a1 / s;
	u2 = a2 / s;
	r__1 = sqrt(u1 * u1 + u2 * u2);
	r__ = r_sign(&r__1, &u1);
	v1 = -(u1 + r__) / r__;
	v2 = -u2 / r__;
	u2 = v2 / v1;

	i__2 = enorn;
	for (j = km1; j <= i__2; ++j) {
	    t = a[k + j * a_dim1] + u2 * a[k1 + j * a_dim1];
	    a[k + j * a_dim1] += t * v1;
	    a[k1 + j * a_dim1] += t * v2;
	    t = b[k + j * b_dim1] + u2 * b[k1 + j * b_dim1];
	    b[k + j * b_dim1] += t * v1;
	    b[k1 + j * b_dim1] += t * v2;
/* L180: */
	}

	if (k != l) {
	    a[k1 + km1 * a_dim1] = 0.f;
	}
	goto L240;
/*     .......... ZERO A(K+1,K-1) AND A(K+2,K-1) .......... */
L190:
	if (k == l) {
	    goto L200;
	}
	a1 = a[k + km1 * a_dim1];
	a2 = a[k1 + km1 * a_dim1];
	a3 = a[k2 + km1 * a_dim1];
L200:
	s = dabs(a1) + dabs(a2) + dabs(a3);
	if (s == 0.f) {
	    goto L260;
	}
	u1 = a1 / s;
	u2 = a2 / s;
	u3 = a3 / s;
	r__1 = sqrt(u1 * u1 + u2 * u2 + u3 * u3);
	r__ = r_sign(&r__1, &u1);
	v1 = -(u1 + r__) / r__;
	v2 = -u2 / r__;
	v3 = -u3 / r__;
	u2 = v2 / v1;
	u3 = v3 / v1;

	i__2 = enorn;
	for (j = km1; j <= i__2; ++j) {
	    t = a[k + j * a_dim1] + u2 * a[k1 + j * a_dim1] + u3 * a[k2 + j * 
		    a_dim1];
	    a[k + j * a_dim1] += t * v1;
	    a[k1 + j * a_dim1] += t * v2;
	    a[k2 + j * a_dim1] += t * v3;
	    t = b[k + j * b_dim1] + u2 * b[k1 + j * b_dim1] + u3 * b[k2 + j * 
		    b_dim1];
	    b[k + j * b_dim1] += t * v1;
	    b[k1 + j * b_dim1] += t * v2;
	    b[k2 + j * b_dim1] += t * v3;
/* L210: */
	}

	if (k == l) {
	    goto L220;
	}
	a[k1 + km1 * a_dim1] = 0.f;
	a[k2 + km1 * a_dim1] = 0.f;
/*     .......... ZERO B(K+2,K+1) AND B(K+2,K) .......... */
L220:
	s = (r__1 = b[k2 + k2 * b_dim1], dabs(r__1)) + (r__2 = b[k2 + k1 * 
		b_dim1], dabs(r__2)) + (r__3 = b[k2 + k * b_dim1], dabs(r__3))
		;
	if (s == 0.f) {
	    goto L240;
	}
	u1 = b[k2 + k2 * b_dim1] / s;
	u2 = b[k2 + k1 * b_dim1] / s;
	u3 = b[k2 + k * b_dim1] / s;
	r__1 = sqrt(u1 * u1 + u2 * u2 + u3 * u3);
	r__ = r_sign(&r__1, &u1);
	v1 = -(u1 + r__) / r__;
	v2 = -u2 / r__;
	v3 = -u3 / r__;
	u2 = v2 / v1;
	u3 = v3 / v1;

	i__2 = ll;
	for (i__ = lor1; i__ <= i__2; ++i__) {
	    t = a[i__ + k2 * a_dim1] + u2 * a[i__ + k1 * a_dim1] + u3 * a[i__ 
		    + k * a_dim1];
	    a[i__ + k2 * a_dim1] += t * v1;
	    a[i__ + k1 * a_dim1] += t * v2;
	    a[i__ + k * a_dim1] += t * v3;
	    t = b[i__ + k2 * b_dim1] + u2 * b[i__ + k1 * b_dim1] + u3 * b[i__ 
		    + k * b_dim1];
	    b[i__ + k2 * b_dim1] += t * v1;
	    b[i__ + k1 * b_dim1] += t * v2;
	    b[i__ + k * b_dim1] += t * v3;
/* L230: */
	}

	b[k2 + k * b_dim1] = 0.f;
	b[k2 + k1 * b_dim1] = 0.f;
	if (! (*matz)) {
	    goto L240;
	}

	i__2 = *n;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    t = z__[i__ + k2 * z_dim1] + u2 * z__[i__ + k1 * z_dim1] + u3 * 
		    z__[i__ + k * z_dim1];
	    z__[i__ + k2 * z_dim1] += t * v1;
	    z__[i__ + k1 * z_dim1] += t * v2;
	    z__[i__ + k * z_dim1] += t * v3;
/* L235: */
	}
/*     .......... ZERO B(K+1,K) .......... */
L240:
	s = (r__1 = b[k1 + k1 * b_dim1], dabs(r__1)) + (r__2 = b[k1 + k * 
		b_dim1], dabs(r__2));
	if (s == 0.f) {
	    goto L260;
	}
	u1 = b[k1 + k1 * b_dim1] / s;
	u2 = b[k1 + k * b_dim1] / s;
	r__1 = sqrt(u1 * u1 + u2 * u2);
	r__ = r_sign(&r__1, &u1);
	v1 = -(u1 + r__) / r__;
	v2 = -u2 / r__;
	u2 = v2 / v1;

	i__2 = ll;
	for (i__ = lor1; i__ <= i__2; ++i__) {
	    t = a[i__ + k1 * a_dim1] + u2 * a[i__ + k * a_dim1];
	    a[i__ + k1 * a_dim1] += t * v1;
	    a[i__ + k * a_dim1] += t * v2;
	    t = b[i__ + k1 * b_dim1] + u2 * b[i__ + k * b_dim1];
	    b[i__ + k1 * b_dim1] += t * v1;
	    b[i__ + k * b_dim1] += t * v2;
/* L250: */
	}

	b[k1 + k * b_dim1] = 0.f;
	if (! (*matz)) {
	    goto L260;
	}

	i__2 = *n;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    t = z__[i__ + k1 * z_dim1] + u2 * z__[i__ + k * z_dim1];
	    z__[i__ + k1 * z_dim1] += t * v1;
	    z__[i__ + k * z_dim1] += t * v2;
/* L255: */
	}

L260:
	;
    }
/*     .......... END QZ STEP .......... */
    goto L70;
/*     .......... SET ERROR -- NEITHER BOTTOM SUBDIAGONAL ELEMENT */
/*                HAS BECOME NEGLIGIBLE AFTER 30*N ITERATIONS .......... */
L1000:
    *ierr = en;
/*     .......... SAVE EPSB FOR USE BY QZVAL AND QZVEC .......... */
L1001:
    if (*n > 1) {
	b[*n + b_dim1] = epsb;
    }
    return 0;
} /* qzit_ */

