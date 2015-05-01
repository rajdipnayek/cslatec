/* qzvec.f -- translated by f2c (version 12.02.01).
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

/* DECK QZVEC */
/* Subroutine */ int qzvec_(integer *nm, integer *n, real *a, real *b, real *
	alfr, real *alfi, real *beta, real *z__)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, z_dim1, z_offset, i__1, i__2, 
	    i__3;
    real r__1, r__2;

    /* Local variables */
    static real d__;
    static integer i__, j, k, m;
    static real q, r__, s, t, w, x, y, t1, t2, w1, x1, z1, di;
    static integer na, ii, en, jj;
    static real ra, dr, sa;
    static integer nn;
    static real ti, rr, tr, zz;
    static integer isw, enm2;
    static real alfm, almi, betm, epsb, almr;

/* ***BEGIN PROLOGUE  QZVEC */
/* ***PURPOSE  The optional fourth step of the QZ algorithm for */
/*            generalized eigenproblems.  Accepts a matrix in */
/*            quasi-triangular form and another in upper triangular */
/*            and computes the eigenvectors of the triangular problem */
/*            and transforms them back to the original coordinates */
/*            Usually preceded by QZHES, QZIT, and QZVAL. */
/* ***LIBRARY   SLATEC (EISPACK) */
/* ***CATEGORY  D4C3 */
/* ***TYPE      SINGLE PRECISION (QZVEC-S) */
/* ***KEYWORDS  EIGENVALUES, EIGENVECTORS, EISPACK */
/* ***AUTHOR  Smith, B. T., et al. */
/* ***DESCRIPTION */

/*     This subroutine is the optional fourth step of the QZ algorithm */
/*     for solving generalized matrix eigenvalue problems, */
/*     SIAM J. NUMER. ANAL. 10, 241-256(1973) by MOLER and STEWART. */

/*     This subroutine accepts a pair of REAL matrices, one of them in */
/*     quasi-triangular form (in which each 2-by-2 block corresponds to */
/*     a pair of complex eigenvalues) and the other in upper triangular */
/*     form.  It computes the eigenvectors of the triangular problem and */
/*     transforms the results back to the original coordinate system. */
/*     It is usually preceded by  QZHES,  QZIT, and  QZVAL. */

/*     On Input */

/*        NM must be set to the row dimension of the two-dimensional */
/*          array parameters, A, B, and Z, as declared in the calling */
/*          program dimension statement.  NM is an INTEGER variable. */

/*        N is the order of the matrices A and B.  N is an INTEGER */
/*          variable.  N must be less than or equal to NM. */

/*        A contains a real upper quasi-triangular matrix.  A is a two- */
/*          dimensional REAL array, dimensioned A(NM,N). */

/*        B contains a real upper triangular matrix.  In addition, */
/*          location B(N,1) contains the tolerance quantity (EPSB) */
/*          computed and saved in  QZIT.  B is a two-dimensional REAL */
/*          array, dimensioned B(NM,N). */

/*        ALFR, ALFI, and BETA are one-dimensional REAL arrays with */
/*          components whose ratios ((ALFR+I*ALFI)/BETA) are the */
/*          generalized eigenvalues.  They are usually obtained from */
/*          QZVAL.  They are dimensioned ALFR(N), ALFI(N), and BETA(N). */

/*        Z contains the transformation matrix produced in the reductions */
/*          by  QZHES,  QZIT, and  QZVAL,  if performed.  If the */
/*          eigenvectors of the triangular problem are desired, Z must */
/*          contain the identity matrix.  Z is a two-dimensional REAL */
/*          array, dimensioned Z(NM,N). */

/*     On Output */

/*        A is unaltered.  Its subdiagonal elements provide information */
/*           about the storage of the complex eigenvectors. */

/*        B has been destroyed. */

/*        ALFR, ALFI, and BETA are unaltered. */

/*        Z contains the real and imaginary parts of the eigenvectors. */
/*          If ALFI(J) .EQ. 0.0, the J-th eigenvalue is real and */
/*            the J-th column of Z contains its eigenvector. */
/*          If ALFI(J) .NE. 0.0, the J-th eigenvalue is complex. */
/*            If ALFI(J) .GT. 0.0, the eigenvalue is the first of */
/*              a complex pair and the J-th and (J+1)-th columns */
/*              of Z contain its eigenvector. */
/*            If ALFI(J) .LT. 0.0, the eigenvalue is the second of */
/*              a complex pair and the (J-1)-th and J-th columns */
/*              of Z contain the conjugate of its eigenvector. */
/*          Each eigenvector is normalized so that the modulus */
/*          of its largest component is 1.0 . */

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
/* ***END PROLOGUE  QZVEC */


/* ***FIRST EXECUTABLE STATEMENT  QZVEC */
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
    --alfr;
    --alfi;
    --beta;

    /* Function Body */
    epsb = b[*n + b_dim1];
    isw = 1;
/*     .......... FOR EN=N STEP -1 UNTIL 1 DO -- .......... */
    i__1 = *n;
    for (nn = 1; nn <= i__1; ++nn) {
	en = *n + 1 - nn;
	na = en - 1;
	if (isw == 2) {
	    goto L795;
	}
	if (alfi[en] != 0.f) {
	    goto L710;
	}
/*     .......... REAL VECTOR .......... */
	m = en;
	b[en + en * b_dim1] = 1.f;
	if (na == 0) {
	    goto L800;
	}
	alfm = alfr[m];
	betm = beta[m];
/*     .......... FOR I=EN-1 STEP -1 UNTIL 1 DO -- .......... */
	i__2 = na;
	for (ii = 1; ii <= i__2; ++ii) {
	    i__ = en - ii;
	    w = betm * a[i__ + i__ * a_dim1] - alfm * b[i__ + i__ * b_dim1];
	    r__ = 0.f;

	    i__3 = en;
	    for (j = m; j <= i__3; ++j) {
/* L610: */
		r__ += (betm * a[i__ + j * a_dim1] - alfm * b[i__ + j * 
			b_dim1]) * b[j + en * b_dim1];
	    }

	    if (i__ == 1 || isw == 2) {
		goto L630;
	    }
	    if (betm * a[i__ + (i__ - 1) * a_dim1] == 0.f) {
		goto L630;
	    }
	    zz = w;
	    s = r__;
	    goto L690;
L630:
	    m = i__;
	    if (isw == 2) {
		goto L640;
	    }
/*     .......... REAL 1-BY-1 BLOCK .......... */
	    t = w;
	    if (w == 0.f) {
		t = epsb;
	    }
	    b[i__ + en * b_dim1] = -r__ / t;
	    goto L700;
/*     .......... REAL 2-BY-2 BLOCK .......... */
L640:
	    x = betm * a[i__ + (i__ + 1) * a_dim1] - alfm * b[i__ + (i__ + 1) 
		    * b_dim1];
	    y = betm * a[i__ + 1 + i__ * a_dim1];
	    q = w * zz - x * y;
	    t = (x * s - zz * r__) / q;
	    b[i__ + en * b_dim1] = t;
	    if (dabs(x) <= dabs(zz)) {
		goto L650;
	    }
	    b[i__ + 1 + en * b_dim1] = (-r__ - w * t) / x;
	    goto L690;
L650:
	    b[i__ + 1 + en * b_dim1] = (-s - y * t) / zz;
L690:
	    isw = 3 - isw;
L700:
	    ;
	}
/*     .......... END REAL VECTOR .......... */
	goto L800;
/*     .......... COMPLEX VECTOR .......... */
L710:
	m = na;
	almr = alfr[m];
	almi = alfi[m];
	betm = beta[m];
/*     .......... LAST VECTOR COMPONENT CHOSEN IMAGINARY SO THAT */
/*                EIGENVECTOR MATRIX IS TRIANGULAR .......... */
	y = betm * a[en + na * a_dim1];
	b[na + na * b_dim1] = -almi * b[en + en * b_dim1] / y;
	b[na + en * b_dim1] = (almr * b[en + en * b_dim1] - betm * a[en + en *
		 a_dim1]) / y;
	b[en + na * b_dim1] = 0.f;
	b[en + en * b_dim1] = 1.f;
	enm2 = na - 1;
	if (enm2 == 0) {
	    goto L795;
	}
/*     .......... FOR I=EN-2 STEP -1 UNTIL 1 DO -- .......... */
	i__2 = enm2;
	for (ii = 1; ii <= i__2; ++ii) {
	    i__ = na - ii;
	    w = betm * a[i__ + i__ * a_dim1] - almr * b[i__ + i__ * b_dim1];
	    w1 = -almi * b[i__ + i__ * b_dim1];
	    ra = 0.f;
	    sa = 0.f;

	    i__3 = en;
	    for (j = m; j <= i__3; ++j) {
		x = betm * a[i__ + j * a_dim1] - almr * b[i__ + j * b_dim1];
		x1 = -almi * b[i__ + j * b_dim1];
		ra = ra + x * b[j + na * b_dim1] - x1 * b[j + en * b_dim1];
		sa = sa + x * b[j + en * b_dim1] + x1 * b[j + na * b_dim1];
/* L760: */
	    }

	    if (i__ == 1 || isw == 2) {
		goto L770;
	    }
	    if (betm * a[i__ + (i__ - 1) * a_dim1] == 0.f) {
		goto L770;
	    }
	    zz = w;
	    z1 = w1;
	    r__ = ra;
	    s = sa;
	    isw = 2;
	    goto L790;
L770:
	    m = i__;
	    if (isw == 2) {
		goto L780;
	    }
/*     .......... COMPLEX 1-BY-1 BLOCK .......... */
	    tr = -ra;
	    ti = -sa;
L773:
	    dr = w;
	    di = w1;
/*     .......... COMPLEX DIVIDE (T1,T2) = (TR,TI) / (DR,DI) .......... */
L775:
	    if (dabs(di) > dabs(dr)) {
		goto L777;
	    }
	    rr = di / dr;
	    d__ = dr + di * rr;
	    t1 = (tr + ti * rr) / d__;
	    t2 = (ti - tr * rr) / d__;
	    switch (isw) {
		case 1:  goto L787;
		case 2:  goto L782;
	    }
L777:
	    rr = dr / di;
	    d__ = dr * rr + di;
	    t1 = (tr * rr + ti) / d__;
	    t2 = (ti * rr - tr) / d__;
	    switch (isw) {
		case 1:  goto L787;
		case 2:  goto L782;
	    }
/*     .......... COMPLEX 2-BY-2 BLOCK .......... */
L780:
	    x = betm * a[i__ + (i__ + 1) * a_dim1] - almr * b[i__ + (i__ + 1) 
		    * b_dim1];
	    x1 = -almi * b[i__ + (i__ + 1) * b_dim1];
	    y = betm * a[i__ + 1 + i__ * a_dim1];
	    tr = y * ra - w * r__ + w1 * s;
	    ti = y * sa - w * s - w1 * r__;
	    dr = w * zz - w1 * z1 - x * y;
	    di = w * z1 + w1 * zz - x1 * y;
	    if (dr == 0.f && di == 0.f) {
		dr = epsb;
	    }
	    goto L775;
L782:
	    b[i__ + 1 + na * b_dim1] = t1;
	    b[i__ + 1 + en * b_dim1] = t2;
	    isw = 1;
	    if (dabs(y) > dabs(w) + dabs(w1)) {
		goto L785;
	    }
	    tr = -ra - x * b[i__ + 1 + na * b_dim1] + x1 * b[i__ + 1 + en * 
		    b_dim1];
	    ti = -sa - x * b[i__ + 1 + en * b_dim1] - x1 * b[i__ + 1 + na * 
		    b_dim1];
	    goto L773;
L785:
	    t1 = (-r__ - zz * b[i__ + 1 + na * b_dim1] + z1 * b[i__ + 1 + en *
		     b_dim1]) / y;
	    t2 = (-s - zz * b[i__ + 1 + en * b_dim1] - z1 * b[i__ + 1 + na * 
		    b_dim1]) / y;
L787:
	    b[i__ + na * b_dim1] = t1;
	    b[i__ + en * b_dim1] = t2;
L790:
	    ;
	}
/*     .......... END COMPLEX VECTOR .......... */
L795:
	isw = 3 - isw;
L800:
	;
    }
/*     .......... END BACK SUBSTITUTION. */
/*                TRANSFORM TO ORIGINAL COORDINATE SYSTEM. */
/*                FOR J=N STEP -1 UNTIL 1 DO -- .......... */
    i__1 = *n;
    for (jj = 1; jj <= i__1; ++jj) {
	j = *n + 1 - jj;

	i__2 = *n;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    zz = 0.f;

	    i__3 = j;
	    for (k = 1; k <= i__3; ++k) {
/* L860: */
		zz += z__[i__ + k * z_dim1] * b[k + j * b_dim1];
	    }

	    z__[i__ + j * z_dim1] = zz;
/* L880: */
	}
    }
/*     .......... NORMALIZE SO THAT MODULUS OF LARGEST */
/*                COMPONENT OF EACH VECTOR IS 1. */
/*                (ISW IS 1 INITIALLY FROM BEFORE) .......... */
    i__2 = *n;
    for (j = 1; j <= i__2; ++j) {
	d__ = 0.f;
	if (isw == 2) {
	    goto L920;
	}
	if (alfi[j] != 0.f) {
	    goto L945;
	}

	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    if ((r__1 = z__[i__ + j * z_dim1], dabs(r__1)) > d__) {
		d__ = (r__2 = z__[i__ + j * z_dim1], dabs(r__2));
	    }
/* L890: */
	}

	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
/* L900: */
	    z__[i__ + j * z_dim1] /= d__;
	}

	goto L950;

L920:
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    r__ = (r__1 = z__[i__ + (j - 1) * z_dim1], dabs(r__1)) + (r__2 = 
		    z__[i__ + j * z_dim1], dabs(r__2));
	    if (r__ != 0.f) {
/* Computing 2nd power */
		r__1 = z__[i__ + (j - 1) * z_dim1] / r__;
/* Computing 2nd power */
		r__2 = z__[i__ + j * z_dim1] / r__;
		r__ *= sqrt(r__1 * r__1 + r__2 * r__2);
	    }
	    if (r__ > d__) {
		d__ = r__;
	    }
/* L930: */
	}

	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    z__[i__ + (j - 1) * z_dim1] /= d__;
	    z__[i__ + j * z_dim1] /= d__;
/* L940: */
	}

L945:
	isw = 3 - isw;
L950:
	;
    }

    return 0;
} /* qzvec_ */

