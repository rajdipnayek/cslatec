/* qzval.f -- translated by f2c (version 12.02.01).
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

/* DECK QZVAL */
/* Subroutine */ int qzval_(integer *nm, integer *n, real *a, real *b, real *
	alfr, real *alfi, real *beta, logical *matz, real *z__)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, z_dim1, z_offset, i__1, i__2;
    real r__1, r__2, r__3, r__4;

    /* Local variables */
    static real c__, d__, e;
    static integer i__, j;
    static real r__, s, t, a1, a2, u1, u2, v1, v2, a11, a12, a21, a22, b11, 
	    b12, b22, di, ei;
    static integer na;
    static real an, bn;
    static integer en;
    static real cq, dr;
    static integer nn;
    static real cz, ti, tr, a1i, a2i, a11i, a12i, a22i, a11r, a12r, a22r, sqi,
	     ssi;
    static integer isw;
    static real sqr, szi, ssr, szr, epsb;

/* ***BEGIN PROLOGUE  QZVAL */
/* ***PURPOSE  The third step of the QZ algorithm for generalized */
/*            eigenproblems.  Accepts a pair of real matrices, one in */
/*            quasi-triangular form and the other in upper triangular */
/*            form and computes the eigenvalues of the associated */
/*            eigenproblem.  Usually preceded by QZHES, QZIT, and */
/*            followed by QZVEC. */
/* ***LIBRARY   SLATEC (EISPACK) */
/* ***CATEGORY  D4C2C */
/* ***TYPE      SINGLE PRECISION (QZVAL-S) */
/* ***KEYWORDS  EIGENVALUES, EIGENVECTORS, EISPACK */
/* ***AUTHOR  Smith, B. T., et al. */
/* ***DESCRIPTION */

/*     This subroutine is the third step of the QZ algorithm */
/*     for solving generalized matrix eigenvalue problems, */
/*     SIAM J. NUMER. ANAL. 10, 241-256(1973) by MOLER and STEWART. */

/*     This subroutine accepts a pair of REAL matrices, one of them */
/*     in quasi-triangular form and the other in upper triangular form. */
/*     It reduces the quasi-triangular matrix further, so that any */
/*     remaining 2-by-2 blocks correspond to pairs of complex */
/*     eigenvalues, and returns quantities whose ratios give the */
/*     generalized eigenvalues.  It is usually preceded by  QZHES */
/*     and  QZIT  and may be followed by  QZVEC. */

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

/*        MATZ should be set to .TRUE. if the right hand transformations */
/*          are to be accumulated for later use in computing */
/*          eigenvectors, and to .FALSE. otherwise.  MATZ is a LOGICAL */
/*          variable. */

/*        Z contains, if MATZ has been set to .TRUE., the transformation */
/*          matrix produced in the reductions by  QZHES  and  QZIT,  if */
/*          performed, or else the identity matrix.  If MATZ has been set */
/*          to .FALSE., Z is not referenced.  Z is a two-dimensional REAL */
/*          array, dimensioned Z(NM,N). */

/*     On Output */

/*        A has been reduced further to a quasi-triangular matrix in */
/*          which all nonzero subdiagonal elements correspond to pairs */
/*          of complex eigenvalues. */

/*        B is still in upper triangular form, although its elements */
/*          have been altered.  B(N,1) is unaltered. */

/*        ALFR and ALFI contain the real and imaginary parts of the */
/*          diagonal elements of the triangular matrix that would be */
/*          obtained if A were reduced completely to triangular form */
/*          by unitary transformations.  Non-zero values of ALFI occur */
/*          in pairs, the first member positive and the second negative. */
/*          ALFR and ALFI are one-dimensional REAL arrays, dimensioned */
/*          ALFR(N) and ALFI(N). */

/*        BETA contains the diagonal elements of the corresponding B, */
/*          normalized to be real and non-negative.  The generalized */
/*          eigenvalues are then the ratios ((ALFR+I*ALFI)/BETA). */
/*          BETA is a one-dimensional REAL array, dimensioned BETA(N). */

/*        Z contains the product of the right hand transformations */
/*          (for all three steps) if MATZ has been set to .TRUE. */

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
/* ***END PROLOGUE  QZVAL */


/* ***FIRST EXECUTABLE STATEMENT  QZVAL */
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
/*     .......... FIND EIGENVALUES OF QUASI-TRIANGULAR MATRICES. */
/*                FOR EN=N STEP -1 UNTIL 1 DO -- .......... */
    i__1 = *n;
    for (nn = 1; nn <= i__1; ++nn) {
	en = *n + 1 - nn;
	na = en - 1;
	if (isw == 2) {
	    goto L505;
	}
	if (en == 1) {
	    goto L410;
	}
	if (a[en + na * a_dim1] != 0.f) {
	    goto L420;
	}
/*     .......... 1-BY-1 BLOCK, ONE REAL ROOT .......... */
L410:
	alfr[en] = a[en + en * a_dim1];
	if (b[en + en * b_dim1] < 0.f) {
	    alfr[en] = -alfr[en];
	}
	beta[en] = (r__1 = b[en + en * b_dim1], dabs(r__1));
	alfi[en] = 0.f;
	goto L510;
/*     .......... 2-BY-2 BLOCK .......... */
L420:
	if ((r__1 = b[na + na * b_dim1], dabs(r__1)) <= epsb) {
	    goto L455;
	}
	if ((r__1 = b[en + en * b_dim1], dabs(r__1)) > epsb) {
	    goto L430;
	}
	a1 = a[en + en * a_dim1];
	a2 = a[en + na * a_dim1];
	bn = 0.f;
	goto L435;
L430:
	an = (r__1 = a[na + na * a_dim1], dabs(r__1)) + (r__2 = a[na + en * 
		a_dim1], dabs(r__2)) + (r__3 = a[en + na * a_dim1], dabs(r__3)
		) + (r__4 = a[en + en * a_dim1], dabs(r__4));
	bn = (r__1 = b[na + na * b_dim1], dabs(r__1)) + (r__2 = b[na + en * 
		b_dim1], dabs(r__2)) + (r__3 = b[en + en * b_dim1], dabs(r__3)
		);
	a11 = a[na + na * a_dim1] / an;
	a12 = a[na + en * a_dim1] / an;
	a21 = a[en + na * a_dim1] / an;
	a22 = a[en + en * a_dim1] / an;
	b11 = b[na + na * b_dim1] / bn;
	b12 = b[na + en * b_dim1] / bn;
	b22 = b[en + en * b_dim1] / bn;
	e = a11 / b11;
	ei = a22 / b22;
	s = a21 / (b11 * b22);
	t = (a22 - e * b22) / b22;
	if (dabs(e) <= dabs(ei)) {
	    goto L431;
	}
	e = ei;
	t = (a11 - e * b11) / b11;
L431:
	c__ = (t - s * b12) * .5f;
	d__ = c__ * c__ + s * (a12 - e * b12);
	if (d__ < 0.f) {
	    goto L480;
	}
/*     .......... TWO REAL ROOTS. */
/*                ZERO BOTH A(EN,NA) AND B(EN,NA) .......... */
	r__1 = sqrt(d__);
	e += c__ + r_sign(&r__1, &c__);
	a11 -= e * b11;
	a12 -= e * b12;
	a22 -= e * b22;
	if (dabs(a11) + dabs(a12) < dabs(a21) + dabs(a22)) {
	    goto L432;
	}
	a1 = a12;
	a2 = a11;
	goto L435;
L432:
	a1 = a22;
	a2 = a21;
/*     .......... CHOOSE AND APPLY REAL Z .......... */
L435:
	s = dabs(a1) + dabs(a2);
	u1 = a1 / s;
	u2 = a2 / s;
	r__1 = sqrt(u1 * u1 + u2 * u2);
	r__ = r_sign(&r__1, &u1);
	v1 = -(u1 + r__) / r__;
	v2 = -u2 / r__;
	u2 = v2 / v1;

	i__2 = en;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    t = a[i__ + en * a_dim1] + u2 * a[i__ + na * a_dim1];
	    a[i__ + en * a_dim1] += t * v1;
	    a[i__ + na * a_dim1] += t * v2;
	    t = b[i__ + en * b_dim1] + u2 * b[i__ + na * b_dim1];
	    b[i__ + en * b_dim1] += t * v1;
	    b[i__ + na * b_dim1] += t * v2;
/* L440: */
	}

	if (! (*matz)) {
	    goto L450;
	}

	i__2 = *n;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    t = z__[i__ + en * z_dim1] + u2 * z__[i__ + na * z_dim1];
	    z__[i__ + en * z_dim1] += t * v1;
	    z__[i__ + na * z_dim1] += t * v2;
/* L445: */
	}

L450:
	if (bn == 0.f) {
	    goto L475;
	}
	if (an < dabs(e) * bn) {
	    goto L455;
	}
	a1 = b[na + na * b_dim1];
	a2 = b[en + na * b_dim1];
	goto L460;
L455:
	a1 = a[na + na * a_dim1];
	a2 = a[en + na * a_dim1];
/*     .......... CHOOSE AND APPLY REAL Q .......... */
L460:
	s = dabs(a1) + dabs(a2);
	if (s == 0.f) {
	    goto L475;
	}
	u1 = a1 / s;
	u2 = a2 / s;
	r__1 = sqrt(u1 * u1 + u2 * u2);
	r__ = r_sign(&r__1, &u1);
	v1 = -(u1 + r__) / r__;
	v2 = -u2 / r__;
	u2 = v2 / v1;

	i__2 = *n;
	for (j = na; j <= i__2; ++j) {
	    t = a[na + j * a_dim1] + u2 * a[en + j * a_dim1];
	    a[na + j * a_dim1] += t * v1;
	    a[en + j * a_dim1] += t * v2;
	    t = b[na + j * b_dim1] + u2 * b[en + j * b_dim1];
	    b[na + j * b_dim1] += t * v1;
	    b[en + j * b_dim1] += t * v2;
/* L470: */
	}

L475:
	a[en + na * a_dim1] = 0.f;
	b[en + na * b_dim1] = 0.f;
	alfr[na] = a[na + na * a_dim1];
	alfr[en] = a[en + en * a_dim1];
	if (b[na + na * b_dim1] < 0.f) {
	    alfr[na] = -alfr[na];
	}
	if (b[en + en * b_dim1] < 0.f) {
	    alfr[en] = -alfr[en];
	}
	beta[na] = (r__1 = b[na + na * b_dim1], dabs(r__1));
	beta[en] = (r__1 = b[en + en * b_dim1], dabs(r__1));
	alfi[en] = 0.f;
	alfi[na] = 0.f;
	goto L505;
/*     .......... TWO COMPLEX ROOTS .......... */
L480:
	e += c__;
	ei = sqrt(-d__);
	a11r = a11 - e * b11;
	a11i = ei * b11;
	a12r = a12 - e * b12;
	a12i = ei * b12;
	a22r = a22 - e * b22;
	a22i = ei * b22;
	if (dabs(a11r) + dabs(a11i) + dabs(a12r) + dabs(a12i) < dabs(a21) + 
		dabs(a22r) + dabs(a22i)) {
	    goto L482;
	}
	a1 = a12r;
	a1i = a12i;
	a2 = -a11r;
	a2i = -a11i;
	goto L485;
L482:
	a1 = a22r;
	a1i = a22i;
	a2 = -a21;
	a2i = 0.f;
/*     .......... CHOOSE COMPLEX Z .......... */
L485:
	cz = sqrt(a1 * a1 + a1i * a1i);
	if (cz == 0.f) {
	    goto L487;
	}
	szr = (a1 * a2 + a1i * a2i) / cz;
	szi = (a1 * a2i - a1i * a2) / cz;
	r__ = sqrt(cz * cz + szr * szr + szi * szi);
	cz /= r__;
	szr /= r__;
	szi /= r__;
	goto L490;
L487:
	szr = 1.f;
	szi = 0.f;
L490:
	if (an < (dabs(e) + ei) * bn) {
	    goto L492;
	}
	a1 = cz * b11 + szr * b12;
	a1i = szi * b12;
	a2 = szr * b22;
	a2i = szi * b22;
	goto L495;
L492:
	a1 = cz * a11 + szr * a12;
	a1i = szi * a12;
	a2 = cz * a21 + szr * a22;
	a2i = szi * a22;
/*     .......... CHOOSE COMPLEX Q .......... */
L495:
	cq = sqrt(a1 * a1 + a1i * a1i);
	if (cq == 0.f) {
	    goto L497;
	}
	sqr = (a1 * a2 + a1i * a2i) / cq;
	sqi = (a1 * a2i - a1i * a2) / cq;
	r__ = sqrt(cq * cq + sqr * sqr + sqi * sqi);
	cq /= r__;
	sqr /= r__;
	sqi /= r__;
	goto L500;
L497:
	sqr = 1.f;
	sqi = 0.f;
/*     .......... COMPUTE DIAGONAL ELEMENTS THAT WOULD RESULT */
/*                IF TRANSFORMATIONS WERE APPLIED .......... */
L500:
	ssr = sqr * szr + sqi * szi;
	ssi = sqr * szi - sqi * szr;
	i__ = 1;
	tr = cq * cz * a11 + cq * szr * a12 + sqr * cz * a21 + ssr * a22;
	ti = cq * szi * a12 - sqi * cz * a21 + ssi * a22;
	dr = cq * cz * b11 + cq * szr * b12 + ssr * b22;
	di = cq * szi * b12 + ssi * b22;
	goto L503;
L502:
	i__ = 2;
	tr = ssr * a11 - sqr * cz * a12 - cq * szr * a21 + cq * cz * a22;
	ti = -ssi * a11 - sqi * cz * a12 + cq * szi * a21;
	dr = ssr * b11 - sqr * cz * b12 + cq * cz * b22;
	di = -ssi * b11 - sqi * cz * b12;
L503:
	t = ti * dr - tr * di;
	j = na;
	if (t < 0.f) {
	    j = en;
	}
	r__ = sqrt(dr * dr + di * di);
	beta[j] = bn * r__;
	alfr[j] = an * (tr * dr + ti * di) / r__;
	alfi[j] = an * t / r__;
	if (i__ == 1) {
	    goto L502;
	}
L505:
	isw = 3 - isw;
L510:
	;
    }

    return 0;
} /* qzval_ */

