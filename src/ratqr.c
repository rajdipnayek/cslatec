/* ratqr.f -- translated by f2c (version 12.02.01).
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

static integer c__4 = 4;

/* DECK RATQR */
/* Subroutine */ int ratqr_(integer *n, real *eps1, real *d__, real *e, real *
	e2, integer *m, real *w, integer *ind, real *bd, logical *type__, 
	integer *idef, integer *ierr)
{
    /* Initialized data */

    static logical first = TRUE_;

    /* System generated locals */
    integer i__1, i__2;
    real r__1, r__2;

    /* Local variables */
    static real f;
    static integer i__, j, k;
    static real p, q, r__, s;
    static integer k1, ii, jj;
    static real ep, qp, err, tot;
    static integer jdef;
    static real delta;
    extern doublereal r1mach_(integer *);
    static real machep;

/* ***BEGIN PROLOGUE  RATQR */
/* ***PURPOSE  Compute the largest or smallest eigenvalues of a symmetric */
/*            tridiagonal matrix using the rational QR method with Newton */
/*            correction. */
/* ***LIBRARY   SLATEC (EISPACK) */
/* ***CATEGORY  D4A5, D4C2A */
/* ***TYPE      SINGLE PRECISION (RATQR-S) */
/* ***KEYWORDS  EIGENVALUES, EIGENVECTORS, EISPACK */
/* ***AUTHOR  Smith, B. T., et al. */
/* ***DESCRIPTION */

/*     This subroutine is a translation of the ALGOL procedure RATQR, */
/*     NUM. MATH. 11, 264-272(1968) by REINSCH and BAUER. */
/*     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 257-265(1971). */

/*     This subroutine finds the algebraically smallest or largest */
/*     eigenvalues of a SYMMETRIC TRIDIAGONAL matrix by the */
/*     rational QR method with Newton corrections. */

/*     On Input */

/*        N is the order of the matrix.  N is an INTEGER variable. */

/*        EPS1 is a theoretical absolute error tolerance for the */
/*          computed eigenvalues.  If the input EPS1 is non-positive, or */
/*          indeed smaller than its default value, it is reset at each */
/*          iteration to the respective default value, namely, the */
/*          product of the relative machine precision and the magnitude */
/*          of the current eigenvalue iterate.  The theoretical absolute */
/*          error in the K-th eigenvalue is usually not greater than */
/*          K times EPS1.  EPS1 is a REAL variable. */

/*        D contains the diagonal elements of the symmetric tridiagonal */
/*          matrix.  D is a one-dimensional REAL array, dimensioned D(N). */

/*        E contains the subdiagonal elements of the symmetric */
/*          tridiagonal matrix in its last N-1 positions.  E(1) is */
/*          arbitrary.  E is a one-dimensional REAL array, dimensioned */
/*          E(N). */

/*        E2 contains the squares of the corresponding elements of E in */
/*          its last N-1 positions.  E2(1) is arbitrary.  E2 is a one- */
/*          dimensional REAL array, dimensioned E2(N). */

/*        M is the number of eigenvalues to be found.  M is an INTEGER */
/*          variable. */

/*        IDEF should be set to 1 if the input matrix is known to be */
/*          positive definite, to -1 if the input matrix is known to */
/*          be negative definite, and to 0 otherwise.  IDEF is an */
/*          INTEGER variable. */

/*        TYPE should be set to .TRUE. if the smallest eigenvalues are */
/*          to be found, and to .FALSE. if the largest eigenvalues are */
/*          to be found.  TYPE is a LOGICAL variable. */

/*     On Output */

/*        EPS1 is unaltered unless it has been reset to its */
/*          (last) default value. */

/*        D and E are unaltered (unless W overwrites D). */

/*        Elements of E2, corresponding to elements of E regarded as */
/*          negligible, have been replaced by zero causing the matrix */
/*          to split into a direct sum of submatrices.  E2(1) is set */
/*          to 0.0e0 if the smallest eigenvalues have been found, and */
/*          to 2.0e0 if the largest eigenvalues have been found.  E2 */
/*          is otherwise unaltered (unless overwritten by BD). */

/*        W contains the M algebraically smallest eigenvalues in */
/*          ascending order, or the M largest eigenvalues in descending */
/*          order.  If an error exit is made because of an incorrect */
/*          specification of IDEF, no eigenvalues are found.  If the */
/*          Newton iterates for a particular eigenvalue are not monotone, */
/*          the best estimate obtained is returned and IERR is set. */
/*          W is a one-dimensional REAL array, dimensioned W(N).  W need */
/*          not be distinct from D. */

/*        IND contains in its first M positions the submatrix indices */
/*          associated with the corresponding eigenvalues in W -- */
/*          1 for eigenvalues belonging to the first submatrix from */
/*          the top, 2 for those belonging to the second submatrix, etc. */
/*          IND is an one-dimensional INTEGER array, dimensioned IND(N). */

/*        BD contains refined bounds for the theoretical errors of the */
/*          corresponding eigenvalues in W.  These bounds are usually */
/*          within the tolerance specified by EPS1.  BD is a one- */
/*          dimensional REAL array, dimensioned BD(N).  BD need not be */
/*          distinct from E2. */

/*        IERR is an INTEGER flag set to */
/*          Zero       for normal return, */
/*          6*N+1      if  IDEF  is set to 1 and  TYPE  to .TRUE. */
/*                     when the matrix is NOT positive definite, or */
/*                     if  IDEF  is set to -1 and  TYPE  to .FALSE. */
/*                     when the matrix is NOT negative definite, */
/*                     no eigenvalues are computed, or */
/*                     M is greater than N, */
/*          5*N+K      if successive iterates to the K-th eigenvalue */
/*                     are NOT monotone increasing, where K refers */
/*                     to the last such occurrence. */

/*     Note that subroutine TRIDIB is generally faster and more */
/*     accurate than RATQR if the eigenvalues are clustered. */

/*     Questions and comments should be directed to B. S. Garbow, */
/*     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY */
/*     ------------------------------------------------------------------ */

/* ***REFERENCES  B. T. Smith, J. M. Boyle, J. J. Dongarra, B. S. Garbow, */
/*                 Y. Ikebe, V. C. Klema and C. B. Moler, Matrix Eigen- */
/*                 system Routines - EISPACK Guide, Springer-Verlag, */
/*                 1976. */
/* ***ROUTINES CALLED  R1MACH */
/* ***REVISION HISTORY  (YYMMDD) */
/*   760101  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890831  Modified array declarations.  (WRB) */
/*   890831  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   920501  Reformatted the REFERENCES section.  (WRB) */
/* ***END PROLOGUE  RATQR */


    /* Parameter adjustments */
    --bd;
    --ind;
    --w;
    --e2;
    --e;
    --d__;

    /* Function Body */
/* ***FIRST EXECUTABLE STATEMENT  RATQR */
    if (first) {
	machep = r1mach_(&c__4);
    }
    first = FALSE_;

    *ierr = 0;
    jdef = *idef;
/*     .......... COPY D ARRAY INTO W .......... */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L20: */
	w[i__] = d__[i__];
    }

    if (*type__) {
	goto L40;
    }
    j = 1;
    goto L400;
L40:
    err = 0.f;
    s = 0.f;
/*     .......... LOOK FOR SMALL SUB-DIAGONAL ENTRIES AND DEFINE */
/*                INITIAL SHIFT FROM LOWER GERSCHGORIN BOUND. */
/*                COPY E2 ARRAY INTO BD .......... */
    tot = w[1];
    q = 0.f;
    j = 0;

    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	p = q;
	if (i__ == 1) {
	    goto L60;
	}
	if (p > machep * ((r__1 = d__[i__], dabs(r__1)) + (r__2 = d__[i__ - 1]
		, dabs(r__2)))) {
	    goto L80;
	}
L60:
	e2[i__] = 0.f;
L80:
	bd[i__] = e2[i__];
/*     .......... COUNT ALSO IF ELEMENT OF E2 HAS UNDERFLOWED .......... */
	if (e2[i__] == 0.f) {
	    ++j;
	}
	ind[i__] = j;
	q = 0.f;
	if (i__ != *n) {
	    q = (r__1 = e[i__ + 1], dabs(r__1));
	}
/* Computing MIN */
	r__1 = w[i__] - p - q;
	tot = dmin(r__1,tot);
/* L100: */
    }

    if (jdef == 1 && tot < 0.f) {
	goto L140;
    }

    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L110: */
	w[i__] -= tot;
    }

    goto L160;
L140:
    tot = 0.f;

L160:
    i__1 = *m;
    for (k = 1; k <= i__1; ++k) {
/*     .......... NEXT QR TRANSFORMATION .......... */
L180:
	tot += s;
	delta = w[*n] - s;
	i__ = *n;
	f = (r__1 = machep * tot, dabs(r__1));
	if (*eps1 < f) {
	    *eps1 = f;
	}
	if (delta > *eps1) {
	    goto L190;
	}
	if (delta < -(*eps1)) {
	    goto L1000;
	}
	goto L300;
/*     .......... REPLACE SMALL SUB-DIAGONAL SQUARES BY ZERO */
/*                TO REDUCE THE INCIDENCE OF UNDERFLOWS .......... */
L190:
	if (k == *n) {
	    goto L210;
	}
	k1 = k + 1;
	i__2 = *n;
	for (j = k1; j <= i__2; ++j) {
/* Computing 2nd power */
	    r__1 = machep * (w[j] + w[j - 1]);
	    if (bd[j] <= r__1 * r__1) {
		bd[j] = 0.f;
	    }
/* L200: */
	}

L210:
	f = bd[*n] / delta;
	qp = delta + f;
	p = 1.f;
	if (k == *n) {
	    goto L260;
	}
	k1 = *n - k;
/*     .......... FOR I=N-1 STEP -1 UNTIL K DO -- .......... */
	i__2 = k1;
	for (ii = 1; ii <= i__2; ++ii) {
	    i__ = *n - ii;
	    q = w[i__] - s - f;
	    r__ = q / qp;
	    p = p * r__ + 1.f;
	    ep = f * r__;
	    w[i__ + 1] = qp + ep;
	    delta = q - ep;
	    if (delta > *eps1) {
		goto L220;
	    }
	    if (delta < -(*eps1)) {
		goto L1000;
	    }
	    goto L300;
L220:
	    f = bd[i__] / q;
	    qp = delta + f;
	    bd[i__ + 1] = qp * ep;
/* L240: */
	}

L260:
	w[k] = qp;
	s = qp / p;
	if (tot + s > tot) {
	    goto L180;
	}
/*     .......... SET ERROR -- IRREGULAR END OF ITERATION. */
/*                DEFLATE MINIMUM DIAGONAL ELEMENT .......... */
	*ierr = *n * 5 + k;
	s = 0.f;
	delta = qp;

	i__2 = *n;
	for (j = k; j <= i__2; ++j) {
	    if (w[j] > delta) {
		goto L280;
	    }
	    i__ = j;
	    delta = w[j];
L280:
	    ;
	}
/*     .......... CONVERGENCE .......... */
L300:
	if (i__ < *n) {
	    bd[i__ + 1] = bd[i__] * f / qp;
	}
	ii = ind[i__];
	if (i__ == k) {
	    goto L340;
	}
	k1 = i__ - k;
/*     .......... FOR J=I-1 STEP -1 UNTIL K DO -- .......... */
	i__2 = k1;
	for (jj = 1; jj <= i__2; ++jj) {
	    j = i__ - jj;
	    w[j + 1] = w[j] - s;
	    bd[j + 1] = bd[j];
	    ind[j + 1] = ind[j];
/* L320: */
	}

L340:
	w[k] = tot;
	err += dabs(delta);
	bd[k] = err;
	ind[k] = ii;
/* L360: */
    }

    if (*type__) {
	goto L1001;
    }
    f = bd[1];
    e2[1] = 2.f;
    bd[1] = f;
    j = 2;
/*     .......... NEGATE ELEMENTS OF W FOR LARGEST VALUES .......... */
L400:
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L500: */
	w[i__] = -w[i__];
    }

    jdef = -jdef;
    switch (j) {
	case 1:  goto L40;
	case 2:  goto L1001;
    }
/*     .......... SET ERROR -- IDEF SPECIFIED INCORRECTLY .......... */
L1000:
    *ierr = *n * 6 + 1;
L1001:
    return 0;
} /* ratqr_ */

