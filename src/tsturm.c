/* tsturm.f -- translated by f2c (version 12.02.01).
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

/* DECK TSTURM */
/* Subroutine */ int tsturm_(integer *nm, integer *n, real *eps1, real *d__, 
	real *e, real *e2, real *lb, real *ub, integer *mm, integer *m, real *
	w, real *z__, integer *ierr, real *rv1, real *rv2, real *rv3, real *
	rv4, real *rv5, real *rv6)
{
    /* Initialized data */

    static logical first = TRUE_;

    /* System generated locals */
    integer z_dim1, z_offset, i__1, i__2, i__3;
    real r__1, r__2, r__3, r__4;

    /* Local variables */
    static integer i__, j, k, p, q, r__, s;
    static real u, v;
    static integer m1, m2;
    static real s1, t1, t2, s2, x0, x1;
    static integer ii, jj, ip;
    static real uk, xu;
    static integer its;
    static real eps2, eps3, eps4, norm;
    static integer group;
    extern doublereal r1mach_(integer *);
    static real machep;
    static integer isturm;

/* ***BEGIN PROLOGUE  TSTURM */
/* ***PURPOSE  Find those eigenvalues of a symmetric tridiagonal matrix */
/*            in a given interval and their associated eigenvectors by */
/*            Sturm sequencing. */
/* ***LIBRARY   SLATEC (EISPACK) */
/* ***CATEGORY  D4A5, D4C2A */
/* ***TYPE      SINGLE PRECISION (TSTURM-S) */
/* ***KEYWORDS  EIGENVALUES, EIGENVECTORS, EISPACK */
/* ***AUTHOR  Smith, B. T., et al. */
/* ***DESCRIPTION */

/*     This subroutine finds those eigenvalues of a TRIDIAGONAL */
/*     SYMMETRIC matrix which lie in a specified interval and their */
/*     associated eigenvectors, using bisection and inverse iteration. */

/*     On Input */

/*        NM must be set to the row dimension of the two-dimensional */
/*          array parameter, Z, as declared in the calling program */
/*          dimension statement.  NM is an INTEGER variable. */

/*        N is the order of the matrix.  N is an INTEGER variable. */
/*          N must be less than or equal to NM. */

/*        EPS1 is an absolute error tolerance for the computed eigen- */
/*          values.  It should be chosen so that the accuracy of these */
/*          eigenvalues is commensurate with relative perturbations of */
/*          the order of the relative machine precision in the matrix */
/*          elements.  If the input EPS1 is non-positive, it is reset */
/*          for each submatrix to a default value, namely, minus the */
/*          product of the relative machine precision and the 1-norm of */
/*          the submatrix.  EPS1 is a REAL variable. */

/*        D contains the diagonal elements of the symmetric tridiagonal */
/*          matrix.  D is a one-dimensional REAL array, dimensioned D(N). */

/*        E contains the subdiagonal elements of the symmetric */
/*          tridiagonal matrix in its last N-1 positions.  E(1) is */
/*          arbitrary.  E is a one-dimensional REAL array, dimensioned */
/*          E(N). */

/*        E2 contains the squares of the corresponding elements of E. */
/*          E2(1) is arbitrary.  E2 is a one-dimensional REAL array, */
/*          dimensioned E2(N). */

/*        LB and UB define the interval to be searched for eigenvalues. */
/*          If LB is not less than UB, no eigenvalues will be found. */
/*          LB and UB are REAL variables. */

/*        MM should be set to an upper bound for the number of */
/*          eigenvalues in the interval.  MM is an INTEGER variable. */
/*          WARNING -  If more than MM eigenvalues are determined to lie */
/*          in the interval, an error return is made with no values or */
/*          vectors found. */

/*     On Output */

/*        EPS1 is unaltered unless it has been reset to its */
/*          (last) default value. */

/*        D and E are unaltered. */

/*        Elements of E2, corresponding to elements of E regarded as */
/*          negligible, have been replaced by zero causing the matrix to */
/*          split into a direct sum of submatrices.  E2(1) is also set */
/*          to zero. */

/*        M is the number of eigenvalues determined to lie in (LB,UB). */
/*          M is an INTEGER variable. */

/*        W contains the M eigenvalues in ascending order if the matrix */
/*          does not split.  If the matrix splits, the eigenvalues are */
/*          in ascending order for each submatrix.  If a vector error */
/*          exit is made, W contains those values already found.  W is a */
/*          one-dimensional REAL array, dimensioned W(MM). */

/*        Z contains the associated set of orthonormal eigenvectors. */
/*          If an error exit is made, Z contains those vectors already */
/*          found.  Z is a one-dimensional REAL array, dimensioned */
/*          Z(NM,MM). */

/*        IERR is an INTEGER flag set to */
/*          Zero       for normal return, */
/*          3*N+1      if M exceeds MM no eigenvalues or eigenvectors */
/*                     are computed, */
/*          4*N+J      if the eigenvector corresponding to the J-th */
/*                     eigenvalue fails to converge in 5 iterations, then */
/*                     the eigenvalues and eigenvectors in W and Z should */
/*                     be correct for indices 1, 2, ..., J-1. */

/*        RV1, RV2, RV3, RV4, RV5, and RV6 are temporary storage arrays, */
/*          dimensioned RV1(N), RV2(N), RV3(N), RV4(N), RV5(N), and */
/*          RV6(N). */

/*     The ALGOL procedure STURMCNT contained in TRISTURM */
/*     appears in TSTURM in-line. */

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
/*   890531  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   920501  Reformatted the REFERENCES section.  (WRB) */
/* ***END PROLOGUE  TSTURM */


    /* Parameter adjustments */
    z_dim1 = *nm;
    z_offset = 1 + z_dim1;
    z__ -= z_offset;
    --d__;
    --e;
    --e2;
    --w;
    --rv1;
    --rv2;
    --rv3;
    --rv4;
    --rv5;
    --rv6;

    /* Function Body */
/* ***FIRST EXECUTABLE STATEMENT  TSTURM */
    if (first) {
	machep = r1mach_(&c__4);
    }
    first = FALSE_;

    *ierr = 0;
    t1 = *lb;
    t2 = *ub;
/*     .......... LOOK FOR SMALL SUB-DIAGONAL ENTRIES .......... */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (i__ == 1) {
	    goto L20;
	}
	s1 = (r__1 = d__[i__], dabs(r__1)) + (r__2 = d__[i__ - 1], dabs(r__2))
		;
	s2 = s1 + (r__1 = e[i__], dabs(r__1));
	if (s2 > s1) {
	    goto L40;
	}
L20:
	e2[i__] = 0.f;
L40:
	;
    }
/*     .......... DETERMINE THE NUMBER OF EIGENVALUES */
/*                IN THE INTERVAL .......... */
    p = 1;
    q = *n;
    x1 = *ub;
    isturm = 1;
    goto L320;
L60:
    *m = s;
    x1 = *lb;
    isturm = 2;
    goto L320;
L80:
    *m -= s;
    if (*m > *mm) {
	goto L980;
    }
    q = 0;
    r__ = 0;
/*     .......... ESTABLISH AND PROCESS NEXT SUBMATRIX, REFINING */
/*                INTERVAL BY THE GERSCHGORIN BOUNDS .......... */
L100:
    if (r__ == *m) {
	goto L1001;
    }
    p = q + 1;
    xu = d__[p];
    x0 = d__[p];
    u = 0.f;

    i__1 = *n;
    for (q = p; q <= i__1; ++q) {
	x1 = u;
	u = 0.f;
	v = 0.f;
	if (q == *n) {
	    goto L110;
	}
	u = (r__1 = e[q + 1], dabs(r__1));
	v = e2[q + 1];
L110:
/* Computing MIN */
	r__1 = d__[q] - (x1 + u);
	xu = dmin(r__1,xu);
/* Computing MAX */
	r__1 = d__[q] + (x1 + u);
	x0 = dmax(r__1,x0);
	if (v == 0.f) {
	    goto L140;
	}
/* L120: */
    }

L140:
/* Computing MAX */
    r__1 = dabs(xu), r__2 = dabs(x0);
    x1 = dmax(r__1,r__2) * machep;
    if (*eps1 <= 0.f) {
	*eps1 = -x1;
    }
    if (p != q) {
	goto L180;
    }
/*     .......... CHECK FOR ISOLATED ROOT WITHIN INTERVAL .......... */
    if (t1 > d__[p] || d__[p] >= t2) {
	goto L940;
    }
    ++r__;

    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L160: */
	z__[i__ + r__ * z_dim1] = 0.f;
    }

    w[r__] = d__[p];
    z__[p + r__ * z_dim1] = 1.f;
    goto L940;
L180:
    x1 *= q - p + 1;
/* Computing MAX */
    r__1 = t1, r__2 = xu - x1;
    *lb = dmax(r__1,r__2);
/* Computing MIN */
    r__1 = t2, r__2 = x0 + x1;
    *ub = dmin(r__1,r__2);
    x1 = *lb;
    isturm = 3;
    goto L320;
L200:
    m1 = s + 1;
    x1 = *ub;
    isturm = 4;
    goto L320;
L220:
    m2 = s;
    if (m1 > m2) {
	goto L940;
    }
/*     .......... FIND ROOTS BY BISECTION .......... */
    x0 = *ub;
    isturm = 5;

    i__1 = m2;
    for (i__ = m1; i__ <= i__1; ++i__) {
	rv5[i__] = *ub;
	rv4[i__] = *lb;
/* L240: */
    }
/*     .......... LOOP FOR K-TH EIGENVALUE */
/*                FOR K=M2 STEP -1 UNTIL M1 DO -- */
/*                (-DO- NOT USED TO LEGALIZE -COMPUTED GO TO-) .......... */
    k = m2;
L250:
    xu = *lb;
/*     .......... FOR I=K STEP -1 UNTIL M1 DO -- .......... */
    i__1 = k;
    for (ii = m1; ii <= i__1; ++ii) {
	i__ = m1 + k - ii;
	if (xu >= rv4[i__]) {
	    goto L260;
	}
	xu = rv4[i__];
	goto L280;
L260:
	;
    }

L280:
    if (x0 > rv5[k]) {
	x0 = rv5[k];
    }
/*     .......... NEXT BISECTION STEP .......... */
L300:
    x1 = (xu + x0) * .5f;
    s1 = (dabs(xu) + dabs(x0) + dabs(*eps1)) * 2.f;
    s2 = s1 + (r__1 = x0 - xu, dabs(r__1));
    if (s2 == s1) {
	goto L420;
    }
/*     .......... IN-LINE PROCEDURE FOR STURM SEQUENCE .......... */
L320:
    s = p - 1;
    u = 1.f;

    i__1 = q;
    for (i__ = p; i__ <= i__1; ++i__) {
	if (u != 0.f) {
	    goto L325;
	}
	v = (r__1 = e[i__], dabs(r__1)) / machep;
	if (e2[i__] == 0.f) {
	    v = 0.f;
	}
	goto L330;
L325:
	v = e2[i__] / u;
L330:
	u = d__[i__] - x1 - v;
	if (u < 0.f) {
	    ++s;
	}
/* L340: */
    }

    switch (isturm) {
	case 1:  goto L60;
	case 2:  goto L80;
	case 3:  goto L200;
	case 4:  goto L220;
	case 5:  goto L360;
    }
/*     .......... REFINE INTERVALS .......... */
L360:
    if (s >= k) {
	goto L400;
    }
    xu = x1;
    if (s >= m1) {
	goto L380;
    }
    rv4[m1] = x1;
    goto L300;
L380:
    rv4[s + 1] = x1;
    if (rv5[s] > x1) {
	rv5[s] = x1;
    }
    goto L300;
L400:
    x0 = x1;
    goto L300;
/*     .......... K-TH EIGENVALUE FOUND .......... */
L420:
    rv5[k] = x1;
    --k;
    if (k >= m1) {
	goto L250;
    }
/*     .......... FIND VECTORS BY INVERSE ITERATION .......... */
    norm = (r__1 = d__[p], dabs(r__1));
    ip = p + 1;

    i__1 = q;
    for (i__ = ip; i__ <= i__1; ++i__) {
/* L500: */
/* Computing MAX */
	r__3 = norm, r__4 = (r__1 = d__[i__], dabs(r__1)) + (r__2 = e[i__], 
		dabs(r__2));
	norm = dmax(r__3,r__4);
    }
/*     .......... EPS2 IS THE CRITERION FOR GROUPING, */
/*                EPS3 REPLACES ZERO PIVOTS AND EQUAL */
/*                ROOTS ARE MODIFIED BY EPS3, */
/*                EPS4 IS TAKEN VERY SMALL TO AVOID OVERFLOW .......... */
    eps2 = norm * .001f;
    uk = sqrt((real) (q - p + 5));
    eps3 = uk * machep * norm;
    eps4 = uk * eps3;
    uk = eps4 / sqrt(uk);
    group = 0;
    s = p;

    i__1 = m2;
    for (k = m1; k <= i__1; ++k) {
	++r__;
	its = 1;
	w[r__] = rv5[k];
	x1 = rv5[k];
/*     .......... LOOK FOR CLOSE OR COINCIDENT ROOTS .......... */
	if (k == m1) {
	    goto L520;
	}
	if (x1 - x0 >= eps2) {
	    group = -1;
	}
	++group;
	if (x1 <= x0) {
	    x1 = x0 + eps3;
	}
/*     .......... ELIMINATION WITH INTERCHANGES AND */
/*                INITIALIZATION OF VECTOR .......... */
L520:
	v = 0.f;

	i__2 = q;
	for (i__ = p; i__ <= i__2; ++i__) {
	    rv6[i__] = uk;
	    if (i__ == p) {
		goto L560;
	    }
	    if ((r__1 = e[i__], dabs(r__1)) < dabs(u)) {
		goto L540;
	    }
	    xu = u / e[i__];
	    rv4[i__] = xu;
	    rv1[i__ - 1] = e[i__];
	    rv2[i__ - 1] = d__[i__] - x1;
	    rv3[i__ - 1] = 0.f;
	    if (i__ != q) {
		rv3[i__ - 1] = e[i__ + 1];
	    }
	    u = v - xu * rv2[i__ - 1];
	    v = -xu * rv3[i__ - 1];
	    goto L580;
L540:
	    xu = e[i__] / u;
	    rv4[i__] = xu;
	    rv1[i__ - 1] = u;
	    rv2[i__ - 1] = v;
	    rv3[i__ - 1] = 0.f;
L560:
	    u = d__[i__] - x1 - xu * v;
	    if (i__ != q) {
		v = e[i__ + 1];
	    }
L580:
	    ;
	}

	if (u == 0.f) {
	    u = eps3;
	}
	rv1[q] = u;
	rv2[q] = 0.f;
	rv3[q] = 0.f;
/*     .......... BACK SUBSTITUTION */
/*                FOR I=Q STEP -1 UNTIL P DO -- .......... */
L600:
	i__2 = q;
	for (ii = p; ii <= i__2; ++ii) {
	    i__ = p + q - ii;
	    rv6[i__] = (rv6[i__] - u * rv2[i__] - v * rv3[i__]) / rv1[i__];
	    v = u;
	    u = rv6[i__];
/* L620: */
	}
/*     .......... ORTHOGONALIZE WITH RESPECT TO PREVIOUS */
/*                MEMBERS OF GROUP .......... */
	if (group == 0) {
	    goto L700;
	}

	i__2 = group;
	for (jj = 1; jj <= i__2; ++jj) {
	    j = r__ - group - 1 + jj;
	    xu = 0.f;

	    i__3 = q;
	    for (i__ = p; i__ <= i__3; ++i__) {
/* L640: */
		xu += rv6[i__] * z__[i__ + j * z_dim1];
	    }

	    i__3 = q;
	    for (i__ = p; i__ <= i__3; ++i__) {
/* L660: */
		rv6[i__] -= xu * z__[i__ + j * z_dim1];
	    }

/* L680: */
	}

L700:
	norm = 0.f;

	i__2 = q;
	for (i__ = p; i__ <= i__2; ++i__) {
/* L720: */
	    norm += (r__1 = rv6[i__], dabs(r__1));
	}

	if (norm >= 1.f) {
	    goto L840;
	}
/*     .......... FORWARD SUBSTITUTION .......... */
	if (its == 5) {
	    goto L960;
	}
	if (norm != 0.f) {
	    goto L740;
	}
	rv6[s] = eps4;
	++s;
	if (s > q) {
	    s = p;
	}
	goto L780;
L740:
	xu = eps4 / norm;

	i__2 = q;
	for (i__ = p; i__ <= i__2; ++i__) {
/* L760: */
	    rv6[i__] *= xu;
	}
/*     .......... ELIMINATION OPERATIONS ON NEXT VECTOR */
/*                ITERATE .......... */
L780:
	i__2 = q;
	for (i__ = ip; i__ <= i__2; ++i__) {
	    u = rv6[i__];
/*     .......... IF RV1(I-1) .EQ. E(I), A ROW INTERCHANGE */
/*                WAS PERFORMED EARLIER IN THE */
/*                TRIANGULARIZATION PROCESS .......... */
	    if (rv1[i__ - 1] != e[i__]) {
		goto L800;
	    }
	    u = rv6[i__ - 1];
	    rv6[i__ - 1] = rv6[i__];
L800:
	    rv6[i__] = u - rv4[i__] * rv6[i__ - 1];
/* L820: */
	}

	++its;
	goto L600;
/*     .......... NORMALIZE SO THAT SUM OF SQUARES IS */
/*                1 AND EXPAND TO FULL ORDER .......... */
L840:
	u = 0.f;

	i__2 = q;
	for (i__ = p; i__ <= i__2; ++i__) {
/* L860: */
/* Computing 2nd power */
	    r__1 = rv6[i__];
	    u += r__1 * r__1;
	}

	xu = 1.f / sqrt(u);

	i__2 = *n;
	for (i__ = 1; i__ <= i__2; ++i__) {
/* L880: */
	    z__[i__ + r__ * z_dim1] = 0.f;
	}

	i__2 = q;
	for (i__ = p; i__ <= i__2; ++i__) {
/* L900: */
	    z__[i__ + r__ * z_dim1] = rv6[i__] * xu;
	}

	x0 = x1;
/* L920: */
    }

L940:
    if (q < *n) {
	goto L100;
    }
    goto L1001;
/*     .......... SET ERROR -- NON-CONVERGED EIGENVECTOR .......... */
L960:
    *ierr = (*n << 2) + r__;
    goto L1001;
/*     .......... SET ERROR -- UNDERESTIMATE OF NUMBER OF */
/*                EIGENVALUES IN INTERVAL .......... */
L980:
    *ierr = *n * 3 + 1;
L1001:
    *lb = t1;
    *ub = t2;
    return 0;
} /* tsturm_ */

