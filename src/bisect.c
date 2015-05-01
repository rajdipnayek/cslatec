/* bisect.f -- translated by f2c (version 12.02.01).
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

/* DECK BISECT */
/* Subroutine */ int bisect_(integer *n, real *eps1, real *d__, real *e, real 
	*e2, real *lb, real *ub, integer *mm, integer *m, real *w, integer *
	ind, integer *ierr, real *rv4, real *rv5)
{
    /* Initialized data */

    static logical first = TRUE_;

    /* System generated locals */
    integer i__1, i__2;
    real r__1, r__2;

    /* Local variables */
    static integer i__, j, k, l, p, q, r__, s;
    static real u, v;
    static integer m1, m2;
    static real s1, t1, t2, s2, x0, x1;
    static integer ii;
    static real xu;
    static integer tag;
    extern doublereal r1mach_(integer *);
    static real machep;
    static integer isturm;

/* ***BEGIN PROLOGUE  BISECT */
/* ***PURPOSE  Compute the eigenvalues of a symmetric tridiagonal matrix */
/*            in a given interval using Sturm sequencing. */
/* ***LIBRARY   SLATEC (EISPACK) */
/* ***CATEGORY  D4A5, D4C2A */
/* ***TYPE      SINGLE PRECISION (BISECT-S) */
/* ***KEYWORDS  EIGENVALUES, EISPACK */
/* ***AUTHOR  Smith, B. T., et al. */
/* ***DESCRIPTION */

/*     This subroutine is a translation of the bisection technique */
/*     in the ALGOL procedure TRISTURM by Peters and Wilkinson. */
/*     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 418-439(1971). */

/*     This subroutine finds those eigenvalues of a TRIDIAGONAL */
/*     SYMMETRIC matrix which lie in a specified interval, */
/*     using bisection. */

/*     On INPUT */

/*        N is the order of the matrix.  N is an INTEGER variable. */

/*        EPS1 is an absolute error tolerance for the computed */
/*          eigenvalues.  If the input EPS1 is non-positive, */
/*          it is reset for each submatrix to a default value, */
/*          namely, minus the product of the relative machine */
/*          precision and the 1-norm of the submatrix. */
/*          EPS1 is a REAL variable. */

/*        D contains the diagonal elements of the input matrix. */
/*          D is a one-dimensional REAL array, dimensioned D(N). */

/*        E contains the subdiagonal elements of the input matrix */
/*          in its last N-1 positions.  E(1) is arbitrary. */
/*          E is a one-dimensional REAL array, dimensioned E(N). */

/*        E2 contains the squares of the corresponding elements of E. */
/*          E2(1) is arbitrary.  E2 is a one-dimensional REAL array, */
/*          dimensioned E2(N). */

/*        LB and UB define the interval to be searched for eigenvalues. */
/*          If LB is not less than UB, no eigenvalues will be found. */
/*          LB and UB are REAL variables. */

/*        MM should be set to an upper bound for the number of */
/*          eigenvalues in the interval.  WARNING - If more than */
/*          MM eigenvalues are determined to lie in the interval, */
/*          an error return is made with no eigenvalues found. */
/*          MM is an INTEGER variable. */

/*     On OUTPUT */

/*        EPS1 is unaltered unless it has been reset to its */
/*          (last) default value. */

/*        D and E are unaltered. */

/*        Elements of E2, corresponding to elements of E regarded */
/*          as negligible, have been replaced by zero causing the */
/*          matrix to split into a direct sum of submatrices. */
/*          E2(1) is also set to zero. */

/*        M is the number of eigenvalues determined to lie in (LB,UB). */
/*          M is an INTEGER variable. */

/*        W contains the M eigenvalues in ascending order. */
/*          W is a one-dimensional REAL array, dimensioned W(MM). */

/*        IND contains in its first M positions the submatrix indices */
/*          associated with the corresponding eigenvalues in W -- */
/*          1 for eigenvalues belonging to the first submatrix from */
/*          the top, 2 for those belonging to the second submatrix, etc. */
/*          IND is an one-dimensional INTEGER array, dimensioned IND(MM). */

/*        IERR is an INTEGER flag set to */
/*          Zero       for normal return, */
/*          3*N+1      if M exceeds MM.  In this case, M contains the */
/*                     number of eigenvalues determined to lie in */
/*                     (LB,UB). */

/*        RV4 and RV5 are one-dimensional REAL arrays used for temporary */
/*          storage, dimensioned RV4(N) and RV5(N). */

/*     The ALGOL procedure STURMCNT contained in TRISTURM */
/*     appears in BISECT in-line. */

/*     Note that subroutine TQL1 or IMTQL1 is generally faster than */
/*     BISECT, if more than N/4 eigenvalues are to be found. */

/*     Questions and comments should be directed to B. S. Garbow, */
/*     Applied Mathematics Division, ARGONNE NATIONAL LABORATORY */
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
/* ***END PROLOGUE  BISECT */


    /* Parameter adjustments */
    --rv5;
    --rv4;
    --ind;
    --w;
    --e2;
    --e;
    --d__;

    /* Function Body */
/* ***FIRST EXECUTABLE STATEMENT  BISECT */
    if (first) {
	machep = r1mach_(&c__4);
    }
    first = FALSE_;

    *ierr = 0;
    tag = 0;
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
    ++tag;
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
    m1 = p;
    m2 = p;
    rv5[p] = d__[p];
    goto L900;
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
/*     .......... ORDER EIGENVALUES TAGGED WITH THEIR */
/*                SUBMATRIX ASSOCIATIONS .......... */
L900:
    s = r__;
    r__ = r__ + m2 - m1 + 1;
    j = 1;
    k = m1;

    i__1 = r__;
    for (l = 1; l <= i__1; ++l) {
	if (j > s) {
	    goto L910;
	}
	if (k > m2) {
	    goto L940;
	}
	if (rv5[k] >= w[l]) {
	    goto L915;
	}

	i__2 = s;
	for (ii = j; ii <= i__2; ++ii) {
	    i__ = l + s - ii;
	    w[i__ + 1] = w[i__];
	    ind[i__ + 1] = ind[i__];
/* L905: */
	}

L910:
	w[l] = rv5[k];
	ind[l] = tag;
	++k;
	goto L920;
L915:
	++j;
L920:
	;
    }

L940:
    if (q < *n) {
	goto L100;
    }
    goto L1001;
/*     .......... SET ERROR -- UNDERESTIMATE OF NUMBER OF */
/*                EIGENVALUES IN INTERVAL .......... */
L980:
    *ierr = *n * 3 + 1;
L1001:
    *lb = t1;
    *ub = t2;
    return 0;
} /* bisect_ */

