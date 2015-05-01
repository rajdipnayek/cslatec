/* tridib.f -- translated by f2c (version 12.02.01).
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

/* DECK TRIDIB */
/* Subroutine */ int tridib_(integer *n, real *eps1, real *d__, real *e, real 
	*e2, real *lb, real *ub, integer *m11, integer *m, real *w, integer *
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
    static integer m22, ii;
    static real xu;
    static integer tag;
    extern doublereal r1mach_(integer *);
    static real machep;
    static integer isturm;

/* ***BEGIN PROLOGUE  TRIDIB */
/* ***PURPOSE  Compute the eigenvalues of a symmetric tridiagonal matrix */
/*            in a given interval using Sturm sequencing. */
/* ***LIBRARY   SLATEC (EISPACK) */
/* ***CATEGORY  D4A5, D4C2A */
/* ***TYPE      SINGLE PRECISION (TRIDIB-S) */
/* ***KEYWORDS  EIGENVALUES OF A REAL SYMMETRIC MATRIX, EISPACK */
/* ***AUTHOR  Smith, B. T., et al. */
/* ***DESCRIPTION */

/*     This subroutine is a translation of the ALGOL procedure BISECT, */
/*     NUM. MATH. 9, 386-393(1967) by Barth, Martin, and Wilkinson. */
/*     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 249-256(1971). */

/*     This subroutine finds those eigenvalues of a TRIDIAGONAL */
/*     SYMMETRIC matrix between specified boundary indices, */
/*     using bisection. */

/*     On Input */

/*        N is the order of the matrix.  N is an INTEGER variable. */

/*        EPS1 is an absolute error tolerance for the computed eigen- */
/*          values.  If the input EPS1 is non-positive, it is reset for */
/*          each submatrix to a default value, namely, minus the product */
/*          of the relative machine precision and the 1-norm of the */
/*          submatrix.  EPS1 is a REAL variable. */

/*        D contains the diagonal elements of the symmetric tridiagonal */
/*          matrix.  D is a one-dimensional REAL array, dimensioned D(N). */

/*        E contains the subdiagonal elements of the symmetric */
/*          tridiagonal matrix in its last N-1 positions.  E(1) is */
/*          arbitrary.  E is a one-dimensional REAL array, dimensioned */
/*          E(N). */

/*        E2 contains the squares of the corresponding elements of E. */
/*          E2(1) is arbitrary.  E2 is a one-dimensional REAL array, */
/*          dimensioned E2(N). */

/*        M11 specifies the lower boundary index for the set of desired */
/*          eigenvalues.  M11 is an INTEGER variable. */

/*        M specifies the number of eigenvalues desired.  The upper */
/*          boundary index M22 is then obtained as M22=M11+M-1. */
/*          M is an INTEGER variable. */

/*     On Output */

/*        EPS1 is unaltered unless it has been reset to its */
/*          (last) default value. */

/*        D and E are unaltered. */

/*        Elements of E2, corresponding to elements of E regarded */
/*          as negligible, have been replaced by zero causing the */
/*          matrix to split into a direct sum of submatrices. */
/*          E2(1) is also set to zero. */

/*        LB and UB define an interval containing exactly the desired */
/*          eigenvalues.  LB and UB are REAL variables. */

/*        W contains, in its first M positions, the eigenvalues */
/*          between indices M11 and M22 in ascending order. */
/*          W is a one-dimensional REAL array, dimensioned W(M). */

/*        IND contains in its first M positions the submatrix indices */
/*          associated with the corresponding eigenvalues in W -- */
/*          1 for eigenvalues belonging to the first submatrix from */
/*          the top, 2 for those belonging to the second submatrix, etc. */
/*          IND is an one-dimensional INTEGER array, dimensioned IND(M). */

/*        IERR is an INTEGER flag set to */
/*          Zero       for normal return, */
/*          3*N+1      if multiple eigenvalues at index M11 make */
/*                     unique selection of LB impossible, */
/*          3*N+2      if multiple eigenvalues at index M22 make */
/*                     unique selection of UB impossible. */

/*        RV4 and RV5 are one-dimensional REAL arrays used for temporary */
/*          storage of the lower and upper bounds for the eigenvalues in */
/*          the bisection process.  RV4 and RV5 are dimensioned RV4(N) */
/*          and RV5(N). */

/*     Note that subroutine TQL1, IMTQL1, or TQLRAT is generally faster */
/*     than TRIDIB, if more than N/4 eigenvalues are to be found. */

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
/* ***END PROLOGUE  TRIDIB */


    /* Parameter adjustments */
    --rv5;
    --rv4;
    --ind;
    --w;
    --e2;
    --e;
    --d__;

    /* Function Body */
/* ***FIRST EXECUTABLE STATEMENT  TRIDIB */
    if (first) {
	machep = r1mach_(&c__4);
    }
    first = FALSE_;

    *ierr = 0;
    tag = 0;
    xu = d__[1];
    x0 = d__[1];
    u = 0.f;
/*     .......... LOOK FOR SMALL SUB-DIAGONAL ENTRIES AND DETERMINE AN */
/*                INTERVAL CONTAINING ALL THE EIGENVALUES .......... */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	x1 = u;
	u = 0.f;
	if (i__ != *n) {
	    u = (r__1 = e[i__ + 1], dabs(r__1));
	}
/* Computing MIN */
	r__1 = d__[i__] - (x1 + u);
	xu = dmin(r__1,xu);
/* Computing MAX */
	r__1 = d__[i__] + (x1 + u);
	x0 = dmax(r__1,x0);
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

/* Computing MAX */
    r__1 = dabs(xu), r__2 = dabs(x0);
    x1 = dmax(r__1,r__2) * machep * *n;
    xu -= x1;
    t1 = xu;
    x0 += x1;
    t2 = x0;
/*     .......... DETERMINE AN INTERVAL CONTAINING EXACTLY */
/*                THE DESIRED EIGENVALUES .......... */
    p = 1;
    q = *n;
    m1 = *m11 - 1;
    if (m1 == 0) {
	goto L75;
    }
    isturm = 1;
L50:
    v = x1;
    x1 = xu + (x0 - xu) * .5f;
    if (x1 == v) {
	goto L980;
    }
    goto L320;
L60:
    if ((i__1 = s - m1) < 0) {
	goto L65;
    } else if (i__1 == 0) {
	goto L73;
    } else {
	goto L70;
    }
L65:
    xu = x1;
    goto L50;
L70:
    x0 = x1;
    goto L50;
L73:
    xu = x1;
    t1 = x1;
L75:
    m22 = m1 + *m;
    if (m22 == *n) {
	goto L90;
    }
    x0 = t2;
    isturm = 2;
    goto L50;
L80:
    if ((i__1 = s - m22) < 0) {
	goto L65;
    } else if (i__1 == 0) {
	goto L85;
    } else {
	goto L70;
    }
L85:
    t2 = x1;
L90:
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
    s1 = dabs(xu) + dabs(x0) + dabs(*eps1);
    s2 = s1 + (r__1 = x0 - xu, dabs(r__1)) / 2.f;
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
/*     .......... SET ERROR -- INTERVAL CANNOT BE FOUND CONTAINING */
/*                EXACTLY THE DESIRED EIGENVALUES .......... */
L980:
    *ierr = *n * 3 + isturm;
L1001:
    *lb = t1;
    *ub = t2;
    return 0;
} /* tridib_ */

