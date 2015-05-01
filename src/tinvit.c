/* tinvit.f -- translated by f2c (version 12.02.01).
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

/* DECK TINVIT */
/* Subroutine */ int tinvit_(integer *nm, integer *n, real *d__, real *e, 
	real *e2, integer *m, real *w, integer *ind, real *z__, integer *ierr,
	 real *rv1, real *rv2, real *rv3, real *rv4, real *rv6)
{
    /* System generated locals */
    integer z_dim1, z_offset, i__1, i__2, i__3;
    real r__1, r__2, r__3, r__4;

    /* Local variables */
    static integer i__, j, p, q, r__, s;
    static real u, v, x0, x1;
    static integer ii, jj, ip;
    static real uk, xu;
    static integer tag, its;
    static real eps2, eps3, eps4, norm, order;
    static integer group;

/* ***BEGIN PROLOGUE  TINVIT */
/* ***PURPOSE  Compute the eigenvectors of symmetric tridiagonal matrix */
/*            corresponding to specified eigenvalues, using inverse */
/*            iteration. */
/* ***LIBRARY   SLATEC (EISPACK) */
/* ***CATEGORY  D4C3 */
/* ***TYPE      SINGLE PRECISION (TINVIT-S) */
/* ***KEYWORDS  EIGENVECTORS, EISPACK */
/* ***AUTHOR  Smith, B. T., et al. */
/* ***DESCRIPTION */

/*     This subroutine is a translation of the inverse iteration tech- */
/*     nique in the ALGOL procedure TRISTURM by Peters and Wilkinson. */
/*     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 418-439(1971). */

/*     This subroutine finds those eigenvectors of a TRIDIAGONAL */
/*     SYMMETRIC matrix corresponding to specified eigenvalues, */
/*     using inverse iteration. */

/*     On Input */

/*        NM must be set to the row dimension of the two-dimensional */
/*          array parameter, Z, as declared in the calling program */
/*          dimension statement.  NM is an INTEGER variable. */

/*        N is the order of the matrix.  N is an INTEGER variable. */
/*          N must be less than or equal to NM. */

/*        D contains the diagonal elements of the symmetric tridiagonal */
/*          matrix.  D is a one-dimensional REAL array, dimensioned D(N). */

/*        E contains the subdiagonal elements of the symmetric */
/*          tridiagonal matrix in its last N-1 positions.  E(1) is */
/*          arbitrary.  E is a one-dimensional REAL array, dimensioned */
/*          E(N). */

/*        E2 contains the squares of the corresponding elements of E, */
/*          with zeros corresponding to negligible elements of E. */
/*          E(I) is considered negligible if it is not larger than */
/*          the product of the relative machine precision and the sum */
/*          of the magnitudes of D(I) and D(I-1).  E2(1) must contain */
/*          0.0e0 if the eigenvalues are in ascending order, or 2.0e0 */
/*          if the eigenvalues are in descending order.  If  BISECT, */
/*          TRIDIB, or  IMTQLV  has been used to find the eigenvalues, */
/*          their output E2 array is exactly what is expected here. */
/*          E2 is a one-dimensional REAL array, dimensioned E2(N). */

/*        M is the number of specified eigenvalues for which eigenvectors */
/*          are to be determined.  M is an INTEGER variable. */

/*        W contains the M eigenvalues in ascending or descending order. */
/*          W is a one-dimensional REAL array, dimensioned W(M). */

/*        IND contains in its first M positions the submatrix indices */
/*          associated with the corresponding eigenvalues in W -- */
/*          1 for eigenvalues belonging to the first submatrix from */
/*          the top, 2 for those belonging to the second submatrix, etc. */
/*          If  BISECT  or  TRIDIB  has been used to determine the */
/*          eigenvalues, their output IND array is suitable for input */
/*          to TINVIT.  IND is a one-dimensional INTEGER array, */
/*          dimensioned IND(M). */

/*     On Output */

/*       ** All input arrays are unaltered.** */

/*        Z contains the associated set of orthonormal eigenvectors. */
/*          Any vector which fails to converge is set to zero. */
/*          Z is a two-dimensional REAL array, dimensioned Z(NM,M). */

/*        IERR is an INTEGER flag set to */
/*          Zero       for normal return, */
/*          -J         if the eigenvector corresponding to the J-th */
/*                     eigenvalue fails to converge in 5 iterations. */

/*        RV1, RV2 and RV3 are one-dimensional REAL arrays used for */
/*          temporary storage.  They are used to store the main diagonal */
/*          and the two adjacent diagonals of the triangular matrix */
/*          produced in the inverse iteration process.  RV1, RV2 and */
/*          RV3 are dimensioned RV1(N), RV2(N) and RV3(N). */

/*        RV4 and RV6 are one-dimensional REAL arrays used for temporary */
/*          storage.  RV4 holds the multipliers of the Gaussian */
/*          elimination process.  RV6 holds the approximate eigenvectors */
/*          in this process.  RV4 and RV6 are dimensioned RV4(N) and */
/*          RV6(N). */

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
/* ***END PROLOGUE  TINVIT */


/* ***FIRST EXECUTABLE STATEMENT  TINVIT */
    /* Parameter adjustments */
    z_dim1 = *nm;
    z_offset = 1 + z_dim1;
    z__ -= z_offset;
    --d__;
    --e;
    --e2;
    --w;
    --ind;
    --rv1;
    --rv2;
    --rv3;
    --rv4;
    --rv6;

    /* Function Body */
    *ierr = 0;
    if (*m == 0) {
	goto L1001;
    }
    tag = 0;
    order = 1.f - e2[1];
    q = 0;
/*     .......... ESTABLISH AND PROCESS NEXT SUBMATRIX .......... */
L100:
    p = q + 1;

    i__1 = *n;
    for (q = p; q <= i__1; ++q) {
	if (q == *n) {
	    goto L140;
	}
	if (e2[q + 1] == 0.f) {
	    goto L140;
	}
/* L120: */
    }
/*     .......... FIND VECTORS BY INVERSE ITERATION .......... */
L140:
    ++tag;
    s = 0;

    i__1 = *m;
    for (r__ = 1; r__ <= i__1; ++r__) {
	if (ind[r__] != tag) {
	    goto L920;
	}
	its = 1;
	x1 = w[r__];
	if (s != 0) {
	    goto L510;
	}
/*     .......... CHECK FOR ISOLATED ROOT .......... */
	xu = 1.f;
	if (p != q) {
	    goto L490;
	}
	rv6[p] = 1.f;
	goto L870;
L490:
	norm = (r__1 = d__[p], dabs(r__1));
	ip = p + 1;

	i__2 = q;
	for (i__ = ip; i__ <= i__2; ++i__) {
/* L500: */
/* Computing MAX */
	    r__3 = norm, r__4 = (r__1 = d__[i__], dabs(r__1)) + (r__2 = e[i__]
		    , dabs(r__2));
	    norm = dmax(r__3,r__4);
	}
/*     .......... EPS2 IS THE CRITERION FOR GROUPING, */
/*                EPS3 REPLACES ZERO PIVOTS AND EQUAL */
/*                ROOTS ARE MODIFIED BY EPS3, */
/*                EPS4 IS TAKEN VERY SMALL TO AVOID OVERFLOW .......... */
	eps2 = norm * .001f;
	eps3 = norm;
L502:
	eps3 *= .5f;
	if (norm + eps3 > norm) {
	    goto L502;
	}
	uk = sqrt((real) (q - p + 5));
	eps3 = uk * eps3;
	eps4 = uk * eps3;
	uk = eps4 / uk;
	s = p;
L505:
	group = 0;
	goto L520;
/*     .......... LOOK FOR CLOSE OR COINCIDENT ROOTS .......... */
L510:
	if ((r__1 = x1 - x0, dabs(r__1)) >= eps2) {
	    goto L505;
	}
	++group;
	if (order * (x1 - x0) <= 0.f) {
	    x1 = x0 + order * eps3;
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
/*     .......... WARNING -- A DIVIDE CHECK MAY OCCUR HERE IF */
/*                E2 ARRAY HAS NOT BEEN SPECIFIED CORRECTLY .......... */
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
	j = r__;

	i__2 = group;
	for (jj = 1; jj <= i__2; ++jj) {
L630:
	    --j;
	    if (ind[j] != tag) {
		goto L630;
	    }
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
	    goto L830;
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
/*     .......... SET ERROR -- NON-CONVERGED EIGENVECTOR .......... */
L830:
	*ierr = -r__;
	xu = 0.f;
	goto L870;
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

L870:
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
L920:
	;
    }

    if (q < *n) {
	goto L100;
    }
L1001:
    return 0;
} /* tinvit_ */

