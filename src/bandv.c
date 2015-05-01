/* bandv.f -- translated by f2c (version 12.02.01).
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

/* DECK BANDV */
/* Subroutine */ int bandv_(integer *nm, integer *n, integer *mbw, real *a, 
	real *e21, integer *m, real *w, real *z__, integer *ierr, integer *nv,
	 real *rv, real *rv6)
{
    /* System generated locals */
    integer a_dim1, a_offset, z_dim1, z_offset, i__1, i__2, i__3, i__4, i__5;
    real r__1;

    /* Local variables */
    static integer i__, j, k, r__;
    static real s, u, v;
    static integer m1;
    static real x0, x1;
    static integer mb, m21, ii, ij, jj, kj;
    static real uk, xu;
    static integer ij1, kj1, its;
    static real eps2, eps3, eps4;
    static integer maxj, maxk;
    static real norm, order;
    static integer group;

/* ***BEGIN PROLOGUE  BANDV */
/* ***PURPOSE  Form the eigenvectors of a real symmetric band matrix */
/*            associated with a set of ordered approximate eigenvalues */
/*            by inverse iteration. */
/* ***LIBRARY   SLATEC (EISPACK) */
/* ***CATEGORY  D4C3 */
/* ***TYPE      SINGLE PRECISION (BANDV-S) */
/* ***KEYWORDS  EIGENVECTORS, EISPACK */
/* ***AUTHOR  Smith, B. T., et al. */
/* ***DESCRIPTION */

/*     This subroutine finds those eigenvectors of a REAL SYMMETRIC */
/*     BAND matrix corresponding to specified eigenvalues, using inverse */
/*     iteration.  The subroutine may also be used to solve systems */
/*     of linear equations with a symmetric or non-symmetric band */
/*     coefficient matrix. */

/*     On INPUT */

/*        NM must be set to the row dimension of the two-dimensional */
/*          array parameters, A and Z, as declared in the calling */
/*          program dimension statement.  NM is an INTEGER variable. */

/*        N is the order of the matrix A.  N is an INTEGER variable. */
/*          N must be less than or equal to NM. */

/*        MBW is the number of columns of the array A used to store the */
/*          band matrix.  If the matrix is symmetric, MBW is its (half) */
/*          band width, denoted MB and defined as the number of adjacent */
/*          diagonals, including the principal diagonal, required to */
/*          specify the non-zero portion of the lower triangle of the */
/*          matrix.  If the subroutine is being used to solve systems */
/*          of linear equations and the coefficient matrix is not */
/*          symmetric, it must however have the same number of adjacent */
/*          diagonals above the main diagonal as below, and in this */
/*          case, MBW=2*MB-1.  MBW is an INTEGER variable.  MB must not */
/*          be greater than N. */

/*        A contains the lower triangle of the symmetric band input */
/*          matrix stored as an N by MB array.  Its lowest subdiagonal */
/*          is stored in the last N+1-MB positions of the first column, */
/*          its next subdiagonal in the last N+2-MB positions of the */
/*          second column, further subdiagonals similarly, and finally */
/*          its principal diagonal in the N positions of column MB. */
/*          If the subroutine is being used to solve systems of linear */
/*          equations and the coefficient matrix is not symmetric, A is */
/*          N by 2*MB-1 instead with lower triangle as above and with */
/*          its first superdiagonal stored in the first N-1 positions of */
/*          column MB+1, its second superdiagonal in the first N-2 */
/*          positions of column MB+2, further superdiagonals similarly, */
/*          and finally its highest superdiagonal in the first N+1-MB */
/*          positions of the last column.  Contents of storage locations */
/*          not part of the matrix are arbitrary.  A is a two-dimensional */
/*          REAL array, dimensioned A(NM,MBW). */

/*        E21 specifies the ordering of the eigenvalues and contains */
/*            0.0E0 if the eigenvalues are in ascending order, or */
/*            2.0E0 if the eigenvalues are in descending order. */
/*          If the subroutine is being used to solve systems of linear */
/*          equations, E21 should be set to 1.0E0 if the coefficient */
/*          matrix is symmetric and to -1.0E0 if not.  E21 is a REAL */
/*          variable. */

/*        M is the number of specified eigenvalues or the number of */
/*          systems of linear equations.  M is an INTEGER variable. */

/*        W contains the M eigenvalues in ascending or descending order. */
/*          If the subroutine is being used to solve systems of linear */
/*          equations (A-W(J)*I)*X(J)=B(J), where I is the identity */
/*          matrix, W(J) should be set accordingly, for J=1,2,...,M. */
/*          W is a one-dimensional REAL array, dimensioned W(M). */

/*        Z contains the constant matrix columns (B(J),J=1,2,...,M), if */
/*          the subroutine is used to solve systems of linear equations. */
/*          Z is a two-dimensional REAL array, dimensioned Z(NM,M). */

/*        NV must be set to the dimension of the array parameter RV */
/*          as declared in the calling program dimension statement. */
/*          NV is an INTEGER variable. */

/*     On OUTPUT */

/*        A and W are unaltered. */

/*        Z contains the associated set of orthogonal eigenvectors. */
/*          Any vector which fails to converge is set to zero.  If the */
/*          subroutine is used to solve systems of linear equations, */
/*          Z contains the solution matrix columns (X(J),J=1,2,...,M). */

/*        IERR is an INTEGER flag set to */
/*          Zero       for normal return, */
/*          -J         if the eigenvector corresponding to the J-th */
/*                     eigenvalue fails to converge, or if the J-th */
/*                     system of linear equations is nearly singular. */

/*        RV and RV6 are temporary storage arrays.  If the subroutine */
/*          is being used to solve systems of linear equations, the */
/*          determinant (up to sign) of A-W(M)*I is available, upon */
/*          return, as the product of the first N elements of RV. */
/*          RV and RV6 are one-dimensional REAL arrays.  Note that RV */
/*          is dimensioned RV(NV), where NV must be at least N*(2*MB-1). */
/*          RV6 is dimensioned RV6(N). */

/*     Questions and comments should be directed to B. S. Garbow, */
/*     Applied Mathematics Division, ARGONNE NATIONAL LABORATORY */
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
/* ***END PROLOGUE  BANDV */


/* ***FIRST EXECUTABLE STATEMENT  BANDV */
    /* Parameter adjustments */
    z_dim1 = *nm;
    z_offset = 1 + z_dim1;
    z__ -= z_offset;
    a_dim1 = *nm;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --w;
    --rv;
    --rv6;

    /* Function Body */
    *ierr = 0;
    if (*m == 0) {
	goto L1001;
    }
    mb = *mbw;
    if (*e21 < 0.f) {
	mb = (*mbw + 1) / 2;
    }
    m1 = mb - 1;
    m21 = m1 + mb;
    order = 1.f - dabs(*e21);
/*     .......... FIND VECTORS BY INVERSE ITERATION .......... */
    i__1 = *m;
    for (r__ = 1; r__ <= i__1; ++r__) {
	its = 1;
	x1 = w[r__];
	if (r__ != 1) {
	    goto L100;
	}
/*     .......... COMPUTE NORM OF MATRIX .......... */
	norm = 0.f;

	i__2 = mb;
	for (j = 1; j <= i__2; ++j) {
	    jj = mb + 1 - j;
	    kj = jj + m1;
	    ij = 1;
	    s = 0.f;

	    i__3 = *n;
	    for (i__ = jj; i__ <= i__3; ++i__) {
		s += (r__1 = a[i__ + j * a_dim1], dabs(r__1));
		if (*e21 >= 0.f) {
		    goto L40;
		}
		s += (r__1 = a[ij + kj * a_dim1], dabs(r__1));
		++ij;
L40:
		;
	    }

	    norm = dmax(norm,s);
/* L60: */
	}

	if (*e21 < 0.f) {
	    norm *= .5f;
	}
/*     .......... EPS2 IS THE CRITERION FOR GROUPING, */
/*                EPS3 REPLACES ZERO PIVOTS AND EQUAL */
/*                ROOTS ARE MODIFIED BY EPS3, */
/*                EPS4 IS TAKEN VERY SMALL TO AVOID OVERFLOW .......... */
	if (norm == 0.f) {
	    norm = 1.f;
	}
	eps2 = norm * .001f * dabs(order);
	eps3 = norm;
L70:
	eps3 *= .5f;
	if (norm + eps3 > norm) {
	    goto L70;
	}
	uk = sqrt((real) (*n));
	eps3 = uk * eps3;
	eps4 = uk * eps3;
L80:
	group = 0;
	goto L120;
/*     .......... LOOK FOR CLOSE OR COINCIDENT ROOTS .......... */
L100:
	if ((r__1 = x1 - x0, dabs(r__1)) >= eps2) {
	    goto L80;
	}
	++group;
	if (order * (x1 - x0) <= 0.f) {
	    x1 = x0 + order * eps3;
	}
/*     .......... EXPAND MATRIX, SUBTRACT EIGENVALUE, */
/*                AND INITIALIZE VECTOR .......... */
L120:
	i__2 = *n;
	for (i__ = 1; i__ <= i__2; ++i__) {
/* Computing MIN */
	    i__3 = 0, i__4 = i__ - m1;
	    ij = i__ + min(i__3,i__4) * *n;
	    kj = ij + mb * *n;
	    ij1 = kj + m1 * *n;
	    if (m1 == 0) {
		goto L180;
	    }

	    i__3 = m1;
	    for (j = 1; j <= i__3; ++j) {
		if (ij > m1) {
		    goto L125;
		}
		if (ij > 0) {
		    goto L130;
		}
		rv[ij1] = 0.f;
		ij1 += *n;
		goto L130;
L125:
		rv[ij] = a[i__ + j * a_dim1];
L130:
		ij += *n;
		ii = i__ + j;
		if (ii > *n) {
		    goto L150;
		}
		jj = mb - j;
		if (*e21 >= 0.f) {
		    goto L140;
		}
		ii = i__;
		jj = mb + j;
L140:
		rv[kj] = a[ii + jj * a_dim1];
		kj += *n;
L150:
		;
	    }

L180:
	    rv[ij] = a[i__ + mb * a_dim1] - x1;
	    rv6[i__] = eps4;
	    if (order == 0.f) {
		rv6[i__] = z__[i__ + r__ * z_dim1];
	    }
/* L200: */
	}

	if (m1 == 0) {
	    goto L600;
	}
/*     .......... ELIMINATION WITH INTERCHANGES .......... */
	i__2 = *n;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    ii = i__ + 1;
/* Computing MIN */
	    i__3 = i__ + m1 - 1;
	    maxk = min(i__3,*n);
/* Computing MIN */
	    i__3 = *n - i__, i__4 = m21 - 2;
	    maxj = min(i__3,i__4) * *n;

	    i__3 = maxk;
	    for (k = i__; k <= i__3; ++k) {
		kj1 = k;
		j = kj1 + *n;
		jj = j + maxj;

		i__4 = jj;
		i__5 = *n;
		for (kj = j; i__5 < 0 ? kj >= i__4 : kj <= i__4; kj += i__5) {
		    rv[kj1] = rv[kj];
		    kj1 = kj;
/* L340: */
		}

		rv[kj1] = 0.f;
/* L360: */
	    }

	    if (i__ == *n) {
		goto L580;
	    }
	    u = 0.f;
/* Computing MIN */
	    i__3 = i__ + m1;
	    maxk = min(i__3,*n);
/* Computing MIN */
	    i__3 = *n - ii, i__5 = m21 - 2;
	    maxj = min(i__3,i__5) * *n;

	    i__3 = maxk;
	    for (j = i__; j <= i__3; ++j) {
		if ((r__1 = rv[j], dabs(r__1)) < dabs(u)) {
		    goto L450;
		}
		u = rv[j];
		k = j;
L450:
		;
	    }

	    j = i__ + *n;
	    jj = j + maxj;
	    if (k == i__) {
		goto L520;
	    }
	    kj = k;

	    i__3 = jj;
	    i__5 = *n;
	    for (ij = i__; i__5 < 0 ? ij >= i__3 : ij <= i__3; ij += i__5) {
		v = rv[ij];
		rv[ij] = rv[kj];
		rv[kj] = v;
		kj += *n;
/* L500: */
	    }

	    if (order != 0.f) {
		goto L520;
	    }
	    v = rv6[i__];
	    rv6[i__] = rv6[k];
	    rv6[k] = v;
L520:
	    if (u == 0.f) {
		goto L580;
	    }

	    i__5 = maxk;
	    for (k = ii; k <= i__5; ++k) {
		v = rv[k] / u;
		kj = k;

		i__3 = jj;
		i__4 = *n;
		for (ij = j; i__4 < 0 ? ij >= i__3 : ij <= i__3; ij += i__4) {
		    kj += *n;
		    rv[kj] -= v * rv[ij];
/* L540: */
		}

		if (order == 0.f) {
		    rv6[k] -= v * rv6[i__];
		}
/* L560: */
	    }

L580:
	    ;
	}
/*     .......... BACK SUBSTITUTION */
/*                FOR I=N STEP -1 UNTIL 1 DO -- .......... */
L600:
	i__2 = *n;
	for (ii = 1; ii <= i__2; ++ii) {
	    i__ = *n + 1 - ii;
	    maxj = min(ii,m21);
	    if (maxj == 1) {
		goto L620;
	    }
	    ij1 = i__;
	    j = ij1 + *n;
	    jj = j + (maxj - 2) * *n;

	    i__5 = jj;
	    i__4 = *n;
	    for (ij = j; i__4 < 0 ? ij >= i__5 : ij <= i__5; ij += i__4) {
		++ij1;
		rv6[i__] -= rv[ij] * rv6[ij1];
/* L610: */
	    }

L620:
	    v = rv[i__];
	    if (dabs(v) >= eps3) {
		goto L625;
	    }
/*     .......... SET ERROR -- NEARLY SINGULAR LINEAR SYSTEM .......... */
	    if (order == 0.f) {
		*ierr = -r__;
	    }
	    v = r_sign(&eps3, &v);
L625:
	    rv6[i__] /= v;
/* L630: */
	}

	xu = 1.f;
	if (order == 0.f) {
	    goto L870;
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

	    i__4 = *n;
	    for (i__ = 1; i__ <= i__4; ++i__) {
/* L640: */
		xu += rv6[i__] * z__[i__ + j * z_dim1];
	    }

	    i__4 = *n;
	    for (i__ = 1; i__ <= i__4; ++i__) {
/* L660: */
		rv6[i__] -= xu * z__[i__ + j * z_dim1];
	    }

/* L680: */
	}

L700:
	norm = 0.f;

	i__2 = *n;
	for (i__ = 1; i__ <= i__2; ++i__) {
/* L720: */
	    norm += (r__1 = rv6[i__], dabs(r__1));
	}

	if (norm >= .1f) {
	    goto L840;
	}
/*     .......... IN-LINE PROCEDURE FOR CHOOSING */
/*                A NEW STARTING VECTOR .......... */
	if (its >= *n) {
	    goto L830;
	}
	++its;
	xu = eps4 / (uk + 1.f);
	rv6[1] = eps4;

	i__2 = *n;
	for (i__ = 2; i__ <= i__2; ++i__) {
/* L760: */
	    rv6[i__] = xu;
	}

	rv6[its] -= eps4 * uk;
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

	i__2 = *n;
	for (i__ = 1; i__ <= i__2; ++i__) {
/* L860: */
/* Computing 2nd power */
	    r__1 = rv6[i__];
	    u += r__1 * r__1;
	}

	xu = 1.f / sqrt(u);

L870:
	i__2 = *n;
	for (i__ = 1; i__ <= i__2; ++i__) {
/* L900: */
	    z__[i__ + r__ * z_dim1] = rv6[i__] * xu;
	}

	x0 = x1;
/* L920: */
    }

L1001:
    return 0;
} /* bandv_ */

