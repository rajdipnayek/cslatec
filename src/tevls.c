/* tevls.f -- translated by f2c (version 12.02.01).
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

/* Common Block Declarations */

struct {
    integer npp, k;
    real machep, cnv;
    integer nm, ncmplx, ik;
} cblkt_;

#define cblkt_1 cblkt_

/* DECK TEVLS */
/* Subroutine */ int tevls_(integer *n, real *d__, real *e2, integer *ierr)
{
    /* System generated locals */
    integer i__1, i__2;
    real r__1, r__2;

    /* Local variables */
    static real b, c__, f, g, h__;
    static integer i__, j, l, m;
    static real p, r__, s;
    static integer l1, ii, mml, ntop, nhalf;
    static real dhold;

/* ***BEGIN PROLOGUE  TEVLS */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to BLKTRI */
/* ***LIBRARY   SLATEC */
/* ***TYPE      SINGLE PRECISION (TEVLS-S) */
/* ***AUTHOR  (UNKNOWN) */
/* ***DESCRIPTION */

/*     This subroutine finds the eigenvalues of a symmetric tridiagonal */
/*     matrix by the rational QL method. */

/*     On Input- */

/*        N is the order of the matrix, */

/*        D contains the diagonal elements of the input matrix, */

/*        E2 contains the subdiagonal elements of the input matrix */
/*           in its last N-1 positions.  E2(1) is arbitrary. */

/*      On Output- */

/*        D contains the eigenvalues in ascending order.  If an */
/*          error exit is made, the eigenvalues are correct and */
/*          ordered for indices 1,2,...IERR-1, but may not be */
/*          the smallest eigenvalues, */

/*        E2 has been destroyed, */

/*        IERR is set to */
/*          ZERO       for normal return, */
/*          J          if the J-th eigenvalue has not been */
/*                     determined after 30 iterations. */

/* ***SEE ALSO  BLKTRI */
/* ***REFERENCES  C. H. Reinsch, Eigenvalues of a real, symmetric, tri- */
/*                 diagonal matrix, Algorithm 464, Communications of the */
/*                 ACM 16, 11 (November 1973), pp. 689. */
/* ***ROUTINES CALLED  (NONE) */
/* ***COMMON BLOCKS    CBLKT */
/* ***REVISION HISTORY  (YYMMDD) */
/*   801001  DATE WRITTEN */
/*   890831  Modified array declarations.  (WRB) */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900402  Added TYPE section.  (WRB) */
/*   920528  DESCRIPTION revised and REFERENCES section added.  (WRB) */
/* ***END PROLOGUE  TEVLS */


/* ***FIRST EXECUTABLE STATEMENT  TEVLS */
    /* Parameter adjustments */
    --e2;
    --d__;

    /* Function Body */
    *ierr = 0;
    if (*n == 1) {
	goto L115;
    }

    i__1 = *n;
    for (i__ = 2; i__ <= i__1; ++i__) {
	e2[i__ - 1] = e2[i__] * e2[i__];
/* L101: */
    }

    f = 0.f;
    b = 0.f;
    e2[*n] = 0.f;

    i__1 = *n;
    for (l = 1; l <= i__1; ++l) {
	j = 0;
	h__ = cblkt_1.machep * ((r__1 = d__[l], dabs(r__1)) + sqrt(e2[l]));
	if (b > h__) {
	    goto L102;
	}
	b = h__;
	c__ = b * b;

/*     ********** LOOK FOR SMALL SQUARED SUB-DIAGONAL ELEMENT ********** */

L102:
	i__2 = *n;
	for (m = l; m <= i__2; ++m) {
	    if (e2[m] <= c__) {
		goto L104;
	    }

/*     ********** E2(N) IS ALWAYS ZERO, SO THERE IS NO EXIT */
/*                THROUGH THE BOTTOM OF THE LOOP ********** */

/* L103: */
	}

L104:
	if (m == l) {
	    goto L108;
	}
L105:
	if (j == 30) {
	    goto L114;
	}
	++j;

/*     ********** FORM SHIFT ********** */

	l1 = l + 1;
	s = sqrt(e2[l]);
	g = d__[l];
	p = (d__[l1] - g) / (s * 2.f);
	r__ = sqrt(p * p + 1.f);
	d__[l] = s / (p + r_sign(&r__, &p));
	h__ = g - d__[l];

	i__2 = *n;
	for (i__ = l1; i__ <= i__2; ++i__) {
	    d__[i__] -= h__;
/* L106: */
	}

	f += h__;

/*     ********** RATIONAL QL TRANSFORMATION ********** */

	g = d__[m];
	if (g == 0.f) {
	    g = b;
	}
	h__ = g;
	s = 0.f;
	mml = m - l;

/*     ********** FOR I=M-1 STEP -1 UNTIL L DO -- ********** */

	i__2 = mml;
	for (ii = 1; ii <= i__2; ++ii) {
	    i__ = m - ii;
	    p = g * h__;
	    r__ = p + e2[i__];
	    e2[i__ + 1] = s * r__;
	    s = e2[i__] / r__;
	    d__[i__ + 1] = h__ + s * (h__ + d__[i__]);
	    g = d__[i__] - e2[i__] / g;
	    if (g == 0.f) {
		g = b;
	    }
	    h__ = g * p / r__;
/* L107: */
	}

	e2[l] = s * g;
	d__[l] = h__;

/*     ********** GUARD AGAINST UNDERFLOWED H ********** */

	if (h__ == 0.f) {
	    goto L108;
	}
	if ((r__1 = e2[l], dabs(r__1)) <= (r__2 = c__ / h__, dabs(r__2))) {
	    goto L108;
	}
	e2[l] = h__ * e2[l];
	if (e2[l] != 0.f) {
	    goto L105;
	}
L108:
	p = d__[l] + f;

/*     ********** ORDER EIGENVALUES ********** */

	if (l == 1) {
	    goto L110;
	}

/*     ********** FOR I=L STEP -1 UNTIL 2 DO -- ********** */

	i__2 = l;
	for (ii = 2; ii <= i__2; ++ii) {
	    i__ = l + 2 - ii;
	    if (p >= d__[i__ - 1]) {
		goto L111;
	    }
	    d__[i__] = d__[i__ - 1];
/* L109: */
	}

L110:
	i__ = 1;
L111:
	d__[i__] = p;
/* L112: */
    }

    if ((r__1 = d__[*n], dabs(r__1)) >= dabs(d__[1])) {
	goto L115;
    }
    nhalf = *n / 2;
    i__1 = nhalf;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ntop = *n - i__;
	dhold = d__[i__];
	d__[i__] = d__[ntop + 1];
	d__[ntop + 1] = dhold;
/* L113: */
    }
    goto L115;

/*     ********** SET ERROR -- NO CONVERGENCE TO AN */
/*                EIGENVALUE AFTER 30 ITERATIONS ********** */

L114:
    *ierr = l;
L115:
    return 0;

/*     ********** LAST CARD OF TQLRAT ********** */

} /* tevls_ */

