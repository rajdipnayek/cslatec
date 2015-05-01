/* dppfa.f -- translated by f2c (version 12.02.01).
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

static integer c__1 = 1;

/* DECK DPPFA */
/* Subroutine */ int dppfa_(doublereal *ap, integer *n, integer *info)
{
    /* System generated locals */
    integer i__1, i__2, i__3;

    /* Local variables */
    static integer j, k;
    static doublereal s, t;
    static integer jj, kj, kk, jm1;
    extern doublereal ddot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *);

/* ***BEGIN PROLOGUE  DPPFA */
/* ***PURPOSE  Factor a real symmetric positive definite matrix stored in */
/*            packed form. */
/* ***LIBRARY   SLATEC (LINPACK) */
/* ***CATEGORY  D2B1B */
/* ***TYPE      DOUBLE PRECISION (SPPFA-S, DPPFA-D, CPPFA-C) */
/* ***KEYWORDS  LINEAR ALGEBRA, LINPACK, MATRIX FACTORIZATION, PACKED, */
/*             POSITIVE DEFINITE */
/* ***AUTHOR  Moler, C. B., (U. of New Mexico) */
/* ***DESCRIPTION */

/*     DPPFA factors a double precision symmetric positive definite */
/*     matrix stored in packed form. */

/*     DPPFA is usually called by DPPCO, but it can be called */
/*     directly with a saving in time if  RCOND  is not needed. */
/*     (time for DPPCO) = (1 + 18/N)*(time for DPPFA) . */

/*     On Entry */

/*        AP      DOUBLE PRECISION (N*(N+1)/2) */
/*                the packed form of a symmetric matrix  A .  The */
/*                columns of the upper triangle are stored sequentially */
/*                in a one-dimensional array of length  N*(N+1)/2 . */
/*                See comments below for details. */

/*        N       INTEGER */
/*                the order of the matrix  A . */

/*     On Return */

/*        AP      an upper triangular matrix  R , stored in packed */
/*                form, so that  A = TRANS(R)*R . */

/*        INFO    INTEGER */
/*                = 0  for normal return. */
/*                = K  if the leading minor of order  K  is not */
/*                     positive definite. */


/*     Packed Storage */

/*          The following program segment will pack the upper */
/*          triangle of a symmetric matrix. */

/*                K = 0 */
/*                DO 20 J = 1, N */
/*                   DO 10 I = 1, J */
/*                      K = K + 1 */
/*                      AP(K) = A(I,J) */
/*             10    CONTINUE */
/*             20 CONTINUE */

/* ***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W. */
/*                 Stewart, LINPACK Users' Guide, SIAM, 1979. */
/* ***ROUTINES CALLED  DDOT */
/* ***REVISION HISTORY  (YYMMDD) */
/*   780814  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890831  Modified array declarations.  (WRB) */
/*   890831  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900326  Removed duplicate information from DESCRIPTION section. */
/*           (WRB) */
/*   920501  Reformatted the REFERENCES section.  (WRB) */
/* ***END PROLOGUE  DPPFA */

/* ***FIRST EXECUTABLE STATEMENT  DPPFA */
    /* Parameter adjustments */
    --ap;

    /* Function Body */
    jj = 0;
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	*info = j;
	s = 0.;
	jm1 = j - 1;
	kj = jj;
	kk = 0;
	if (jm1 < 1) {
	    goto L20;
	}
	i__2 = jm1;
	for (k = 1; k <= i__2; ++k) {
	    ++kj;
	    i__3 = k - 1;
	    t = ap[kj] - ddot_(&i__3, &ap[kk + 1], &c__1, &ap[jj + 1], &c__1);
	    kk += k;
	    t /= ap[kk];
	    ap[kj] = t;
	    s += t * t;
/* L10: */
	}
L20:
	jj += j;
	s = ap[jj] - s;
	if (s <= 0.) {
	    goto L40;
	}
	ap[jj] = sqrt(s);
/* L30: */
    }
    *info = 0;
L40:
    return 0;
} /* dppfa_ */

