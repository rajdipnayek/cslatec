/* dhfti.f -- translated by f2c (version 12.02.01).
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
static integer c__1 = 1;
static integer c__2 = 2;

/* DECK DHFTI */
/* Subroutine */ int dhfti_(doublereal *a, integer *mda, integer *m, integer *
	n, doublereal *b, integer *mdb, integer *nb, doublereal *tau, integer 
	*krank, doublereal *rnorm, doublereal *h__, doublereal *g, integer *
	ip)
{
    /* Initialized data */

    static doublereal releps = 0.;

    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, i__1, i__2, i__3;
    doublereal d__1;

    /* Local variables */
    static integer i__, j, k, l, jb, ii, jj;
    static doublereal sm;
    static integer ip1, kp1;
    static doublereal sm1;
    extern /* Subroutine */ int dh12_(integer *, integer *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    integer *, integer *, integer *);
    static doublereal tmp, hmax;
    static integer lmax, nerr, iopt, ldiag;
    static doublereal dzero;
    extern doublereal d1mach_(integer *);
    static doublereal szero, factor;
    extern /* Subroutine */ int xermsg_(char *, char *, char *, integer *, 
	    integer *, ftnlen, ftnlen, ftnlen);

/* ***BEGIN PROLOGUE  DHFTI */
/* ***PURPOSE  Solve a least squares problem for banded matrices using */
/*            sequential accumulation of rows of the data matrix. */
/*            Exactly one right-hand side vector is permitted. */
/* ***LIBRARY   SLATEC */
/* ***CATEGORY  D9 */
/* ***TYPE      DOUBLE PRECISION (HFTI-S, DHFTI-D) */
/* ***KEYWORDS  CURVE FITTING, LEAST SQUARES */
/* ***AUTHOR  Lawson, C. L., (JPL) */
/*           Hanson, R. J., (SNLA) */
/* ***DESCRIPTION */

/*     DIMENSION A(MDA,N),(B(MDB,NB) or B(M)),RNORM(NB),H(N),G(N),IP(N) */

/*     This subroutine solves a linear least squares problem or a set of */
/*     linear least squares problems having the same matrix but different */
/*     right-side vectors.  The problem data consists of an M by N matrix */
/*     A, an M by NB matrix B, and an absolute tolerance parameter TAU */
/*     whose usage is described below.  The NB column vectors of B */
/*     represent right-side vectors for NB distinct linear least squares */
/*     problems. */

/*     This set of problems can also be written as the matrix least */
/*     squares problem */

/*                       AX = B, */

/*     where X is the N by NB solution matrix. */

/*     Note that if B is the M by M identity matrix, then X will be the */
/*     pseudo-inverse of A. */

/*     This subroutine first transforms the augmented matrix (A B) to a */
/*     matrix (R C) using premultiplying Householder transformations with */
/*     column interchanges.  All subdiagonal elements in the matrix R are */
/*     zero and its diagonal elements satisfy */

/*                       ABS(R(I,I)).GE.ABS(R(I+1,I+1)), */

/*                       I = 1,...,L-1, where */

/*                       L = MIN(M,N). */

/*     The subroutine will compute an integer, KRANK, equal to the number */
/*     of diagonal terms of R that exceed TAU in magnitude. Then a */
/*     solution of minimum Euclidean length is computed using the first */
/*     KRANK rows of (R C). */

/*     To be specific we suggest that the user consider an easily */
/*     computable matrix norm, such as, the maximum of all column sums of */
/*     magnitudes. */

/*     Now if the relative uncertainty of B is EPS, (norm of uncertainty/ */
/*     norm of B), it is suggested that TAU be set approximately equal to */
/*     EPS*(norm of A). */

/*     The user must dimension all arrays appearing in the call list.. */
/*     A(MDA,N),(B(MDB,NB) or B(M)),RNORM(NB),H(N),G(N),IP(N).  This */
/*     permits the solution of a range of problems in the same array */
/*     space. */

/*     The entire set of parameters for DHFTI are */

/*     INPUT.. All TYPE REAL variables are DOUBLE PRECISION */

/*     A(*,*),MDA,M,N    The array A(*,*) initially contains the M by N */
/*                       matrix A of the least squares problem AX = B. */
/*                       The first dimensioning parameter of the array */
/*                       A(*,*) is MDA, which must satisfy MDA.GE.M */
/*                       Either M.GE.N or M.LT.N is permitted.  There */
/*                       is no restriction on the rank of A.  The */
/*                       condition MDA.LT.M is considered an error. */

/*     B(*),MDB,NB       If NB = 0 the subroutine will perform the */
/*                       orthogonal decomposition but will make no */
/*                       references to the array B(*).  If NB.GT.0 */
/*                       the array B(*) must initially contain the M by */
/*                       NB matrix B of the least squares problem AX = */
/*                       B.  If NB.GE.2 the array B(*) must be doubly */
/*                       subscripted with first dimensioning parameter */
/*                       MDB.GE.MAX(M,N).  If NB = 1 the array B(*) may */
/*                       be either doubly or singly subscripted.  In */
/*                       the latter case the value of MDB is arbitrary */
/*                       but it should be set to some valid integer */
/*                       value such as MDB = M. */

/*                       The condition of NB.GT.1.AND.MDB.LT. MAX(M,N) */
/*                       is considered an error. */

/*     TAU               Absolute tolerance parameter provided by user */
/*                       for pseudorank determination. */

/*     H(*),G(*),IP(*)   Arrays of working space used by DHFTI. */

/*     OUTPUT.. All TYPE REAL variables are DOUBLE PRECISION */

/*     A(*,*)            The contents of the array A(*,*) will be */
/*                       modified by the subroutine. These contents */
/*                       are not generally required by the user. */

/*     B(*)              On return the array B(*) will contain the N by */
/*                       NB solution matrix X. */

/*     KRANK             Set by the subroutine to indicate the */
/*                       pseudorank of A. */

/*     RNORM(*)          On return, RNORM(J) will contain the Euclidean */
/*                       norm of the residual vector for the problem */
/*                       defined by the J-th column vector of the array */
/*                       B(*,*) for J = 1,...,NB. */

/*     H(*),G(*)         On return these arrays respectively contain */
/*                       elements of the pre- and post-multiplying */
/*                       Householder transformations used to compute */
/*                       the minimum Euclidean length solution. */

/*     IP(*)             Array in which the subroutine records indices */
/*                       describing the permutation of column vectors. */
/*                       The contents of arrays H(*),G(*) and IP(*) */
/*                       are not generally required by the user. */

/* ***REFERENCES  C. L. Lawson and R. J. Hanson, Solving Least Squares */
/*                 Problems, Prentice-Hall, Inc., 1974, Chapter 14. */
/* ***ROUTINES CALLED  D1MACH, DH12, XERMSG */
/* ***REVISION HISTORY  (YYMMDD) */
/*   790101  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   891006  Cosmetic changes to prologue.  (WRB) */
/*   891006  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ) */
/*   901005  Replace usage of DDIFF with usage of D1MACH.  (RWC) */
/*   920501  Reformatted the REFERENCES section.  (WRB) */
/* ***END PROLOGUE  DHFTI */
    /* Parameter adjustments */
    a_dim1 = *mda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    b_dim1 = *mdb;
    b_offset = 1 + b_dim1;
    b -= b_offset;
    --rnorm;
    --h__;
    --g;
    --ip;

    /* Function Body */
/*     BEGIN BLOCK PERMITTING ...EXITS TO 360 */
/* ***FIRST EXECUTABLE STATEMENT  DHFTI */
    if (releps == 0.) {
	releps = d1mach_(&c__4);
    }
    szero = 0.;
    dzero = 0.;
    factor = .001;

    k = 0;
    ldiag = min(*m,*n);
    if (ldiag <= 0) {
	goto L350;
    }
/*           BEGIN BLOCK PERMITTING ...EXITS TO 130 */
/*              BEGIN BLOCK PERMITTING ...EXITS TO 120 */
    if (*mda >= *m) {
	goto L10;
    }
    nerr = 1;
    iopt = 2;
    xermsg_("SLATEC", "DHFTI", "MDA.LT.M, PROBABLE ERROR.", &nerr, &iopt, (
	    ftnlen)6, (ftnlen)5, (ftnlen)25);
/*     ...............EXIT */
    goto L360;
L10:

    if (*nb <= 1 || max(*m,*n) <= *mdb) {
	goto L20;
    }
    nerr = 2;
    iopt = 2;
    xermsg_("SLATEC", "DHFTI", "MDB.LT.MAX(M,N).AND.NB.GT.1. PROBABLE ERROR.",
	     &nerr, &iopt, (ftnlen)6, (ftnlen)5, (ftnlen)44);
/*     ...............EXIT */
    goto L360;
L20:

    i__1 = ldiag;
    for (j = 1; j <= i__1; ++j) {
/*                    BEGIN BLOCK PERMITTING ...EXITS TO 70 */
	if (j == 1) {
	    goto L40;
	}

/*                           UPDATE SQUARED COLUMN LENGTHS AND FIND LMAX */
/*                          .. */
	lmax = j;
	i__2 = *n;
	for (l = j; l <= i__2; ++l) {
/* Computing 2nd power */
	    d__1 = a[j - 1 + l * a_dim1];
	    h__[l] -= d__1 * d__1;
	    if (h__[l] > h__[lmax]) {
		lmax = l;
	    }
/* L30: */
	}
/*                    ......EXIT */
	if (factor * h__[lmax] > hmax * releps) {
	    goto L70;
	}
L40:

/*                        COMPUTE SQUARED COLUMN LENGTHS AND FIND LMAX */
/*                       .. */
	lmax = j;
	i__2 = *n;
	for (l = j; l <= i__2; ++l) {
	    h__[l] = 0.;
	    i__3 = *m;
	    for (i__ = j; i__ <= i__3; ++i__) {
/* Computing 2nd power */
		d__1 = a[i__ + l * a_dim1];
		h__[l] += d__1 * d__1;
/* L50: */
	    }
	    if (h__[l] > h__[lmax]) {
		lmax = l;
	    }
/* L60: */
	}
	hmax = h__[lmax];
L70:
/*                    .. */
/*                     LMAX HAS BEEN DETERMINED */

/*                     DO COLUMN INTERCHANGES IF NEEDED. */
/*                    .. */
	ip[j] = lmax;
	if (ip[j] == j) {
	    goto L90;
	}
	i__2 = *m;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    tmp = a[i__ + j * a_dim1];
	    a[i__ + j * a_dim1] = a[i__ + lmax * a_dim1];
	    a[i__ + lmax * a_dim1] = tmp;
/* L80: */
	}
	h__[lmax] = h__[j];
L90:

/*                     COMPUTE THE J-TH TRANSFORMATION AND APPLY IT TO A */
/*                     AND B. */
/*                    .. */
	i__2 = j + 1;
	i__3 = *n - j;
	dh12_(&c__1, &j, &i__2, m, &a[j * a_dim1 + 1], &c__1, &h__[j], &a[(j 
		+ 1) * a_dim1 + 1], &c__1, mda, &i__3);
	i__2 = j + 1;
	dh12_(&c__2, &j, &i__2, m, &a[j * a_dim1 + 1], &c__1, &h__[j], &b[
		b_offset], &c__1, mdb, nb);
/* L100: */
    }

/*                  DETERMINE THE PSEUDORANK, K, USING THE TOLERANCE, */
/*                  TAU. */
/*                 .. */
    i__1 = ldiag;
    for (j = 1; j <= i__1; ++j) {
/*              ......EXIT */
	if ((d__1 = a[j + j * a_dim1], abs(d__1)) <= *tau) {
	    goto L120;
	}
/* L110: */
    }
    k = ldiag;
/*           ......EXIT */
    goto L130;
L120:
    k = j - 1;
L130:
    kp1 = k + 1;

/*           COMPUTE THE NORMS OF THE RESIDUAL VECTORS. */

    if (*nb < 1) {
	goto L170;
    }
    i__1 = *nb;
    for (jb = 1; jb <= i__1; ++jb) {
	tmp = szero;
	if (*m < kp1) {
	    goto L150;
	}
	i__2 = *m;
	for (i__ = kp1; i__ <= i__2; ++i__) {
/* Computing 2nd power */
	    d__1 = b[i__ + jb * b_dim1];
	    tmp += d__1 * d__1;
/* L140: */
	}
L150:
	rnorm[jb] = sqrt(tmp);
/* L160: */
    }
L170:
/*           SPECIAL FOR PSEUDORANK = 0 */
    if (k > 0) {
	goto L210;
    }
    if (*nb < 1) {
	goto L200;
    }
    i__1 = *nb;
    for (jb = 1; jb <= i__1; ++jb) {
	i__2 = *n;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    b[i__ + jb * b_dim1] = szero;
/* L180: */
	}
/* L190: */
    }
L200:
    goto L340;
L210:

/*               IF THE PSEUDORANK IS LESS THAN N COMPUTE HOUSEHOLDER */
/*               DECOMPOSITION OF FIRST K ROWS. */
/*              .. */
    if (k == *n) {
	goto L230;
    }
    i__1 = k;
    for (ii = 1; ii <= i__1; ++ii) {
	i__ = kp1 - ii;
	i__2 = i__ - 1;
	dh12_(&c__1, &i__, &kp1, n, &a[i__ + a_dim1], mda, &g[i__], &a[
		a_offset], mda, &c__1, &i__2);
/* L220: */
    }
L230:


    if (*nb < 1) {
	goto L330;
    }
    i__1 = *nb;
    for (jb = 1; jb <= i__1; ++jb) {

/*                  SOLVE THE K BY K TRIANGULAR SYSTEM. */
/*                 .. */
	i__2 = k;
	for (l = 1; l <= i__2; ++l) {
	    sm = dzero;
	    i__ = kp1 - l;
	    ip1 = i__ + 1;
	    if (k < ip1) {
		goto L250;
	    }
	    i__3 = k;
	    for (j = ip1; j <= i__3; ++j) {
		sm += a[i__ + j * a_dim1] * b[j + jb * b_dim1];
/* L240: */
	    }
L250:
	    sm1 = sm;
	    b[i__ + jb * b_dim1] = (b[i__ + jb * b_dim1] - sm1) / a[i__ + i__ 
		    * a_dim1];
/* L260: */
	}

/*                  COMPLETE COMPUTATION OF SOLUTION VECTOR. */
/*                 .. */
	if (k == *n) {
	    goto L290;
	}
	i__2 = *n;
	for (j = kp1; j <= i__2; ++j) {
	    b[j + jb * b_dim1] = szero;
/* L270: */
	}
	i__2 = k;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    dh12_(&c__2, &i__, &kp1, n, &a[i__ + a_dim1], mda, &g[i__], &b[jb 
		    * b_dim1 + 1], &c__1, mdb, &c__1);
/* L280: */
	}
L290:

/*                   RE-ORDER THE SOLUTION VECTOR TO COMPENSATE FOR THE */
/*                   COLUMN INTERCHANGES. */
/*                 .. */
	i__2 = ldiag;
	for (jj = 1; jj <= i__2; ++jj) {
	    j = ldiag + 1 - jj;
	    if (ip[j] == j) {
		goto L300;
	    }
	    l = ip[j];
	    tmp = b[l + jb * b_dim1];
	    b[l + jb * b_dim1] = b[j + jb * b_dim1];
	    b[j + jb * b_dim1] = tmp;
L300:
/* L310: */
	    ;
	}
/* L320: */
    }
L330:
L340:
L350:
/*        .. */
/*         THE SOLUTION VECTORS, X, ARE NOW */
/*         IN THE FIRST  N  ROWS OF THE ARRAY B(,). */

    *krank = k;
L360:
    return 0;
} /* dhfti_ */

