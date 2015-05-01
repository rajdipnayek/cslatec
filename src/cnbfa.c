/* cnbfa.f -- translated by f2c (version 12.02.01).
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

static complex c_b9 = {-1.f,-0.f};

/* DECK CNBFA */
/* Subroutine */ int cnbfa_(complex *abe, integer *lda, integer *n, integer *
	ml, integer *mu, integer *ipvt, integer *info)
{
    /* System generated locals */
    integer abe_dim1, abe_offset, i__1, i__2, i__3;
    real r__1, r__2;
    complex q__1;

    /* Local variables */
    static integer i__, j, k, l, m;
    static complex t;
    static integer n1, mb, lm, mp, ml1, lm1, lm2, ldb;
    extern /* Subroutine */ int cscal_(integer *, complex *, complex *, 
	    integer *), cswap_(integer *, complex *, integer *, complex *, 
	    integer *), caxpy_(integer *, complex *, complex *, integer *, 
	    complex *, integer *);
    extern integer icamax_(integer *, complex *, integer *);

/* ***BEGIN PROLOGUE  CNBFA */
/* ***PURPOSE  Factor a band matrix by elimination. */
/* ***LIBRARY   SLATEC */
/* ***CATEGORY  D2C2 */
/* ***TYPE      COMPLEX (SNBFA-S, DNBFA-D, CNBFA-C) */
/* ***KEYWORDS  BANDED, LINEAR EQUATIONS, MATRIX FACTORIZATION, */
/*             NONSYMMETRIC */
/* ***AUTHOR  Voorhees, E. A., (LANL) */
/* ***DESCRIPTION */

/*    CNBFA factors a complex band matrix by elimination. */

/*     CNBFA is usually called by CNBCO, but it can be called */
/*     directly with a saving in time if RCOND is not needed. */

/*     On Entry */

/*        ABE     COMPLEX(LDA, NC) */
/*                contains the matrix in band storage.  The rows */
/*                of the original matrix are stored in the rows */
/*                of ABE and the diagonals of the original matrix */
/*                are stored in columns 1 through ML+MU+1 of ABE. */
/*                NC must be .GE. 2*ML+MU+1 . */
/*                See the comments below for details. */

/*        LDA     INTEGER */
/*                the leading dimension of the array ABE. */
/*                LDA must be .GE. N . */

/*        N       INTEGER */
/*                the order of the original matrix. */

/*        ML      INTEGER */
/*                number of diagonals below the main diagonal. */
/*                0 .LE. ML .LT. N . */

/*        MU      INTEGER */
/*                number of diagonals above the main diagonal. */
/*                0 .LE. MU .LT. N . */
/*                More efficient if ML .LE. MU . */

/*     On Return */

/*        ABE     an upper triangular matrix in band storage */
/*                and the multipliers which were used to obtain it. */
/*                the factorization can be written  A = L*U  where */
/*                L is a product of permutation and unit lower */
/*                triangular matrices and  U  is upper triangular. */

/*        IPVT    INTEGER(N) */
/*                an integer vector of pivot indices. */

/*        INFO    INTEGER */
/*                =0  normal value */
/*                =K  if  U(K,K) .EQ. 0.0 .  This is not an error */
/*                condition for this subroutine, but it does */
/*                indicate that CNBSL will divide by zero if */
/*                called.  Use RCOND in CNBCO for a reliable */
/*                indication of singularity. */

/*     Band Storage */

/*           If  A  is a band matrix, the following program segment */
/*           will set up the input. */

/*                   ML = (band width below the diagonal) */
/*                   MU = (band width above the diagonal) */
/*                   DO 20 I = 1, N */
/*                      J1 = MAX(1, I-ML) */
/*                      J2 = MIN(N, I+MU) */
/*                      DO 10 J = J1, J2 */
/*                         K = J - I + ML + 1 */
/*                         ABE(I,K) = A(I,J) */
/*                10    CONTINUE */
/*                20 CONTINUE */

/*           This uses columns  1  through  ML+MU+1  of ABE . */
/*           Furthermore,  ML  additional columns are needed in */
/*           ABE  starting with column  ML+MU+2  for elements */
/*           generated during the triangularization.  The total */
/*           number of columns needed in  ABE  is  2*ML+MU+1 . */

/*     Example:  If the original matrix is */

/*           11 12 13  0  0  0 */
/*           21 22 23 24  0  0 */
/*            0 32 33 34 35  0 */
/*            0  0 43 44 45 46 */
/*            0  0  0 54 55 56 */
/*            0  0  0  0 65 66 */

/*      then  N = 6, ML = 1, MU = 2, LDA .GE. 5  and ABE should contain */

/*            * 11 12 13  +     , * = not used */
/*           21 22 23 24  +     , + = used for pivoting */
/*           32 33 34 35  + */
/*           43 44 45 46  + */
/*           54 55 56  *  + */
/*           65 66  *  *  + */

/* ***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W. */
/*                 Stewart, LINPACK Users' Guide, SIAM, 1979. */
/* ***ROUTINES CALLED  CAXPY, CSCAL, CSWAP, ICAMAX */
/* ***REVISION HISTORY  (YYMMDD) */
/*   800730  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890831  Modified array declarations.  (WRB) */
/*   890831  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   920501  Reformatted the REFERENCES section.  (WRB) */
/* ***END PROLOGUE  CNBFA */


/* ***FIRST EXECUTABLE STATEMENT  CNBFA */
    /* Parameter adjustments */
    abe_dim1 = *lda;
    abe_offset = 1 + abe_dim1;
    abe -= abe_offset;
    --ipvt;

    /* Function Body */
    ml1 = *ml + 1;
    mb = *ml + *mu;
    m = *ml + *mu + 1;
    n1 = *n - 1;
    ldb = *lda - 1;
    *info = 0;

/*     SET FILL-IN COLUMNS TO ZERO */

    if (*n <= 1) {
	goto L50;
    }
    if (*ml <= 0) {
	goto L7;
    }
    i__1 = *ml;
    for (j = 1; j <= i__1; ++j) {
	i__2 = *n;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    i__3 = i__ + (m + j) * abe_dim1;
	    abe[i__3].r = 0.f, abe[i__3].i = 0.f;
/* L5: */
	}
/* L6: */
    }
L7:

/*     GAUSSIAN ELIMINATION WITH PARTIAL ELIMINATION */

    i__1 = n1;
    for (k = 1; k <= i__1; ++k) {
/* Computing MIN */
	i__2 = *n - k;
	lm = min(i__2,*ml);
	lm1 = lm + 1;
	lm2 = ml1 - lm;

/*     SEARCH FOR PIVOT INDEX */

	l = -icamax_(&lm1, &abe[lm + k + lm2 * abe_dim1], &ldb) + lm1 + k;
	ipvt[k] = l;
/* Computing MIN */
	i__2 = mb, i__3 = *n - k;
	mp = min(i__2,i__3);

/*     SWAP ROWS IF NECESSARY */

	if (l != k) {
	    i__2 = mp + 1;
	    cswap_(&i__2, &abe[k + ml1 * abe_dim1], lda, &abe[l + (ml1 + k - 
		    l) * abe_dim1], lda);
	}

/*     SKIP COLUMN REDUCTION IF PIVOT IS ZERO */

	i__2 = k + ml1 * abe_dim1;
	if ((r__1 = abe[i__2].r, dabs(r__1)) + (r__2 = r_imag(&abe[k + ml1 * 
		abe_dim1]), dabs(r__2)) == 0.f) {
	    goto L20;
	}

/*     COMPUTE MULTIPLIERS */

	c_div(&q__1, &c_b9, &abe[k + ml1 * abe_dim1]);
	t.r = q__1.r, t.i = q__1.i;
	cscal_(&lm, &t, &abe[lm + k + lm2 * abe_dim1], &ldb);

/*     ROW ELIMINATION WITH COLUMN INDEXING */

	i__2 = mp;
	for (j = 1; j <= i__2; ++j) {
	    caxpy_(&lm, &abe[k + (ml1 + j) * abe_dim1], &abe[lm + k + lm2 * 
		    abe_dim1], &ldb, &abe[lm + k + (lm2 + j) * abe_dim1], &
		    ldb);
/* L10: */
	}
	goto L30;
L20:
	*info = k;
L30:
/* L40: */
	;
    }
L50:
    ipvt[*n] = *n;
    i__1 = *n + ml1 * abe_dim1;
    if ((r__1 = abe[i__1].r, dabs(r__1)) + (r__2 = r_imag(&abe[*n + ml1 * 
	    abe_dim1]), dabs(r__2)) == 0.f) {
	*info = *n;
    }
    return 0;
} /* cnbfa_ */

