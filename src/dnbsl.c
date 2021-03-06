/* dnbsl.f -- translated by f2c (version 12.02.01).
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

/* DECK DNBSL */
/* Subroutine */ int dnbsl_(doublereal *abe, integer *lda, integer *n, 
	integer *ml, integer *mu, integer *ipvt, doublereal *b, integer *job)
{
    /* System generated locals */
    integer abe_dim1, abe_offset, i__1, i__2, i__3;

    /* Local variables */
    static integer k, l, m;
    static doublereal t;
    static integer kb, lb, lm, nm1, ldb, mlm;
    extern doublereal ddot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    extern /* Subroutine */ int daxpy_(integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *);

/* ***BEGIN PROLOGUE  DNBSL */
/* ***PURPOSE  Solve a real band system using the factors computed by */
/*            DNBCO or DNBFA. */
/* ***LIBRARY   SLATEC */
/* ***CATEGORY  D2A2 */
/* ***TYPE      DOUBLE PRECISION (SNBSL-S, DNBSL-D, CNBSL-C) */
/* ***KEYWORDS  BANDED, LINEAR EQUATIONS, NONSYMMETRIC, SOLVE */
/* ***AUTHOR  Voorhees, E. A., (LANL) */
/* ***DESCRIPTION */

/*     DNBSL solves the double precision band system */
/*     A * X = B  or  TRANS(A) * X = B */
/*     using the factors computed by DNBCO or DNBFA. */

/*     On Entry */

/*        ABE     DOUBLE PRECISION(LDA, NC) */
/*                the output from DNBCO or DNBFA. */
/*                NC must be .GE. 2*ML+MU+1 . */

/*        LDA     INTEGER */
/*                the leading dimension of the array  ABE . */

/*        N       INTEGER */
/*                the order of the original matrix. */

/*        ML      INTEGER */
/*                number of diagonals below the main diagonal. */

/*        MU      INTEGER */
/*                number of diagonals above the main diagonal. */

/*        IPVT    INTEGER(N) */
/*                the pivot vector from DNBCO or DNBFA. */

/*        B       DOUBLE PRECISION(N) */
/*                the right hand side vector. */

/*        JOB     INTEGER */
/*                = 0         to solve  A*X = B . */
/*                = nonzero   to solve  TRANS(A)*X = B , where */
/*                            TRANS(A)  is the transpose. */

/*     On Return */

/*        B       the solution vector  X . */

/*     Error Condition */

/*        A division by zero will occur if the input factor contains a */
/*        zero on the diagonal.  Technically this indicates singularity */
/*        but it is often caused by improper arguments or improper */
/*        setting of LDA.  It will not occur if the subroutines are */
/*        called correctly and if DNBCO has set RCOND .GT. 0.0 */
/*        or DNBFA has set INFO .EQ. 0 . */

/*     To compute  INVERSE(A) * C  where  C  is a matrix */
/*     with  P  columns */
/*           CALL DNBCO(ABE,LDA,N,ML,MU,IPVT,RCOND,Z) */
/*           IF (RCOND is too small) GO TO ... */
/*           DO 10 J = 1, P */
/*             CALL DNBSL(ABE,LDA,N,ML,MU,IPVT,C(1,J),0) */
/*        10 CONTINUE */

/* ***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W. */
/*                 Stewart, LINPACK Users' Guide, SIAM, 1979. */
/* ***ROUTINES CALLED  DAXPY, DDOT */
/* ***REVISION HISTORY  (YYMMDD) */
/*   800728  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890831  Modified array declarations.  (WRB) */
/*   890831  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   920501  Reformatted the REFERENCES section.  (WRB) */
/* ***END PROLOGUE  DNBSL */

/* ***FIRST EXECUTABLE STATEMENT  DNBSL */
    /* Parameter adjustments */
    abe_dim1 = *lda;
    abe_offset = 1 + abe_dim1;
    abe -= abe_offset;
    --ipvt;
    --b;

    /* Function Body */
    m = *mu + *ml + 1;
    nm1 = *n - 1;
    ldb = 1 - *lda;
    if (*job != 0) {
	goto L50;
    }

/*       JOB = 0 , SOLVE  A * X = B */
/*       FIRST SOLVE L*Y = B */

    if (*ml == 0) {
	goto L30;
    }
    if (nm1 < 1) {
	goto L30;
    }
    i__1 = nm1;
    for (k = 1; k <= i__1; ++k) {
/* Computing MIN */
	i__2 = *ml, i__3 = *n - k;
	lm = min(i__2,i__3);
	l = ipvt[k];
	t = b[l];
	if (l == k) {
	    goto L10;
	}
	b[l] = b[k];
	b[k] = t;
L10:
	mlm = *ml - (lm - 1);
	daxpy_(&lm, &t, &abe[k + lm + mlm * abe_dim1], &ldb, &b[k + 1], &c__1)
		;
/* L20: */
    }
L30:

/*       NOW SOLVE  U*X = Y */

    i__1 = *n;
    for (kb = 1; kb <= i__1; ++kb) {
	k = *n + 1 - kb;
	b[k] /= abe[k + (*ml + 1) * abe_dim1];
	lm = min(k,m) - 1;
	lb = k - lm;
	t = -b[k];
	daxpy_(&lm, &t, &abe[k - 1 + (*ml + 2) * abe_dim1], &ldb, &b[lb], &
		c__1);
/* L40: */
    }
    goto L100;
L50:

/*       JOB = NONZERO, SOLVE TRANS(A) * X = B */
/*       FIRST SOLVE  TRANS(U)*Y = B */

    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
	lm = min(k,m) - 1;
	lb = k - lm;
	t = ddot_(&lm, &abe[k - 1 + (*ml + 2) * abe_dim1], &ldb, &b[lb], &
		c__1);
	b[k] = (b[k] - t) / abe[k + (*ml + 1) * abe_dim1];
/* L60: */
    }

/*       NOW SOLVE TRANS(L)*X = Y */

    if (*ml == 0) {
	goto L90;
    }
    if (nm1 < 1) {
	goto L90;
    }
    i__1 = nm1;
    for (kb = 1; kb <= i__1; ++kb) {
	k = *n - kb;
/* Computing MIN */
	i__2 = *ml, i__3 = *n - k;
	lm = min(i__2,i__3);
	mlm = *ml - (lm - 1);
	b[k] += ddot_(&lm, &abe[k + lm + mlm * abe_dim1], &ldb, &b[k + 1], &
		c__1);
	l = ipvt[k];
	if (l == k) {
	    goto L70;
	}
	t = b[l];
	b[l] = b[k];
	b[k] = t;
L70:
/* L80: */
	;
    }
L90:
L100:
    return 0;
} /* dnbsl_ */

