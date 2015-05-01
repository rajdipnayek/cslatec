/* spbsl.f -- translated by f2c (version 12.02.01).
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

/* DECK SPBSL */
/* Subroutine */ int spbsl_(real *abd, integer *lda, integer *n, integer *m, 
	real *b)
{
    /* System generated locals */
    integer abd_dim1, abd_offset, i__1, i__2;

    /* Local variables */
    static integer k;
    static real t;
    static integer kb, la, lb, lm;
    extern doublereal sdot_(integer *, real *, integer *, real *, integer *);
    extern /* Subroutine */ int saxpy_(integer *, real *, real *, integer *, 
	    real *, integer *);

/* ***BEGIN PROLOGUE  SPBSL */
/* ***PURPOSE  Solve a real symmetric positive definite band system */
/*            using the factors computed by SPBCO or SPBFA. */
/* ***LIBRARY   SLATEC (LINPACK) */
/* ***CATEGORY  D2B2 */
/* ***TYPE      SINGLE PRECISION (SPBSL-S, DPBSL-D, CPBSL-C) */
/* ***KEYWORDS  BANDED, LINEAR ALGEBRA, LINPACK, MATRIX, */
/*             POSITIVE DEFINITE, SOLVE */
/* ***AUTHOR  Moler, C. B., (U. of New Mexico) */
/* ***DESCRIPTION */

/*     SPBSL solves the real symmetric positive definite band */
/*     system  A*X = B */
/*     using the factors computed by SPBCO or SPBFA. */

/*     On Entry */

/*        ABD     REAL(LDA, N) */
/*                the output from SPBCO or SPBFA. */

/*        LDA     INTEGER */
/*                the leading dimension of the array  ABD . */

/*        N       INTEGER */
/*                the order of the matrix  A . */

/*        M       INTEGER */
/*                the number of diagonals above the main diagonal. */

/*        B       REAL(N) */
/*                the right hand side vector. */

/*     On Return */

/*        B       the solution vector  X . */

/*     Error Condition */

/*        A division by zero will occur if the input factor contains */
/*        a zero on the diagonal.  Technically, this indicates */
/*        singularity, but it is usually caused by improper subroutine */
/*        arguments.  It will not occur if the subroutines are called */
/*        correctly and  INFO .EQ. 0 . */

/*     To compute  INVERSE(A) * C  where  C  is a matrix */
/*     with  P  columns */
/*           CALL SPBCO(ABD,LDA,N,RCOND,Z,INFO) */
/*           IF (RCOND is too small .OR. INFO .NE. 0) GO TO ... */
/*           DO 10 J = 1, P */
/*              CALL SPBSL(ABD,LDA,N,C(1,J)) */
/*        10 CONTINUE */

/* ***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W. */
/*                 Stewart, LINPACK Users' Guide, SIAM, 1979. */
/* ***ROUTINES CALLED  SAXPY, SDOT */
/* ***REVISION HISTORY  (YYMMDD) */
/*   780814  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890831  Modified array declarations.  (WRB) */
/*   890831  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900326  Removed duplicate information from DESCRIPTION section. */
/*           (WRB) */
/*   920501  Reformatted the REFERENCES section.  (WRB) */
/* ***END PROLOGUE  SPBSL */


/*     SOLVE TRANS(R)*Y = B */

/* ***FIRST EXECUTABLE STATEMENT  SPBSL */
    /* Parameter adjustments */
    abd_dim1 = *lda;
    abd_offset = 1 + abd_dim1;
    abd -= abd_offset;
    --b;

    /* Function Body */
    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
/* Computing MIN */
	i__2 = k - 1;
	lm = min(i__2,*m);
	la = *m + 1 - lm;
	lb = k - lm;
	t = sdot_(&lm, &abd[la + k * abd_dim1], &c__1, &b[lb], &c__1);
	b[k] = (b[k] - t) / abd[*m + 1 + k * abd_dim1];
/* L10: */
    }

/*     SOLVE R*X = Y */

    i__1 = *n;
    for (kb = 1; kb <= i__1; ++kb) {
	k = *n + 1 - kb;
/* Computing MIN */
	i__2 = k - 1;
	lm = min(i__2,*m);
	la = *m + 1 - lm;
	lb = k - lm;
	b[k] /= abd[*m + 1 + k * abd_dim1];
	t = -b[k];
	saxpy_(&lm, &t, &abd[la + k * abd_dim1], &c__1, &b[lb], &c__1);
/* L20: */
    }
    return 0;
} /* spbsl_ */

