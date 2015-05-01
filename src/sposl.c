/* sposl.f -- translated by f2c (version 12.02.01).
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

/* DECK SPOSL */
/* Subroutine */ int sposl_(real *a, integer *lda, integer *n, real *b)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;

    /* Local variables */
    static integer k;
    static real t;
    static integer kb;
    extern doublereal sdot_(integer *, real *, integer *, real *, integer *);
    extern /* Subroutine */ int saxpy_(integer *, real *, real *, integer *, 
	    real *, integer *);

/* ***BEGIN PROLOGUE  SPOSL */
/* ***PURPOSE  Solve the real symmetric positive definite linear system */
/*            using the factors computed by SPOCO or SPOFA. */
/* ***LIBRARY   SLATEC (LINPACK) */
/* ***CATEGORY  D2B1B */
/* ***TYPE      SINGLE PRECISION (SPOSL-S, DPOSL-D, CPOSL-C) */
/* ***KEYWORDS  LINEAR ALGEBRA, LINPACK, MATRIX, POSITIVE DEFINITE, SOLVE */
/* ***AUTHOR  Moler, C. B., (U. of New Mexico) */
/* ***DESCRIPTION */

/*     SPOSL solves the real symmetric positive definite system */
/*     A * X = B */
/*     using the factors computed by SPOCO or SPOFA. */

/*     On Entry */

/*        A       REAL(LDA, N) */
/*                the output from SPOCO or SPOFA. */

/*        LDA     INTEGER */
/*                the leading dimension of the array  A . */

/*        N       INTEGER */
/*                the order of the matrix  A . */

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
/*           CALL SPOCO(A,LDA,N,RCOND,Z,INFO) */
/*           IF (RCOND is too small .OR. INFO .NE. 0) GO TO ... */
/*           DO 10 J = 1, P */
/*              CALL SPOSL(A,LDA,N,C(1,J)) */
/*        10 CONTINUE */

/* ***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W. */
/*                 Stewart, LINPACK Users' Guide, SIAM, 1979. */
/* ***ROUTINES CALLED  SAXPY, SDOT */
/* ***REVISION HISTORY  (YYMMDD) */
/*   780814  DATE WRITTEN */
/*   890831  Modified array declarations.  (WRB) */
/*   890831  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900326  Removed duplicate information from DESCRIPTION section. */
/*           (WRB) */
/*   920501  Reformatted the REFERENCES section.  (WRB) */
/* ***END PROLOGUE  SPOSL */


/*     SOLVE TRANS(R)*Y = B */

/* ***FIRST EXECUTABLE STATEMENT  SPOSL */
    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --b;

    /* Function Body */
    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
	i__2 = k - 1;
	t = sdot_(&i__2, &a[k * a_dim1 + 1], &c__1, &b[1], &c__1);
	b[k] = (b[k] - t) / a[k + k * a_dim1];
/* L10: */
    }

/*     SOLVE R*X = Y */

    i__1 = *n;
    for (kb = 1; kb <= i__1; ++kb) {
	k = *n + 1 - kb;
	b[k] /= a[k + k * a_dim1];
	t = -b[k];
	i__2 = k - 1;
	saxpy_(&i__2, &t, &a[k * a_dim1 + 1], &c__1, &b[1], &c__1);
/* L20: */
    }
    return 0;
} /* sposl_ */

