/* dposl.f -- translated by f2c (version 12.02.01).
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

/* DECK DPOSL */
/* Subroutine */ int dposl_(doublereal *a, integer *lda, integer *n, 
	doublereal *b)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;

    /* Local variables */
    static integer k;
    static doublereal t;
    static integer kb;
    extern doublereal ddot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    extern /* Subroutine */ int daxpy_(integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *);

/* ***BEGIN PROLOGUE  DPOSL */
/* ***PURPOSE  Solve the real symmetric positive definite linear system */
/*            using the factors computed by DPOCO or DPOFA. */
/* ***LIBRARY   SLATEC (LINPACK) */
/* ***CATEGORY  D2B1B */
/* ***TYPE      DOUBLE PRECISION (SPOSL-S, DPOSL-D, CPOSL-C) */
/* ***KEYWORDS  LINEAR ALGEBRA, LINPACK, MATRIX, POSITIVE DEFINITE, SOLVE */
/* ***AUTHOR  Moler, C. B., (U. of New Mexico) */
/* ***DESCRIPTION */

/*     DPOSL solves the double precision symmetric positive definite */
/*     system A * X = B */
/*     using the factors computed by DPOCO or DPOFA. */

/*     On Entry */

/*        A       DOUBLE PRECISION(LDA, N) */
/*                the output from DPOCO or DPOFA. */

/*        LDA     INTEGER */
/*                the leading dimension of the array  A . */

/*        N       INTEGER */
/*                the order of the matrix  A . */

/*        B       DOUBLE PRECISION(N) */
/*                the right hand side vector. */

/*     On Return */

/*        B       the solution vector  X . */

/*     Error Condition */

/*        A division by zero will occur if the input factor contains */
/*        a zero on the diagonal.  Technically this indicates */
/*        singularity, but it is usually caused by improper subroutine */
/*        arguments.  It will not occur if the subroutines are called */
/*        correctly and  INFO .EQ. 0 . */

/*     To compute  INVERSE(A) * C  where  C  is a matrix */
/*     with  P  columns */
/*           CALL DPOCO(A,LDA,N,RCOND,Z,INFO) */
/*           IF (RCOND is too small .OR. INFO .NE. 0) GO TO ... */
/*           DO 10 J = 1, P */
/*              CALL DPOSL(A,LDA,N,C(1,J)) */
/*        10 CONTINUE */

/* ***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W. */
/*                 Stewart, LINPACK Users' Guide, SIAM, 1979. */
/* ***ROUTINES CALLED  DAXPY, DDOT */
/* ***REVISION HISTORY  (YYMMDD) */
/*   780814  DATE WRITTEN */
/*   890831  Modified array declarations.  (WRB) */
/*   890831  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900326  Removed duplicate information from DESCRIPTION section. */
/*           (WRB) */
/*   920501  Reformatted the REFERENCES section.  (WRB) */
/* ***END PROLOGUE  DPOSL */


/*     SOLVE TRANS(R)*Y = B */

/* ***FIRST EXECUTABLE STATEMENT  DPOSL */
    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --b;

    /* Function Body */
    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
	i__2 = k - 1;
	t = ddot_(&i__2, &a[k * a_dim1 + 1], &c__1, &b[1], &c__1);
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
	daxpy_(&i__2, &t, &a[k * a_dim1 + 1], &c__1, &b[1], &c__1);
/* L20: */
    }
    return 0;
} /* dposl_ */

