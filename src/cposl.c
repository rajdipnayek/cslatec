/* cposl.f -- translated by f2c (version 12.02.01).
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

/* DECK CPOSL */
/* Subroutine */ int cposl_(complex *a, integer *lda, integer *n, complex *b)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3;
    complex q__1, q__2;

    /* Local variables */
    static integer k;
    static complex t;
    static integer kb;
    extern /* Complex */ void cdotc_(complex *, integer *, complex *, integer 
	    *, complex *, integer *);
    extern /* Subroutine */ int caxpy_(integer *, complex *, complex *, 
	    integer *, complex *, integer *);

/* ***BEGIN PROLOGUE  CPOSL */
/* ***PURPOSE  Solve the complex Hermitian positive definite linear system */
/*            using the factors computed by CPOCO or CPOFA. */
/* ***LIBRARY   SLATEC (LINPACK) */
/* ***CATEGORY  D2D1B */
/* ***TYPE      COMPLEX (SPOSL-S, DPOSL-D, CPOSL-C) */
/* ***KEYWORDS  LINEAR ALGEBRA, LINPACK, MATRIX, POSITIVE DEFINITE, SOLVE */
/* ***AUTHOR  Moler, C. B., (U. of New Mexico) */
/* ***DESCRIPTION */

/*     CPOSL solves the COMPLEX Hermitian positive definite system */
/*     A * X = B */
/*     using the factors computed by CPOCO or CPOFA. */

/*     On Entry */

/*        A       COMPLEX(LDA, N) */
/*                the output from CPOCO or CPOFA. */

/*        LDA     INTEGER */
/*                the leading dimension of the array  A . */

/*        N       INTEGER */
/*                the order of the matrix  A . */

/*        B       COMPLEX(N) */
/*                the right hand side vector. */

/*     On Return */

/*        B       the solution vector  X . */

/*     Error Condition */

/*        A division by zero will occur if the input factor contains */
/*        a zero on the diagonal.  Technically this indicates */
/*        singularity but it is usually caused by improper subroutine */
/*        arguments.  It will not occur if the subroutines are called */
/*        correctly and  INFO .EQ. 0 . */

/*     To compute  INVERSE(A) * C  where  C  is a matrix */
/*     with  P  columns */
/*           CALL CPOCO(A,LDA,N,RCOND,Z,INFO) */
/*           IF (RCOND is too small .OR. INFO .NE. 0) GO TO ... */
/*           DO 10 J = 1, P */
/*              CALL CPOSL(A,LDA,N,C(1,J)) */
/*        10 CONTINUE */

/* ***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W. */
/*                 Stewart, LINPACK Users' Guide, SIAM, 1979. */
/* ***ROUTINES CALLED  CAXPY, CDOTC */
/* ***REVISION HISTORY  (YYMMDD) */
/*   780814  DATE WRITTEN */
/*   890831  Modified array declarations.  (WRB) */
/*   890831  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900326  Removed duplicate information from DESCRIPTION section. */
/*           (WRB) */
/*   920501  Reformatted the REFERENCES section.  (WRB) */
/* ***END PROLOGUE  CPOSL */


/*     SOLVE CTRANS(R)*Y = B */

/* ***FIRST EXECUTABLE STATEMENT  CPOSL */
    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --b;

    /* Function Body */
    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
	i__2 = k - 1;
	cdotc_(&q__1, &i__2, &a[k * a_dim1 + 1], &c__1, &b[1], &c__1);
	t.r = q__1.r, t.i = q__1.i;
	i__2 = k;
	i__3 = k;
	q__2.r = b[i__3].r - t.r, q__2.i = b[i__3].i - t.i;
	c_div(&q__1, &q__2, &a[k + k * a_dim1]);
	b[i__2].r = q__1.r, b[i__2].i = q__1.i;
/* L10: */
    }

/*     SOLVE R*X = Y */

    i__1 = *n;
    for (kb = 1; kb <= i__1; ++kb) {
	k = *n + 1 - kb;
	i__2 = k;
	c_div(&q__1, &b[k], &a[k + k * a_dim1]);
	b[i__2].r = q__1.r, b[i__2].i = q__1.i;
	i__2 = k;
	q__1.r = -b[i__2].r, q__1.i = -b[i__2].i;
	t.r = q__1.r, t.i = q__1.i;
	i__2 = k - 1;
	caxpy_(&i__2, &t, &a[k * a_dim1 + 1], &c__1, &b[1], &c__1);
/* L20: */
    }
    return 0;
} /* cposl_ */

