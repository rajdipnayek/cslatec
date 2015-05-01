/* dppsl.f -- translated by f2c (version 12.02.01).
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

/* DECK DPPSL */
/* Subroutine */ int dppsl_(doublereal *ap, integer *n, doublereal *b)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer k;
    static doublereal t;
    static integer kb, kk;
    extern doublereal ddot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    extern /* Subroutine */ int daxpy_(integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *);

/* ***BEGIN PROLOGUE  DPPSL */
/* ***PURPOSE  Solve the real symmetric positive definite system using */
/*            the factors computed by DPPCO or DPPFA. */
/* ***LIBRARY   SLATEC (LINPACK) */
/* ***CATEGORY  D2B1B */
/* ***TYPE      DOUBLE PRECISION (SPPSL-S, DPPSL-D, CPPSL-C) */
/* ***KEYWORDS  LINEAR ALGEBRA, LINPACK, MATRIX, PACKED, */
/*             POSITIVE DEFINITE, SOLVE */
/* ***AUTHOR  Moler, C. B., (U. of New Mexico) */
/* ***DESCRIPTION */

/*     DPPSL solves the double precision symmetric positive definite */
/*     system A * X = B */
/*     using the factors computed by DPPCO or DPPFA. */

/*     On Entry */

/*        AP      DOUBLE PRECISION (N*(N+1)/2) */
/*                the output from DPPCO or DPPFA. */

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
/*           CALL DPPCO(AP,N,RCOND,Z,INFO) */
/*           IF (RCOND is too small .OR. INFO .NE. 0) GO TO ... */
/*           DO 10 J = 1, P */
/*              CALL DPPSL(AP,N,C(1,J)) */
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
/* ***END PROLOGUE  DPPSL */

/* ***FIRST EXECUTABLE STATEMENT  DPPSL */
    /* Parameter adjustments */
    --b;
    --ap;

    /* Function Body */
    kk = 0;
    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
	i__2 = k - 1;
	t = ddot_(&i__2, &ap[kk + 1], &c__1, &b[1], &c__1);
	kk += k;
	b[k] = (b[k] - t) / ap[kk];
/* L10: */
    }
    i__1 = *n;
    for (kb = 1; kb <= i__1; ++kb) {
	k = *n + 1 - kb;
	b[k] /= ap[kk];
	kk -= k;
	t = -b[k];
	i__2 = k - 1;
	daxpy_(&i__2, &t, &ap[kk + 1], &c__1, &b[1], &c__1);
/* L20: */
    }
    return 0;
} /* dppsl_ */

