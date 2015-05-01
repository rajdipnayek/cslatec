/* rpzero.f -- translated by f2c (version 12.02.01).
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

/* DECK RPZERO */
/* Subroutine */ int rpzero_(integer *n, real *a, complex *r__, complex *t, 
	integer *iflg, real *s)
{
    /* System generated locals */
    integer i__1, i__2, i__3;
    complex q__1;

    /* Local variables */
    static integer i__, n1;
    extern /* Subroutine */ int cpzero_(integer *, complex *, complex *, 
	    complex *, integer *, real *);

/* ***BEGIN PROLOGUE  RPZERO */
/* ***PURPOSE  Find the zeros of a polynomial with real coefficients. */
/* ***LIBRARY   SLATEC */
/* ***CATEGORY  F1A1A */
/* ***TYPE      SINGLE PRECISION (RPZERO-S, CPZERO-C) */
/* ***KEYWORDS  POLYNOMIAL ROOTS, POLYNOMIAL ZEROS, REAL ROOTS */
/* ***AUTHOR  Kahaner, D. K., (NBS) */
/* ***DESCRIPTION */

/*      Find the zeros of the real polynomial */
/*         P(X)= A(1)*X**N + A(2)*X**(N-1) +...+ A(N+1) */

/*    Input... */
/*       N = degree of P(X) */
/*       A = real vector containing coefficients of P(X), */
/*            A(I) = coefficient of X**(N+1-I) */
/*       R = N word complex vector containing initial estimates for zeros */
/*            if these are known. */
/*       T = 6(N+1) word array used for temporary storage */
/*       IFLG = flag to indicate if initial estimates of */
/*              zeros are input. */
/*            If IFLG .EQ. 0, no estimates are input. */
/*            If IFLG .NE. 0, the vector R contains estimates of */
/*               the zeros */
/*       ** Warning ****** If estimates are input, they must */
/*                         be separated; that is, distinct or */
/*                         not repeated. */
/*       S = an N word array */

/*    Output... */
/*       R(I) = ith zero, */
/*       S(I) = bound for R(I) . */
/*       IFLG = error diagnostic */
/*    Error Diagnostics... */
/*       If IFLG .EQ. 0 on return, all is well. */
/*       If IFLG .EQ. 1 on return, A(1)=0.0 or N=0 on input. */
/*       If IFLG .EQ. 2 on return, the program failed to converge */
/*                after 25*N iterations.  Best current estimates of the */
/*                zeros are in R(I).  Error bounds are not calculated. */

/* ***REFERENCES  (NONE) */
/* ***ROUTINES CALLED  CPZERO */
/* ***REVISION HISTORY  (YYMMDD) */
/*   810223  DATE WRITTEN */
/*   890206  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/* ***END PROLOGUE  RPZERO */

/* ***FIRST EXECUTABLE STATEMENT  RPZERO */
    /* Parameter adjustments */
    --s;
    --t;
    --r__;
    --a;

    /* Function Body */
    n1 = *n + 1;
    i__1 = n1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = i__;
	i__3 = i__;
	q__1.r = a[i__3], q__1.i = 0.f;
	t[i__2].r = q__1.r, t[i__2].i = q__1.i;
/* L1: */
    }
    cpzero_(n, &t[1], &r__[1], &t[*n + 2], iflg, &s[1]);
    return 0;
} /* rpzero_ */

