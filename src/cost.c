/* cost.f -- translated by f2c (version 12.02.01).
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

/* DECK COST */
/* Subroutine */ int cost_(integer *n, real *x, real *wsave)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__, k;
    static real c1, t1, t2;
    static integer kc;
    static real xi;
    static integer nm1, np1;
    static real x1h;
    static integer ns2;
    static real tx2, x1p3, xim2;
    static integer modn;
    extern /* Subroutine */ int rfftf_(integer *, real *, real *);

/* ***BEGIN PROLOGUE  COST */
/* ***PURPOSE  Compute the cosine transform of a real, even sequence. */
/* ***LIBRARY   SLATEC (FFTPACK) */
/* ***CATEGORY  J1A3 */
/* ***TYPE      SINGLE PRECISION (COST-S) */
/* ***KEYWORDS  COSINE FOURIER TRANSFORM, FFTPACK */
/* ***AUTHOR  Swarztrauber, P. N., (NCAR) */
/* ***DESCRIPTION */

/*  Subroutine COST computes the discrete Fourier cosine transform */
/*  of an even sequence X(I).  The transform is defined below at output */
/*  parameter X. */

/*  COST is the unnormalized inverse of itself since a call of COST */
/*  followed by another call of COST will multiply the input sequence */
/*  X by 2*(N-1).  The transform is defined below at output parameter X. */

/*  The array WSAVE which is used by subroutine COST must be */
/*  initialized by calling subroutine COSTI(N,WSAVE). */

/*  Input Parameters */

/*  N       the length of the sequence X.  N must be greater than 1. */
/*          The method is most efficient when N-1 is a product of */
/*          small primes. */

/*  X       an array which contains the sequence to be transformed */

/*  WSAVE   a work array which must be dimensioned at least 3*N+15 */
/*          in the program that calls COST.  The WSAVE array must be */
/*          initialized by calling subroutine COSTI(N,WSAVE), and a */
/*          different WSAVE array must be used for each different */
/*          value of N.  This initialization does not have to be */
/*          repeated so long as N remains unchanged.  Thus subsequent */
/*          transforms can be obtained faster than the first. */

/*  Output Parameters */

/*  X       For I=1,...,N */

/*             X(I) = X(1)+(-1)**(I-1)*X(N) */

/*               + the sum from K=2 to K=N-1 */

/*                 2*X(K)*COS((K-1)*(I-1)*PI/(N-1)) */

/*               A call of COST followed by another call of */
/*               COST will multiply the sequence X by 2*(N-1). */
/*               Hence COST is the unnormalized inverse */
/*               of itself. */

/*  WSAVE   contains initialization calculations which must not be */
/*          destroyed between calls of COST. */

/* ***REFERENCES  P. N. Swarztrauber, Vectorizing the FFTs, in Parallel */
/*                 Computations (G. Rodrigue, ed.), Academic Press, */
/*                 1982, pp. 51-83. */
/* ***ROUTINES CALLED  RFFTF */
/* ***REVISION HISTORY  (YYMMDD) */
/*   790601  DATE WRITTEN */
/*   830401  Modified to use SLATEC library source file format. */
/*   860115  Modified by Ron Boisvert to adhere to Fortran 77 by */
/*           changing dummy array size declarations (1) to (*) */
/*   861211  REVISION DATE from Version 3.2 */
/*   881128  Modified by Dick Valent to meet prologue standards. */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   920501  Reformatted the REFERENCES section.  (WRB) */
/* ***END PROLOGUE  COST */
/* ***FIRST EXECUTABLE STATEMENT  COST */
    /* Parameter adjustments */
    --wsave;
    --x;

    /* Function Body */
    nm1 = *n - 1;
    np1 = *n + 1;
    ns2 = *n / 2;
    if ((i__1 = *n - 2) < 0) {
	goto L106;
    } else if (i__1 == 0) {
	goto L101;
    } else {
	goto L102;
    }
L101:
    x1h = x[1] + x[2];
    x[2] = x[1] - x[2];
    x[1] = x1h;
    return 0;
L102:
    if (*n > 3) {
	goto L103;
    }
    x1p3 = x[1] + x[3];
    tx2 = x[2] + x[2];
    x[2] = x[1] - x[3];
    x[1] = x1p3 + tx2;
    x[3] = x1p3 - tx2;
    return 0;
L103:
    c1 = x[1] - x[*n];
    x[1] += x[*n];
    i__1 = ns2;
    for (k = 2; k <= i__1; ++k) {
	kc = np1 - k;
	t1 = x[k] + x[kc];
	t2 = x[k] - x[kc];
	c1 += wsave[kc] * t2;
	t2 = wsave[k] * t2;
	x[k] = t1 - t2;
	x[kc] = t1 + t2;
/* L104: */
    }
    modn = *n % 2;
    if (modn != 0) {
	x[ns2 + 1] += x[ns2 + 1];
    }
    rfftf_(&nm1, &x[1], &wsave[*n + 1]);
    xim2 = x[2];
    x[2] = c1;
    i__1 = *n;
    for (i__ = 4; i__ <= i__1; i__ += 2) {
	xi = x[i__];
	x[i__] = x[i__ - 2] - x[i__ - 1];
	x[i__ - 1] = xim2;
	xim2 = xi;
/* L105: */
    }
    if (modn != 0) {
	x[*n] = xim2;
    }
L106:
    return 0;
} /* cost_ */

