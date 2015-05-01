/* sinqf.f -- translated by f2c (version 12.02.01).
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

/* DECK SINQF */
/* Subroutine */ int sinqf_(integer *n, real *x, real *wsave)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer k, kc, ns2;
    extern /* Subroutine */ int cosqf_(integer *, real *, real *);
    static real xhold;

/* ***BEGIN PROLOGUE  SINQF */
/* ***PURPOSE  Compute the forward sine transform with odd wave numbers. */
/* ***LIBRARY   SLATEC (FFTPACK) */
/* ***CATEGORY  J1A3 */
/* ***TYPE      SINGLE PRECISION (SINQF-S) */
/* ***KEYWORDS  FFTPACK, FOURIER TRANSFORM */
/* ***AUTHOR  Swarztrauber, P. N., (NCAR) */
/* ***DESCRIPTION */

/*  Subroutine SINQF computes the fast Fourier transform of quarter */
/*  wave data.  That is, SINQF computes the coefficients in a sine */
/*  series representation with only odd wave numbers.  The transform */
/*  is defined below at output parameter X. */

/*  SINQB is the unnormalized inverse of SINQF since a call of SINQF */
/*  followed by a call of SINQB will multiply the input sequence X */
/*  by 4*N. */

/*  The array WSAVE which is used by subroutine SINQF must be */
/*  initialized by calling subroutine SINQI(N,WSAVE). */

/*  Input Parameters */

/*  N       the length of the array X to be transformed.  The method */
/*          is most efficient when N is a product of small primes. */

/*  X       an array which contains the sequence to be transformed */

/*  WSAVE   a work array which must be dimensioned at least 3*N+15 */
/*          in the program that calls SINQF.  The WSAVE array must be */
/*          initialized by calling subroutine SINQI(N,WSAVE), and a */
/*          different WSAVE array must be used for each different */
/*          value of N.  This initialization does not have to be */
/*          repeated so long as N remains unchanged.  Thus subsequent */
/*          transforms can be obtained faster than the first. */

/*  Output Parameters */

/*  X       For I=1,...,N */

/*               X(I) = (-1)**(I-1)*X(N) */

/*                  + the sum from K=1 to K=N-1 of */

/*                  2*X(K)*SIN((2*I-1)*K*PI/(2*N)) */

/*               A call of SINQF followed by a call of */
/*               SINQB will multiply the sequence X by 4*N. */
/*               Therefore SINQB is the unnormalized inverse */
/*               of SINQF. */

/*  WSAVE   contains initialization calculations which must not */
/*          be destroyed between calls of SINQF or SINQB. */

/* ***REFERENCES  P. N. Swarztrauber, Vectorizing the FFTs, in Parallel */
/*                 Computations (G. Rodrigue, ed.), Academic Press, */
/*                 1982, pp. 51-83. */
/* ***ROUTINES CALLED  COSQF */
/* ***REVISION HISTORY  (YYMMDD) */
/*   790601  DATE WRITTEN */
/*   830401  Modified to use SLATEC library source file format. */
/*   860115  Modified by Ron Boisvert to adhere to Fortran 77 by */
/*           changing dummy array size declarations (1) to (*) */
/*   861211  REVISION DATE from Version 3.2 */
/*   881128  Modified by Dick Valent to meet prologue standards. */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   920501  Reformatted the REFERENCES section.  (WRB) */
/* ***END PROLOGUE  SINQF */
/* ***FIRST EXECUTABLE STATEMENT  SINQF */
    /* Parameter adjustments */
    --wsave;
    --x;

    /* Function Body */
    if (*n == 1) {
	return 0;
    }
    ns2 = *n / 2;
    i__1 = ns2;
    for (k = 1; k <= i__1; ++k) {
	kc = *n - k;
	xhold = x[k];
	x[k] = x[kc + 1];
	x[kc + 1] = xhold;
/* L101: */
    }
    cosqf_(n, &x[1], &wsave[1]);
    i__1 = *n;
    for (k = 2; k <= i__1; k += 2) {
	x[k] = -x[k];
/* L102: */
    }
    return 0;
} /* sinqf_ */

