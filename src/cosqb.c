/* cosqb.f -- translated by f2c (version 12.02.01).
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

/* DECK COSQB */
/* Subroutine */ int cosqb_(integer *n, real *x, real *wsave)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static real x1;
    extern /* Subroutine */ int cosqb1_(integer *, real *, real *, real *);
    static real tsqrt2;

/* ***BEGIN PROLOGUE  COSQB */
/* ***PURPOSE  Compute the unnormalized inverse cosine transform. */
/* ***LIBRARY   SLATEC (FFTPACK) */
/* ***CATEGORY  J1A3 */
/* ***TYPE      SINGLE PRECISION (COSQB-S) */
/* ***KEYWORDS  FFTPACK, INVERSE COSINE FOURIER TRANSFORM */
/* ***AUTHOR  Swarztrauber, P. N., (NCAR) */
/* ***DESCRIPTION */

/*  Subroutine COSQB computes the fast Fourier transform of quarter */
/*  wave data. That is, COSQB computes a sequence from its */
/*  representation in terms of a cosine series with odd wave numbers. */
/*  The transform is defined below at output parameter X. */

/*  COSQB is the unnormalized inverse of COSQF since a call of COSQB */
/*  followed by a call of COSQF will multiply the input sequence X */
/*  by 4*N. */

/*  The array WSAVE which is used by subroutine COSQB must be */
/*  initialized by calling subroutine COSQI(N,WSAVE). */


/*  Input Parameters */

/*  N       the length of the array X to be transformed.  The method */
/*          is most efficient when N is a product of small primes. */

/*  X       an array which contains the sequence to be transformed */

/*  WSAVE   a work array which must be dimensioned at least 3*N+15 */
/*          in the program that calls COSQB.  The WSAVE array must be */
/*          initialized by calling subroutine COSQI(N,WSAVE), and a */
/*          different WSAVE array must be used for each different */
/*          value of N.  This initialization does not have to be */
/*          repeated so long as N remains unchanged.  Thus subsequent */
/*          transforms can be obtained faster than the first. */

/*  Output Parameters */

/*  X       For I=1,...,N */

/*               X(I)= the sum from K=1 to K=N of */

/*                  2*X(K)*COS((2*K-1)*(I-1)*PI/(2*N)) */

/*               A call of COSQB followed by a call of */
/*               COSQF will multiply the sequence X by 4*N. */
/*               Therefore COSQF is the unnormalized inverse */
/*               of COSQB. */

/*  WSAVE   contains initialization calculations which must not */
/*          be destroyed between calls of COSQB or COSQF. */

/* ***REFERENCES  P. N. Swarztrauber, Vectorizing the FFTs, in Parallel */
/*                 Computations (G. Rodrigue, ed.), Academic Press, */
/*                 1982, pp. 51-83. */
/* ***ROUTINES CALLED  COSQB1 */
/* ***REVISION HISTORY  (YYMMDD) */
/*   790601  DATE WRITTEN */
/*   830401  Modified to use SLATEC library source file format. */
/*   860115  Modified by Ron Boisvert to adhere to Fortran 77 by */
/*           (a) changing dummy array size declarations (1) to (*), */
/*           (b) changing definition of variable TSQRT2 by using */
/*               FORTRAN intrinsic function SQRT instead of a DATA */
/*               statement. */
/*   861211  REVISION DATE from Version 3.2 */
/*   881128  Modified by Dick Valent to meet prologue standards. */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   920501  Reformatted the REFERENCES section.  (WRB) */
/* ***END PROLOGUE  COSQB */
/* ***FIRST EXECUTABLE STATEMENT  COSQB */
    /* Parameter adjustments */
    --wsave;
    --x;

    /* Function Body */
    tsqrt2 = sqrt(2.f) * 2.f;
    if ((i__1 = *n - 2) < 0) {
	goto L101;
    } else if (i__1 == 0) {
	goto L102;
    } else {
	goto L103;
    }
L101:
    x[1] *= 4.f;
    return 0;
L102:
    x1 = (x[1] + x[2]) * 4.f;
    x[2] = tsqrt2 * (x[1] - x[2]);
    x[1] = x1;
    return 0;
L103:
    cosqb1_(n, &x[1], &wsave[1], &wsave[*n + 1]);
    return 0;
} /* cosqb_ */

