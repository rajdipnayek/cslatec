/* cosqi.f -- translated by f2c (version 12.02.01).
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

/* DECK COSQI */
/* Subroutine */ int cosqi_(integer *n, real *wsave)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer k;
    static real fk, dt, pih;
    extern /* Subroutine */ int rffti_(integer *, real *);

/* ***BEGIN PROLOGUE  COSQI */
/* ***PURPOSE  Initialize a work array for COSQF and COSQB. */
/* ***LIBRARY   SLATEC (FFTPACK) */
/* ***CATEGORY  J1A3 */
/* ***TYPE      SINGLE PRECISION (COSQI-S) */
/* ***KEYWORDS  COSINE FOURIER TRANSFORM, FFTPACK */
/* ***AUTHOR  Swarztrauber, P. N., (NCAR) */
/* ***DESCRIPTION */

/*  Subroutine COSQI initializes the work array WSAVE which is used in */
/*  both COSQF1 and COSQB1.  The prime factorization of N together with */
/*  a tabulation of the trigonometric functions are computed and */
/*  stored in WSAVE. */

/*  Input Parameter */

/*  N       the length of the array to be transformed.  The method */
/*          is most efficient when N is a product of small primes. */

/*  Output Parameter */

/*  WSAVE   a work array which must be dimensioned at least 3*N+15. */
/*          The same work array can be used for both COSQF1 and COSQB1 */
/*          as long as N remains unchanged.  Different WSAVE arrays */
/*          are required for different values of N.  The contents of */
/*          WSAVE must not be changed between calls of COSQF1 or COSQB1. */

/* ***REFERENCES  P. N. Swarztrauber, Vectorizing the FFTs, in Parallel */
/*                 Computations (G. Rodrigue, ed.), Academic Press, */
/*                 1982, pp. 51-83. */
/* ***ROUTINES CALLED  RFFTI */
/* ***REVISION HISTORY  (YYMMDD) */
/*   790601  DATE WRITTEN */
/*   830401  Modified to use SLATEC library source file format. */
/*   860115  Modified by Ron Boisvert to adhere to Fortran 77 by */
/*           (a) changing dummy array size declarations (1) to (*), */
/*           (b) changing references to intrinsic function FLOAT */
/*               to REAL, and */
/*           (c) changing definition of variable PIH by using */
/*               FORTRAN intrinsic function ATAN instead of a DATA */
/*               statement. */
/*   881128  Modified by Dick Valent to meet prologue standards. */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890531  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   920501  Reformatted the REFERENCES section.  (WRB) */
/* ***END PROLOGUE  COSQI */
/* ***FIRST EXECUTABLE STATEMENT  COSQI */
    /* Parameter adjustments */
    --wsave;

    /* Function Body */
    pih = atan(1.f) * 2.f;
    dt = pih / *n;
    fk = 0.f;
    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
	fk += 1.f;
	wsave[k] = cos(fk * dt);
/* L101: */
    }
    rffti_(n, &wsave[*n + 1]);
    return 0;
} /* cosqi_ */

