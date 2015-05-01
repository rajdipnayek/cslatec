/* sinqi.f -- translated by f2c (version 12.02.01).
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

/* DECK SINQI */
/* Subroutine */ int sinqi_(integer *n, real *wsave)
{
    extern /* Subroutine */ int cosqi_(integer *, real *);

/* ***BEGIN PROLOGUE  SINQI */
/* ***PURPOSE  Initialize a work array for SINQF and SINQB. */
/* ***LIBRARY   SLATEC (FFTPACK) */
/* ***CATEGORY  J1A3 */
/* ***TYPE      SINGLE PRECISION (SINQI-S) */
/* ***KEYWORDS  FFTPACK, FOURIER TRANSFORM */
/* ***AUTHOR  Swarztrauber, P. N., (NCAR) */
/* ***DESCRIPTION */

/*  Subroutine SINQI initializes the array WSAVE which is used in */
/*  both SINQF and SINQB.  The prime factorization of N together with */
/*  a tabulation of the trigonometric functions are computed and */
/*  stored in WSAVE. */

/*  Input Parameter */

/*  N       the length of the sequence to be transformed.  The method */
/*          is most efficient when N is a product of small primes. */

/*  Output Parameter */

/*  WSAVE   a work array which must be dimensioned at least 3*N+15. */
/*          The same work array can be used for both SINQF and SINQB */
/*          as long as N remains unchanged.  Different WSAVE arrays */
/*          are required for different values of N.  The contents of */
/*          WSAVE must not be changed between calls of SINQF or SINQB. */

/* ***REFERENCES  P. N. Swarztrauber, Vectorizing the FFTs, in Parallel */
/*                 Computations (G. Rodrigue, ed.), Academic Press, */
/*                 1982, pp. 51-83. */
/* ***ROUTINES CALLED  COSQI */
/* ***REVISION HISTORY  (YYMMDD) */
/*   790601  DATE WRITTEN */
/*   830401  Modified to use SLATEC library source file format. */
/*   860115  Modified by Ron Boisvert to adhere to Fortran 77 by */
/*           changing dummy array size declarations (1) to (*) */
/*   861211  REVISION DATE from Version 3.2 */
/*   881128  Modified by Dick Valent to meet prologue standards. */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   920501  Reformatted the REFERENCES section.  (WRB) */
/* ***END PROLOGUE  SINQI */
/* ***FIRST EXECUTABLE STATEMENT  SINQI */
    /* Parameter adjustments */
    --wsave;

    /* Function Body */
    cosqi_(n, &wsave[1]);
    return 0;
} /* sinqi_ */

