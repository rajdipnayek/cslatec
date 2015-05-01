/* sinti.f -- translated by f2c (version 12.02.01).
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

/* DECK SINTI */
/* Subroutine */ int sinti_(integer *n, real *wsave)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer k, kf;
    static real fk, dt, pi;
    static integer ks, np1, ns2;
    extern /* Subroutine */ int rffti_(integer *, real *);

/* ***BEGIN PROLOGUE  SINTI */
/* ***PURPOSE  Initialize a work array for SINT. */
/* ***LIBRARY   SLATEC (FFTPACK) */
/* ***CATEGORY  J1A3 */
/* ***TYPE      SINGLE PRECISION (SINTI-S) */
/* ***KEYWORDS  FFTPACK, FOURIER TRANSFORM */
/* ***AUTHOR  Swarztrauber, P. N., (NCAR) */
/* ***DESCRIPTION */

/*  Subroutine SINTI initializes the array WSAVE which is used in */
/*  subroutine SINT.  The prime factorization of N together with */
/*  a tabulation of the trigonometric functions are computed and */
/*  stored in WSAVE. */

/*  Input Parameter */

/*  N       the length of the sequence to be transformed.  The method */
/*          is most efficient when N+1 is a product of small primes. */

/*  Output Parameter */

/*  WSAVE   a work array with at least INT(3.5*N+16) locations. */
/*          Different WSAVE arrays are required for different values */
/*          of N.  The contents of WSAVE must not be changed between */
/*          calls of SINT. */

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
/*           (c) changing definition of variable PI by using */
/*               FORTRAN intrinsic function ATAN instead of a DATA */
/*               statement. */
/*   881128  Modified by Dick Valent to meet prologue standards. */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890531  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   920501  Reformatted the REFERENCES section.  (WRB) */
/* ***END PROLOGUE  SINTI */
/* ***FIRST EXECUTABLE STATEMENT  SINTI */
    /* Parameter adjustments */
    --wsave;

    /* Function Body */
    if (*n <= 1) {
	return 0;
    }
    pi = atan(1.f) * 4.f;
    np1 = *n + 1;
    ns2 = *n / 2;
    dt = pi / np1;
    ks = *n + 2;
    kf = ks + ns2 - 1;
    fk = 0.f;
    i__1 = kf;
    for (k = ks; k <= i__1; ++k) {
	fk += 1.f;
	wsave[k] = sin(fk * dt) * 2.f;
/* L101: */
    }
    rffti_(&np1, &wsave[kf + 1]);
    return 0;
} /* sinti_ */

