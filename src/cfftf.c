/* cfftf.f -- translated by f2c (version 12.02.01).
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

/* DECK CFFTF */
/* Subroutine */ int cfftf_(integer *n, complex *c__, real *wsave)
{
    static integer iw1, iw2;
    extern /* Subroutine */ int cfftf1_(integer *, complex *, real *, real *, 
	    real *);

/* ***BEGIN PROLOGUE  CFFTF */
/* ***SUBSIDIARY */
/* ***PURPOSE  Compute the forward transform of a complex, periodic */
/*            sequence. */
/* ***LIBRARY   SLATEC (FFTPACK) */
/* ***CATEGORY  J1A2 */
/* ***TYPE      COMPLEX (RFFTF-S, CFFTF-C) */
/* ***KEYWORDS  FFTPACK, FOURIER TRANSFORM */
/* ***AUTHOR  Swarztrauber, P. N., (NCAR) */
/* ***DESCRIPTION */

/*  ******************************************************************** */
/*  *   NOTICE   NOTICE   NOTICE   NOTICE   NOTICE   NOTICE   NOTICE   * */
/*  ******************************************************************** */
/*  *                                                                  * */
/*  *   This routine uses non-standard Fortran 77 constructs and will  * */
/*  *   be removed from the library at a future date.  You are         * */
/*  *   requested to use CFFTF1.                                       * */
/*  *                                                                  * */
/*  ******************************************************************** */

/*  Subroutine CFFTF computes the forward complex discrete Fourier */
/*  transform (the Fourier analysis).  Equivalently, CFFTF computes */
/*  the Fourier coefficients of a complex periodic sequence. */
/*  The transform is defined below at output parameter C. */

/*  The transform is not normalized.  To obtain a normalized transform */
/*  the output must be divided by N.  Otherwise a call of CFFTF */
/*  followed by a call of CFFTB will multiply the sequence by N. */

/*  The array WSAVE which is used by subroutine CFFTF must be */
/*  initialized by calling subroutine CFFTI(N,WSAVE). */

/*  Input Parameters */

/*  N       the length of the complex sequence C.  The method is */
/*          more efficient when N is the product of small primes. */

/*  C       a complex array of length N which contains the sequence */

/*  WSAVE   a real work array which must be dimensioned at least 4*N+15 */
/*          in the program that calls CFFTF.  The WSAVE array must be */
/*          initialized by calling subroutine CFFTI(N,WSAVE), and a */
/*          different WSAVE array must be used for each different */
/*          value of N.  This initialization does not have to be */
/*          repeated so long as N remains unchanged.  Thus subsequent */
/*          transforms can be obtained faster than the first. */
/*          The same WSAVE array can be used by CFFTF and CFFTB. */

/*  Output Parameters */

/*  C       For J=1,...,N */

/*              C(J)=the sum from K=1,...,N of */

/*                 C(K)*EXP(-I*(J-1)*(K-1)*2*PI/N) */

/*                         where I=SQRT(-1) */

/*  WSAVE   contains initialization calculations which must not be */
/*          destroyed between calls of subroutine CFFTF or CFFTB */

/* ***REFERENCES  P. N. Swarztrauber, Vectorizing the FFTs, in Parallel */
/*                 Computations (G. Rodrigue, ed.), Academic Press, */
/*                 1982, pp. 51-83. */
/* ***ROUTINES CALLED  CFFTF1 */
/* ***REVISION HISTORY  (YYMMDD) */
/*   790601  DATE WRITTEN */
/*   830401  Modified to use SLATEC library source file format. */
/*   860115  Modified by Ron Boisvert to adhere to Fortran 77 by */
/*           changing dummy array size declarations (1) to (*). */
/*   861211  REVISION DATE from Version 3.2 */
/*   881128  Modified by Dick Valent to meet prologue standards. */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900131  Routine changed from user-callable to subsidiary */
/*           because of non-standard Fortran 77 arguments in the */
/*           call to CFFTB1.  (WRB) */
/*   920501  Reformatted the REFERENCES section.  (WRB) */
/* ***END PROLOGUE  CFFTF */
/* ***FIRST EXECUTABLE STATEMENT  CFFTF */
    /* Parameter adjustments */
    --wsave;
    --c__;

    /* Function Body */
    if (*n == 1) {
	return 0;
    }
    iw1 = *n + *n + 1;
    iw2 = iw1 + *n + *n;
    cfftf1_(n, &c__[1], &wsave[1], &wsave[iw1], &wsave[iw2]);
    return 0;
} /* cfftf_ */

