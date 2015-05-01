/* sint.f -- translated by f2c (version 12.02.01).
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

/* DECK SINT */
/* Subroutine */ int sint_(integer *n, real *x, real *wsave)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__, k;
    static real t1, t2;
    static integer kc, nf;
    static real xh;
    static integer kw, np1, ns2, modn;
    static real sqrt3;
    extern /* Subroutine */ int rfftf_(integer *, real *, real *);

/* ***BEGIN PROLOGUE  SINT */
/* ***PURPOSE  Compute the sine transform of a real, odd sequence. */
/* ***LIBRARY   SLATEC (FFTPACK) */
/* ***CATEGORY  J1A3 */
/* ***TYPE      SINGLE PRECISION (SINT-S) */
/* ***KEYWORDS  FFTPACK, FOURIER TRANSFORM */
/* ***AUTHOR  Swarztrauber, P. N., (NCAR) */
/* ***DESCRIPTION */

/*  Subroutine SINT computes the discrete Fourier sine transform */
/*  of an odd sequence X(I).  The transform is defined below at */
/*  output parameter X. */

/*  SINT is the unnormalized inverse of itself since a call of SINT */
/*  followed by another call of SINT will multiply the input sequence */
/*  X by 2*(N+1). */

/*  The array WSAVE which is used by subroutine SINT must be */
/*  initialized by calling subroutine SINTI(N,WSAVE). */

/*  Input Parameters */

/*  N       the length of the sequence to be transformed.  The method */
/*          is most efficient when N+1 is the product of small primes. */

/*  X       an array which contains the sequence to be transformed */


/*  WSAVE   a work array with dimension at least INT(3.5*N+16) */
/*          in the program that calls SINT.  The WSAVE array must be */
/*          initialized by calling subroutine SINTI(N,WSAVE), and a */
/*          different WSAVE array must be used for each different */
/*          value of N.  This initialization does not have to be */
/*          repeated so long as N remains unchanged.  Thus subsequent */
/*          transforms can be obtained faster than the first. */

/*  Output Parameters */

/*  X       For I=1,...,N */

/*               X(I)= the sum from K=1 to K=N */

/*                    2*X(K)*SIN(K*I*PI/(N+1)) */

/*               A call of SINT followed by another call of */
/*               SINT will multiply the sequence X by 2*(N+1). */
/*               Hence SINT is the unnormalized inverse */
/*               of itself. */

/*  WSAVE   contains initialization calculations which must not be */
/*          destroyed between calls of SINT. */

/* ***REFERENCES  P. N. Swarztrauber, Vectorizing the FFTs, in Parallel */
/*                 Computations (G. Rodrigue, ed.), Academic Press, */
/*                 1982, pp. 51-83. */
/* ***ROUTINES CALLED  RFFTF */
/* ***REVISION HISTORY  (YYMMDD) */
/*   790601  DATE WRITTEN */
/*   830401  Modified to use SLATEC library source file format. */
/*   860115  Modified by Ron Boisvert to adhere to Fortran 77 by */
/*           (a) changing dummy array size declarations (1) to (*), */
/*           (b) changing definition of variable SQRT3 by using */
/*               FORTRAN intrinsic function SQRT instead of a DATA */
/*               statement. */
/*   881128  Modified by Dick Valent to meet prologue standards. */
/*   891009  Removed unreferenced statement label.  (WRB) */
/*   891009  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   920501  Reformatted the REFERENCES section.  (WRB) */
/* ***END PROLOGUE  SINT */
/* ***FIRST EXECUTABLE STATEMENT  SINT */
    /* Parameter adjustments */
    --wsave;
    --x;

    /* Function Body */
    sqrt3 = sqrt(3.f);
    if ((i__1 = *n - 2) < 0) {
	goto L101;
    } else if (i__1 == 0) {
	goto L102;
    } else {
	goto L103;
    }
L101:
    x[1] += x[1];
    return 0;
L102:
    xh = sqrt3 * (x[1] + x[2]);
    x[2] = sqrt3 * (x[1] - x[2]);
    x[1] = xh;
    return 0;
L103:
    np1 = *n + 1;
    ns2 = *n / 2;
    wsave[1] = 0.f;
    kw = np1;
    i__1 = ns2;
    for (k = 1; k <= i__1; ++k) {
	++kw;
	kc = np1 - k;
	t1 = x[k] - x[kc];
	t2 = wsave[kw] * (x[k] + x[kc]);
	wsave[k + 1] = t1 + t2;
	wsave[kc + 1] = t2 - t1;
/* L104: */
    }
    modn = *n % 2;
    if (modn != 0) {
	wsave[ns2 + 2] = x[ns2 + 1] * 4.f;
    }
    nf = np1 + ns2 + 1;
    rfftf_(&np1, &wsave[1], &wsave[nf]);
    x[1] = wsave[1] * .5f;
    i__1 = *n;
    for (i__ = 3; i__ <= i__1; i__ += 2) {
	x[i__ - 1] = -wsave[i__];
	x[i__] = x[i__ - 2] + wsave[i__ - 1];
/* L105: */
    }
    if (modn != 0) {
	return 0;
    }
    x[*n] = -wsave[*n + 1];
    return 0;
} /* sint_ */

