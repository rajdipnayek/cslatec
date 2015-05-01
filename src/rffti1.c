/* rffti1.f -- translated by f2c (version 12.02.01).
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

/* DECK RFFTI1 */
/* Subroutine */ int rffti1_(integer *n, real *wa, integer *ifac)
{
    /* Initialized data */

    static integer ntryh[4] = { 4,2,3,5 };

    /* System generated locals */
    integer i__1, i__2, i__3;

    /* Local variables */
    static integer i__, j, k1, l1, l2, ib;
    static real fi;
    static integer ld, ii, nf, ip, nl, is, nq, nr;
    static real arg;
    static integer ido, ipm;
    static real tpi;
    static integer nfm1;
    static real argh;
    static integer ntry;
    static real argld;

/* ***BEGIN PROLOGUE  RFFTI1 */
/* ***PURPOSE  Initialize a real and an integer work array for RFFTF1 and */
/*            RFFTB1. */
/* ***LIBRARY   SLATEC (FFTPACK) */
/* ***CATEGORY  J1A1 */
/* ***TYPE      SINGLE PRECISION (RFFTI1-S, CFFTI1-C) */
/* ***KEYWORDS  FFTPACK, FOURIER TRANSFORM */
/* ***AUTHOR  Swarztrauber, P. N., (NCAR) */
/* ***DESCRIPTION */

/*   Subroutine RFFTI1 initializes the work arrays WA and IFAC which are */
/*   used in both RFFTF1 and RFFTB1.  The prime factorization of N and a */
/*   tabulation of the trigonometric functions are computed and stored in */
/*   IFAC and WA, respectively. */

/*   Input Argument */

/*   N       the length of the sequence to be transformed. */

/*   Output Arguments */

/*   WA      a real work array which must be dimensioned at least N. */

/*   IFAC    an integer work array which must be dimensioned at least 15. */

/*   The same work arrays can be used for both RFFTF1 and RFFTB1 as long */
/*   as N remains unchanged.  Different WA and IFAC arrays are required */
/*   for different values of N.  The contents of WA and IFAC must not be */
/*   changed between calls of RFFTF1 or RFFTB1. */

/* ***REFERENCES  P. N. Swarztrauber, Vectorizing the FFTs, in Parallel */
/*                 Computations (G. Rodrigue, ed.), Academic Press, */
/*                 1982, pp. 51-83. */
/* ***ROUTINES CALLED  (NONE) */
/* ***REVISION HISTORY  (YYMMDD) */
/*   790601  DATE WRITTEN */
/*   830401  Modified to use SLATEC library source file format. */
/*   860115  Modified by Ron Boisvert to adhere to Fortran 77 by */
/*           (a) changing dummy array size declarations (1) to (*), */
/*           (b) changing references to intrinsic function FLOAT */
/*               to REAL, and */
/*           (c) changing definition of variable TPI by using */
/*               FORTRAN intrinsic functions instead of DATA */
/*               statements. */
/*   881128  Modified by Dick Valent to meet prologue standards. */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900131  Routine changed from subsidiary to user-callable.  (WRB) */
/*   920501  Reformatted the REFERENCES section.  (WRB) */
/* ***END PROLOGUE  RFFTI1 */
    /* Parameter adjustments */
    --ifac;
    --wa;

    /* Function Body */
/* ***FIRST EXECUTABLE STATEMENT  RFFTI1 */
    nl = *n;
    nf = 0;
    j = 0;
L101:
    ++j;
    if (j - 4 <= 0) {
	goto L102;
    } else {
	goto L103;
    }
L102:
    ntry = ntryh[j - 1];
    goto L104;
L103:
    ntry += 2;
L104:
    nq = nl / ntry;
    nr = nl - ntry * nq;
    if (nr != 0) {
	goto L101;
    } else {
	goto L105;
    }
L105:
    ++nf;
    ifac[nf + 2] = ntry;
    nl = nq;
    if (ntry != 2) {
	goto L107;
    }
    if (nf == 1) {
	goto L107;
    }
    i__1 = nf;
    for (i__ = 2; i__ <= i__1; ++i__) {
	ib = nf - i__ + 2;
	ifac[ib + 2] = ifac[ib + 1];
/* L106: */
    }
    ifac[3] = 2;
L107:
    if (nl != 1) {
	goto L104;
    }
    ifac[1] = *n;
    ifac[2] = nf;
    tpi = atan(1.f) * 8.f;
    argh = tpi / *n;
    is = 0;
    nfm1 = nf - 1;
    l1 = 1;
    if (nfm1 == 0) {
	return 0;
    }
    i__1 = nfm1;
    for (k1 = 1; k1 <= i__1; ++k1) {
	ip = ifac[k1 + 2];
	ld = 0;
	l2 = l1 * ip;
	ido = *n / l2;
	ipm = ip - 1;
	i__2 = ipm;
	for (j = 1; j <= i__2; ++j) {
	    ld += l1;
	    i__ = is;
	    argld = ld * argh;
	    fi = 0.f;
	    i__3 = ido;
	    for (ii = 3; ii <= i__3; ii += 2) {
		i__ += 2;
		fi += 1.f;
		arg = fi * argld;
		wa[i__ - 1] = cos(arg);
		wa[i__] = sin(arg);
/* L108: */
	    }
	    is += ido;
/* L109: */
	}
	l1 = l2;
/* L110: */
    }
    return 0;
} /* rffti1_ */

