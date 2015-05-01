/* cfftf1.f -- translated by f2c (version 12.02.01).
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

/* DECK CFFTF1 */
/* Subroutine */ int cfftf1_(integer *n, real *c__, real *ch, real *wa, 
	integer *ifac)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__, k1, l1, l2, n2, na, nf, ip, iw, ix2, ix3, ix4, nac, 
	    ido, idl1, idot;
    extern /* Subroutine */ int passf_(integer *, integer *, integer *, 
	    integer *, integer *, real *, real *, real *, real *, real *, 
	    real *), passf2_(integer *, integer *, real *, real *, real *), 
	    passf3_(integer *, integer *, real *, real *, real *, real *), 
	    passf4_(integer *, integer *, real *, real *, real *, real *, 
	    real *), passf5_(integer *, integer *, real *, real *, real *, 
	    real *, real *, real *);

/* ***BEGIN PROLOGUE  CFFTF1 */
/* ***PURPOSE  Compute the forward transform of a complex, periodic */
/*            sequence. */
/* ***LIBRARY   SLATEC (FFTPACK) */
/* ***CATEGORY  J1A2 */
/* ***TYPE      COMPLEX (RFFTF1-S, CFFTF1-C) */
/* ***KEYWORDS  FFTPACK, FOURIER TRANSFORM */
/* ***AUTHOR  Swarztrauber, P. N., (NCAR) */
/* ***DESCRIPTION */

/*  Subroutine CFFTF1 computes the forward complex discrete Fourier */
/*  transform (the Fourier analysis).  Equivalently, CFFTF1 computes */
/*  the Fourier coefficients of a complex periodic sequence. */
/*  The transform is defined below at output parameter C. */

/*  The transform is not normalized.  To obtain a normalized transform */
/*  the output must be divided by N.  Otherwise a call of CFFTF1 */
/*  followed by a call of CFFTB1 will multiply the sequence by N. */

/*  The arrays WA and IFAC which are used by subroutine CFFTB1 must be */
/*  initialized by calling subroutine CFFTI1 (N, WA, IFAC). */

/*  Input Parameters */

/*  N       the length of the complex sequence C.  The method is */
/*          more efficient when N is the product of small primes. */

/*  C       a complex array of length N which contains the sequence */

/*  CH      a real work array of length at least 2*N */

/*  WA      a real work array which must be dimensioned at least 2*N. */

/*  IFAC    an integer work array which must be dimensioned at least 15. */

/*          The WA and IFAC arrays must be initialized by calling */
/*          subroutine CFFTI1 (N, WA, IFAC), and different WA and IFAC */
/*          arrays must be used for each different value of N.  This */
/*          initialization does not have to be repeated so long as N */
/*          remains unchanged.  Thus subsequent transforms can be */
/*          obtained faster than the first.  The same WA and IFAC arrays */
/*          can be used by CFFTF1 and CFFTB1. */

/*  Output Parameters */

/*  C       For J=1,...,N */

/*              C(J)=the sum from K=1,...,N of */

/*                 C(K)*EXP(-I*(J-1)*(K-1)*2*PI/N) */

/*                         where I=SQRT(-1) */

/*  NOTE:   WA and IFAC contain initialization calculations which must */
/*          not be destroyed between calls of subroutine CFFTF1 or CFFTB1 */

/* ***REFERENCES  P. N. Swarztrauber, Vectorizing the FFTs, in Parallel */
/*                 Computations (G. Rodrigue, ed.), Academic Press, */
/*                 1982, pp. 51-83. */
/* ***ROUTINES CALLED  PASSF, PASSF2, PASSF3, PASSF4, PASSF5 */
/* ***REVISION HISTORY  (YYMMDD) */
/*   790601  DATE WRITTEN */
/*   830401  Modified to use SLATEC library source file format. */
/*   860115  Modified by Ron Boisvert to adhere to Fortran 77 by */
/*           changing dummy array size declarations (1) to (*). */
/*   881128  Modified by Dick Valent to meet prologue standards. */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900131  Routine changed from subsidiary to user-callable.  (WRB) */
/*   920501  Reformatted the REFERENCES section.  (WRB) */
/* ***END PROLOGUE  CFFTF1 */
/* ***FIRST EXECUTABLE STATEMENT  CFFTF1 */
    /* Parameter adjustments */
    --ifac;
    --wa;
    --ch;
    --c__;

    /* Function Body */
    nf = ifac[2];
    na = 0;
    l1 = 1;
    iw = 1;
    i__1 = nf;
    for (k1 = 1; k1 <= i__1; ++k1) {
	ip = ifac[k1 + 2];
	l2 = ip * l1;
	ido = *n / l2;
	idot = ido + ido;
	idl1 = idot * l1;
	if (ip != 4) {
	    goto L103;
	}
	ix2 = iw + idot;
	ix3 = ix2 + idot;
	if (na != 0) {
	    goto L101;
	}
	passf4_(&idot, &l1, &c__[1], &ch[1], &wa[iw], &wa[ix2], &wa[ix3]);
	goto L102;
L101:
	passf4_(&idot, &l1, &ch[1], &c__[1], &wa[iw], &wa[ix2], &wa[ix3]);
L102:
	na = 1 - na;
	goto L115;
L103:
	if (ip != 2) {
	    goto L106;
	}
	if (na != 0) {
	    goto L104;
	}
	passf2_(&idot, &l1, &c__[1], &ch[1], &wa[iw]);
	goto L105;
L104:
	passf2_(&idot, &l1, &ch[1], &c__[1], &wa[iw]);
L105:
	na = 1 - na;
	goto L115;
L106:
	if (ip != 3) {
	    goto L109;
	}
	ix2 = iw + idot;
	if (na != 0) {
	    goto L107;
	}
	passf3_(&idot, &l1, &c__[1], &ch[1], &wa[iw], &wa[ix2]);
	goto L108;
L107:
	passf3_(&idot, &l1, &ch[1], &c__[1], &wa[iw], &wa[ix2]);
L108:
	na = 1 - na;
	goto L115;
L109:
	if (ip != 5) {
	    goto L112;
	}
	ix2 = iw + idot;
	ix3 = ix2 + idot;
	ix4 = ix3 + idot;
	if (na != 0) {
	    goto L110;
	}
	passf5_(&idot, &l1, &c__[1], &ch[1], &wa[iw], &wa[ix2], &wa[ix3], &wa[
		ix4]);
	goto L111;
L110:
	passf5_(&idot, &l1, &ch[1], &c__[1], &wa[iw], &wa[ix2], &wa[ix3], &wa[
		ix4]);
L111:
	na = 1 - na;
	goto L115;
L112:
	if (na != 0) {
	    goto L113;
	}
	passf_(&nac, &idot, &ip, &l1, &idl1, &c__[1], &c__[1], &c__[1], &ch[1]
		, &ch[1], &wa[iw]);
	goto L114;
L113:
	passf_(&nac, &idot, &ip, &l1, &idl1, &ch[1], &ch[1], &ch[1], &c__[1], 
		&c__[1], &wa[iw]);
L114:
	if (nac != 0) {
	    na = 1 - na;
	}
L115:
	l1 = l2;
	iw += (ip - 1) * idot;
/* L116: */
    }
    if (na == 0) {
	return 0;
    }
    n2 = *n + *n;
    i__1 = n2;
    for (i__ = 1; i__ <= i__1; ++i__) {
	c__[i__] = ch[i__];
/* L117: */
    }
    return 0;
} /* cfftf1_ */

