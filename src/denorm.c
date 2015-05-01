/* denorm.f -- translated by f2c (version 12.02.01).
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

/* DECK DENORM */
doublereal denorm_(integer *n, doublereal *x)
{
    /* Initialized data */

    static doublereal one = 1.;
    static doublereal zero = 0.;
    static doublereal rdwarf = 3.834e-20;
    static doublereal rgiant = 1.304e19;

    /* System generated locals */
    integer i__1;
    doublereal ret_val, d__1;

    /* Local variables */
    static integer i__;
    static doublereal s1, s2, s3, xabs, x1max, x3max, agiant, floatn;

/* ***BEGIN PROLOGUE  DENORM */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to DNSQ and DNSQE */
/* ***LIBRARY   SLATEC */
/* ***TYPE      DOUBLE PRECISION (ENORM-S, DENORM-D) */
/* ***AUTHOR  (UNKNOWN) */
/* ***DESCRIPTION */

/*     Given an N-vector X, this function calculates the */
/*     Euclidean norm of X. */

/*     The Euclidean norm is computed by accumulating the sum of */
/*     squares in three different sums. The sums of squares for the */
/*     small and large components are scaled so that no overflows */
/*     occur. Non-destructive underflows are permitted. Underflows */
/*     and overflows do not occur in the computation of the unscaled */
/*     sum of squares for the intermediate components. */
/*     The definitions of small, intermediate and large components */
/*     depend on two constants, RDWARF and RGIANT. The main */
/*     restrictions on these constants are that RDWARF**2 not */
/*     underflow and RGIANT**2 not overflow. The constants */
/*     given here are suitable for every known computer. */

/*     The function statement is */

/*       DOUBLE PRECISION FUNCTION DENORM(N,X) */

/*     where */

/*       N is a positive integer input variable. */

/*       X is an input array of length N. */

/* ***SEE ALSO  DNSQ, DNSQE */
/* ***ROUTINES CALLED  (NONE) */
/* ***REVISION HISTORY  (YYMMDD) */
/*   800301  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890831  Modified array declarations.  (WRB) */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900326  Removed duplicate information from DESCRIPTION section. */
/*           (WRB) */
/*   900328  Added TYPE section.  (WRB) */
/* ***END PROLOGUE  DENORM */
    /* Parameter adjustments */
    --x;

    /* Function Body */
/* ***FIRST EXECUTABLE STATEMENT  DENORM */
    s1 = zero;
    s2 = zero;
    s3 = zero;
    x1max = zero;
    x3max = zero;
    floatn = (doublereal) (*n);
    agiant = rgiant / floatn;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	xabs = (d__1 = x[i__], abs(d__1));
	if (xabs > rdwarf && xabs < agiant) {
	    goto L70;
	}
	if (xabs <= rdwarf) {
	    goto L30;
	}

/*              SUM FOR LARGE COMPONENTS. */

	if (xabs <= x1max) {
	    goto L10;
	}
/* Computing 2nd power */
	d__1 = x1max / xabs;
	s1 = one + s1 * (d__1 * d__1);
	x1max = xabs;
	goto L20;
L10:
/* Computing 2nd power */
	d__1 = xabs / x1max;
	s1 += d__1 * d__1;
L20:
	goto L60;
L30:

/*              SUM FOR SMALL COMPONENTS. */

	if (xabs <= x3max) {
	    goto L40;
	}
/* Computing 2nd power */
	d__1 = x3max / xabs;
	s3 = one + s3 * (d__1 * d__1);
	x3max = xabs;
	goto L50;
L40:
	if (xabs != zero) {
/* Computing 2nd power */
	    d__1 = xabs / x3max;
	    s3 += d__1 * d__1;
	}
L50:
L60:
	goto L80;
L70:

/*           SUM FOR INTERMEDIATE COMPONENTS. */

/* Computing 2nd power */
	d__1 = xabs;
	s2 += d__1 * d__1;
L80:
/* L90: */
	;
    }

/*     CALCULATION OF NORM. */

    if (s1 == zero) {
	goto L100;
    }
    ret_val = x1max * sqrt(s1 + s2 / x1max / x1max);
    goto L130;
L100:
    if (s2 == zero) {
	goto L110;
    }
    if (s2 >= x3max) {
	ret_val = sqrt(s2 * (one + x3max / s2 * (x3max * s3)));
    }
    if (s2 < x3max) {
	ret_val = sqrt(x3max * (s2 / x3max + x3max * s3));
    }
    goto L120;
L110:
    ret_val = x3max * sqrt(s3);
L120:
L130:
    return ret_val;

/*     LAST CARD OF FUNCTION DENORM. */

} /* denorm_ */

