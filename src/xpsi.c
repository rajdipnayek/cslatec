/* xpsi.f -- translated by f2c (version 12.02.01).
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

/* DECK XPSI */
doublereal xpsi_(real *a, integer *ipsik, integer *ipsix)
{
    /* Initialized data */

    static real cnum[12] = { 1.f,-1.f,1.f,-1.f,1.f,-691.f,1.f,-3617.f,43867.f,
	    -174611.f,77683.f,-236364091.f };
    static real cdenom[12] = { 12.f,120.f,252.f,240.f,132.f,32760.f,12.f,
	    8160.f,14364.f,6600.f,276.f,65520.f };

    /* System generated locals */
    integer i__1, i__2;
    real ret_val, r__1;

    /* Local variables */
    static real b, c__;
    static integer i__, k, m, n, k1;

/* ***BEGIN PROLOGUE  XPSI */
/* ***SUBSIDIARY */
/* ***PURPOSE  To compute values of the Psi function for XLEGF. */
/* ***LIBRARY   SLATEC */
/* ***CATEGORY  C7C */
/* ***TYPE      SINGLE PRECISION (XPSI-S, DXPSI-D) */
/* ***KEYWORDS  PSI FUNCTION */
/* ***AUTHOR  Smith, John M., (NBS and George Mason University) */
/* ***ROUTINES CALLED  (NONE) */
/* ***REVISION HISTORY  (YYMMDD) */
/*   820728  DATE WRITTEN */
/*   890126  Revised to meet SLATEC CML recommendations.  (DWL and JMS) */
/*   901019  Revisions to prologue.  (DWL and WRB) */
/*   901106  Changed all specific intrinsics to generic.  (WRB) */
/*           Corrected order of sections in prologue and added TYPE */
/*           section.  (WRB) */
/*   920127  Revised PURPOSE section of prologue.  (DWL) */
/* ***END PROLOGUE  XPSI */

/*        CNUM(I) AND CDENOM(I) ARE THE ( REDUCED ) NUMERATOR */
/*        AND 2*I*DENOMINATOR RESPECTIVELY OF THE 2*I TH BERNOULLI */
/*        NUMBER. */

/* ***FIRST EXECUTABLE STATEMENT  XPSI */
/* Computing MAX */
    i__1 = 0, i__2 = *ipsix - (integer) (*a);
    n = max(i__1,i__2);
    b = n + *a;
    k1 = *ipsik - 1;

/*        SERIES EXPANSION FOR A .GT. IPSIX USING IPSIK-1 TERMS. */

    c__ = 0.f;
    i__1 = k1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	k = *ipsik - i__;
/* L12: */
/* Computing 2nd power */
	r__1 = b;
	c__ = (c__ + cnum[k - 1] / cdenom[k - 1]) / (r__1 * r__1);
    }
    ret_val = log(b) - (c__ + .5f / b);
    if (n == 0) {
	goto L20;
    }
    b = 0.f;

/*        RECURRENCE FOR A .LE. IPSIX. */

    i__1 = n;
    for (m = 1; m <= i__1; ++m) {
/* L15: */
	b += 1.f / (n - m + *a);
    }
    ret_val -= b;
L20:
    return ret_val;
} /* xpsi_ */

