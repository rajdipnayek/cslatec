/* cdcst.f -- translated by f2c (version 12.02.01).
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

/* DECK CDCST */
/* Subroutine */ int cdcst_(integer *maxord, integer *mint, integer *iswflg, 
	real *el, real *tq)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer i__, j;
    static real sum;
    static integer mxrd;
    static real gamma[14], factrl[12];

/* ***BEGIN PROLOGUE  CDCST */
/* ***SUBSIDIARY */
/* ***PURPOSE  CDCST sets coefficients used by the core integrator CDSTP. */
/* ***LIBRARY   SLATEC (SDRIVE) */
/* ***TYPE      COMPLEX (SDCST-S, DDCST-D, CDCST-C) */
/* ***AUTHOR  Kahaner, D. K., (NIST) */
/*             National Institute of Standards and Technology */
/*             Gaithersburg, MD  20899 */
/*           Sutherland, C. D., (LANL) */
/*             Mail Stop D466 */
/*             Los Alamos National Laboratory */
/*             Los Alamos, NM  87545 */
/* ***DESCRIPTION */

/*  CDCST is called by CDNTL.  The array EL determines the basic method. */
/*  The array TQ is involved in adjusting the step size in relation */
/*  to truncation error.  EL and TQ depend upon MINT, and are calculated */
/*  for orders 1 to MAXORD(.LE. 12).  For each order NQ, the coefficients */
/*  EL are calculated from the generating polynomial: */
/*    L(T) = EL(1,NQ) + EL(2,NQ)*T + ... + EL(NQ+1,NQ)*T**NQ. */
/*  For the implicit Adams methods, L(T) is given by */
/*    dL/dT = (1+T)*(2+T)* ... *(NQ-1+T)/K,   L(-1) = 0, */
/*    where      K = factorial(NQ-1). */
/*  For the Gear methods, */
/*    L(T) = (1+T)*(2+T)* ... *(NQ+T)/K, */
/*    where      K = factorial(NQ)*(1 + 1/2 + ... + 1/NQ). */
/*  For each order NQ, there are three components of TQ. */

/* ***ROUTINES CALLED  (NONE) */
/* ***REVISION HISTORY  (YYMMDD) */
/*   790601  DATE WRITTEN */
/*   900329  Initial submission to SLATEC. */
/* ***END PROLOGUE  CDCST */
/* ***FIRST EXECUTABLE STATEMENT  CDCST */
    /* Parameter adjustments */
    tq -= 4;
    el -= 14;

    /* Function Body */
    factrl[0] = 1.f;
    i__1 = *maxord;
    for (i__ = 2; i__ <= i__1; ++i__) {
/* L10: */
	factrl[i__ - 1] = i__ * factrl[i__ - 2];
    }
/*                                             Compute Adams coefficients */
    if (*mint == 1) {
	gamma[0] = 1.f;
	i__1 = *maxord + 1;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    sum = 0.f;
	    i__2 = i__;
	    for (j = 1; j <= i__2; ++j) {
/* L30: */
		sum -= gamma[j - 1] / (i__ - j + 2);
	    }
/* L40: */
	    gamma[i__] = sum;
	}
	el[14] = 1.f;
	el[15] = 1.f;
	el[28] = 1.f;
	el[29] = 1.f;
	i__1 = *maxord;
	for (j = 3; j <= i__1; ++j) {
	    el[j * 13 + 2] = factrl[j - 2];
	    i__2 = j;
	    for (i__ = 3; i__ <= i__2; ++i__) {
/* L50: */
		el[i__ + j * 13] = (j - 1) * el[i__ + (j - 1) * 13] + el[i__ 
			- 1 + (j - 1) * 13];
	    }
/* L60: */
	    el[j + 1 + j * 13] = 1.f;
	}
	i__1 = *maxord;
	for (j = 2; j <= i__1; ++j) {
	    el[j * 13 + 1] = el[(j - 1) * 13 + 1] + gamma[j - 1];
	    el[j * 13 + 2] = 1.f;
	    i__2 = j + 1;
	    for (i__ = 3; i__ <= i__2; ++i__) {
/* L80: */
		el[i__ + j * 13] /= (i__ - 1) * factrl[j - 2];
	    }
	}
	i__2 = *maxord;
	for (j = 1; j <= i__2; ++j) {
	    tq[j * 3 + 1] = -1.f / (factrl[j - 1] * gamma[j - 1]);
	    tq[j * 3 + 2] = -1.f / gamma[j];
/* L100: */
	    tq[j * 3 + 3] = -1.f / gamma[j + 1];
	}
/*                                              Compute Gear coefficients */
    } else if (*mint == 2) {
	el[14] = 1.f;
	el[15] = 1.f;
	i__2 = *maxord;
	for (j = 2; j <= i__2; ++j) {
	    el[j * 13 + 1] = factrl[j - 1];
	    i__1 = j;
	    for (i__ = 2; i__ <= i__1; ++i__) {
/* L120: */
		el[i__ + j * 13] = j * el[i__ + (j - 1) * 13] + el[i__ - 1 + (
			j - 1) * 13];
	    }
/* L130: */
	    el[j + 1 + j * 13] = 1.f;
	}
	sum = 1.f;
	i__2 = *maxord;
	for (j = 2; j <= i__2; ++j) {
	    sum += 1.f / j;
	    i__1 = j + 1;
	    for (i__ = 1; i__ <= i__1; ++i__) {
/* L150: */
		el[i__ + j * 13] /= factrl[j - 1] * sum;
	    }
	}
	i__1 = *maxord;
	for (j = 1; j <= i__1; ++j) {
	    if (j > 1) {
		tq[j * 3 + 1] = 1.f / factrl[j - 2];
	    }
	    tq[j * 3 + 2] = (j + 1) / el[j * 13 + 1];
/* L170: */
	    tq[j * 3 + 3] = (j + 2) / el[j * 13 + 1];
	}
    }
/*                          Compute constants used in the stiffness test. */
/*                          These are the ratio of TQ(2,NQ) for the Gear */
/*                          methods to those for the Adams methods. */
    if (*iswflg == 3) {
	mxrd = min(*maxord,5);
	if (*mint == 2) {
	    gamma[0] = 1.f;
	    i__1 = mxrd;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		sum = 0.f;
		i__2 = i__;
		for (j = 1; j <= i__2; ++j) {
/* L180: */
		    sum -= gamma[j - 1] / (i__ - j + 2);
		}
/* L190: */
		gamma[i__] = sum;
	    }
	}
	sum = 1.f;
	i__1 = mxrd;
	for (i__ = 2; i__ <= i__1; ++i__) {
	    sum += 1.f / i__;
/* L200: */
	    el[i__ + 14] = -(i__ + 1) * sum * gamma[i__];
	}
    }
    return 0;
} /* cdcst_ */

