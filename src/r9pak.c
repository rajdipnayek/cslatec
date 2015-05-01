/* r9pak.f -- translated by f2c (version 12.02.01).
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

/* Table of constant values */

static integer c__10 = 10;
static integer c__5 = 5;
static integer c__12 = 12;
static integer c__13 = 13;
static integer c__2 = 2;
static integer c__1 = 1;

/* DECK R9PAK */
doublereal r9pak_(real *y, integer *n)
{
    /* Initialized data */

    static real a1n210 = 3.321928094887362f;
    static logical first = TRUE_;

    /* System generated locals */
    real ret_val;

    /* Local variables */
    static integer ny;
    static real a1n2b;
    static integer nmin, nmax, nsum;
    extern integer i1mach_(integer *);
    extern doublereal r1mach_(integer *);
    extern /* Subroutine */ int r9upak_(real *, real *, integer *), xermsg_(
	    char *, char *, char *, integer *, integer *, ftnlen, ftnlen, 
	    ftnlen);

/* ***BEGIN PROLOGUE  R9PAK */
/* ***PURPOSE  Pack a base 2 exponent into a floating point number. */
/* ***LIBRARY   SLATEC (FNLIB) */
/* ***CATEGORY  A6B */
/* ***TYPE      SINGLE PRECISION (R9PAK-S, D9PAK-D) */
/* ***KEYWORDS  FNLIB, PACK */
/* ***AUTHOR  Fullerton, W., (LANL) */
/* ***DESCRIPTION */

/* Pack a base 2 exponent into floating point number Y.  This */
/* routine is almost the inverse of R9UPAK.  It is not exactly */
/* the inverse, because ABS(X) need not be between 0.5 and */
/* 1.0.  If both R9PAK and 2.0**N were known to be in range, we */
/* could compute */
/*       R9PAK = Y * 2.0**N . */

/* ***REFERENCES  (NONE) */
/* ***ROUTINES CALLED  I1MACH, R1MACH, R9UPAK, XERMSG */
/* ***REVISION HISTORY  (YYMMDD) */
/*   790801  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890531  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ) */
/*   901009  Routine used I1MACH(7) where it should use I1MACH(10), */
/*           Corrected (RWC) */
/* ***END PROLOGUE  R9PAK */
/* ***FIRST EXECUTABLE STATEMENT  R9PAK */
    if (first) {
	a1n2b = 1.f;
	if (i1mach_(&c__10) != 2) {
	    a1n2b = r1mach_(&c__5) * a1n210;
	}
	nmin = a1n2b * i1mach_(&c__12);
	nmax = a1n2b * i1mach_(&c__13);
    }
    first = FALSE_;

    r9upak_(y, &ret_val, &ny);

    nsum = *n + ny;
    if (nsum < nmin) {
	goto L40;
    }
    if (nsum > nmax) {
	xermsg_("SLATEC", "R9PAK", "PACKED NUMBER OVERFLOWS", &c__2, &c__2, (
		ftnlen)6, (ftnlen)5, (ftnlen)23);
    }

    if (nsum == 0) {
	return ret_val;
    }
    if (nsum > 0) {
	goto L30;
    }

L20:
    ret_val *= .5f;
    ++nsum;
    if (nsum != 0) {
	goto L20;
    }
    return ret_val;

L30:
    ret_val *= 2.f;
    --nsum;
    if (nsum != 0) {
	goto L30;
    }
    return ret_val;

L40:
    xermsg_("SLATEC", "R9PAK", "PACKED NUMBER UNDERFLOWS", &c__1, &c__1, (
	    ftnlen)6, (ftnlen)5, (ftnlen)24);
    ret_val = 0.f;
    return ret_val;

} /* r9pak_ */

