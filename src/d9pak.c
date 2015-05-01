/* d9pak.f -- translated by f2c (version 12.02.01).
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
static integer c__15 = 15;
static integer c__16 = 16;
static integer c__1 = 1;
static integer c__2 = 2;

/* DECK D9PAK */
doublereal d9pak_(doublereal *y, integer *n)
{
    /* Initialized data */

    static doublereal a1n210 = 3.321928094887362347870319429489;
    static logical first = TRUE_;

    /* System generated locals */
    doublereal ret_val;

    /* Local variables */
    static integer ny;
    static doublereal a1n2b;
    static integer nmin, nmax, nsum;
    extern doublereal d1mach_(integer *);
    extern integer i1mach_(integer *);
    extern /* Subroutine */ int d9upak_(doublereal *, doublereal *, integer *)
	    , xermsg_(char *, char *, char *, integer *, integer *, ftnlen, 
	    ftnlen, ftnlen);

/* ***BEGIN PROLOGUE  D9PAK */
/* ***PURPOSE  Pack a base 2 exponent into a floating point number. */
/* ***LIBRARY   SLATEC (FNLIB) */
/* ***CATEGORY  A6B */
/* ***TYPE      DOUBLE PRECISION (R9PAK-S, D9PAK-D) */
/* ***KEYWORDS  FNLIB, PACK */
/* ***AUTHOR  Fullerton, W., (LANL) */
/* ***DESCRIPTION */

/* Pack a base 2 exponent into floating point number X.  This routine is */
/* almost the inverse of D9UPAK.  It is not exactly the inverse, because */
/* ABS(X) need not be between 0.5 and 1.0.  If both D9PAK and 2.d0**N */
/* were known to be in range we could compute */
/*               D9PAK = X *2.0d0**N */

/* ***REFERENCES  (NONE) */
/* ***ROUTINES CALLED  D1MACH, D9UPAK, I1MACH, XERMSG */
/* ***REVISION HISTORY  (YYMMDD) */
/*   790801  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890911  Removed unnecessary intrinsics.  (WRB) */
/*   891009  Corrected error when XERROR called.  (WRB) */
/*   891009  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ) */
/*   901009  Routine used I1MACH(7) where it should use I1MACH(10), */
/*           Corrected (RWC) */
/* ***END PROLOGUE  D9PAK */
/* ***FIRST EXECUTABLE STATEMENT  D9PAK */
    if (first) {
	a1n2b = 1.;
	if (i1mach_(&c__10) != 2) {
	    a1n2b = d1mach_(&c__5) * a1n210;
	}
	nmin = (integer) (a1n2b * i1mach_(&c__15));
	nmax = (integer) (a1n2b * i1mach_(&c__16));
    }
    first = FALSE_;

    d9upak_(y, &ret_val, &ny);

    nsum = *n + ny;
    if (nsum < nmin) {
	goto L40;
    }
    if (nsum > nmax) {
	xermsg_("SLATEC", "D9PAK", "PACKED NUMBER OVERFLOWS", &c__1, &c__2, (
		ftnlen)6, (ftnlen)5, (ftnlen)23);
    }

    if (nsum == 0) {
	return ret_val;
    }
    if (nsum > 0) {
	goto L30;
    }

L20:
    ret_val *= .5;
    ++nsum;
    if (nsum != 0) {
	goto L20;
    }
    return ret_val;

L30:
    ret_val *= 2.;
    --nsum;
    if (nsum != 0) {
	goto L30;
    }
    return ret_val;

L40:
    xermsg_("SLATEC", "D9PAK", "PACKED NUMBER UNDERFLOWS", &c__1, &c__1, (
	    ftnlen)6, (ftnlen)5, (ftnlen)24);
    ret_val = 0.;
    return ret_val;

} /* d9pak_ */

