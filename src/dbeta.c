/* dbeta.f -- translated by f2c (version 12.02.01).
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

static integer c__1 = 1;
static integer c__2 = 2;

/* DECK DBETA */
doublereal dbeta_(doublereal *a, doublereal *b)
{
    /* Initialized data */

    static logical first = TRUE_;

    /* System generated locals */
    doublereal ret_val, d__1;

    /* Local variables */
    static doublereal xmin, xmax;
    extern doublereal d1mach_(integer *), dgamma_(doublereal *), dlbeta_(
	    doublereal *, doublereal *);
    extern /* Subroutine */ int dgamlm_(doublereal *, doublereal *);
    static doublereal alnsml;
    extern /* Subroutine */ int xermsg_(char *, char *, char *, integer *, 
	    integer *, ftnlen, ftnlen, ftnlen);

/* ***BEGIN PROLOGUE  DBETA */
/* ***PURPOSE  Compute the complete Beta function. */
/* ***LIBRARY   SLATEC (FNLIB) */
/* ***CATEGORY  C7B */
/* ***TYPE      DOUBLE PRECISION (BETA-S, DBETA-D, CBETA-C) */
/* ***KEYWORDS  COMPLETE BETA FUNCTION, FNLIB, SPECIAL FUNCTIONS */
/* ***AUTHOR  Fullerton, W., (LANL) */
/* ***DESCRIPTION */

/* DBETA(A,B) calculates the double precision complete beta function */
/* for double precision arguments A and B. */

/* ***REFERENCES  (NONE) */
/* ***ROUTINES CALLED  D1MACH, DGAMLM, DGAMMA, DLBETA, XERMSG */
/* ***REVISION HISTORY  (YYMMDD) */
/*   770601  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890531  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ) */
/*   900727  Added EXTERNAL statement.  (WRB) */
/* ***END PROLOGUE  DBETA */
/* ***FIRST EXECUTABLE STATEMENT  DBETA */
    if (first) {
	dgamlm_(&xmin, &xmax);
	alnsml = log(d1mach_(&c__1));
    }
    first = FALSE_;

    if (*a <= 0. || *b <= 0.) {
	xermsg_("SLATEC", "DBETA", "BOTH ARGUMENTS MUST BE GT 0", &c__2, &
		c__2, (ftnlen)6, (ftnlen)5, (ftnlen)27);
    }

    if (*a + *b < xmax) {
	d__1 = *a + *b;
	ret_val = dgamma_(a) * dgamma_(b) / dgamma_(&d__1);
    }
    if (*a + *b < xmax) {
	return ret_val;
    }

    ret_val = dlbeta_(a, b);
    if (ret_val < alnsml) {
	goto L20;
    }
    ret_val = exp(ret_val);
    return ret_val;

L20:
    ret_val = 0.;
    xermsg_("SLATEC", "DBETA", "A AND/OR B SO BIG BETA UNDERFLOWS", &c__1, &
	    c__1, (ftnlen)6, (ftnlen)5, (ftnlen)33);
    return ret_val;

} /* dbeta_ */

