/* ali.f -- translated by f2c (version 12.02.01).
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

/* DECK ALI */
doublereal ali_(real *x)
{
    /* System generated locals */
    real ret_val, r__1;

    /* Local variables */
    extern doublereal ei_(real *);
    extern /* Subroutine */ int xermsg_(char *, char *, char *, integer *, 
	    integer *, ftnlen, ftnlen, ftnlen);

/* ***BEGIN PROLOGUE  ALI */
/* ***PURPOSE  Compute the logarithmic integral. */
/* ***LIBRARY   SLATEC (FNLIB) */
/* ***CATEGORY  C5 */
/* ***TYPE      SINGLE PRECISION (ALI-S, DLI-D) */
/* ***KEYWORDS  FNLIB, LOGARITHMIC INTEGRAL, SPECIAL FUNCTIONS */
/* ***AUTHOR  Fullerton, W., (LANL) */
/* ***DESCRIPTION */

/* ALI(X) computes the logarithmic integral; i.e., the */
/* integral from 0.0 to X of (1.0/ln(t))dt. */

/* ***REFERENCES  (NONE) */
/* ***ROUTINES CALLED  EI, XERMSG */
/* ***REVISION HISTORY  (YYMMDD) */
/*   770601  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890531  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ) */
/*   900326  Removed duplicate information from DESCRIPTION section. */
/*           (WRB) */
/* ***END PROLOGUE  ALI */
/* ***FIRST EXECUTABLE STATEMENT  ALI */
    if (*x <= 0.f) {
	xermsg_("SLATEC", "ALI", "LOG INTEGRAL UNDEFINED FOR X LE 0", &c__1, &
		c__2, (ftnlen)6, (ftnlen)3, (ftnlen)33);
    }
    if (*x == 1.f) {
	xermsg_("SLATEC", "ALI", "LOG INTEGRAL UNDEFINED FOR X = 1", &c__2, &
		c__2, (ftnlen)6, (ftnlen)3, (ftnlen)32);
    }

    r__1 = log(*x);
    ret_val = ei_(&r__1);

    return ret_val;
} /* ali_ */

