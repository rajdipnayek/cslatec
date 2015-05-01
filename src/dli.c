/* dli.f -- translated by f2c (version 12.02.01).
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

/* DECK DLI */
doublereal dli_(doublereal *x)
{
    /* System generated locals */
    doublereal ret_val, d__1;

    /* Local variables */
    extern doublereal dei_(doublereal *);
    extern /* Subroutine */ int xermsg_(char *, char *, char *, integer *, 
	    integer *, ftnlen, ftnlen, ftnlen);

/* ***BEGIN PROLOGUE  DLI */
/* ***PURPOSE  Compute the logarithmic integral. */
/* ***LIBRARY   SLATEC (FNLIB) */
/* ***CATEGORY  C5 */
/* ***TYPE      DOUBLE PRECISION (ALI-S, DLI-D) */
/* ***KEYWORDS  FNLIB, LOGARITHMIC INTEGRAL, SPECIAL FUNCTIONS */
/* ***AUTHOR  Fullerton, W., (LANL) */
/* ***DESCRIPTION */

/* DLI(X) calculates the double precision logarithmic integral */
/* for double precision argument X. */

/* ***REFERENCES  (NONE) */
/* ***ROUTINES CALLED  DEI, XERMSG */
/* ***REVISION HISTORY  (YYMMDD) */
/*   770701  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890531  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ) */
/* ***END PROLOGUE  DLI */
/* ***FIRST EXECUTABLE STATEMENT  DLI */
    if (*x <= 0.) {
	xermsg_("SLATEC", "DLI", "LOG INTEGRAL UNDEFINED FOR X LE 0", &c__1, &
		c__2, (ftnlen)6, (ftnlen)3, (ftnlen)33);
    }
    if (*x == 1.) {
	xermsg_("SLATEC", "DLI", "LOG INTEGRAL UNDEFINED FOR X = 0", &c__2, &
		c__2, (ftnlen)6, (ftnlen)3, (ftnlen)32);
    }

    d__1 = log(*x);
    ret_val = dei_(&d__1);

    return ret_val;
} /* dli_ */

