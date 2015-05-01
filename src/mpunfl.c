/* mpunfl.f -- translated by f2c (version 12.02.01).
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
static integer c__4 = 4;

/* DECK MPUNFL */
/* Subroutine */ int mpunfl_(integer *x)
{
    extern /* Subroutine */ int mpchk_(integer *, integer *);

/* ***BEGIN PROLOGUE  MPUNFL */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to DQDOTA and DQDOTI */
/* ***LIBRARY   SLATEC */
/* ***TYPE      ALL (MPUNFL-A) */
/* ***AUTHOR  (UNKNOWN) */
/* ***DESCRIPTION */

/* Called on multiple-precision underflow, i.e.  when the */
/* exponent of 'mp' number X would be less than -M. */

/* ***SEE ALSO  DQDOTA, DQDOTI */
/* ***ROUTINES CALLED  MPCHK */
/* ***REVISION HISTORY  (YYMMDD) */
/*   791001  DATE WRITTEN */
/*   890831  Modified array declarations.  (WRB) */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900402  Added TYPE section.  (WRB) */
/* ***END PROLOGUE  MPUNFL */
/* ***FIRST EXECUTABLE STATEMENT  MPUNFL */
    /* Parameter adjustments */
    --x;

    /* Function Body */
    mpchk_(&c__1, &c__4);
/* THE UNDERFLOWING NUMBER IS SET TO ZERO */
/* AN ALTERNATIVE WOULD BE TO CALL MPMINR (X) AND RETURN, */
/* POSSIBLY UPDATING A COUNTER AND TERMINATING EXECUTION */
/* AFTER A PRESET NUMBER OF UNDERFLOWS.  ACTION COULD EASILY */
/* BE DETERMINED BY A FLAG IN LABELLED COMMON. */
    x[1] = 0;
    return 0;
} /* mpunfl_ */

