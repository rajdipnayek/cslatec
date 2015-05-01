/* mpmlp.f -- translated by f2c (version 12.02.01).
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

/* DECK MPMLP */
/* Subroutine */ int mpmlp_(integer *u, integer *v, integer *w, integer *j)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__;

/* ***BEGIN PROLOGUE  MPMLP */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to DQDOTA and DQDOTI */
/* ***LIBRARY   SLATEC */
/* ***TYPE      ALL (MPMLP-A) */
/* ***AUTHOR  (UNKNOWN) */
/* ***DESCRIPTION */

/* Performs inner multiplication loop for MPMUL. Carries are not pro- */
/* pagated in inner loop, which saves time at the expense of space. */

/* ***SEE ALSO  DQDOTA, DQDOTI */
/* ***ROUTINES CALLED  (NONE) */
/* ***REVISION HISTORY  (YYMMDD) */
/*   791001  DATE WRITTEN */
/*   890831  Modified array declarations.  (WRB) */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900402  Added TYPE section.  (WRB) */
/* ***END PROLOGUE  MPMLP */
/* ***FIRST EXECUTABLE STATEMENT  MPMLP */
    /* Parameter adjustments */
    --v;
    --u;

    /* Function Body */
    i__1 = *j;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L10: */
	u[i__] += *w * v[i__];
    }
    return 0;
} /* mpmlp_ */

