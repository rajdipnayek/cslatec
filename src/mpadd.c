/* mpadd.f -- translated by f2c (version 12.02.01).
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

static integer c__0 = 0;

/* DECK MPADD */
/* Subroutine */ int mpadd_(integer *x, integer *y, integer *z__)
{
    extern /* Subroutine */ int mpadd2_(integer *, integer *, integer *, 
	    integer *, integer *);

/* ***BEGIN PROLOGUE  MPADD */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to DQDOTA and DQDOTI */
/* ***LIBRARY   SLATEC */
/* ***TYPE      ALL (MPADD-A) */
/* ***AUTHOR  (UNKNOWN) */
/* ***DESCRIPTION */

/* Adds X and Y, forming result in Z, where X, Y and Z are 'mp' */
/*  (multiple precision) numbers.  Four guard digits are used, */
/*  and then R*-rounding. */

/* ***SEE ALSO  DQDOTA, DQDOTI */
/* ***ROUTINES CALLED  MPADD2 */
/* ***REVISION HISTORY  (YYMMDD) */
/*   791001  DATE WRITTEN */
/*   890831  Modified array declarations.  (WRB) */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900402  Added TYPE section.  (WRB) */
/* ***END PROLOGUE  MPADD */
/* ***FIRST EXECUTABLE STATEMENT  MPADD */
    /* Parameter adjustments */
    --z__;
    --y;
    --x;

    /* Function Body */
    mpadd2_(&x[1], &y[1], &z__[1], &y[1], &c__0);
    return 0;
} /* mpadd_ */

