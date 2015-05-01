/* mpmuli.f -- translated by f2c (version 12.02.01).
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

/* DECK MPMULI */
/* Subroutine */ int mpmuli_(integer *x, integer *iy, integer *z__)
{
    extern /* Subroutine */ int mpmul2_(integer *, integer *, integer *, 
	    integer *);

/* ***BEGIN PROLOGUE  MPMULI */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to DQDOTA and DQDOTI */
/* ***LIBRARY   SLATEC */
/* ***TYPE      ALL (MPMULI-A) */
/* ***AUTHOR  (UNKNOWN) */
/* ***DESCRIPTION */

/* Multiplies 'mp' X by single-precision integer IY giving 'mp' Z. */
/* This is faster than using MPMUL.  Result is ROUNDED. */
/* Multiplication by 1 may be used to normalize a number */
/* even if the last digit is B. */

/* ***SEE ALSO  DQDOTA, DQDOTI */
/* ***ROUTINES CALLED  MPMUL2 */
/* ***REVISION HISTORY  (YYMMDD) */
/*   791001  DATE WRITTEN */
/*   890831  Modified array declarations.  (WRB) */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900402  Added TYPE section.  (WRB) */
/* ***END PROLOGUE  MPMULI */
/* ***FIRST EXECUTABLE STATEMENT  MPMULI */
    /* Parameter adjustments */
    --z__;
    --x;

    /* Function Body */
    mpmul2_(&x[1], iy, &z__[1], &c__0);
    return 0;
} /* mpmuli_ */

