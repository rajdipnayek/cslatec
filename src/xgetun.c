/* xgetun.f -- translated by f2c (version 12.02.01).
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

static integer c__3 = 3;
static integer c__0 = 0;
static logical c_false = FALSE_;

/* DECK XGETUN */
/* Subroutine */ int xgetun_(integer *iunit)
{
    extern integer j4save_(integer *, integer *, logical *);

/* ***BEGIN PROLOGUE  XGETUN */
/* ***PURPOSE  Return the (first) output file to which error messages */
/*            are being sent. */
/* ***LIBRARY   SLATEC (XERROR) */
/* ***CATEGORY  R3C */
/* ***TYPE      ALL (XGETUN-A) */
/* ***KEYWORDS  ERROR, XERROR */
/* ***AUTHOR  Jones, R. E., (SNLA) */
/* ***DESCRIPTION */

/*     Abstract */
/*        XGETUN gets the (first) output file to which error messages */
/*        are being sent.  To find out if more than one file is being */
/*        used, one must use the XGETUA routine. */

/*     Description of Parameter */
/*      --Output-- */
/*        IUNIT - the logical unit number of the  (first) unit to */
/*                which error messages are being sent. */
/*                A value of zero means that the default file, as */
/*                defined by the I1MACH routine, is being used. */

/* ***REFERENCES  R. E. Jones and D. K. Kahaner, XERROR, the SLATEC */
/*                 Error-handling Package, SAND82-0800, Sandia */
/*                 Laboratories, 1982. */
/* ***ROUTINES CALLED  J4SAVE */
/* ***REVISION HISTORY  (YYMMDD) */
/*   790801  DATE WRITTEN */
/*   861211  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   920501  Reformatted the REFERENCES section.  (WRB) */
/* ***END PROLOGUE  XGETUN */
/* ***FIRST EXECUTABLE STATEMENT  XGETUN */
    *iunit = j4save_(&c__3, &c__0, &c_false);
    return 0;
} /* xgetun_ */

