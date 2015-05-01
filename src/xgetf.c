/* xgetf.f -- translated by f2c (version 12.02.01).
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

static integer c__2 = 2;
static integer c__0 = 0;
static logical c_false = FALSE_;

/* DECK XGETF */
/* Subroutine */ int xgetf_(integer *kontrl)
{
    extern integer j4save_(integer *, integer *, logical *);

/* ***BEGIN PROLOGUE  XGETF */
/* ***PURPOSE  Return the current value of the error control flag. */
/* ***LIBRARY   SLATEC (XERROR) */
/* ***CATEGORY  R3C */
/* ***TYPE      ALL (XGETF-A) */
/* ***KEYWORDS  ERROR, XERROR */
/* ***AUTHOR  Jones, R. E., (SNLA) */
/* ***DESCRIPTION */

/*   Abstract */
/*        XGETF returns the current value of the error control flag */
/*        in KONTRL.  See subroutine XSETF for flag value meanings. */
/*        (KONTRL is an output parameter only.) */

/* ***REFERENCES  R. E. Jones and D. K. Kahaner, XERROR, the SLATEC */
/*                 Error-handling Package, SAND82-0800, Sandia */
/*                 Laboratories, 1982. */
/* ***ROUTINES CALLED  J4SAVE */
/* ***REVISION HISTORY  (YYMMDD) */
/*   790801  DATE WRITTEN */
/*   861211  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   920501  Reformatted the REFERENCES section.  (WRB) */
/* ***END PROLOGUE  XGETF */
/* ***FIRST EXECUTABLE STATEMENT  XGETF */
    *kontrl = j4save_(&c__2, &c__0, &c_false);
    return 0;
} /* xgetf_ */

