/* xsetun.f -- translated by f2c (version 12.02.01).
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
static logical c_true = TRUE_;
static integer c__5 = 5;
static integer c__1 = 1;

/* DECK XSETUN */
/* Subroutine */ int xsetun_(integer *iunit)
{
    static integer junk;
    extern integer j4save_(integer *, integer *, logical *);

/* ***BEGIN PROLOGUE  XSETUN */
/* ***PURPOSE  Set output file to which error messages are to be sent. */
/* ***LIBRARY   SLATEC (XERROR) */
/* ***CATEGORY  R3B */
/* ***TYPE      ALL (XSETUN-A) */
/* ***KEYWORDS  ERROR, XERROR */
/* ***AUTHOR  Jones, R. E., (SNLA) */
/* ***DESCRIPTION */

/*     Abstract */
/*        XSETUN sets the output file to which error messages are to */
/*        be sent.  Only one file will be used.  See XSETUA for */
/*        how to declare more than one file. */

/*     Description of Parameter */
/*      --Input-- */
/*        IUNIT - an input parameter giving the logical unit number */
/*                to which error messages are to be sent. */

/* ***REFERENCES  R. E. Jones and D. K. Kahaner, XERROR, the SLATEC */
/*                 Error-handling Package, SAND82-0800, Sandia */
/*                 Laboratories, 1982. */
/* ***ROUTINES CALLED  J4SAVE */
/* ***REVISION HISTORY  (YYMMDD) */
/*   790801  DATE WRITTEN */
/*   861211  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   920501  Reformatted the REFERENCES section.  (WRB) */
/* ***END PROLOGUE  XSETUN */
/* ***FIRST EXECUTABLE STATEMENT  XSETUN */
    junk = j4save_(&c__3, iunit, &c_true);
    junk = j4save_(&c__5, &c__1, &c_true);
    return 0;
} /* xsetun_ */

