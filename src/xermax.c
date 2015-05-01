/* xermax.f -- translated by f2c (version 12.02.01).
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

static integer c__4 = 4;
static logical c_true = TRUE_;

/* DECK XERMAX */
/* Subroutine */ int xermax_(integer *max__)
{
    static integer junk;
    extern integer j4save_(integer *, integer *, logical *);

/* ***BEGIN PROLOGUE  XERMAX */
/* ***PURPOSE  Set maximum number of times any error message is to be */
/*            printed. */
/* ***LIBRARY   SLATEC (XERROR) */
/* ***CATEGORY  R3C */
/* ***TYPE      ALL (XERMAX-A) */
/* ***KEYWORDS  ERROR, XERROR */
/* ***AUTHOR  Jones, R. E., (SNLA) */
/* ***DESCRIPTION */

/*     Abstract */
/*        XERMAX sets the maximum number of times any message */
/*        is to be printed.  That is, non-fatal messages are */
/*        not to be printed after they have occurred MAX times. */
/*        Such non-fatal messages may be printed less than */
/*        MAX times even if they occur MAX times, if error */
/*        suppression mode (KONTRL=0) is ever in effect. */

/*     Description of Parameter */
/*      --Input-- */
/*        MAX - the maximum number of times any one message */
/*              is to be printed. */

/* ***REFERENCES  R. E. Jones and D. K. Kahaner, XERROR, the SLATEC */
/*                 Error-handling Package, SAND82-0800, Sandia */
/*                 Laboratories, 1982. */
/* ***ROUTINES CALLED  J4SAVE */
/* ***REVISION HISTORY  (YYMMDD) */
/*   790801  DATE WRITTEN */
/*   861211  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   920501  Reformatted the REFERENCES section.  (WRB) */
/* ***END PROLOGUE  XERMAX */
/* ***FIRST EXECUTABLE STATEMENT  XERMAX */
    junk = j4save_(&c__4, max__, &c_true);
    return 0;
} /* xermax_ */

