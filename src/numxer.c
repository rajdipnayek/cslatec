/* numxer.f -- translated by f2c (version 12.02.01).
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
static integer c__0 = 0;
static logical c_false = FALSE_;

/* DECK NUMXER */
integer numxer_(integer *nerr)
{
    /* System generated locals */
    integer ret_val;

    /* Local variables */
    extern integer j4save_(integer *, integer *, logical *);

/* ***BEGIN PROLOGUE  NUMXER */
/* ***PURPOSE  Return the most recent error number. */
/* ***LIBRARY   SLATEC (XERROR) */
/* ***CATEGORY  R3C */
/* ***TYPE      INTEGER (NUMXER-I) */
/* ***KEYWORDS  ERROR NUMBER, XERROR */
/* ***AUTHOR  Jones, R. E., (SNLA) */
/* ***DESCRIPTION */

/*     Abstract */
/*        NUMXER returns the most recent error number, */
/*        in both NUMXER and the parameter NERR. */

/* ***REFERENCES  R. E. Jones and D. K. Kahaner, XERROR, the SLATEC */
/*                 Error-handling Package, SAND82-0800, Sandia */
/*                 Laboratories, 1982. */
/* ***ROUTINES CALLED  J4SAVE */
/* ***REVISION HISTORY  (YYMMDD) */
/*   790801  DATE WRITTEN */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900402  Added TYPE section.  (WRB) */
/*   910411  Made user-callable and added KEYWORDS section.  (WRB) */
/*   920501  Reformatted the REFERENCES section.  (WRB) */
/* ***END PROLOGUE  NUMXER */
/* ***FIRST EXECUTABLE STATEMENT  NUMXER */
    *nerr = j4save_(&c__1, &c__0, &c_false);
    ret_val = *nerr;
    return ret_val;
} /* numxer_ */

