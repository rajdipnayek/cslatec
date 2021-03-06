/* xerdmp.f -- translated by f2c (version 12.02.01).
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

/* DECK XERDMP */
/* Subroutine */ int xerdmp_(void)
{
    static integer kount;
    extern /* Subroutine */ int xersve_(char *, char *, char *, integer *, 
	    integer *, integer *, integer *, ftnlen, ftnlen, ftnlen);

/* ***BEGIN PROLOGUE  XERDMP */
/* ***PURPOSE  Print the error tables and then clear them. */
/* ***LIBRARY   SLATEC (XERROR) */
/* ***CATEGORY  R3C */
/* ***TYPE      ALL (XERDMP-A) */
/* ***KEYWORDS  ERROR, XERROR */
/* ***AUTHOR  Jones, R. E., (SNLA) */
/* ***DESCRIPTION */

/*     Abstract */
/*        XERDMP prints the error tables, then clears them. */

/* ***REFERENCES  R. E. Jones and D. K. Kahaner, XERROR, the SLATEC */
/*                 Error-handling Package, SAND82-0800, Sandia */
/*                 Laboratories, 1982. */
/* ***ROUTINES CALLED  XERSVE */
/* ***REVISION HISTORY  (YYMMDD) */
/*   790801  DATE WRITTEN */
/*   861211  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900510  Changed call of XERSAV to XERSVE.  (RWC) */
/*   920501  Reformatted the REFERENCES section.  (WRB) */
/* ***END PROLOGUE  XERDMP */
/* ***FIRST EXECUTABLE STATEMENT  XERDMP */
    xersve_(" ", " ", " ", &c__0, &c__0, &c__0, &kount, (ftnlen)1, (ftnlen)1, 
	    (ftnlen)1);
    return 0;
} /* xerdmp_ */

