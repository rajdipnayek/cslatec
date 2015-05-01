/* svco.f -- translated by f2c (version 12.02.01).
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

/* Common Block Declarations */

struct {
    real rls[218];
    integer ils[33];
} debdf1_;

#define debdf1_1 debdf1_

/* DECK SVCO */
/* Subroutine */ int svco_(real *rsav, integer *isav)
{
    /* Initialized data */

    static integer lenrls = 218;
    static integer lenils = 33;

    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__;

/* ***BEGIN PROLOGUE  SVCO */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to DEBDF */
/* ***LIBRARY   SLATEC */
/* ***TYPE      SINGLE PRECISION (SVCO-S, DSVCO-D) */
/* ***AUTHOR  (UNKNOWN) */
/* ***DESCRIPTION */

/*   SVCO transfers data from a common block to arrays within the */
/*   integrator package DEBDF. */

/* ***SEE ALSO  DEBDF */
/* ***ROUTINES CALLED  (NONE) */
/* ***COMMON BLOCKS    DEBDF1 */
/* ***REVISION HISTORY  (YYMMDD) */
/*   800901  DATE WRITTEN */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900328  Added TYPE section.  (WRB) */
/* ***END PROLOGUE  SVCO */


/* ----------------------------------------------------------------------- */
/* THIS ROUTINE STORES IN RSAV AND ISAV THE CONTENTS OF COMMON BLOCK */
/* DEBDF1  , WHICH IS USED INTERNALLY IN THE DEBDF PACKAGE. */

/* RSAV = REAL ARRAY OF LENGTH 218 OR MORE. */
/* ISAV = INTEGER ARRAY OF LENGTH 33 OR MORE. */
/* ----------------------------------------------------------------------- */
    /* Parameter adjustments */
    --isav;
    --rsav;

    /* Function Body */

/* ***FIRST EXECUTABLE STATEMENT  SVCO */
    i__1 = lenrls;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L10: */
	rsav[i__] = debdf1_1.rls[i__ - 1];
    }
    i__1 = lenils;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L20: */
	isav[i__] = debdf1_1.ils[i__ - 1];
    }
    return 0;
/* ----------------------- END OF SUBROUTINE SVCO ----------------------- */
} /* svco_ */

