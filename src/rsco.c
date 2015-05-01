/* rsco.f -- translated by f2c (version 12.02.01).
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

/* DECK RSCO */
/* Subroutine */ int rsco_(real *rsav, integer *isav)
{
    /* Initialized data */

    static integer lenrls = 218;
    static integer lenils = 33;

    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__;

/* ***BEGIN PROLOGUE  RSCO */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to DEBDF */
/* ***LIBRARY   SLATEC */
/* ***TYPE      SINGLE PRECISION (RSCO-S, DRSCO-D) */
/* ***AUTHOR  Watts, H. A., (SNLA) */
/* ***DESCRIPTION */

/*   RSCO transfers data from arrays to a common block within the */
/*   integrator package DEBDF. */

/* ***SEE ALSO  DEBDF */
/* ***ROUTINES CALLED  (NONE) */
/* ***COMMON BLOCKS    DEBDF1 */
/* ***REVISION HISTORY  (YYMMDD) */
/*   800901  DATE WRITTEN */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900328  Added TYPE section.  (WRB) */
/*   910722  Updated AUTHOR section.  (ALS) */
/* ***END PROLOGUE  RSCO */


/* ----------------------------------------------------------------------- */
/* THIS ROUTINE RESTORES FROM RSAV AND ISAV THE CONTENTS OF COMMON */
/* BLOCK DEBDF1  , WHICH IS USED INTERNALLY IN THE DEBDF */
/* PACKAGE.  THIS PRESUMES THAT RSAV AND ISAV WERE LOADED BY MEANS */
/* OF SUBROUTINE SVCO OR THE EQUIVALENT. */
/* ----------------------------------------------------------------------- */
    /* Parameter adjustments */
    --isav;
    --rsav;

    /* Function Body */

/* ***FIRST EXECUTABLE STATEMENT  RSCO */
    i__1 = lenrls;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L10: */
	debdf1_1.rls[i__ - 1] = rsav[i__];
    }
    i__1 = lenils;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L20: */
	debdf1_1.ils[i__ - 1] = isav[i__];
    }
    return 0;
/* ----------------------- END OF SUBROUTINE RSCO ----------------------- */
} /* rsco_ */

