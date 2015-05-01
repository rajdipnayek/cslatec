/* drsco.f -- translated by f2c (version 12.02.01).
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
    doublereal rls[218];
    integer ils[33];
} ddebd1_;

#define ddebd1_1 ddebd1_

/* DECK DRSCO */
/* Subroutine */ int drsco_(doublereal *rsav, integer *isav)
{
    /* Initialized data */

    static integer lenrls = 218;
    static integer lenils = 33;

    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__;

/* ***BEGIN PROLOGUE  DRSCO */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to DDEBDF */
/* ***LIBRARY   SLATEC */
/* ***TYPE      DOUBLE PRECISION (RSCO-S, DRSCO-D) */
/* ***AUTHOR  Watts, H. A., (SNLA) */
/* ***DESCRIPTION */

/*   DRSCO transfers data from arrays to a common block within the */
/*   integrator package DDEBDF. */

/* ***SEE ALSO  DDEBDF */
/* ***ROUTINES CALLED  (NONE) */
/* ***COMMON BLOCKS    DDEBD1 */
/* ***REVISION HISTORY  (YYMMDD) */
/*   820301  DATE WRITTEN */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900328  Added TYPE section.  (WRB) */
/*   910722  Updated AUTHOR section.  (ALS) */
/* ***END PROLOGUE  DRSCO */
/* ----------------------------------------------------------------------- */
/* THIS ROUTINE RESTORES FROM RSAV AND ISAV THE CONTENTS OF COMMON */
/* BLOCK DDEBD1  , WHICH IS USED INTERNALLY IN THE DDEBDF */
/* PACKAGE.  THIS PRESUMES THAT RSAV AND ISAV WERE LOADED BY MEANS */
/* OF SUBROUTINE DSVCO OR THE EQUIVALENT. */
/* ----------------------------------------------------------------------- */

    /* Parameter adjustments */
    --isav;
    --rsav;

    /* Function Body */

/* ***FIRST EXECUTABLE STATEMENT  DRSCO */
    i__1 = lenrls;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ddebd1_1.rls[i__ - 1] = rsav[i__];
/* L10: */
    }
    i__1 = lenils;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ddebd1_1.ils[i__ - 1] = isav[i__];
/* L20: */
    }
    return 0;
/*     ----------------------- END OF SUBROUTINE DRSCO */
/*     ----------------------- */
} /* drsco_ */

