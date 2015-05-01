/* dsvco.f -- translated by f2c (version 12.02.01).
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

/* DECK DSVCO */
/* Subroutine */ int dsvco_(doublereal *rsav, integer *isav)
{
    /* Initialized data */

    static integer lenrls = 218;
    static integer lenils = 33;

    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__;

/* ***BEGIN PROLOGUE  DSVCO */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to DDEBDF */
/* ***LIBRARY   SLATEC */
/* ***TYPE      DOUBLE PRECISION (SVCO-S, DSVCO-D) */
/* ***AUTHOR  (UNKNOWN) */
/* ***DESCRIPTION */

/*   DSVCO transfers data from a common block to arrays within the */
/*   integrator package DDEBDF. */

/* ***SEE ALSO  DDEBDF */
/* ***ROUTINES CALLED  (NONE) */
/* ***COMMON BLOCKS    DDEBD1 */
/* ***REVISION HISTORY  (YYMMDD) */
/*   820301  DATE WRITTEN */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900328  Added TYPE section.  (WRB) */
/* ***END PROLOGUE  DSVCO */
/* ----------------------------------------------------------------------- */
/* THIS ROUTINE STORES IN RSAV AND ISAV THE CONTENTS OF COMMON BLOCK */
/* DDEBD1  , WHICH IS USED INTERNALLY IN THE DDEBDF PACKAGE. */

/* RSAV = DOUBLE PRECISION ARRAY OF LENGTH 218 OR MORE. */
/* ISAV = INTEGER ARRAY OF LENGTH 33 OR MORE. */
/* ----------------------------------------------------------------------- */
    /* Parameter adjustments */
    --isav;
    --rsav;

    /* Function Body */

/* ***FIRST EXECUTABLE STATEMENT  DSVCO */
    i__1 = lenrls;
    for (i__ = 1; i__ <= i__1; ++i__) {
	rsav[i__] = ddebd1_1.rls[i__ - 1];
/* L10: */
    }
    i__1 = lenils;
    for (i__ = 1; i__ <= i__1; ++i__) {
	isav[i__] = ddebd1_1.ils[i__ - 1];
/* L20: */
    }
    return 0;
/*     ----------------------- END OF SUBROUTINE DSVCO */
/*     ----------------------- */
} /* dsvco_ */

