/* macon.f -- translated by f2c (version 12.02.01).
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
    real uro, sru, eps, sqovfl, twou, fouru;
    integer lpar;
} ml5mco_;

#define ml5mco_1 ml5mco_

/* Table of constant values */

static integer c__4 = 4;
static real c_b3 = 10.f;
static integer c__2 = 2;

/* DECK MACON */
/* Subroutine */ int macon_(void)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static real dd;
    static integer ke;
    extern doublereal r1mach_(integer *);

/* ***BEGIN PROLOGUE  MACON */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to BVSUP */
/* ***LIBRARY   SLATEC */
/* ***TYPE      SINGLE PRECISION (MACON-S, DMACON-D) */
/* ***AUTHOR  (UNKNOWN) */
/* ***DESCRIPTION */

/*    Sets up machine constants using R1MACH */

/* ***SEE ALSO  BVSUP */
/* ***ROUTINES CALLED  R1MACH */
/* ***COMMON BLOCKS    ML5MCO */
/* ***REVISION HISTORY  (YYMMDD) */
/*   750601  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890921  Realigned order of variables in certain COMMON blocks. */
/*           (WRB) */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900328  Added TYPE section.  (WRB) */
/* ***END PROLOGUE  MACON */
/* ***FIRST EXECUTABLE STATEMENT  MACON */
    ml5mco_1.uro = r1mach_(&c__4);
    ml5mco_1.sru = sqrt(ml5mco_1.uro);
    dd = -r_lg10(&ml5mco_1.uro);
    ml5mco_1.lpar = dd * .5f;
    ke = dd * .75f + .5f;
    i__1 = ke * -2;
    ml5mco_1.eps = pow_ri(&c_b3, &i__1);
    ml5mco_1.sqovfl = sqrt(r1mach_(&c__2));
    ml5mco_1.twou = ml5mco_1.uro * 2.f;
    ml5mco_1.fouru = ml5mco_1.uro * 4.f;
    return 0;
} /* macon_ */

