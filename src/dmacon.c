/* dmacon.f -- translated by f2c (version 12.02.01).
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
    doublereal uro, sru, eps, sqovfl, twou, fouru;
    integer lpar;
} dml5mc_;

#define dml5mc_1 dml5mc_

/* Table of constant values */

static integer c__4 = 4;
static doublereal c_b3 = 10.;
static integer c__2 = 2;

/* DECK DMACON */
/* Subroutine */ int dmacon_(void)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static doublereal dd;
    static integer ke;
    extern doublereal d1mach_(integer *);

/* ***BEGIN PROLOGUE  DMACON */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to DBVSUP */
/* ***LIBRARY   SLATEC */
/* ***TYPE      DOUBLE PRECISION (MACON-S, DMACON-D) */
/* ***AUTHOR  (UNKNOWN) */
/* ***SEE ALSO  DBVSUP */
/* ***ROUTINES CALLED  D1MACH */
/* ***COMMON BLOCKS    DML5MC */
/* ***REVISION HISTORY  (YYMMDD) */
/*   750601  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890921  Realigned order of variables in certain COMMON blocks. */
/*           (WRB) */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900328  Added TYPE section.  (WRB) */
/* ***END PROLOGUE  DMACON */
/* ***FIRST EXECUTABLE STATEMENT  DMACON */
    dml5mc_1.uro = d1mach_(&c__4);
    dml5mc_1.sru = sqrt(dml5mc_1.uro);
    dd = -d_lg10(&dml5mc_1.uro);
    dml5mc_1.lpar = (integer) (dd * .5);
    ke = (integer) (dd * .75 + .5);
    i__1 = ke * -2;
    dml5mc_1.eps = pow_di(&c_b3, &i__1);
    dml5mc_1.sqovfl = sqrt(d1mach_(&c__2));
    dml5mc_1.twou = dml5mc_1.uro * 2.;
    dml5mc_1.fouru = dml5mc_1.uro * 4.;
    return 0;
} /* dmacon_ */

