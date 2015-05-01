/* stway.f -- translated by f2c (version 12.02.01).
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
    real c__, xsav;
    integer igofx, inhomo, ivp, ncomp, nfc;
} ml8sz_;

#define ml8sz_1 ml8sz_

struct {
    real px, pwcnd, tnd, x, xbeg, xend, xot, xop;
    integer info[15], istkop, knswot, kop, lotjp, mnswot, nswot;
} ml15to_;

#define ml15to_1 ml15to_

struct {
    real ae, re, tol;
    integer nxpts, nic, nopg, mxnon, ndisk, ntape, neq, indpvt, integ, nps, 
	    ntp, neqivp, numort, nfcc, icoco;
} ml18jr_;

#define ml18jr_1 ml18jr_

/* Table of constant values */

static integer c__1 = 1;
static integer c__0 = 0;

/* DECK STWAY */
/* Subroutine */ int stway_(real *u, real *v, real *yhp, integer *inout, real 
	*stowa)
{
    /* System generated locals */
    integer i__1;
    alist al__1;

    /* Local variables */
    static integer j, k, ko, ks, ksj;
    extern /* Subroutine */ int stor1_(real *, real *, real *, real *, 
	    integer *, integer *, integer *);

/* ***BEGIN PROLOGUE  STWAY */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to BVSUP */
/* ***LIBRARY   SLATEC */
/* ***TYPE      SINGLE PRECISION (STWAY-S, DSTWAY-D) */
/* ***AUTHOR  Watts, H. A., (SNLA) */
/* ***DESCRIPTION */

/*  This subroutine stores (recalls) integration data in the event */
/*  that a restart is needed (the homogeneous solution vectors become */
/*  too dependent to continue) */

/* ***SEE ALSO  BVSUP */
/* ***ROUTINES CALLED  STOR1 */
/* ***COMMON BLOCKS    ML15TO, ML18JR, ML8SZ */
/* ***REVISION HISTORY  (YYMMDD) */
/*   750601  DATE WRITTEN */
/*   890921  Realigned order of variables in certain COMMON blocks. */
/*           (WRB) */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900328  Added TYPE section.  (WRB) */
/*   910722  Updated AUTHOR section.  (ALS) */
/* ***END PROLOGUE  STWAY */



/* ***FIRST EXECUTABLE STATEMENT  STWAY */
    /* Parameter adjustments */
    --stowa;
    --yhp;
    --v;
    --u;

    /* Function Body */
    if (*inout == 1) {
	goto L100;
    }

/*     SAVE IN STOWA ARRAY AND ISTKOP */

    ks = ml8sz_1.nfc * ml8sz_1.ncomp;
    stor1_(&stowa[1], &u[1], &stowa[ks + 1], &v[1], &c__1, &c__0, &c__0);
    ks += ml8sz_1.ncomp;
    if (ml18jr_1.neqivp == 0) {
	goto L50;
    }
    i__1 = ml18jr_1.neqivp;
    for (j = 1; j <= i__1; ++j) {
	ksj = ks + j;
/* L25: */
	stowa[ksj] = yhp[ksj];
    }
L50:
    ks += ml18jr_1.neqivp;
    stowa[ks + 1] = ml15to_1.x;
    ml15to_1.istkop = ml15to_1.kop;
    if (ml15to_1.xop == ml15to_1.x) {
	ml15to_1.istkop = ml15to_1.kop + 1;
    }
    return 0;

/*     RECALL FROM STOWA ARRAY AND ISTKOP */

L100:
    ks = ml8sz_1.nfc * ml8sz_1.ncomp;
    stor1_(&yhp[1], &stowa[1], &yhp[ks + 1], &stowa[ks + 1], &c__1, &c__0, &
	    c__0);
    ks += ml8sz_1.ncomp;
    if (ml18jr_1.neqivp == 0) {
	goto L150;
    }
    i__1 = ml18jr_1.neqivp;
    for (j = 1; j <= i__1; ++j) {
	ksj = ks + j;
/* L125: */
	yhp[ksj] = stowa[ksj];
    }
L150:
    ks += ml18jr_1.neqivp;
    ml15to_1.x = stowa[ks + 1];
    ml15to_1.info[0] = 0;
    ko = ml15to_1.kop - ml15to_1.istkop;
    ml15to_1.kop = ml15to_1.istkop;
    if (ml18jr_1.ndisk == 0 || ko == 0) {
	return 0;
    }
    i__1 = ko;
    for (k = 1; k <= i__1; ++k) {
/* L175: */
	al__1.aerr = 0;
	al__1.aunit = ml18jr_1.ntape;
	f_back(&al__1);
    }
    return 0;
} /* stway_ */

