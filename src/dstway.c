/* dstway.f -- translated by f2c (version 12.02.01).
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
    doublereal c__, xsav;
    integer igofx, inhomo, ivp, ncomp, nfc;
} dml8sz_;

#define dml8sz_1 dml8sz_

struct {
    doublereal px, pwcnd, tnd, x, xbeg, xend, xot, xop;
    integer info[15], istkop, knswot, kop, lotjp, mnswot, nswot;
} dml15t_;

#define dml15t_1 dml15t_

struct {
    doublereal ae, re, tol;
    integer nxpts, nic, nopg, mxnon, ndisk, ntape, neq, indpvt, integ, nps, 
	    ntp, neqivp, numort, nfcc, icoco;
} dml18j_;

#define dml18j_1 dml18j_

/* Table of constant values */

static integer c__1 = 1;
static integer c__0 = 0;

/* DECK DSTWAY */
/* Subroutine */ int dstway_(doublereal *u, doublereal *v, doublereal *yhp, 
	integer *inout, doublereal *stowa)
{
    /* System generated locals */
    integer i__1;
    alist al__1;

    /* Local variables */
    static integer j, k, ko, ks, ksj;
    extern /* Subroutine */ int dstor1_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, integer *, integer *, integer *);

/* ***BEGIN PROLOGUE  DSTWAY */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to DBVSUP */
/* ***LIBRARY   SLATEC */
/* ***TYPE      DOUBLE PRECISION (STWAY-S, DSTWAY-D) */
/* ***AUTHOR  Watts, H. A., (SNLA) */
/* ***DESCRIPTION */

/*  This subroutine stores (recalls) integration data in the event */
/*  that a restart is needed (the homogeneous solution vectors become */
/*  too dependent to continue). */

/* ***SEE ALSO  DBVSUP */
/* ***ROUTINES CALLED  DSTOR1 */
/* ***COMMON BLOCKS    DML15T, DML18J, DML8SZ */
/* ***REVISION HISTORY  (YYMMDD) */
/*   750601  DATE WRITTEN */
/*   890831  Modified array declarations.  (WRB) */
/*   890921  Realigned order of variables in certain COMMON blocks. */
/*           (WRB) */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900328  Added TYPE section.  (WRB) */
/*   910722  Updated AUTHOR section.  (ALS) */
/* ***END PROLOGUE  DSTWAY */



/* ***FIRST EXECUTABLE STATEMENT  DSTWAY */
    /* Parameter adjustments */
    --stowa;
    --yhp;
    --v;
    --u;

    /* Function Body */
    if (*inout == 1) {
	goto L30;
    }

/*        SAVE IN STOWA ARRAY AND ISTKOP */

    ks = dml8sz_1.nfc * dml8sz_1.ncomp;
    dstor1_(&stowa[1], &u[1], &stowa[ks + 1], &v[1], &c__1, &c__0, &c__0);
    ks += dml8sz_1.ncomp;
    if (dml18j_1.neqivp < 1) {
	goto L20;
    }
    i__1 = dml18j_1.neqivp;
    for (j = 1; j <= i__1; ++j) {
	ksj = ks + j;
	stowa[ksj] = yhp[ksj];
/* L10: */
    }
L20:
    ks += dml18j_1.neqivp;
    stowa[ks + 1] = dml15t_1.x;
    dml15t_1.istkop = dml15t_1.kop;
    if (dml15t_1.xop == dml15t_1.x) {
	dml15t_1.istkop = dml15t_1.kop + 1;
    }
    goto L80;
L30:

/*        RECALL FROM STOWA ARRAY AND ISTKOP */

    ks = dml8sz_1.nfc * dml8sz_1.ncomp;
    dstor1_(&yhp[1], &stowa[1], &yhp[ks + 1], &stowa[ks + 1], &c__1, &c__0, &
	    c__0);
    ks += dml8sz_1.ncomp;
    if (dml18j_1.neqivp < 1) {
	goto L50;
    }
    i__1 = dml18j_1.neqivp;
    for (j = 1; j <= i__1; ++j) {
	ksj = ks + j;
	yhp[ksj] = stowa[ksj];
/* L40: */
    }
L50:
    ks += dml18j_1.neqivp;
    dml15t_1.x = stowa[ks + 1];
    dml15t_1.info[0] = 0;
    ko = dml15t_1.kop - dml15t_1.istkop;
    dml15t_1.kop = dml15t_1.istkop;
    if (dml18j_1.ndisk == 0 || ko == 0) {
	goto L70;
    }
    i__1 = ko;
    for (k = 1; k <= i__1; ++k) {
	al__1.aerr = 0;
	al__1.aunit = dml18j_1.ntape;
	f_back(&al__1);
/* L60: */
    }
L70:
L80:
    return 0;
} /* dstway_ */

