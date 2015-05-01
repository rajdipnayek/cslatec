/* dbvder.f -- translated by f2c (version 12.02.01).
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
    integer nofst;
} dmlivp_;

#define dmlivp_1 dmlivp_

/* DECK DBVDER */
/* Subroutine */ int dbvder_(doublereal *x, doublereal *y, doublereal *yp, 
	doublereal *g, integer *ipar)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer j, k, l, na;
    extern /* Subroutine */ int dgvec_(doublereal *, doublereal *), dfmat_(
	    doublereal *, doublereal *, doublereal *), duvec_(doublereal *, 
	    doublereal *, doublereal *), duivp_(doublereal *, doublereal *, 
	    doublereal *);

/* ***BEGIN PROLOGUE  DBVDER */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to DBVSUP */
/* ***LIBRARY   SLATEC */
/* ***TYPE      DOUBLE PRECISION (BVDER-S, DBVDER-D) */
/* ***AUTHOR  Watts, H. A., (SNLA) */
/* ***DESCRIPTION */

/* ********************************************************************** */
/*     NFC = Number of base solution vectors */

/*     NCOMP = Number of components per solution vector */

/*              1 -- Nonzero particular solution */
/*     INHOMO = */
/*              2 or 3 -- Zero particular solution */

/*             0 -- Inhomogeneous vector term G(X) identically zero */
/*     IGOFX = */
/*             1 -- Inhomogeneous vector term G(X) not identically zero */

/*     G = Inhomogeneous vector term G(X) */

/*     XSAV = Previous value of X */

/*     C = Normalization factor for the particular solution */

/*           0   ( if  NEQIVP = 0 ) */
/*     IVP = */
/*           Number of differential equations integrated due to */
/*           the original boundary value problem   ( if  NEQIVP .GT. 0 ) */

/*     NOFST - For problems with auxiliary initial value equations, */
/*             NOFST communicates to the routine DFMAT how to access */
/*             the dependent variables corresponding to this initial */
/*             value problem.  For example, during any call to DFMAT, */
/*             the first dependent variable for the initial value */
/*             problem is in position  Y(NOFST + 1). */
/*             See example in SAND77-1328. */
/* ********************************************************************** */

/* ***SEE ALSO  DBVSUP */
/* ***ROUTINES CALLED  (NONE) */
/* ***COMMON BLOCKS    DML8SZ, DMLIVP */
/* ***REVISION HISTORY  (YYMMDD) */
/*   750601  DATE WRITTEN */
/*   890831  Modified array declarations.  (WRB) */
/*   890921  Realigned order of variables in certain COMMON blocks. */
/*           (WRB) */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900328  Added TYPE section.  (WRB) */
/*   910701  Corrected ROUTINES CALLED section.  (WRB) */
/*   910722  Updated AUTHOR section.  (ALS) */
/*   920618  Minor restructuring of code.  (RWC, WRB) */
/* ***END PROLOGUE  DBVDER */

/* ********************************************************************** */


/* ********************************************************************** */
/*     The COMMON block below is used to communicate with the user */
/*     supplied subroutine DFMAT.  The user should not alter this */
/*     COMMON block. */

/* ********************************************************************** */

/* ***FIRST EXECUTABLE STATEMENT  DBVDER */
    /* Parameter adjustments */
    --g;
    --yp;
    --y;

    /* Function Body */
    if (dml8sz_1.ivp > 0) {
	duivp_(x, &y[dml8sz_1.ivp + 1], &yp[dml8sz_1.ivp + 1]);
    }
    dmlivp_1.nofst = dml8sz_1.ivp;
    na = 1;
    i__1 = dml8sz_1.nfc;
    for (k = 1; k <= i__1; ++k) {
	dfmat_(x, &y[na], &yp[na]);
	dmlivp_1.nofst -= dml8sz_1.ncomp;
	na += dml8sz_1.ncomp;
/* L10: */
    }

    if (dml8sz_1.inhomo != 1) {
	return 0;
    }
    dfmat_(x, &y[na], &yp[na]);

    if (dml8sz_1.igofx == 0) {
	return 0;
    }
    if (*x != dml8sz_1.xsav) {
	if (dml8sz_1.ivp == 0) {
	    dgvec_(x, &g[1]);
	}
	if (dml8sz_1.ivp > 0) {
	    duvec_(x, &y[dml8sz_1.ivp + 1], &g[1]);
	}
	dml8sz_1.xsav = *x;
    }

/*     If the user has chosen not to normalize the particular */
/*     solution, then C is defined in DBVPOR to be 1.0 */

/*     The following loop is just */
/*     CALL DAXPY (NCOMP, 1.0D0/C, G, 1, YP(NA), 1) */

    i__1 = dml8sz_1.ncomp;
    for (j = 1; j <= i__1; ++j) {
	l = na + j - 1;
	yp[l] += g[j] / dml8sz_1.c__;
/* L20: */
    }
    return 0;
} /* dbvder_ */

