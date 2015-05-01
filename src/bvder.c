/* bvder.f -- translated by f2c (version 12.02.01).
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
    integer nofst;
} mlivp_;

#define mlivp_1 mlivp_

/* DECK BVDER */
/* Subroutine */ int bvder_(real *x, real *y, real *yp, real *g, integer *
	ipar)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer j, k, l, na;
    extern /* Subroutine */ int gvec_(real *, real *), fmat_(real *, real *, 
	    real *), uvec_(real *, real *, real *), uivp_(real *, real *, 
	    real *);

/* ***BEGIN PROLOGUE  BVDER */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to BVSUP */
/* ***LIBRARY   SLATEC */
/* ***TYPE      SINGLE PRECISION (BVDER-S, DBVDER-D) */
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
/*             NOFST communicates to the routine FMAT how to access */
/*             the dependent variables corresponding to this initial */
/*             value problem.  For example, during any call to FMAT, */
/*             the first dependent variable for the initial value */
/*             problem is in position  Y(NOFST + 1). */
/*             See example in SAND77-1328. */
/* ********************************************************************** */

/* ***SEE ALSO  BVSUP */
/* ***ROUTINES CALLED  (NONE) */
/* ***COMMON BLOCKS    ML8SZ, MLIVP */
/* ***REVISION HISTORY  (YYMMDD) */
/*   750601  DATE WRITTEN */
/*   890921  Realigned order of variables in certain COMMON blocks. */
/*           (WRB) */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900328  Added TYPE section.  (WRB) */
/*   910701  Corrected ROUTINES CALLED section.  (WRB) */
/*   910722  Updated AUTHOR section.  (ALS) */
/*   920618  Minor restructuring of code.  (RWC, WRB) */
/* ***END PROLOGUE  BVDER */

/* ********************************************************************** */


/* ********************************************************************** */
/*     The COMMON block below is used to communicate with the user */
/*     supplied subroutine FMAT.  The user should not alter this */
/*     COMMON block. */

/* ********************************************************************** */

/* ***FIRST EXECUTABLE STATEMENT  BVDER */
    /* Parameter adjustments */
    --g;
    --yp;
    --y;

    /* Function Body */
    if (ml8sz_1.ivp > 0) {
	uivp_(x, &y[ml8sz_1.ivp + 1], &yp[ml8sz_1.ivp + 1]);
    }
    mlivp_1.nofst = ml8sz_1.ivp;
    na = 1;
    i__1 = ml8sz_1.nfc;
    for (k = 1; k <= i__1; ++k) {
	fmat_(x, &y[na], &yp[na]);
	mlivp_1.nofst -= ml8sz_1.ncomp;
	na += ml8sz_1.ncomp;
/* L10: */
    }

    if (ml8sz_1.inhomo != 1) {
	return 0;
    }
    fmat_(x, &y[na], &yp[na]);

    if (ml8sz_1.igofx == 0) {
	return 0;
    }
    if (*x != ml8sz_1.xsav) {
	if (ml8sz_1.ivp == 0) {
	    gvec_(x, &g[1]);
	}
	if (ml8sz_1.ivp > 0) {
	    uvec_(x, &y[ml8sz_1.ivp + 1], &g[1]);
	}
	ml8sz_1.xsav = *x;
    }

/*     If the user has chosen not to normalize the particular */
/*     solution, then C is defined in BVPOR to be 1.0 */

/*     The following loop is just */
/*     CALL SAXPY (NCOMP, 1.0E0/C, G, 1, YP(NA), 1) */

    i__1 = ml8sz_1.ncomp;
    for (j = 1; j <= i__1; ++j) {
	l = na + j - 1;
	yp[l] += g[j] / ml8sz_1.c__;
/* L20: */
    }
    return 0;
} /* bvder_ */

