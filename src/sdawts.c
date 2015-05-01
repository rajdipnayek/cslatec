/* sdawts.f -- translated by f2c (version 12.02.01).
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

/* DECK SDAWTS */
/* Subroutine */ int sdawts_(integer *neq, integer *iwt, real *rtol, real *
	atol, real *y, real *wt, real *rpar, integer *ipar)
{
    /* System generated locals */
    integer i__1;
    real r__1;

    /* Local variables */
    static integer i__;
    static real atoli, rtoli;

/* ***BEGIN PROLOGUE  SDAWTS */
/* ***SUBSIDIARY */
/* ***PURPOSE  Set error weight vector for SDASSL. */
/* ***LIBRARY   SLATEC (DASSL) */
/* ***TYPE      SINGLE PRECISION (SDAWTS-S, DDAWTS-D) */
/* ***AUTHOR  Petzold, Linda R., (LLNL) */
/* ***DESCRIPTION */
/* ----------------------------------------------------------------------- */
/*     THIS SUBROUTINE SETS THE ERROR WEIGHT VECTOR */
/*     WT ACCORDING TO WT(I)=RTOL(I)*ABS(Y(I))+ATOL(I), */
/*     I=1,-,N. */
/*     RTOL AND ATOL ARE SCALARS IF IWT = 0, */
/*     AND VECTORS IF IWT = 1. */
/* ----------------------------------------------------------------------- */
/* ***ROUTINES CALLED  (NONE) */
/* ***REVISION HISTORY  (YYMMDD) */
/*   830315  DATE WRITTEN */
/*   901009  Finished conversion to SLATEC 4.0 format (F.N.Fritsch) */
/*   901019  Merged changes made by C. Ulrich with SLATEC 4.0 format. */
/*   901026  Added explicit declarations for all variables and minor */
/*           cosmetic changes to prologue.  (FNF) */
/* ***END PROLOGUE  SDAWTS */



/* ***FIRST EXECUTABLE STATEMENT  SDAWTS */
    /* Parameter adjustments */
    --ipar;
    --rpar;
    --wt;
    --y;
    --atol;
    --rtol;

    /* Function Body */
    rtoli = rtol[1];
    atoli = atol[1];
    i__1 = *neq;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (*iwt == 0) {
	    goto L10;
	}
	rtoli = rtol[i__];
	atoli = atol[i__];
L10:
	wt[i__] = rtoli * (r__1 = y[i__], dabs(r__1)) + atoli;
/* L20: */
    }
    return 0;
/* -----------END OF SUBROUTINE SDAWTS------------------------------------ */
} /* sdawts_ */

