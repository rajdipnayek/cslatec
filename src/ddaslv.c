/* ddaslv.f -- translated by f2c (version 12.02.01).
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

/* Table of constant values */

static integer c__0 = 0;

/* DECK DDASLV */
/* Subroutine */ int ddaslv_(integer *neq, doublereal *delta, doublereal *wm, 
	integer *iwm)
{
    extern /* Subroutine */ int dgbsl_(doublereal *, integer *, integer *, 
	    integer *, integer *, integer *, doublereal *, integer *), dgesl_(
	    doublereal *, integer *, integer *, integer *, doublereal *, 
	    integer *);
    static integer mtype, meband;

/* ***BEGIN PROLOGUE  DDASLV */
/* ***SUBSIDIARY */
/* ***PURPOSE  Linear system solver for DDASSL. */
/* ***LIBRARY   SLATEC (DASSL) */
/* ***TYPE      DOUBLE PRECISION (SDASLV-S, DDASLV-D) */
/* ***AUTHOR  Petzold, Linda R., (LLNL) */
/* ***DESCRIPTION */
/* ----------------------------------------------------------------------- */
/*     THIS ROUTINE MANAGES THE SOLUTION OF THE LINEAR */
/*     SYSTEM ARISING IN THE NEWTON ITERATION. */
/*     MATRICES AND REAL TEMPORARY STORAGE AND */
/*     REAL INFORMATION ARE STORED IN THE ARRAY WM. */
/*     INTEGER MATRIX INFORMATION IS STORED IN */
/*     THE ARRAY IWM. */
/*     FOR A DENSE MATRIX, THE LINPACK ROUTINE */
/*     DGESL IS CALLED. */
/*     FOR A BANDED MATRIX,THE LINPACK ROUTINE */
/*     DGBSL IS CALLED. */
/* ----------------------------------------------------------------------- */
/* ***ROUTINES CALLED  DGBSL, DGESL */
/* ***REVISION HISTORY  (YYMMDD) */
/*   830315  DATE WRITTEN */
/*   901009  Finished conversion to SLATEC 4.0 format (F.N.Fritsch) */
/*   901019  Merged changes made by C. Ulrich with SLATEC 4.0 format. */
/*   901026  Added explicit declarations for all variables and minor */
/*           cosmetic changes to prologue.  (FNF) */
/* ***END PROLOGUE  DDASLV */




/* ***FIRST EXECUTABLE STATEMENT  DDASLV */
    /* Parameter adjustments */
    --iwm;
    --wm;
    --delta;

    /* Function Body */
    mtype = iwm[4];
    switch (mtype) {
	case 1:  goto L100;
	case 2:  goto L100;
	case 3:  goto L300;
	case 4:  goto L400;
	case 5:  goto L400;
    }

/*     DENSE MATRIX */
L100:
    dgesl_(&wm[1], neq, neq, &iwm[21], &delta[1], &c__0);
    return 0;

/*     DUMMY SECTION FOR MTYPE=3 */
L300:
    return 0;

/*     BANDED MATRIX */
L400:
    meband = (iwm[1] << 1) + iwm[2] + 1;
    dgbsl_(&wm[1], &meband, neq, &iwm[1], &iwm[2], &iwm[21], &delta[1], &c__0)
	    ;
    return 0;
/* ------END OF SUBROUTINE DDASLV------ */
} /* ddaslv_ */

