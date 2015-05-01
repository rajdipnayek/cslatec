/* dstor1.f -- translated by f2c (version 12.02.01).
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

/* Table of constant values */

static integer c__1 = 1;

/* DECK DSTOR1 */
/* Subroutine */ int dstor1_(doublereal *u, doublereal *yh, doublereal *v, 
	doublereal *yp, integer *ntemp, integer *ndisk, integer *ntape)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer j, nctnf;

    /* Fortran I/O blocks */
    static cilist io___3 = { 0, 0, 0, 0, 0 };


/* ***BEGIN PROLOGUE  DSTOR1 */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to DBVSUP */
/* ***LIBRARY   SLATEC */
/* ***TYPE      DOUBLE PRECISION (STOR1-S, DSTOR1-D) */
/* ***AUTHOR  Watts, H. A., (SNLA) */
/* ***DESCRIPTION */

/* ********************************************************************** */
/*             0 -- storage at output points. */
/*     NTEMP = */
/*             1 -- temporary storage */
/* ********************************************************************** */

/* ***SEE ALSO  DBVSUP */
/* ***ROUTINES CALLED  (NONE) */
/* ***COMMON BLOCKS    DML8SZ */
/* ***REVISION HISTORY  (YYMMDD) */
/*   750601  DATE WRITTEN */
/*   890831  Modified array declarations.  (WRB) */
/*   890921  Realigned order of variables in certain COMMON blocks. */
/*           (WRB) */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900328  Added TYPE section.  (WRB) */
/*   910722  Updated AUTHOR section.  (ALS) */
/* ***END PROLOGUE  DSTOR1 */

/*     ****************************************************************** */


/*      ***************************************************************** */

/*     BEGIN BLOCK PERMITTING ...EXITS TO 80 */
/* ***FIRST EXECUTABLE STATEMENT  DSTOR1 */
    /* Parameter adjustments */
    --yp;
    --v;
    --yh;
    --u;

    /* Function Body */
    nctnf = dml8sz_1.ncomp * dml8sz_1.nfc;
    i__1 = nctnf;
    for (j = 1; j <= i__1; ++j) {
	u[j] = yh[j];
/* L10: */
    }
    if (dml8sz_1.inhomo == 1) {
	goto L30;
    }

/*           ZERO PARTICULAR SOLUTION */

/*     ......EXIT */
    if (*ntemp == 1) {
	goto L80;
    }
    i__1 = dml8sz_1.ncomp;
    for (j = 1; j <= i__1; ++j) {
	v[j] = 0.;
/* L20: */
    }
    goto L70;
L30:

/*           NONZERO PARTICULAR SOLUTION */

    if (*ntemp == 0) {
	goto L50;
    }

    i__1 = dml8sz_1.ncomp;
    for (j = 1; j <= i__1; ++j) {
	v[j] = yp[j];
/* L40: */
    }
/*     .........EXIT */
    goto L80;
L50:

    i__1 = dml8sz_1.ncomp;
    for (j = 1; j <= i__1; ++j) {
	v[j] = dml8sz_1.c__ * yp[j];
/* L60: */
    }
L70:

/*        IS OUTPUT INFORMATION TO BE WRITTEN TO DISK */

    if (*ndisk == 1) {
	io___3.ciunit = *ntape;
	s_wsue(&io___3);
	i__1 = dml8sz_1.ncomp;
	for (j = 1; j <= i__1; ++j) {
	    do_uio(&c__1, (char *)&v[j], (ftnlen)sizeof(doublereal));
	}
	i__2 = nctnf;
	for (j = 1; j <= i__2; ++j) {
	    do_uio(&c__1, (char *)&u[j], (ftnlen)sizeof(doublereal));
	}
	e_wsue();
    }
L80:

    return 0;
} /* dstor1_ */

