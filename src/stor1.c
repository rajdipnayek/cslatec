/* stor1.f -- translated by f2c (version 12.02.01).
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

/* Table of constant values */

static integer c__1 = 1;

/* DECK STOR1 */
/* Subroutine */ int stor1_(real *u, real *yh, real *v, real *yp, integer *
	ntemp, integer *ndisk, integer *ntape)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer j, nctnf;

    /* Fortran I/O blocks */
    static cilist io___3 = { 0, 0, 0, 0, 0 };


/* ***BEGIN PROLOGUE  STOR1 */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to BVSUP */
/* ***LIBRARY   SLATEC */
/* ***TYPE      SINGLE PRECISION (STOR1-S, DSTOR1-D) */
/* ***AUTHOR  Watts, H. A., (SNLA) */
/* ***DESCRIPTION */

/* ********************************************************************** */
/*             0 -- Storage at output points. */
/*     NTEMP = */
/*             1 -- Temporary storage */
/* ********************************************************************** */

/* ***SEE ALSO  BVSUP */
/* ***ROUTINES CALLED  (NONE) */
/* ***COMMON BLOCKS    ML8SZ */
/* ***REVISION HISTORY  (YYMMDD) */
/*   750601  DATE WRITTEN */
/*   890921  Realigned order of variables in certain COMMON blocks. */
/*           (WRB) */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900328  Added TYPE section.  (WRB) */
/*   910722  Updated AUTHOR section.  (ALS) */
/* ***END PROLOGUE  STOR1 */

/* ********************************************************************** */


/* ********************************************************************** */

/* ***FIRST EXECUTABLE STATEMENT  STOR1 */
    /* Parameter adjustments */
    --yp;
    --v;
    --yh;
    --u;

    /* Function Body */
    nctnf = ml8sz_1.ncomp * ml8sz_1.nfc;
    i__1 = nctnf;
    for (j = 1; j <= i__1; ++j) {
/* L10: */
	u[j] = yh[j];
    }
    if (ml8sz_1.inhomo == 1) {
	goto L30;
    }

/*   ZERO PARTICULAR SOLUTION */

    if (*ntemp == 1) {
	return 0;
    }
    i__1 = ml8sz_1.ncomp;
    for (j = 1; j <= i__1; ++j) {
/* L20: */
	v[j] = 0.f;
    }
    goto L70;

/*   NONZERO PARTICULAR SOLUTION */

L30:
    if (*ntemp == 0) {
	goto L50;
    }

    i__1 = ml8sz_1.ncomp;
    for (j = 1; j <= i__1; ++j) {
/* L40: */
	v[j] = yp[j];
    }
    return 0;

L50:
    i__1 = ml8sz_1.ncomp;
    for (j = 1; j <= i__1; ++j) {
/* L60: */
	v[j] = ml8sz_1.c__ * yp[j];
    }

/*  IS OUTPUT INFORMATION TO BE WRITTEN TO DISK */

L70:
    if (*ndisk == 1) {
	io___3.ciunit = *ntape;
	s_wsue(&io___3);
	i__1 = ml8sz_1.ncomp;
	for (j = 1; j <= i__1; ++j) {
	    do_uio(&c__1, (char *)&v[j], (ftnlen)sizeof(real));
	}
	i__2 = nctnf;
	for (j = 1; j <= i__2; ++j) {
	    do_uio(&c__1, (char *)&u[j], (ftnlen)sizeof(real));
	}
	e_wsue();
    }

    return 0;
} /* stor1_ */

