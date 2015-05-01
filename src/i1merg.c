/* i1merg.f -- translated by f2c (version 12.02.01).
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

static integer c__1 = 1;

/* DECK I1MERG */
/* Subroutine */ int i1merg_(real *icos, integer *i1, integer *m1, integer *
	i2, integer *m2, integer *i3)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer j1, j2, j3;
    extern /* Subroutine */ int icopy_(integer *, real *, integer *, real *, 
	    integer *);

/* ***BEGIN PROLOGUE  I1MERG */
/* ***SUBSIDIARY */
/* ***PURPOSE  Merge two strings of ascending integers. */
/* ***LIBRARY   SLATEC */
/* ***TYPE      INTEGER (S1MERG-S, D1MERG-D, C1MERG-C, I1MERG-I) */
/* ***AUTHOR  Boland, W. Robert, (LANL) */
/*           Clemens, Reginald, (PLK) */
/* ***DESCRIPTION */

/*   This subroutine merges two ascending strings of integers in the */
/*   array ICOS.  The first string is of length M1 and starts at */
/*   ICOS(I1+1).  The second string is of length M2 and starts at */
/*   ICOS(I2+1).  The merged string goes into ICOS(I3+1). */

/* ***ROUTINES CALLED  ICOPY */
/* ***REVISION HISTORY  (YYMMDD) */
/*   920202  DATE WRITTEN */
/* ***END PROLOGUE  I1MERG */


/* ***FIRST EXECUTABLE STATEMENT  I1MERG */
    /* Parameter adjustments */
    --icos;

    /* Function Body */
    if (*m1 == 0 && *m2 == 0) {
	return 0;
    }

    if (*m1 == 0 && *m2 != 0) {
	icopy_(m2, &icos[*i2 + 1], &c__1, &icos[*i3 + 1], &c__1);
	return 0;
    }

    if (*m1 != 0 && *m2 == 0) {
	icopy_(m1, &icos[*i1 + 1], &c__1, &icos[*i3 + 1], &c__1);
	return 0;
    }

    j1 = 1;
    j2 = 1;
    j3 = 1;

L10:
    if (icos[*i1 + j1] <= icos[*i2 + j2]) {
	icos[*i3 + j3] = icos[*i1 + j1];
	++j1;
	if (j1 > *m1) {
	    i__1 = *m2 - j2 + 1;
	    icopy_(&i__1, &icos[*i2 + j2], &c__1, &icos[*i3 + j3 + 1], &c__1);
	    return 0;
	}
    } else {
	icos[*i3 + j3] = icos[*i2 + j2];
	++j2;
	if (j2 > *m2) {
	    i__1 = *m1 - j1 + 1;
	    icopy_(&i__1, &icos[*i1 + j1], &c__1, &icos[*i3 + j3 + 1], &c__1);
	    return 0;
	}
    }
    ++j3;
    goto L10;
} /* i1merg_ */

