/* s1merg.f -- translated by f2c (version 12.02.01).
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

/* DECK S1MERG */
/* Subroutine */ int s1merg_(real *tcos, integer *i1, integer *m1, integer *
	i2, integer *m2, integer *i3)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer j1, j2, j3;
    extern /* Subroutine */ int scopy_(integer *, real *, integer *, real *, 
	    integer *);

/* ***BEGIN PROLOGUE  S1MERG */
/* ***SUBSIDIARY */
/* ***PURPOSE  Merge two strings of ascending real numbers. */
/* ***LIBRARY   SLATEC */
/* ***TYPE      SINGLE PRECISION (S1MERG-S, D1MERG-D, C1MERG-C, I1MERG-I) */
/* ***AUTHOR  (UNKNOWN) */
/* ***DESCRIPTION */

/*   This subroutine merges two ascending strings of numbers in the */
/*   array TCOS.  The first string is of length M1 and starts at */
/*   TCOS(I1+1).  The second string is of length M2 and starts at */
/*   TCOS(I2+1).  The merged string goes into TCOS(I3+1). */

/* ***SEE ALSO  GENBUN */
/* ***ROUTINES CALLED  SCOPY */
/* ***REVISION HISTORY  (YYMMDD) */
/*   801001  DATE WRITTEN */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900402  Added TYPE section.  (WRB) */
/*   901120  Modified to use IF-THEN-ELSE.  Previous spaghetti code did */
/*           not compile correctly with optimization on the IBM RS6000. */
/*           (RWC) */
/*   920130  Code name changed from MERGE to S1MERG.  (WRB) */
/* ***END PROLOGUE  S1MERG */


/* ***FIRST EXECUTABLE STATEMENT  S1MERG */
    /* Parameter adjustments */
    --tcos;

    /* Function Body */
    if (*m1 == 0 && *m2 == 0) {
	return 0;
    }

    if (*m1 == 0 && *m2 != 0) {
	scopy_(m2, &tcos[*i2 + 1], &c__1, &tcos[*i3 + 1], &c__1);
	return 0;
    }

    if (*m1 != 0 && *m2 == 0) {
	scopy_(m1, &tcos[*i1 + 1], &c__1, &tcos[*i3 + 1], &c__1);
	return 0;
    }

    j1 = 1;
    j2 = 1;
    j3 = 1;

L10:
    if (tcos[*i1 + j1] <= tcos[*i2 + j2]) {
	tcos[*i3 + j3] = tcos[*i1 + j1];
	++j1;
	if (j1 > *m1) {
	    i__1 = *m2 - j2 + 1;
	    scopy_(&i__1, &tcos[*i2 + j2], &c__1, &tcos[*i3 + j3 + 1], &c__1);
	    return 0;
	}
    } else {
	tcos[*i3 + j3] = tcos[*i2 + j2];
	++j2;
	if (j2 > *m2) {
	    i__1 = *m1 - j1 + 1;
	    scopy_(&i__1, &tcos[*i1 + j1], &c__1, &tcos[*i3 + j3 + 1], &c__1);
	    return 0;
	}
    }
    ++j3;
    goto L10;
} /* s1merg_ */

