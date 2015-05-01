/* c1merg.f -- translated by f2c (version 12.02.01).
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

/* DECK C1MERG */
/* Subroutine */ int c1merg_(complex *tcos, integer *i1, integer *m1, integer 
	*i2, integer *m2, integer *i3)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer j1, j2, j3;
    extern /* Subroutine */ int ccopy_(integer *, complex *, integer *, 
	    complex *, integer *);

/* ***BEGIN PROLOGUE  C1MERG */
/* ***SUBSIDIARY */
/* ***PURPOSE  Merge two strings of complex numbers.  Each string is */
/*            ascending by the real part. */
/* ***LIBRARY   SLATEC */
/* ***TYPE      COMPLEX (S1MERG-S, D1MERG-D, C1MERG-C, I1MERG-I) */
/* ***AUTHOR  (UNKNOWN) */
/* ***DESCRIPTION */

/*   This subroutine merges two ascending strings of numbers in the */
/*   array TCOS.  The first string is of length M1 and starts at */
/*   TCOS(I1+1).  The second string is of length M2 and starts at */
/*   TCOS(I2+1).  The merged string goes into TCOS(I3+1).  The ordering */
/*   is on the real part. */

/* ***SEE ALSO  CMGNBN */
/* ***ROUTINES CALLED  CCOPY */
/* ***REVISION HISTORY  (YYMMDD) */
/*   801001  DATE WRITTEN */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900402  Added TYPE section.  (WRB) */
/*   910408  Modified to use IF-THEN-ELSE.  Make it look like MERGE */
/*           which was modified earlier due to compiler problems on */
/*           the IBM RS6000.  (RWC) */
/*   920130  Code name changed from CMPMRG to C1MERG.  (WRB) */
/* ***END PROLOGUE  C1MERG */


/* ***FIRST EXECUTABLE STATEMENT  C1MERG */
    /* Parameter adjustments */
    --tcos;

    /* Function Body */
    if (*m1 == 0 && *m2 == 0) {
	return 0;
    }

    if (*m1 == 0 && *m2 != 0) {
	ccopy_(m2, &tcos[*i2 + 1], &c__1, &tcos[*i3 + 1], &c__1);
	return 0;
    }

    if (*m1 != 0 && *m2 == 0) {
	ccopy_(m1, &tcos[*i1 + 1], &c__1, &tcos[*i3 + 1], &c__1);
	return 0;
    }

    j1 = 1;
    j2 = 1;
    j3 = 1;

L10:
    i__1 = j1 + *i1;
    i__2 = *i2 + j2;
    if (tcos[i__1].r <= tcos[i__2].r) {
	i__1 = *i3 + j3;
	i__2 = *i1 + j1;
	tcos[i__1].r = tcos[i__2].r, tcos[i__1].i = tcos[i__2].i;
	++j1;
	if (j1 > *m1) {
	    i__1 = *m2 - j2 + 1;
	    ccopy_(&i__1, &tcos[*i2 + j2], &c__1, &tcos[*i3 + j3 + 1], &c__1);
	    return 0;
	}
    } else {
	i__1 = *i3 + j3;
	i__2 = *i2 + j2;
	tcos[i__1].r = tcos[i__2].r, tcos[i__1].i = tcos[i__2].i;
	++j2;
	if (j2 > *m2) {
	    i__1 = *m1 - j1 + 1;
	    ccopy_(&i__1, &tcos[*i1 + j1], &c__1, &tcos[*i3 + j3 + 1], &c__1);
	    return 0;
	}
    }
    ++j3;
    goto L10;
} /* c1merg_ */

