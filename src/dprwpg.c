/* dprwpg.f -- translated by f2c (version 12.02.01).
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

static integer c__55 = 55;
static integer c__1 = 1;

/* DECK DPRWPG */
/* Subroutine */ int dprwpg_(integer *key, integer *ipage, integer *lpg, 
	doublereal *sx, integer *ix)
{
    extern /* Subroutine */ int xermsg_(char *, char *, char *, integer *, 
	    integer *, ftnlen, ftnlen, ftnlen), dprwvr_(integer *, integer *, 
	    integer *, doublereal *, integer *);

/* ***BEGIN PROLOGUE  DPRWPG */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to DSPLP */
/* ***LIBRARY   SLATEC */
/* ***TYPE      DOUBLE PRECISION (PRWPGE-S, DPRWPG-D) */
/* ***AUTHOR  Hanson, R. J., (SNLA) */
/*           Wisniewski, J. A., (SNLA) */
/* ***DESCRIPTION */

/*     DPRWPG LIMITS THE TYPE OF STORAGE TO A SEQUENTIAL SCHEME. */
/*     VIRTUAL MEMORY PAGE READ/WRITE SUBROUTINE. */

/*     DEPENDING ON THE VALUE OF KEY, SUBROUTINE DPRWPG() PERFORMS A PAGE */
/*     READ OR WRITE OF PAGE IPAGE. THE PAGE HAS LENGTH LPG. */

/*     KEY       IS A FLAG INDICATING WHETHER A PAGE READ OR WRITE IS */
/*               TO BE PERFORMED. */
/*               IF KEY = 1 DATA IS READ. */
/*               IF KEY = 2 DATA IS WRITTEN. */
/*     IPAGE     IS THE PAGE NUMBER OF THE MATRIX TO BE ACCESSED. */
/*     LPG       IS THE LENGTH OF THE PAGE OF THE MATRIX TO BE ACCESSED. */
/*   SX(*),IX(*) IS THE MATRIX TO BE ACCESSED. */

/*     THIS SUBROUTINE IS A MODIFICATION OF THE SUBROUTINE LRWPGE, */
/*     SANDIA LABS. REPT. SAND78-0785. */
/*     MODIFICATIONS BY K.L. HIEBERT AND R.J. HANSON */
/*     REVISED 811130-1000 */
/*     REVISED YYMMDD-HHMM */

/* ***SEE ALSO  DSPLP */
/* ***ROUTINES CALLED  DPRWVR, XERMSG */
/* ***REVISION HISTORY  (YYMMDD) */
/*   811215  DATE WRITTEN */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ) */
/*   900328  Added TYPE section.  (WRB) */
/*   900510  Fixed error messages and replaced GOTOs with */
/*           IF-THEN-ELSE.  (RWC) */
/*   910403  Updated AUTHOR and DESCRIPTION sections.  (WRB) */
/* ***END PROLOGUE  DPRWPG */
/* ***FIRST EXECUTABLE STATEMENT  DPRWPG */

/*     CHECK IF IPAGE IS IN RANGE. */

    /* Parameter adjustments */
    --ix;
    --sx;

    /* Function Body */
    if (*ipage < 1) {
	xermsg_("SLATEC", "DPRWPG", "THE VALUE OF IPAGE (PAGE NUMBER) WAS NO"
		"T IN THE RANGE1.LE.IPAGE.LE.MAXPGE.", &c__55, &c__1, (ftnlen)
		6, (ftnlen)6, (ftnlen)74);
    }

/*     CHECK IF LPG IS POSITIVE. */

    if (*lpg <= 0) {
	xermsg_("SLATEC", "DPRWPG", "THE VALUE OF LPG (PAGE LENGTH) WAS NONP"
		"OSITIVE.", &c__55, &c__1, (ftnlen)6, (ftnlen)6, (ftnlen)47);
    }

/*     DECIDE IF WE ARE READING OR WRITING. */

    if (*key == 1) {

/*        CODE TO DO A PAGE READ. */

	dprwvr_(key, ipage, lpg, &sx[1], &ix[1]);
    } else if (*key == 2) {

/*        CODE TO DO A PAGE WRITE. */

	dprwvr_(key, ipage, lpg, &sx[1], &ix[1]);
    } else {
	xermsg_("SLATEC", "DPRWPG", "THE VALUE OF KEY (READ-WRITE FLAG) WAS "
		"NOT 1 OR 2.", &c__55, &c__1, (ftnlen)6, (ftnlen)6, (ftnlen)50)
		;
    }
    return 0;
} /* dprwpg_ */

