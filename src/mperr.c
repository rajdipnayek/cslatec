/* mperr.f -- translated by f2c (version 12.02.01).
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
    integer b, t, m, lun, mxr, r__[30];
} mpcom_;

#define mpcom_1 mpcom_

/* Table of constant values */

static integer c__1 = 1;
static integer c__2 = 2;

/* DECK MPERR */
/* Subroutine */ int mperr_(void)
{
    /* Local variables */
    extern /* Subroutine */ int xermsg_(char *, char *, char *, integer *, 
	    integer *, ftnlen, ftnlen, ftnlen);

/* ***BEGIN PROLOGUE  MPERR */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to DQDOTA and DQDOTI */
/* ***LIBRARY   SLATEC */
/* ***TYPE      ALL (MPERR-A) */
/* ***AUTHOR  (UNKNOWN) */
/* ***DESCRIPTION */

/*  This routine is called when a fatal error condition is */
/*  encountered, and after a message has been written on */
/*  logical unit LUN. */

/* ***SEE ALSO  DQDOTA, DQDOTI, MPBLAS */
/* ***ROUTINES CALLED  (NONE) */
/* ***COMMON BLOCKS    MPCOM */
/* ***REVISION HISTORY  (YYMMDD) */
/*   791001  DATE WRITTEN */
/*   ??????  Modified for use with BLAS.  Blank COMMON changed to named */
/*           COMMON.  R given dimension 12. */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900402  Added TYPE section.  (WRB) */
/*   930124  Increased Array size in MPCON for SUN -r8.  (RWC) */
/* ***END PROLOGUE  MPERR */
/* ***FIRST EXECUTABLE STATEMENT  MPERR */
    xermsg_("SLATEC", "MPERR", " *** EXECUTION TERMINATED BY CALL TO MPERR I"
	    "N MP VERSION 770217 ***", &c__1, &c__2, (ftnlen)6, (ftnlen)5, (
	    ftnlen)67);

/* AT PRESENT JUST STOP, BUT COULD DUMP B, T, ETC. HERE. */
/* ACTION COULD EASILY BE CONTROLLED BY A FLAG IN LABELLED COMMON. */
/* ANSI VERSION USES STOP, UNIVAC 1108 VERSION USES */
/* RETURN 0 IN ORDER TO GIVE A TRACE-BACK. */
/* FOR DEBUGGING PURPOSES IT MAY BE USEFUL SIMPLY TO */
/* RETURN HERE.  MOST MP ROUTINES RETURN WITH RESULT */
/* ZERO AFTER CALLING MPERR. */
    s_stop("", (ftnlen)0);
    return 0;
} /* mperr_ */

