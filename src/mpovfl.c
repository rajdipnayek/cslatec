/* mpovfl.f -- translated by f2c (version 12.02.01).
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
static integer c__4 = 4;

/* DECK MPOVFL */
/* Subroutine */ int mpovfl_(integer *x)
{
    /* Format strings */
    static char fmt_10[] = "(\002 *** CALL TO MPOVFL, MP OVERFLOW OCCURRED *"
	    "**\002)";

    /* Local variables */
    extern /* Subroutine */ int mpchk_(integer *, integer *), mperr_(void), 
	    mpmaxr_(integer *);

    /* Fortran I/O blocks */
    static cilist io___1 = { 0, 0, 0, fmt_10, 0 };


/* ***BEGIN PROLOGUE  MPOVFL */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to DQDOTA and DQDOTI */
/* ***LIBRARY   SLATEC */
/* ***TYPE      ALL (MPOVFL-A) */
/* ***AUTHOR  (UNKNOWN) */
/* ***DESCRIPTION */

/*  Called on multiple-precision overflow, i.e. when the */
/*  exponent of 'mp' number X would exceed M.  At present execution is */
/*  terminated with an error message after calling MPMAXR(X), but it */
/*  would be possible to return, possibly updating a counter and */
/*  terminating execution after a preset number of overflows.  Action */
/*  could easily be determined by a flag in labelled common. */

/*  The argument X(*) is an INTEGER array of size 30.  See the comments */
/*  in the routine MPBLAS for the reason for this choice. */

/* ***SEE ALSO  DQDOTA, DQDOTI, MPBLAS */
/* ***ROUTINES CALLED  MPCHK, MPERR, MPMAXR */
/* ***COMMON BLOCKS    MPCOM */
/* ***REVISION HISTORY  (YYMMDD) */
/*   791001  DATE WRITTEN */
/*   ??????  Modified for use with BLAS.  Blank COMMON changed to named */
/*           COMMON.  R given dimension 12. */
/*   890831  Modified array declarations.  (WRB) */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900402  Added TYPE section.  (WRB) */
/*   930124  Increased Array size in MPCON for SUN -r8.  (RWC) */
/* ***END PROLOGUE  MPOVFL */
/* ***FIRST EXECUTABLE STATEMENT  MPOVFL */
    /* Parameter adjustments */
    --x;

    /* Function Body */
    mpchk_(&c__1, &c__4);
/* SET X TO LARGEST POSSIBLE POSITIVE NUMBER */
    mpmaxr_(&x[1]);
    io___1.ciunit = mpcom_1.lun;
    s_wsfe(&io___1);
    e_wsfe();
/* TERMINATE EXECUTION BY CALLING MPERR */
    mperr_();
    return 0;
} /* mpovfl_ */

