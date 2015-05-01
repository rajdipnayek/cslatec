/* prvec.f -- translated by f2c (version 12.02.01).
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

/* DECK PRVEC */
doublereal prvec_(integer *m, real *u, real *v)
{
    /* System generated locals */
    real ret_val;

    /* Local variables */
    static integer n, np;
    static real vp;
    extern doublereal sdot_(integer *, real *, integer *, real *, integer *);

/* ***BEGIN PROLOGUE  PRVEC */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to BVSUP */
/* ***LIBRARY   SLATEC */
/* ***TYPE      SINGLE PRECISION (PRVEC-S, DPRVEC-D) */
/* ***AUTHOR  Watts, H. A., (SNLA) */
/* ***DESCRIPTION */

/*  This subroutine computes the inner product of a vector U */
/*  with the imaginary product or mate vector corresponding to V */

/* ***SEE ALSO  BVSUP */
/* ***ROUTINES CALLED  SDOT */
/* ***REVISION HISTORY  (YYMMDD) */
/*   750601  DATE WRITTEN */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900328  Added TYPE section.  (WRB) */
/*   910722  Updated AUTHOR section.  (ALS) */
/* ***END PROLOGUE  PRVEC */

/* ***FIRST EXECUTABLE STATEMENT  PRVEC */
    /* Parameter adjustments */
    --v;
    --u;

    /* Function Body */
    n = *m / 2;
    np = n + 1;
    vp = sdot_(&n, &u[1], &c__1, &v[np], &c__1);
    ret_val = sdot_(&n, &u[np], &c__1, &v[1], &c__1) - vp;
    return ret_val;
} /* prvec_ */

