/* dfspvn.f -- translated by f2c (version 12.02.01).
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

/* DECK DFSPVN */
/* Subroutine */ int dfspvn_(doublereal *t, integer *jhigh, integer *index, 
	doublereal *x, integer *ileft, doublereal *vnikx)
{
    /* Initialized data */

    static integer j = 1;
    static doublereal deltam[20] = { 0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
	    0.,0.,0.,0.,0.,0.,0. };
    static doublereal deltap[20] = { 0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
	    0.,0.,0.,0.,0.,0.,0. };

    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer l;
    static doublereal vm;
    static integer jp1, ipj, imjp1, jp1ml;
    static doublereal vmprev;

/* ***BEGIN PROLOGUE  DFSPVN */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to DFC */
/* ***LIBRARY   SLATEC */
/* ***TYPE      DOUBLE PRECISION (BSPLVN-S, DFSPVN-D) */
/* ***AUTHOR  (UNKNOWN) */
/* ***DESCRIPTION */

/*  **** Double Precision version of BSPLVN **** */

/* Calculates the value of all possibly nonzero B-splines at *X* of */
/*  order MAX(JHIGH,(J+1)(INDEX-1)) on *T*. */

/* ***SEE ALSO  DFC */
/* ***ROUTINES CALLED  (NONE) */
/* ***REVISION HISTORY  (YYMMDD) */
/*   780801  DATE WRITTEN */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900328  Added TYPE section.  (WRB) */
/* ***END PROLOGUE  DFSPVN */
    /* Parameter adjustments */
    --vnikx;
    --t;

    /* Function Body */
/* ***FIRST EXECUTABLE STATEMENT  DFSPVN */
    switch (*index) {
	case 1:  goto L10;
	case 2:  goto L20;
    }
L10:
    j = 1;
    vnikx[1] = 1.;
    if (j >= *jhigh) {
	goto L99;
    }

L20:
    ipj = *ileft + j;
    deltap[j - 1] = t[ipj] - *x;
    imjp1 = *ileft - j + 1;
    deltam[j - 1] = *x - t[imjp1];
    vmprev = 0.;
    jp1 = j + 1;
    i__1 = j;
    for (l = 1; l <= i__1; ++l) {
	jp1ml = jp1 - l;
	vm = vnikx[l] / (deltap[l - 1] + deltam[jp1ml - 1]);
	vnikx[l] = vm * deltap[l - 1] + vmprev;
/* L26: */
	vmprev = vm * deltam[jp1ml - 1];
    }
    vnikx[jp1] = vmprev;
    j = jp1;
    if (j < *jhigh) {
	goto L20;
    }

L99:
    return 0;
} /* dfspvn_ */

