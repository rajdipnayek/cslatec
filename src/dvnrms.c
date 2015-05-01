/* dvnrms.f -- translated by f2c (version 12.02.01).
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

/* DECK DVNRMS */
doublereal dvnrms_(integer *n, doublereal *v, doublereal *w)
{
    /* System generated locals */
    integer i__1;
    doublereal ret_val, d__1;

    /* Local variables */
    static integer i__;
    static doublereal sum;

/* ***BEGIN PROLOGUE  DVNRMS */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to DDEBDF */
/* ***LIBRARY   SLATEC */
/* ***TYPE      DOUBLE PRECISION (VNWRMS-S, DVNRMS-D) */
/* ***AUTHOR  (UNKNOWN) */
/* ***DESCRIPTION */

/*   DVNRMS computes a weighted root-mean-square vector norm for the */
/*   integrator package DDEBDF. */

/* ***SEE ALSO  DDEBDF */
/* ***ROUTINES CALLED  (NONE) */
/* ***REVISION HISTORY  (YYMMDD) */
/*   820301  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890831  Modified array declarations.  (WRB) */
/*   890911  Removed unnecessary intrinsics.  (WRB) */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900328  Added TYPE section.  (WRB) */
/* ***END PROLOGUE  DVNRMS */
/* ***FIRST EXECUTABLE STATEMENT  DVNRMS */
    /* Parameter adjustments */
    --w;
    --v;

    /* Function Body */
    sum = 0.;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing 2nd power */
	d__1 = v[i__] / w[i__];
	sum += d__1 * d__1;
/* L10: */
    }
    ret_val = sqrt(sum / *n);
    return ret_val;
/*     ----------------------- END OF FUNCTION DVNRMS */
/*     ------------------------ */
} /* dvnrms_ */

