/* chksn4.f -- translated by f2c (version 12.02.01).
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
    integer kswx, kswy, k, l;
    real ait, bit, cit, dit;
    integer mit, nit, is, ms, js, ns;
    real dlx, dly, tdlx3, tdly3, dlx4, dly4;
} spl4_;

#define spl4_1 spl4_

/* DECK CHKSN4 */
/* Subroutine */ int chksn4_(integer *mbdcnd, integer *nbdcnd, real *alpha, 
	real *beta, S_fp cofx, logical *singlr)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__;
    static real ai, bi, ci, xi;

/* ***BEGIN PROLOGUE  CHKSN4 */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to SEPX4 */
/* ***LIBRARY   SLATEC */
/* ***TYPE      SINGLE PRECISION (CHKSN4-S) */
/* ***AUTHOR  (UNKNOWN) */
/* ***DESCRIPTION */

/*     This subroutine checks if the PDE SEPX4 */
/*     must solve is a singular operator. */

/* ***SEE ALSO  SEPX4 */
/* ***ROUTINES CALLED  (NONE) */
/* ***COMMON BLOCKS    SPL4 */
/* ***REVISION HISTORY  (YYMMDD) */
/*   801001  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900402  Added TYPE section.  (WRB) */
/* ***END PROLOGUE  CHKSN4 */

/* ***FIRST EXECUTABLE STATEMENT  CHKSN4 */
    *singlr = FALSE_;

/*     CHECK IF THE BOUNDARY CONDITIONS ARE */
/*     ENTIRELY PERIODIC AND/OR MIXED */

    if (*mbdcnd != 0 && *mbdcnd != 3 || *nbdcnd != 0 && *nbdcnd != 3) {
	return 0;
    }

/*     CHECK THAT MIXED CONDITIONS ARE PURE NEUMAN */

    if (*mbdcnd != 3) {
	goto L10;
    }
    if (*alpha != 0.f || *beta != 0.f) {
	return 0;
    }
L10:

/*     CHECK THAT NON-DERIVATIVE COEFFICIENT FUNCTIONS */
/*     ARE ZERO */

    i__1 = spl4_1.ms;
    for (i__ = spl4_1.is; i__ <= i__1; ++i__) {
	xi = spl4_1.ait + (i__ - 1) * spl4_1.dlx;
	(*cofx)(&xi, &ai, &bi, &ci);
	if (ci != 0.f) {
	    return 0;
	}
/* L30: */
    }

/*     THE OPERATOR MUST BE SINGULAR IF THIS POINT IS REACHED */

    *singlr = TRUE_;
    return 0;
} /* chksn4_ */

