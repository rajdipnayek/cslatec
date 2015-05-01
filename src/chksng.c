/* chksng.f -- translated by f2c (version 12.02.01).
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
} splpcm_;

#define splpcm_1 splpcm_

/* DECK CHKSNG */
/* Subroutine */ int chksng_(integer *mbdcnd, integer *nbdcnd, real *alpha, 
	real *beta, real *gama, real *xnu, S_fp cofx, S_fp cofy, logical *
	singlr)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__, j;
    static real ai, bi, ci, dj, ej, fj, xi, yj;

/* ***BEGIN PROLOGUE  CHKSNG */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to SEPELI */
/* ***LIBRARY   SLATEC */
/* ***TYPE      SINGLE PRECISION (CHKSNG-S) */
/* ***AUTHOR  (UNKNOWN) */
/* ***DESCRIPTION */

/*     This subroutine checks if the PDE SEPELI */
/*     must solve is a singular operator. */

/* ***SEE ALSO  SEPELI */
/* ***ROUTINES CALLED  (NONE) */
/* ***COMMON BLOCKS    SPLPCM */
/* ***REVISION HISTORY  (YYMMDD) */
/*   801001  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900402  Added TYPE section.  (WRB) */
/* ***END PROLOGUE  CHKSNG */

/* ***FIRST EXECUTABLE STATEMENT  CHKSNG */
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
    if (*nbdcnd != 3) {
	goto L20;
    }
    if (*gama != 0.f || *xnu != 0.f) {
	return 0;
    }
L20:

/*     CHECK THAT NON-DERIVATIVE COEFFICIENT FUNCTIONS */
/*     ARE ZERO */

    i__1 = splpcm_1.ms;
    for (i__ = splpcm_1.is; i__ <= i__1; ++i__) {
	xi = splpcm_1.ait + (i__ - 1) * splpcm_1.dlx;
	(*cofx)(&xi, &ai, &bi, &ci);
	if (ci != 0.f) {
	    return 0;
	}
/* L30: */
    }
    i__1 = splpcm_1.ns;
    for (j = splpcm_1.js; j <= i__1; ++j) {
	yj = splpcm_1.cit + (j - 1) * splpcm_1.dly;
	(*cofy)(&yj, &dj, &ej, &fj);
	if (fj != 0.f) {
	    return 0;
	}
/* L40: */
    }

/*     THE OPERATOR MUST BE SINGULAR IF THIS POINT IS REACHED */

    *singlr = TRUE_;
    return 0;
} /* chksng_ */

