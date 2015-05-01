/* orthog.f -- translated by f2c (version 12.02.01).
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

/* DECK ORTHOG */
/* Subroutine */ int orthog_(real *usol, integer *idmn, real *zn, real *zm, 
	real *pertrb)
{
    /* System generated locals */
    integer usol_dim1, usol_offset, i__1, i__2;

    /* Local variables */
    static integer i__, j, ii, jj;
    static real ete, ute;
    static integer ifnl, jfnl, istr, jstr;

/* ***BEGIN PROLOGUE  ORTHOG */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to SEPELI */
/* ***LIBRARY   SLATEC */
/* ***TYPE      SINGLE PRECISION (ORTHOG-S) */
/* ***AUTHOR  (UNKNOWN) */
/* ***DESCRIPTION */

/*     This subroutine orthogonalizes the array USOL with respect to */
/*     the constant array in a weighted least squares norm. */

/* ***SEE ALSO  SEPELI */
/* ***ROUTINES CALLED  (NONE) */
/* ***COMMON BLOCKS    SPLPCM */
/* ***REVISION HISTORY  (YYMMDD) */
/*   801001  DATE WRITTEN */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900402  Added TYPE section.  (WRB) */
/* ***END PROLOGUE  ORTHOG */

/* ***FIRST EXECUTABLE STATEMENT  ORTHOG */
    /* Parameter adjustments */
    usol_dim1 = *idmn;
    usol_offset = 1 + usol_dim1;
    usol -= usol_offset;
    --zn;
    --zm;

    /* Function Body */
    istr = splpcm_1.is;
    ifnl = splpcm_1.ms;
    jstr = splpcm_1.js;
    jfnl = splpcm_1.ns;

/*     COMPUTE WEIGHTED INNER PRODUCTS */

    ute = 0.f;
    ete = 0.f;
    i__1 = splpcm_1.ms;
    for (i__ = splpcm_1.is; i__ <= i__1; ++i__) {
	ii = i__ - splpcm_1.is + 1;
	i__2 = splpcm_1.ns;
	for (j = splpcm_1.js; j <= i__2; ++j) {
	    jj = j - splpcm_1.js + 1;
	    ete += zm[ii] * zn[jj];
	    ute += usol[i__ + j * usol_dim1] * zm[ii] * zn[jj];
/* L10: */
	}
/* L20: */
    }

/*     SET PERTURBATION PARAMETER */

    *pertrb = ute / ete;

/*     SUBTRACT OFF CONSTANT PERTRB */

    i__1 = ifnl;
    for (i__ = istr; i__ <= i__1; ++i__) {
	i__2 = jfnl;
	for (j = jstr; j <= i__2; ++j) {
	    usol[i__ + j * usol_dim1] -= *pertrb;
/* L30: */
	}
/* L40: */
    }
    return 0;
} /* orthog_ */

