/* minso4.f -- translated by f2c (version 12.02.01).
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

/* DECK MINSO4 */
/* Subroutine */ int minso4_(real *usol, integer *idmn, real *zn, real *zm, 
	real *pertb)
{
    /* System generated locals */
    integer usol_dim1, usol_offset, i__1, i__2;

    /* Local variables */
    static integer i__, j, ii, jj;
    static real ete, ute;
    static integer ifnl, jfnl, istr, jstr;
    static real pertrb;

/* ***BEGIN PROLOGUE  MINSO4 */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to SEPX4 */
/* ***LIBRARY   SLATEC */
/* ***TYPE      SINGLE PRECISION (MINSO4-S) */
/* ***AUTHOR  (UNKNOWN) */
/* ***DESCRIPTION */

/*     This subroutine orthogonalizes the array USOL with respect to */
/*     the constant array in a weighted least squares norm. */

/*     Entry at MINSO4 occurs when the final solution is */
/*     to be minimized with respect to the weighted */
/*     least squares norm. */

/* ***SEE ALSO  SEPX4 */
/* ***ROUTINES CALLED  (NONE) */
/* ***COMMON BLOCKS    SPL4 */
/* ***REVISION HISTORY  (YYMMDD) */
/*   801001  DATE WRITTEN */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900402  Added TYPE section.  (WRB) */
/* ***END PROLOGUE  MINSO4 */

/* ***FIRST EXECUTABLE STATEMENT  MINSO4 */
    /* Parameter adjustments */
    usol_dim1 = *idmn;
    usol_offset = 1 + usol_dim1;
    usol -= usol_offset;
    --zn;
    --zm;

    /* Function Body */
    istr = 1;
    ifnl = spl4_1.k;
    jstr = 1;
    jfnl = spl4_1.l;

/*     COMPUTE WEIGHTED INNER PRODUCTS */

    ute = 0.f;
    ete = 0.f;
    i__1 = spl4_1.ms;
    for (i__ = spl4_1.is; i__ <= i__1; ++i__) {
	ii = i__ - spl4_1.is + 1;
	i__2 = spl4_1.ns;
	for (j = spl4_1.js; j <= i__2; ++j) {
	    jj = j - spl4_1.js + 1;
	    ete += zm[ii] * zn[jj];
	    ute += usol[i__ + j * usol_dim1] * zm[ii] * zn[jj];
/* L10: */
	}
/* L20: */
    }

/*     SET PERTURBATION PARAMETER */

    pertrb = ute / ete;

/*     SUBTRACT OFF CONSTANT PERTRB */

    i__1 = ifnl;
    for (i__ = istr; i__ <= i__1; ++i__) {
	i__2 = jfnl;
	for (j = jstr; j <= i__2; ++j) {
	    usol[i__ + j * usol_dim1] -= pertrb;
/* L30: */
	}
/* L40: */
    }
    return 0;
} /* minso4_ */

