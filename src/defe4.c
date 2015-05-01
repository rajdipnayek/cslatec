/* defe4.f -- translated by f2c (version 12.02.01).
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

/* DECK DEFE4 */
/* Subroutine */ int defe4_(S_fp cofx, integer *idmn, real *usol, real *grhs)
{
    /* System generated locals */
    integer grhs_dim1, grhs_offset, usol_dim1, usol_offset, i__1, i__2;
    real r__1, r__2, r__3;

    /* Local variables */
    static integer i__, j;
    static real ai, bi, ci, xi, tx, ty;
    extern /* Subroutine */ int dx4_(real *, integer *, integer *, integer *, 
	    real *, real *), dy4_(real *, integer *, integer *, integer *, 
	    real *, real *);
    static real uxxx, uyyy, uxxxx, uyyyy;

/* ***BEGIN PROLOGUE  DEFE4 */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to SEPX4 */
/* ***LIBRARY   SLATEC */
/* ***TYPE      SINGLE PRECISION (DEFE4-S) */
/* ***AUTHOR  (UNKNOWN) */
/* ***DESCRIPTION */

/*     This subroutine first approximates the truncation error given by */
/*     TRUN1(X,Y)=DLX**2*TX+DLY**2*TY where */
/*     TX=AFUN(X)*UXXXX/12.0+BFUN(X)*UXXX/6.0 on the interior and */
/*     at the boundaries if periodic (here UXXX,UXXXX are the third */
/*     and fourth partial derivatives of U with respect to X). */
/*     TX is of the form AFUN(X)/3.0*(UXXXX/4.0+UXXX/DLX) */
/*     at X=A or X=B if the boundary condition there is mixed. */
/*     TX=0.0 along specified boundaries.  TY has symmetric form */
/*     in Y with X,AFUN(X),BFUN(X) replaced by Y,DFUN(Y),EFUN(Y). */
/*     The second order solution in USOL is used to approximate */
/*     (via second order finite differencing) the truncation error */
/*     and the result is added to the right hand side in GRHS */
/*     and then transferred to USOL to be used as a new right */
/*     hand side when calling BLKTRI for a fourth order solution. */

/* ***SEE ALSO  SEPX4 */
/* ***ROUTINES CALLED  DX4, DY4 */
/* ***COMMON BLOCKS    SPL4 */
/* ***REVISION HISTORY  (YYMMDD) */
/*   801001  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900402  Added TYPE section.  (WRB) */
/* ***END PROLOGUE  DEFE4 */

/* ***FIRST EXECUTABLE STATEMENT  DEFE4 */
    /* Parameter adjustments */
    grhs_dim1 = *idmn;
    grhs_offset = 1 + grhs_dim1;
    grhs -= grhs_offset;
    usol_dim1 = *idmn;
    usol_offset = 1 + usol_dim1;
    usol -= usol_offset;

    /* Function Body */
    i__1 = spl4_1.ms;
    for (i__ = spl4_1.is; i__ <= i__1; ++i__) {
	xi = spl4_1.ait + (i__ - 1) * spl4_1.dlx;
	(*cofx)(&xi, &ai, &bi, &ci);
	i__2 = spl4_1.ns;
	for (j = spl4_1.js; j <= i__2; ++j) {

/*     COMPUTE PARTIAL DERIVATIVE APPROXIMATIONS AT (XI,YJ) */

	    dx4_(&usol[usol_offset], idmn, &i__, &j, &uxxx, &uxxxx);
	    dy4_(&usol[usol_offset], idmn, &i__, &j, &uyyy, &uyyyy);
	    tx = ai * uxxxx / 12.f + bi * uxxx / 6.f;
	    ty = uyyyy / 12.f;

/*     RESET FORM OF TRUNCATION IF AT BOUNDARY WHICH IS NON-PERIODIC */

	    if (spl4_1.kswx == 1 || i__ > 1 && i__ < spl4_1.k) {
		goto L10;
	    }
	    tx = ai / 3.f * (uxxxx / 4.f + uxxx / spl4_1.dlx);
L10:
	    if (spl4_1.kswy == 1 || j > 1 && j < spl4_1.l) {
		goto L20;
	    }
	    ty = (uyyyy / 4.f + uyyy / spl4_1.dly) / 3.f;
L20:
/* Computing 2nd power */
	    r__1 = spl4_1.dly;
/* Computing 2nd power */
	    r__2 = spl4_1.dlx;
/* Computing 2nd power */
	    r__3 = spl4_1.dly;
	    grhs[i__ + j * grhs_dim1] += r__1 * r__1 * (r__2 * r__2 * tx + 
		    r__3 * r__3 * ty);
/* L30: */
	}
    }

/*     RESET THE RIGHT HAND SIDE IN USOL */

    i__2 = spl4_1.ms;
    for (i__ = spl4_1.is; i__ <= i__2; ++i__) {
	i__1 = spl4_1.ns;
	for (j = spl4_1.js; j <= i__1; ++j) {
	    usol[i__ + j * usol_dim1] = grhs[i__ + j * grhs_dim1];
/* L50: */
	}
/* L60: */
    }
    return 0;
} /* defe4_ */

