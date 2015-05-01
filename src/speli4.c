/* speli4.f -- translated by f2c (version 12.02.01).
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

/* DECK SPELI4 */
/* Subroutine */ int speli4_(integer *iorder, real *a, real *b, integer *m, 
	integer *mbdcnd, real *bda, real *alpha, real *bdb, real *beta, real *
	c__, real *d__, integer *n, integer *nbdcnd, real *bdc, real *bdd, 
	S_fp cofx, real *an, real *bn, real *cn, real *dn, real *un, real *zn,
	 real *am, real *bm, real *cm, real *dm, real *um, real *zm, real *
	grhs, real *usol, integer *idmn, real *w, real *pertrb, integer *
	ierror)
{
    /* System generated locals */
    integer grhs_dim1, grhs_offset, usol_dim1, usol_offset, i__1, i__2;
    real r__1;

    /* Local variables */
    static integer i__, j;
    static real ai, bi, ci;
    static integer mp, np;
    static real xi, ax1, dy1, axi, bxi, cxi, dyj, eyj, fyj, cxm, fyn, xnu, 
	    gama;
    static integer iord;
    extern /* Subroutine */ int defe4_(S_fp, integer *, real *, real *), 
	    tris4_(integer *, real *, real *, real *, real *, real *, real *);
    static integer ieror;
    static real prtrb;
    extern /* Subroutine */ int chksn4_(integer *, integer *, real *, real *, 
	    S_fp, logical *), minso4_(real *, integer *, real *, real *, real 
	    *), ortho4_(real *, integer *, real *, real *, real *), genbun_(
	    integer *, integer *, integer *, integer *, real *, real *, real *
	    , integer *, real *, integer *, real *);
    static logical singlr;

/* ***BEGIN PROLOGUE  SPELI4 */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to SEPX4 */
/* ***LIBRARY   SLATEC */
/* ***TYPE      SINGLE PRECISION (SPELI4-S) */
/* ***AUTHOR  (UNKNOWN) */
/* ***DESCRIPTION */

/*     SPELI4 sets up vectors and arrays for input to BLKTRI */
/*     and computes a second order solution in USOL.  A return jump to */
/*     SEPX4 occurs if IORDER=2.  If IORDER=4 a fourth order */
/*     solution is generated in USOL. */

/* ***SEE ALSO  SEPX4 */
/* ***ROUTINES CALLED  CHKSN4, DEFE4, GENBUN, MINSO4, ORTHO4, TRIS4 */
/* ***COMMON BLOCKS    SPL4 */
/* ***REVISION HISTORY  (YYMMDD) */
/*   801001  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   891009  Removed unreferenced variable.  (WRB) */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900402  Added TYPE section.  (WRB) */
/* ***END PROLOGUE  SPELI4 */

/* ***FIRST EXECUTABLE STATEMENT  SPELI4 */
    /* Parameter adjustments */
    --bda;
    --bdb;
    --bdc;
    --bdd;
    --an;
    --bn;
    --cn;
    --dn;
    --un;
    --zn;
    --am;
    --bm;
    --cm;
    --dm;
    --um;
    --zm;
    usol_dim1 = *idmn;
    usol_offset = 1 + usol_dim1;
    usol -= usol_offset;
    grhs_dim1 = *idmn;
    grhs_offset = 1 + grhs_dim1;
    grhs -= grhs_offset;
    --w;

    /* Function Body */
    spl4_1.kswx = *mbdcnd + 1;
    spl4_1.kswy = *nbdcnd + 1;
    spl4_1.k = *m + 1;
    spl4_1.l = *n + 1;
    spl4_1.ait = *a;
    spl4_1.bit = *b;
    spl4_1.cit = *c__;
    spl4_1.dit = *d__;
    spl4_1.dly = (spl4_1.dit - spl4_1.cit) / *n;

/*     SET RIGHT HAND SIDE VALUES FROM GRHS IN USOL ON THE INTERIOR */
/*     AND NON-SPECIFIED BOUNDARIES. */

    i__1 = *m;
    for (i__ = 2; i__ <= i__1; ++i__) {
	i__2 = *n;
	for (j = 2; j <= i__2; ++j) {
/* Computing 2nd power */
	    r__1 = spl4_1.dly;
	    usol[i__ + j * usol_dim1] = r__1 * r__1 * grhs[i__ + j * 
		    grhs_dim1];
/* L10: */
	}
/* L20: */
    }
    if (spl4_1.kswx == 2 || spl4_1.kswx == 3) {
	goto L40;
    }
    i__1 = *n;
    for (j = 2; j <= i__1; ++j) {
/* Computing 2nd power */
	r__1 = spl4_1.dly;
	usol[j * usol_dim1 + 1] = r__1 * r__1 * grhs[j * grhs_dim1 + 1];
/* L30: */
    }
L40:
    if (spl4_1.kswx == 2 || spl4_1.kswx == 5) {
	goto L60;
    }
    i__1 = *n;
    for (j = 2; j <= i__1; ++j) {
/* Computing 2nd power */
	r__1 = spl4_1.dly;
	usol[spl4_1.k + j * usol_dim1] = r__1 * r__1 * grhs[spl4_1.k + j * 
		grhs_dim1];
/* L50: */
    }
L60:
    if (spl4_1.kswy == 2 || spl4_1.kswy == 3) {
	goto L80;
    }
    i__1 = *m;
    for (i__ = 2; i__ <= i__1; ++i__) {
/* Computing 2nd power */
	r__1 = spl4_1.dly;
	usol[i__ + usol_dim1] = r__1 * r__1 * grhs[i__ + grhs_dim1];
/* L70: */
    }
L80:
    if (spl4_1.kswy == 2 || spl4_1.kswy == 5) {
	goto L100;
    }
    i__1 = *m;
    for (i__ = 2; i__ <= i__1; ++i__) {
/* Computing 2nd power */
	r__1 = spl4_1.dly;
	usol[i__ + spl4_1.l * usol_dim1] = r__1 * r__1 * grhs[i__ + spl4_1.l *
		 grhs_dim1];
/* L90: */
    }
L100:
    if (spl4_1.kswx != 2 && spl4_1.kswx != 3 && spl4_1.kswy != 2 && 
	    spl4_1.kswy != 3) {
/* Computing 2nd power */
	r__1 = spl4_1.dly;
	usol[usol_dim1 + 1] = r__1 * r__1 * grhs[grhs_dim1 + 1];
    }
    if (spl4_1.kswx != 2 && spl4_1.kswx != 5 && spl4_1.kswy != 2 && 
	    spl4_1.kswy != 3) {
/* Computing 2nd power */
	r__1 = spl4_1.dly;
	usol[spl4_1.k + usol_dim1] = r__1 * r__1 * grhs[spl4_1.k + grhs_dim1];
    }
    if (spl4_1.kswx != 2 && spl4_1.kswx != 3 && spl4_1.kswy != 2 && 
	    spl4_1.kswy != 5) {
/* Computing 2nd power */
	r__1 = spl4_1.dly;
	usol[spl4_1.l * usol_dim1 + 1] = r__1 * r__1 * grhs[spl4_1.l * 
		grhs_dim1 + 1];
    }
    if (spl4_1.kswx != 2 && spl4_1.kswx != 5 && spl4_1.kswy != 2 && 
	    spl4_1.kswy != 5) {
/* Computing 2nd power */
	r__1 = spl4_1.dly;
	usol[spl4_1.k + spl4_1.l * usol_dim1] = r__1 * r__1 * grhs[spl4_1.k + 
		spl4_1.l * grhs_dim1];
    }

/*     SET SWITCHES FOR PERIODIC OR NON-PERIODIC BOUNDARIES */

    mp = 1;
    if (spl4_1.kswx == 1) {
	mp = 0;
    }
    np = *nbdcnd;

/*     SET DLX,DLY AND SIZE OF BLOCK TRI-DIAGONAL SYSTEM GENERATED */
/*     IN NINT,MINT */

    spl4_1.dlx = (spl4_1.bit - spl4_1.ait) / *m;
    spl4_1.mit = spl4_1.k - 1;
    if (spl4_1.kswx == 2) {
	spl4_1.mit = spl4_1.k - 2;
    }
    if (spl4_1.kswx == 4) {
	spl4_1.mit = spl4_1.k;
    }
    spl4_1.dly = (spl4_1.dit - spl4_1.cit) / *n;
    spl4_1.nit = spl4_1.l - 1;
    if (spl4_1.kswy == 2) {
	spl4_1.nit = spl4_1.l - 2;
    }
    if (spl4_1.kswy == 4) {
	spl4_1.nit = spl4_1.l;
    }
/* Computing 3rd power */
    r__1 = spl4_1.dlx;
    spl4_1.tdlx3 = r__1 * (r__1 * r__1) * 2.f;
/* Computing 4th power */
    r__1 = spl4_1.dlx, r__1 *= r__1;
    spl4_1.dlx4 = r__1 * r__1;
/* Computing 3rd power */
    r__1 = spl4_1.dly;
    spl4_1.tdly3 = r__1 * (r__1 * r__1) * 2.f;
/* Computing 4th power */
    r__1 = spl4_1.dly, r__1 *= r__1;
    spl4_1.dly4 = r__1 * r__1;

/*     SET SUBSCRIPT LIMITS FOR PORTION OF ARRAY TO INPUT TO BLKTRI */

    spl4_1.is = 1;
    spl4_1.js = 1;
    if (spl4_1.kswx == 2 || spl4_1.kswx == 3) {
	spl4_1.is = 2;
    }
    if (spl4_1.kswy == 2 || spl4_1.kswy == 3) {
	spl4_1.js = 2;
    }
    spl4_1.ns = spl4_1.nit + spl4_1.js - 1;
    spl4_1.ms = spl4_1.mit + spl4_1.is - 1;

/*     SET X - DIRECTION */

    i__1 = spl4_1.mit;
    for (i__ = 1; i__ <= i__1; ++i__) {
	xi = spl4_1.ait + (spl4_1.is + i__ - 2) * spl4_1.dlx;
	(*cofx)(&xi, &ai, &bi, &ci);
	axi = (ai / spl4_1.dlx - bi * .5f) / spl4_1.dlx;
/* Computing 2nd power */
	r__1 = spl4_1.dlx;
	bxi = ai * -2.f / (r__1 * r__1) + ci;
	cxi = (ai / spl4_1.dlx + bi * .5f) / spl4_1.dlx;
/* Computing 2nd power */
	r__1 = spl4_1.dly;
	am[i__] = r__1 * r__1 * axi;
/* Computing 2nd power */
	r__1 = spl4_1.dly;
	bm[i__] = r__1 * r__1 * bxi;
/* Computing 2nd power */
	r__1 = spl4_1.dly;
	cm[i__] = r__1 * r__1 * cxi;
/* L110: */
    }

/*     SET Y DIRECTION */

    i__1 = spl4_1.nit;
    for (j = 1; j <= i__1; ++j) {
	dyj = 1.f;
	eyj = -2.f;
	fyj = 1.f;
	an[j] = dyj;
	bn[j] = eyj;
	cn[j] = fyj;
/* L120: */
    }

/*     ADJUST EDGES IN X DIRECTION UNLESS PERIODIC */

    ax1 = am[1];
    cxm = cm[spl4_1.mit];
    switch (spl4_1.kswx) {
	case 1:  goto L170;
	case 2:  goto L130;
	case 3:  goto L150;
	case 4:  goto L160;
	case 5:  goto L140;
    }

/*     DIRICHLET-DIRICHLET IN X DIRECTION */

L130:
    am[1] = 0.f;
    cm[spl4_1.mit] = 0.f;
    goto L170;

/*     MIXED-DIRICHLET IN X DIRECTION */

L140:
    am[1] = 0.f;
    bm[1] += *alpha * 2.f * spl4_1.dlx * ax1;
    cm[1] += ax1;
    cm[spl4_1.mit] = 0.f;
    goto L170;

/*     DIRICHLET-MIXED IN X DIRECTION */

L150:
    am[1] = 0.f;
    am[spl4_1.mit] += cxm;
    bm[spl4_1.mit] -= *beta * 2.f * spl4_1.dlx * cxm;
    cm[spl4_1.mit] = 0.f;
    goto L170;

/*     MIXED - MIXED IN X DIRECTION */

L160:
    am[1] = 0.f;
    bm[1] += spl4_1.dlx * 2.f * *alpha * ax1;
    cm[1] += ax1;
    am[spl4_1.mit] += cxm;
    bm[spl4_1.mit] -= spl4_1.dlx * 2.f * *beta * cxm;
    cm[spl4_1.mit] = 0.f;
L170:

/*     ADJUST IN Y DIRECTION UNLESS PERIODIC */

    dy1 = an[1];
    fyn = cn[spl4_1.nit];
    gama = 0.f;
    xnu = 0.f;
    switch (spl4_1.kswy) {
	case 1:  goto L220;
	case 2:  goto L180;
	case 3:  goto L200;
	case 4:  goto L210;
	case 5:  goto L190;
    }

/*     DIRICHLET-DIRICHLET IN Y DIRECTION */

L180:
    an[1] = 0.f;
    cn[spl4_1.nit] = 0.f;
    goto L220;

/*     MIXED-DIRICHLET IN Y DIRECTION */

L190:
    an[1] = 0.f;
    bn[1] += spl4_1.dly * 2.f * gama * dy1;
    cn[1] += dy1;
    cn[spl4_1.nit] = 0.f;
    goto L220;

/*     DIRICHLET-MIXED IN Y DIRECTION */

L200:
    an[1] = 0.f;
    an[spl4_1.nit] += fyn;
    bn[spl4_1.nit] -= spl4_1.dly * 2.f * xnu * fyn;
    cn[spl4_1.nit] = 0.f;
    goto L220;

/*     MIXED - MIXED DIRECTION IN Y DIRECTION */

L210:
    an[1] = 0.f;
    bn[1] += spl4_1.dly * 2.f * gama * dy1;
    cn[1] += dy1;
    an[spl4_1.nit] += fyn;
    bn[spl4_1.nit] -= spl4_1.dly * 2.f * xnu * fyn;
    cn[spl4_1.nit] = 0.f;
L220:
    if (spl4_1.kswx == 1) {
	goto L270;
    }

/*     ADJUST USOL ALONG X EDGE */

    i__1 = spl4_1.ns;
    for (j = spl4_1.js; j <= i__1; ++j) {
	if (spl4_1.kswx != 2 && spl4_1.kswx != 3) {
	    goto L230;
	}
	usol[spl4_1.is + j * usol_dim1] -= ax1 * usol[j * usol_dim1 + 1];
	goto L240;
L230:
	usol[spl4_1.is + j * usol_dim1] += spl4_1.dlx * 2.f * ax1 * bda[j];
L240:
	if (spl4_1.kswx != 2 && spl4_1.kswx != 5) {
	    goto L250;
	}
	usol[spl4_1.ms + j * usol_dim1] -= cxm * usol[spl4_1.k + j * 
		usol_dim1];
	goto L260;
L250:
	usol[spl4_1.ms + j * usol_dim1] -= spl4_1.dlx * 2.f * cxm * bdb[j];
L260:
	;
    }
L270:
    if (spl4_1.kswy == 1) {
	goto L320;
    }

/*     ADJUST USOL ALONG Y EDGE */

    i__1 = spl4_1.ms;
    for (i__ = spl4_1.is; i__ <= i__1; ++i__) {
	if (spl4_1.kswy != 2 && spl4_1.kswy != 3) {
	    goto L280;
	}
	usol[i__ + spl4_1.js * usol_dim1] -= dy1 * usol[i__ + usol_dim1];
	goto L290;
L280:
	usol[i__ + spl4_1.js * usol_dim1] += spl4_1.dly * 2.f * dy1 * bdc[i__]
		;
L290:
	if (spl4_1.kswy != 2 && spl4_1.kswy != 5) {
	    goto L300;
	}
	usol[i__ + spl4_1.ns * usol_dim1] -= fyn * usol[i__ + spl4_1.l * 
		usol_dim1];
	goto L310;
L300:
	usol[i__ + spl4_1.ns * usol_dim1] -= spl4_1.dly * 2.f * fyn * bdd[i__]
		;
L310:
	;
    }
L320:

/*     SAVE ADJUSTED EDGES IN GRHS IF IORDER=4 */

    if (*iorder != 4) {
	goto L350;
    }
    i__1 = spl4_1.ns;
    for (j = spl4_1.js; j <= i__1; ++j) {
	grhs[spl4_1.is + j * grhs_dim1] = usol[spl4_1.is + j * usol_dim1];
	grhs[spl4_1.ms + j * grhs_dim1] = usol[spl4_1.ms + j * usol_dim1];
/* L330: */
    }
    i__1 = spl4_1.ms;
    for (i__ = spl4_1.is; i__ <= i__1; ++i__) {
	grhs[i__ + spl4_1.js * grhs_dim1] = usol[i__ + spl4_1.js * usol_dim1];
	grhs[i__ + spl4_1.ns * grhs_dim1] = usol[i__ + spl4_1.ns * usol_dim1];
/* L340: */
    }
L350:
    iord = *iorder;
    *pertrb = 0.f;

/*     CHECK IF OPERATOR IS SINGULAR */

    chksn4_(mbdcnd, nbdcnd, alpha, beta, (S_fp)cofx, &singlr);

/*     COMPUTE NON-ZERO EIGENVECTOR IN NULL SPACE OF TRANSPOSE */
/*     IF SINGULAR */

    if (singlr) {
	tris4_(&spl4_1.mit, &am[1], &bm[1], &cm[1], &dm[1], &um[1], &zm[1]);
    }
    if (singlr) {
	tris4_(&spl4_1.nit, &an[1], &bn[1], &cn[1], &dn[1], &un[1], &zn[1]);
    }

/*     ADJUST RIGHT HAND SIDE IF NECESSARY */

L360:
    if (singlr) {
	ortho4_(&usol[usol_offset], idmn, &zn[1], &zm[1], pertrb);
    }

/*     COMPUTE SOLUTION */

/*     SAVE ADJUSTED RIGHT HAND SIDE IN GRHS */
    i__1 = spl4_1.ns;
    for (j = spl4_1.js; j <= i__1; ++j) {
	i__2 = spl4_1.ms;
	for (i__ = spl4_1.is; i__ <= i__2; ++i__) {
	    grhs[i__ + j * grhs_dim1] = usol[i__ + j * usol_dim1];
/* L444: */
	}
    }
    genbun_(&np, &spl4_1.nit, &mp, &spl4_1.mit, &am[1], &bm[1], &cm[1], idmn, 
	    &usol[spl4_1.is + spl4_1.js * usol_dim1], &ieror, &w[1]);
/*     CHECK IF ERROR DETECTED IN POIS */
/*     THIS CAN ONLY CORRESPOND TO IERROR=12 */
    if (ieror == 0) {
	goto L224;
    }
/*     SET ERROR FLAG IF IMPROPER COEFFICIENTS INPUT TO POIS */
    *ierror = 12;
    return 0;
L224:
    if (*ierror != 0) {
	return 0;
    }

/*     SET PERIODIC BOUNDARIES IF NECESSARY */

    if (spl4_1.kswx != 1) {
	goto L380;
    }
    i__2 = spl4_1.l;
    for (j = 1; j <= i__2; ++j) {
	usol[spl4_1.k + j * usol_dim1] = usol[j * usol_dim1 + 1];
/* L370: */
    }
L380:
    if (spl4_1.kswy != 1) {
	goto L400;
    }
    i__2 = spl4_1.k;
    for (i__ = 1; i__ <= i__2; ++i__) {
	usol[i__ + spl4_1.l * usol_dim1] = usol[i__ + usol_dim1];
/* L390: */
    }
L400:

/*     MINIMIZE SOLUTION WITH RESPECT TO WEIGHTED LEAST SQUARES */
/*     NORM IF OPERATOR IS SINGULAR */

    if (singlr) {
	minso4_(&usol[usol_offset], idmn, &zn[1], &zm[1], &prtrb);
    }

/*     RETURN IF DEFERRED CORRECTIONS AND A FOURTH ORDER SOLUTION ARE */
/*     NOT FLAGGED */

    if (iord == 2) {
	return 0;
    }
    iord = 2;

/*     COMPUTE NEW RIGHT HAND SIDE FOR FOURTH ORDER SOLUTION */

    defe4_((S_fp)cofx, idmn, &usol[usol_offset], &grhs[grhs_offset]);
    goto L360;
} /* speli4_ */

