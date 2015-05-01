/* spelip.f -- translated by f2c (version 12.02.01).
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

/* DECK SPELIP */
/* Subroutine */ int spelip_(integer *intl, integer *iorder, real *a, real *b,
	 integer *m, integer *mbdcnd, real *bda, real *alpha, real *bdb, real 
	*beta, real *c__, real *d__, integer *n, integer *nbdcnd, real *bdc, 
	real *gama, real *bdd, real *xnu, S_fp cofx, S_fp cofy, real *an, 
	real *bn, real *cn, real *dn, real *un, real *zn, real *am, real *bm, 
	real *cm, real *dm, real *um, real *zm, real *grhs, real *usol, 
	integer *idmn, real *w, real *pertrb, integer *ierror)
{
    /* System generated locals */
    integer grhs_dim1, grhs_offset, usol_dim1, usol_offset, i__1, i__2;
    real r__1;

    /* Local variables */
    static integer i__, j, i1;
    static real ai, bi, ci, dj, ej, fj;
    static integer mp, np;
    static real xi, yj, ax1, dy1, axi, bxi, cxi, dyj, eyj, fyj, cxm, fyn;
    static integer iord;
    extern /* Subroutine */ int defer_(S_fp, S_fp, integer *, real *, real *);
    static real prtrb;
    extern /* Subroutine */ int trisp_(integer *, real *, real *, real *, 
	    real *, real *, real *), chksng_(integer *, integer *, real *, 
	    real *, real *, real *, S_fp, S_fp, logical *), blktri_(integer *,
	     integer *, integer *, real *, real *, real *, integer *, integer 
	    *, real *, real *, real *, integer *, real *, integer *, real *);
    static logical singlr;
    extern /* Subroutine */ int minsol_(real *, integer *, real *, real *, 
	    real *), orthog_(real *, integer *, real *, real *, real *);

/* ***BEGIN PROLOGUE  SPELIP */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to SEPELI */
/* ***LIBRARY   SLATEC */
/* ***TYPE      SINGLE PRECISION (SPELIP-S) */
/* ***AUTHOR  (UNKNOWN) */
/* ***DESCRIPTION */

/*     SPELIP sets up vectors and arrays for input to BLKTRI */
/*     and computes a second order solution in USOL.  A return jump to */
/*     SEPELI occurs if IORDER=2.  If IORDER=4 a fourth order */
/*     solution is generated in USOL. */

/* ***SEE ALSO  SEPELI */
/* ***ROUTINES CALLED  BLKTRI, CHKSNG, DEFER, MINSOL, ORTHOG, TRISP */
/* ***COMMON BLOCKS    SPLPCM */
/* ***REVISION HISTORY  (YYMMDD) */
/*   801001  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900402  Added TYPE section.  (WRB) */
/* ***END PROLOGUE  SPELIP */

/* ***FIRST EXECUTABLE STATEMENT  SPELIP */
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
    splpcm_1.kswx = *mbdcnd + 1;
    splpcm_1.kswy = *nbdcnd + 1;
    splpcm_1.k = *m + 1;
    splpcm_1.l = *n + 1;
    splpcm_1.ait = *a;
    splpcm_1.bit = *b;
    splpcm_1.cit = *c__;
    splpcm_1.dit = *d__;

/*     SET RIGHT HAND SIDE VALUES FROM GRHS IN USOL ON THE INTERIOR */
/*     AND NON-SPECIFIED BOUNDARIES. */

    i__1 = *m;
    for (i__ = 2; i__ <= i__1; ++i__) {
	i__2 = *n;
	for (j = 2; j <= i__2; ++j) {
	    usol[i__ + j * usol_dim1] = grhs[i__ + j * grhs_dim1];
/* L10: */
	}
/* L20: */
    }
    if (splpcm_1.kswx == 2 || splpcm_1.kswx == 3) {
	goto L40;
    }
    i__1 = *n;
    for (j = 2; j <= i__1; ++j) {
	usol[j * usol_dim1 + 1] = grhs[j * grhs_dim1 + 1];
/* L30: */
    }
L40:
    if (splpcm_1.kswx == 2 || splpcm_1.kswx == 5) {
	goto L60;
    }
    i__1 = *n;
    for (j = 2; j <= i__1; ++j) {
	usol[splpcm_1.k + j * usol_dim1] = grhs[splpcm_1.k + j * grhs_dim1];
/* L50: */
    }
L60:
    if (splpcm_1.kswy == 2 || splpcm_1.kswy == 3) {
	goto L80;
    }
    i__1 = *m;
    for (i__ = 2; i__ <= i__1; ++i__) {
	usol[i__ + usol_dim1] = grhs[i__ + grhs_dim1];
/* L70: */
    }
L80:
    if (splpcm_1.kswy == 2 || splpcm_1.kswy == 5) {
	goto L100;
    }
    i__1 = *m;
    for (i__ = 2; i__ <= i__1; ++i__) {
	usol[i__ + splpcm_1.l * usol_dim1] = grhs[i__ + splpcm_1.l * 
		grhs_dim1];
/* L90: */
    }
L100:
    if (splpcm_1.kswx != 2 && splpcm_1.kswx != 3 && splpcm_1.kswy != 2 && 
	    splpcm_1.kswy != 3) {
	usol[usol_dim1 + 1] = grhs[grhs_dim1 + 1];
    }
    if (splpcm_1.kswx != 2 && splpcm_1.kswx != 5 && splpcm_1.kswy != 2 && 
	    splpcm_1.kswy != 3) {
	usol[splpcm_1.k + usol_dim1] = grhs[splpcm_1.k + grhs_dim1];
    }
    if (splpcm_1.kswx != 2 && splpcm_1.kswx != 3 && splpcm_1.kswy != 2 && 
	    splpcm_1.kswy != 5) {
	usol[splpcm_1.l * usol_dim1 + 1] = grhs[splpcm_1.l * grhs_dim1 + 1];
    }
    if (splpcm_1.kswx != 2 && splpcm_1.kswx != 5 && splpcm_1.kswy != 2 && 
	    splpcm_1.kswy != 5) {
	usol[splpcm_1.k + splpcm_1.l * usol_dim1] = grhs[splpcm_1.k + 
		splpcm_1.l * grhs_dim1];
    }
    i1 = 1;

/*     SET SWITCHES FOR PERIODIC OR NON-PERIODIC BOUNDARIES */

    mp = 1;
    np = 1;
    if (splpcm_1.kswx == 1) {
	mp = 0;
    }
    if (splpcm_1.kswy == 1) {
	np = 0;
    }

/*     SET DLX,DLY AND SIZE OF BLOCK TRI-DIAGONAL SYSTEM GENERATED */
/*     IN NINT,MINT */

    splpcm_1.dlx = (splpcm_1.bit - splpcm_1.ait) / *m;
    splpcm_1.mit = splpcm_1.k - 1;
    if (splpcm_1.kswx == 2) {
	splpcm_1.mit = splpcm_1.k - 2;
    }
    if (splpcm_1.kswx == 4) {
	splpcm_1.mit = splpcm_1.k;
    }
    splpcm_1.dly = (splpcm_1.dit - splpcm_1.cit) / *n;
    splpcm_1.nit = splpcm_1.l - 1;
    if (splpcm_1.kswy == 2) {
	splpcm_1.nit = splpcm_1.l - 2;
    }
    if (splpcm_1.kswy == 4) {
	splpcm_1.nit = splpcm_1.l;
    }
/* Computing 3rd power */
    r__1 = splpcm_1.dlx;
    splpcm_1.tdlx3 = r__1 * (r__1 * r__1) * 2.f;
/* Computing 4th power */
    r__1 = splpcm_1.dlx, r__1 *= r__1;
    splpcm_1.dlx4 = r__1 * r__1;
/* Computing 3rd power */
    r__1 = splpcm_1.dly;
    splpcm_1.tdly3 = r__1 * (r__1 * r__1) * 2.f;
/* Computing 4th power */
    r__1 = splpcm_1.dly, r__1 *= r__1;
    splpcm_1.dly4 = r__1 * r__1;

/*     SET SUBSCRIPT LIMITS FOR PORTION OF ARRAY TO INPUT TO BLKTRI */

    splpcm_1.is = 1;
    splpcm_1.js = 1;
    if (splpcm_1.kswx == 2 || splpcm_1.kswx == 3) {
	splpcm_1.is = 2;
    }
    if (splpcm_1.kswy == 2 || splpcm_1.kswy == 3) {
	splpcm_1.js = 2;
    }
    splpcm_1.ns = splpcm_1.nit + splpcm_1.js - 1;
    splpcm_1.ms = splpcm_1.mit + splpcm_1.is - 1;

/*     SET X - DIRECTION */

    i__1 = splpcm_1.mit;
    for (i__ = 1; i__ <= i__1; ++i__) {
	xi = splpcm_1.ait + (splpcm_1.is + i__ - 2) * splpcm_1.dlx;
	(*cofx)(&xi, &ai, &bi, &ci);
	axi = (ai / splpcm_1.dlx - bi * .5f) / splpcm_1.dlx;
/* Computing 2nd power */
	r__1 = splpcm_1.dlx;
	bxi = ai * -2.f / (r__1 * r__1) + ci;
	cxi = (ai / splpcm_1.dlx + bi * .5f) / splpcm_1.dlx;
	am[i__] = axi;
	bm[i__] = bxi;
	cm[i__] = cxi;
/* L110: */
    }

/*     SET Y DIRECTION */

    i__1 = splpcm_1.nit;
    for (j = 1; j <= i__1; ++j) {
	yj = splpcm_1.cit + (splpcm_1.js + j - 2) * splpcm_1.dly;
	(*cofy)(&yj, &dj, &ej, &fj);
	dyj = (dj / splpcm_1.dly - ej * .5f) / splpcm_1.dly;
/* Computing 2nd power */
	r__1 = splpcm_1.dly;
	eyj = dj * -2.f / (r__1 * r__1) + fj;
	fyj = (dj / splpcm_1.dly + ej * .5f) / splpcm_1.dly;
	an[j] = dyj;
	bn[j] = eyj;
	cn[j] = fyj;
/* L120: */
    }

/*     ADJUST EDGES IN X DIRECTION UNLESS PERIODIC */

    ax1 = am[1];
    cxm = cm[splpcm_1.mit];
    switch (splpcm_1.kswx) {
	case 1:  goto L170;
	case 2:  goto L130;
	case 3:  goto L150;
	case 4:  goto L160;
	case 5:  goto L140;
    }

/*     DIRICHLET-DIRICHLET IN X DIRECTION */

L130:
    am[1] = 0.f;
    cm[splpcm_1.mit] = 0.f;
    goto L170;

/*     MIXED-DIRICHLET IN X DIRECTION */

L140:
    am[1] = 0.f;
    bm[1] += *alpha * 2.f * splpcm_1.dlx * ax1;
    cm[1] += ax1;
    cm[splpcm_1.mit] = 0.f;
    goto L170;

/*     DIRICHLET-MIXED IN X DIRECTION */

L150:
    am[1] = 0.f;
    am[splpcm_1.mit] += cxm;
    bm[splpcm_1.mit] -= *beta * 2.f * splpcm_1.dlx * cxm;
    cm[splpcm_1.mit] = 0.f;
    goto L170;

/*     MIXED - MIXED IN X DIRECTION */

L160:
    am[1] = 0.f;
    bm[1] += splpcm_1.dlx * 2.f * *alpha * ax1;
    cm[1] += ax1;
    am[splpcm_1.mit] += cxm;
    bm[splpcm_1.mit] -= splpcm_1.dlx * 2.f * *beta * cxm;
    cm[splpcm_1.mit] = 0.f;
L170:

/*     ADJUST IN Y DIRECTION UNLESS PERIODIC */

    dy1 = an[1];
    fyn = cn[splpcm_1.nit];
    switch (splpcm_1.kswy) {
	case 1:  goto L220;
	case 2:  goto L180;
	case 3:  goto L200;
	case 4:  goto L210;
	case 5:  goto L190;
    }

/*     DIRICHLET-DIRICHLET IN Y DIRECTION */

L180:
    an[1] = 0.f;
    cn[splpcm_1.nit] = 0.f;
    goto L220;

/*     MIXED-DIRICHLET IN Y DIRECTION */

L190:
    an[1] = 0.f;
    bn[1] += splpcm_1.dly * 2.f * *gama * dy1;
    cn[1] += dy1;
    cn[splpcm_1.nit] = 0.f;
    goto L220;

/*     DIRICHLET-MIXED IN Y DIRECTION */

L200:
    an[1] = 0.f;
    an[splpcm_1.nit] += fyn;
    bn[splpcm_1.nit] -= splpcm_1.dly * 2.f * *xnu * fyn;
    cn[splpcm_1.nit] = 0.f;
    goto L220;

/*     MIXED - MIXED DIRECTION IN Y DIRECTION */

L210:
    an[1] = 0.f;
    bn[1] += splpcm_1.dly * 2.f * *gama * dy1;
    cn[1] += dy1;
    an[splpcm_1.nit] += fyn;
    bn[splpcm_1.nit] -= splpcm_1.dly * 2.f * *xnu * fyn;
    cn[splpcm_1.nit] = 0.f;
L220:
    if (splpcm_1.kswx == 1) {
	goto L270;
    }

/*     ADJUST USOL ALONG X EDGE */

    i__1 = splpcm_1.ns;
    for (j = splpcm_1.js; j <= i__1; ++j) {
	if (splpcm_1.kswx != 2 && splpcm_1.kswx != 3) {
	    goto L230;
	}
	usol[splpcm_1.is + j * usol_dim1] -= ax1 * usol[j * usol_dim1 + 1];
	goto L240;
L230:
	usol[splpcm_1.is + j * usol_dim1] += splpcm_1.dlx * 2.f * ax1 * bda[j]
		;
L240:
	if (splpcm_1.kswx != 2 && splpcm_1.kswx != 5) {
	    goto L250;
	}
	usol[splpcm_1.ms + j * usol_dim1] -= cxm * usol[splpcm_1.k + j * 
		usol_dim1];
	goto L260;
L250:
	usol[splpcm_1.ms + j * usol_dim1] -= splpcm_1.dlx * 2.f * cxm * bdb[j]
		;
L260:
	;
    }
L270:
    if (splpcm_1.kswy == 1) {
	goto L320;
    }

/*     ADJUST USOL ALONG Y EDGE */

    i__1 = splpcm_1.ms;
    for (i__ = splpcm_1.is; i__ <= i__1; ++i__) {
	if (splpcm_1.kswy != 2 && splpcm_1.kswy != 3) {
	    goto L280;
	}
	usol[i__ + splpcm_1.js * usol_dim1] -= dy1 * usol[i__ + usol_dim1];
	goto L290;
L280:
	usol[i__ + splpcm_1.js * usol_dim1] += splpcm_1.dly * 2.f * dy1 * bdc[
		i__];
L290:
	if (splpcm_1.kswy != 2 && splpcm_1.kswy != 5) {
	    goto L300;
	}
	usol[i__ + splpcm_1.ns * usol_dim1] -= fyn * usol[i__ + splpcm_1.l * 
		usol_dim1];
	goto L310;
L300:
	usol[i__ + splpcm_1.ns * usol_dim1] -= splpcm_1.dly * 2.f * fyn * bdd[
		i__];
L310:
	;
    }
L320:

/*     SAVE ADJUSTED EDGES IN GRHS IF IORDER=4 */

    if (*iorder != 4) {
	goto L350;
    }
    i__1 = splpcm_1.ns;
    for (j = splpcm_1.js; j <= i__1; ++j) {
	grhs[splpcm_1.is + j * grhs_dim1] = usol[splpcm_1.is + j * usol_dim1];
	grhs[splpcm_1.ms + j * grhs_dim1] = usol[splpcm_1.ms + j * usol_dim1];
/* L330: */
    }
    i__1 = splpcm_1.ms;
    for (i__ = splpcm_1.is; i__ <= i__1; ++i__) {
	grhs[i__ + splpcm_1.js * grhs_dim1] = usol[i__ + splpcm_1.js * 
		usol_dim1];
	grhs[i__ + splpcm_1.ns * grhs_dim1] = usol[i__ + splpcm_1.ns * 
		usol_dim1];
/* L340: */
    }
L350:
    iord = *iorder;
    *pertrb = 0.f;

/*     CHECK IF OPERATOR IS SINGULAR */

    chksng_(mbdcnd, nbdcnd, alpha, beta, gama, xnu, (S_fp)cofx, (S_fp)cofy, &
	    singlr);

/*     COMPUTE NON-ZERO EIGENVECTOR IN NULL SPACE OF TRANSPOSE */
/*     IF SINGULAR */

    if (singlr) {
	trisp_(&splpcm_1.mit, &am[1], &bm[1], &cm[1], &dm[1], &um[1], &zm[1]);
    }
    if (singlr) {
	trisp_(&splpcm_1.nit, &an[1], &bn[1], &cn[1], &dn[1], &un[1], &zn[1]);
    }

/*     MAKE INITIALIZATION CALL TO BLKTRI */

    if (*intl == 0) {
	blktri_(intl, &np, &splpcm_1.nit, &an[1], &bn[1], &cn[1], &mp, &
		splpcm_1.mit, &am[1], &bm[1], &cm[1], idmn, &usol[splpcm_1.is 
		+ splpcm_1.js * usol_dim1], ierror, &w[1]);
    }
    if (*ierror != 0) {
	return 0;
    }

/*     ADJUST RIGHT HAND SIDE IF NECESSARY */

L360:
    if (singlr) {
	orthog_(&usol[usol_offset], idmn, &zn[1], &zm[1], pertrb);
    }

/*     COMPUTE SOLUTION */

    blktri_(&i1, &np, &splpcm_1.nit, &an[1], &bn[1], &cn[1], &mp, &
	    splpcm_1.mit, &am[1], &bm[1], &cm[1], idmn, &usol[splpcm_1.is + 
	    splpcm_1.js * usol_dim1], ierror, &w[1]);
    if (*ierror != 0) {
	return 0;
    }

/*     SET PERIODIC BOUNDARIES IF NECESSARY */

    if (splpcm_1.kswx != 1) {
	goto L380;
    }
    i__1 = splpcm_1.l;
    for (j = 1; j <= i__1; ++j) {
	usol[splpcm_1.k + j * usol_dim1] = usol[j * usol_dim1 + 1];
/* L370: */
    }
L380:
    if (splpcm_1.kswy != 1) {
	goto L400;
    }
    i__1 = splpcm_1.k;
    for (i__ = 1; i__ <= i__1; ++i__) {
	usol[i__ + splpcm_1.l * usol_dim1] = usol[i__ + usol_dim1];
/* L390: */
    }
L400:

/*     MINIMIZE SOLUTION WITH RESPECT TO WEIGHTED LEAST SQUARES */
/*     NORM IF OPERATOR IS SINGULAR */

    if (singlr) {
	minsol_(&usol[usol_offset], idmn, &zn[1], &zm[1], &prtrb);
    }

/*     RETURN IF DEFERRED CORRECTIONS AND A FOURTH ORDER SOLUTION ARE */
/*     NOT FLAGGED */

    if (iord == 2) {
	return 0;
    }
    iord = 2;

/*     COMPUTE NEW RIGHT HAND SIDE FOR FOURTH ORDER SOLUTION */

    defer_((S_fp)cofx, (S_fp)cofy, idmn, &usol[usol_offset], &grhs[
	    grhs_offset]);
    goto L360;
} /* spelip_ */

