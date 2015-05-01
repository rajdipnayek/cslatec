/* dplpce.f -- translated by f2c (version 12.02.01).
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
static integer c__0 = 0;

/* DECK DPLPCE */
/* Subroutine */ int dplpce_(integer *mrelas, integer *nvars__, integer *lmx, 
	integer *lbm, integer *itlp, integer *itbrc, integer *ibasis, integer 
	*imat, integer *ibrc, integer *ipr, integer *iwr, integer *ind, 
	integer *ibb, doublereal *erdnrm, doublereal *eps, doublereal *tune, 
	doublereal *gg, doublereal *amat, doublereal *basmat, doublereal *csc,
	 doublereal *wr, doublereal *ww, doublereal *primal, doublereal *erd, 
	doublereal *erp, logical *singlr, logical *redbas)
{
    /* System generated locals */
    integer ibrc_dim1, ibrc_offset, i__1;
    doublereal d__1, d__2;

    /* Local variables */
    static integer i__, j, l, il1, iu1, ihi;
    static doublereal one;
    static integer lpg;
    static doublereal ten;
    static integer key, n20002, n20012, n20023, n20016, n20061, n20047, 
	    n20057, ilow;
    static doublereal zero;
    extern /* Subroutine */ int la05bd_(doublereal *, integer *, integer *, 
	    integer *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, logical *);
    static integer ipage;
    extern integer idloc_(integer *, doublereal *, integer *);
    extern doublereal dasum_(integer *, doublereal *, integer *);
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *);
    static logical trans, pagepl;
    static doublereal factor;
    extern /* Subroutine */ int dprwpg_(integer *, integer *, integer *, 
	    doublereal *, integer *);

/* ***BEGIN PROLOGUE  DPLPCE */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to DSPLP */
/* ***LIBRARY   SLATEC */
/* ***TYPE      DOUBLE PRECISION (SPLPCE-S, DPLPCE-D) */
/* ***AUTHOR  (UNKNOWN) */
/* ***DESCRIPTION */

/*     THE EDITING REQUIRED TO CONVERT THIS SUBROUTINE FROM SINGLE TO */
/*     DOUBLE PRECISION INVOLVES THE FOLLOWING CHARACTER STRING CHANGES. */

/*     USE AN EDITING COMMAND (CHANGE) /STRING-1/(TO)STRING-2/. */
/*     /REAL (12 BLANKS)/DOUBLE PRECISION/, */
/*     /SASUM/DASUM/,/DCOPY/,DCOPY/. */

/*     REVISED 811219-1630 */
/*     REVISED YYMMDD-HHMM */

/*     THIS SUBPROGRAM IS FROM THE DSPLP( ) PACKAGE.  IT CALCULATES */
/*     THE APPROXIMATE ERROR IN THE PRIMAL AND DUAL SYSTEMS.  IT IS */
/*     THE MAIN PART OF THE PROCEDURE (COMPUTE ERROR IN DUAL AND PRIMAL */
/*     SYSTEMS). */

/* ***SEE ALSO  DSPLP */
/* ***ROUTINES CALLED  DASUM, DCOPY, DPRWPG, IDLOC, LA05BD */
/* ***REVISION HISTORY  (YYMMDD) */
/*   811215  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890605  Removed unreferenced labels.  (WRB) */
/*   890606  Changed references from IPLOC to IDLOC.  (WRB) */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900328  Added TYPE section.  (WRB) */
/* ***END PROLOGUE  DPLPCE */
/* ***FIRST EXECUTABLE STATEMENT  DPLPCE */
    /* Parameter adjustments */
    ibrc_dim1 = *lbm;
    ibrc_offset = 1 + ibrc_dim1;
    ibrc -= ibrc_offset;
    --ibasis;
    --imat;
    --ipr;
    --iwr;
    --ind;
    --ibb;
    --amat;
    --basmat;
    --csc;
    --wr;
    --ww;
    --primal;
    --erd;
    --erp;

    /* Function Body */
    zero = 0.;
    one = 1.;
    ten = 10.;
    lpg = *lmx - (*nvars__ + 4);
    *singlr = FALSE_;
    factor = .01f;

/*     COPY COLSUMS IN WW(*), AND SOLVE TRANSPOSED SYSTEM. */
    i__ = 1;
    n20002 = *mrelas;
    goto L20003;
L20002:
    ++i__;
L20003:
    if (n20002 - i__ < 0) {
	goto L20004;
    }
    j = ibasis[i__];
    if (! (j <= *nvars__)) {
	goto L20006;
    }
    ww[i__] = primal[j];
    goto L20007;
L20006:
    if (! (ind[j] == 2)) {
	goto L20009;
    }
    ww[i__] = one;
    goto L20010;
L20009:
    ww[i__] = -one;
L20010:
L20007:
    goto L20002;

/*     PERTURB RIGHT-SIDE IN UNITS OF LAST BITS TO BETTER REFLECT */
/*     ERRORS IN THE CHECK SUM SOLNS. */
L20004:
    i__ = 1;
    n20012 = *mrelas;
    goto L20013;
L20012:
    ++i__;
L20013:
    if (n20012 - i__ < 0) {
	goto L20014;
    }
    ww[i__] += ten * *eps * ww[i__];
    goto L20012;
L20014:
    trans = TRUE_;
    la05bd_(&basmat[1], &ibrc[ibrc_offset], lbm, mrelas, &ipr[1], &iwr[1], &
	    wr[1], gg, &ww[1], &trans);
    i__ = 1;
    n20016 = *mrelas;
    goto L20017;
L20016:
    ++i__;
L20017:
    if (n20016 - i__ < 0) {
	goto L20018;
    }
/* Computing MAX */
    d__2 = (d__1 = ww[i__] - one, abs(d__1));
    erd[i__] = max(d__2,*eps) * *tune;

/*     SYSTEM BECOMES SINGULAR WHEN ACCURACY OF SOLUTION IS .GT. FACTOR. */
/*     THIS VALUE (FACTOR) MIGHT NEED TO BE CHANGED. */
    *singlr = *singlr || erd[i__] >= factor;
    goto L20016;
L20018:
    *erdnrm = dasum_(mrelas, &erd[1], &c__1);

/*     RECALCULATE ROW CHECK SUMS EVERY ITBRC ITERATIONS OR WHEN */
/*     A REDECOMPOSITION HAS OCCURRED. */
    if (! (*itlp % *itbrc == 0 || *redbas)) {
	goto L20020;
    }

/*     COMPUTE ROW SUMS, STORE IN WW(*), SOLVE PRIMAL SYSTEM. */
    ww[1] = zero;
    dcopy_(mrelas, &ww[1], &c__0, &ww[1], &c__1);
    pagepl = TRUE_;
    j = 1;
    n20023 = *nvars__;
    goto L20024;
L20023:
    ++j;
L20024:
    if (n20023 - j < 0) {
	goto L20025;
    }
    if (! ((doublereal) ibb[j] >= zero)) {
	goto L20027;
    }

/*     THE VARIABLE IS NON-BASIC. */
    pagepl = TRUE_;
    goto L20023;
L20027:
    if (! (j == 1)) {
	goto L20030;
    }
    ilow = *nvars__ + 5;
    goto L20031;
L20030:
    ilow = imat[j + 3] + 1;
L20031:
    if (! pagepl) {
	goto L20033;
    }
    il1 = idloc_(&ilow, &amat[1], &imat[1]);
    if (! (il1 >= *lmx - 1)) {
	goto L20036;
    }
    ilow += 2;
    il1 = idloc_(&ilow, &amat[1], &imat[1]);
L20036:
    ipage = (i__1 = imat[*lmx - 1], abs(i__1));
    goto L20034;
L20033:
    il1 = ihi + 1;
L20034:
    ihi = imat[j + 4] - (ilow - il1);
L20039:
/* Computing MIN */
    i__1 = *lmx - 2;
    iu1 = min(i__1,ihi);
    if (! (il1 > iu1)) {
	goto L20041;
    }
    goto L20040;
L20041:
    i__1 = iu1;
    for (i__ = il1; i__ <= i__1; ++i__) {
	ww[imat[i__]] += amat[i__] * csc[j];
/* L20: */
    }
    if (! (ihi <= *lmx - 2)) {
	goto L20044;
    }
    goto L20040;
L20044:
    ++ipage;
    key = 1;
    dprwpg_(&key, &ipage, &lpg, &amat[1], &imat[1]);
    il1 = *nvars__ + 5;
    ihi -= lpg;
    goto L20039;
L20040:
    pagepl = ihi == *lmx - 2;
    goto L20023;
L20025:
    l = 1;
    n20047 = *mrelas;
    goto L20048;
L20047:
    ++l;
L20048:
    if (n20047 - l < 0) {
	goto L20049;
    }
    j = ibasis[l];
    if (! (j > *nvars__)) {
	goto L20051;
    }
    i__ = j - *nvars__;
    if (! (ind[j] == 2)) {
	goto L20054;
    }
    ww[i__] += one;
    goto L20055;
L20054:
    ww[i__] -= one;
L20055:
L20051:
    goto L20047;

/*     PERTURB RIGHT-SIDE IN UNITS OF LAST BIT POSITIONS. */
L20049:
    i__ = 1;
    n20057 = *mrelas;
    goto L20058;
L20057:
    ++i__;
L20058:
    if (n20057 - i__ < 0) {
	goto L20059;
    }
    ww[i__] += ten * *eps * ww[i__];
    goto L20057;
L20059:
    trans = FALSE_;
    la05bd_(&basmat[1], &ibrc[ibrc_offset], lbm, mrelas, &ipr[1], &iwr[1], &
	    wr[1], gg, &ww[1], &trans);
    i__ = 1;
    n20061 = *mrelas;
    goto L20062;
L20061:
    ++i__;
L20062:
    if (n20061 - i__ < 0) {
	goto L20063;
    }
/* Computing MAX */
    d__2 = (d__1 = ww[i__] - one, abs(d__1));
    erp[i__] = max(d__2,*eps) * *tune;

/*     SYSTEM BECOMES SINGULAR WHEN ACCURACY OF SOLUTION IS .GT. FACTOR. */
/*     THIS VALUE (FACTOR) MIGHT NEED TO BE CHANGED. */
    *singlr = *singlr || erp[i__] >= factor;
    goto L20061;
L20063:

L20020:
    return 0;
} /* dplpce_ */

