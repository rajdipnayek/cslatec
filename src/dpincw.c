/* dpincw.f -- translated by f2c (version 12.02.01).
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

static integer c__0 = 0;
static integer c__1 = 1;

/* DECK DPINCW */
/* Subroutine */ int dpincw_(integer *mrelas, integer *nvars__, integer *lmx, 
	integer *lbm, integer *npp, integer *jstrt, integer *ibasis, integer *
	imat, integer *ibrc, integer *ipr, integer *iwr, integer *ind, 
	integer *ibb, doublereal *costsc, doublereal *gg, doublereal *erdnrm, 
	doublereal *dulnrm, doublereal *amat, doublereal *basmat, doublereal *
	csc, doublereal *wr, doublereal *ww, doublereal *rz, doublereal *rg, 
	doublereal *costs, doublereal *colnrm, doublereal *duals, logical *
	stpedg)
{
    /* System generated locals */
    integer ibrc_dim1, ibrc_offset, i__1;

    /* Local variables */
    static integer i__, j, il1, iu1, ihi;
    static doublereal one;
    static integer lpg, key;
    static doublereal rzj;
    extern doublereal ddot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    static integer ilow;
    static doublereal zero;
    extern /* Subroutine */ int la05bd_(doublereal *, integer *, integer *, 
	    integer *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, logical *);
    static integer ipage;
    extern integer idloc_(integer *, doublereal *, integer *);
    static doublereal scalr, cnorm;
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *);
    static logical trans;
    static doublereal rcost;
    static logical pagepl;
    static integer nnegrc;
    extern /* Subroutine */ int dprwpg_(integer *, integer *, integer *, 
	    doublereal *, integer *);

/* ***BEGIN PROLOGUE  DPINCW */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to DSPLP */
/* ***LIBRARY   SLATEC */
/* ***TYPE      DOUBLE PRECISION (SPINCW-S, DPINCW-D) */
/* ***AUTHOR  (UNKNOWN) */
/* ***DESCRIPTION */

/*     THE EDITING REQUIRED TO CONVERT THIS SUBROUTINE FROM SINGLE TO */
/*     DOUBLE PRECISION INVOLVES THE FOLLOWING CHARACTER STRING CHANGES. */

/*     USE AN EDITING COMMAND (CHANGE) /STRING-1/(TO)STRING-2/, */
/*     REAL (12 BLANKS)/DOUBLE PRECISION/,/SCOPY/DCOPY/,/SDOT/DDOT/. */

/*     THIS SUBPROGRAM IS PART OF THE DSPLP( ) PACKAGE. */
/*     IT IMPLEMENTS THE PROCEDURE (INITIALIZE REDUCED COSTS AND */
/*     STEEPEST EDGE WEIGHTS). */

/* ***SEE ALSO  DSPLP */
/* ***ROUTINES CALLED  DCOPY, DDOT, DPRWPG, IDLOC, LA05BD */
/* ***REVISION HISTORY  (YYMMDD) */
/*   811215  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890605  Removed unreferenced labels.  (WRB) */
/*   890606  Changed references from IPLOC to IDLOC.  (WRB) */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900328  Added TYPE section.  (WRB) */
/* ***END PROLOGUE  DPINCW */
/* ***FIRST EXECUTABLE STATEMENT  DPINCW */
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
    --rz;
    --rg;
    --costs;
    --colnrm;
    --duals;

    /* Function Body */
    lpg = *lmx - (*nvars__ + 4);
    zero = 0.;
    one = 1.;

/*     FORM REDUCED COSTS, RZ(*), AND STEEPEST EDGE WEIGHTS, RG(*). */
    pagepl = TRUE_;
    rz[1] = zero;
    i__1 = *nvars__ + *mrelas;
    dcopy_(&i__1, &rz[1], &c__0, &rz[1], &c__1);
    rg[1] = one;
    i__1 = *nvars__ + *mrelas;
    dcopy_(&i__1, &rg[1], &c__0, &rg[1], &c__1);
    nnegrc = 0;
    j = *jstrt;
L20002:
    if (! (ibb[j] <= 0)) {
	goto L20004;
    }
    pagepl = TRUE_;
    goto L20005;

/*     THESE ARE NONBASIC INDEPENDENT VARIABLES. THE COLS. ARE IN SPARSE */
/*     MATRIX FORMAT. */
L20004:
    if (! (j <= *nvars__)) {
	goto L20007;
    }
    rzj = *costsc * costs[j];
    ww[1] = zero;
    dcopy_(mrelas, &ww[1], &c__0, &ww[1], &c__1);
    if (! (j == 1)) {
	goto L20010;
    }
    ilow = *nvars__ + 5;
    goto L20011;
L20010:
    ilow = imat[j + 3] + 1;
L20011:
    if (! pagepl) {
	goto L20013;
    }
    il1 = idloc_(&ilow, &amat[1], &imat[1]);
    if (! (il1 >= *lmx - 1)) {
	goto L20016;
    }
    ilow += 2;
    il1 = idloc_(&ilow, &amat[1], &imat[1]);
L20016:
    ipage = (i__1 = imat[*lmx - 1], abs(i__1));
    goto L20014;
L20013:
    il1 = ihi + 1;
L20014:
    ihi = imat[j + 4] - (ilow - il1);
L20019:
/* Computing MIN */
    i__1 = *lmx - 2;
    iu1 = min(i__1,ihi);
    if (! (il1 > iu1)) {
	goto L20021;
    }
    goto L20020;
L20021:
    i__1 = iu1;
    for (i__ = il1; i__ <= i__1; ++i__) {
	rzj -= amat[i__] * duals[imat[i__]];
	ww[imat[i__]] = amat[i__] * csc[j];
/* L60: */
    }
    if (! (ihi <= *lmx - 2)) {
	goto L20024;
    }
    goto L20020;
L20024:
    ++ipage;
    key = 1;
    dprwpg_(&key, &ipage, &lpg, &amat[1], &imat[1]);
    il1 = *nvars__ + 5;
    ihi -= lpg;
    goto L20019;
L20020:
    pagepl = ihi == *lmx - 2;
    rz[j] = rzj * csc[j];
    if (! (*stpedg)) {
	goto L20027;
    }
    trans = FALSE_;
    la05bd_(&basmat[1], &ibrc[ibrc_offset], lbm, mrelas, &ipr[1], &iwr[1], &
	    wr[1], gg, &ww[1], &trans);
    rg[j] = ddot_(mrelas, &ww[1], &c__1, &ww[1], &c__1) + one;
L20027:

/*     THESE ARE NONBASIC DEPENDENT VARIABLES. THE COLS. ARE IMPLICITLY */
/*     DEFINED. */
    goto L20008;
L20007:
    pagepl = TRUE_;
    ww[1] = zero;
    dcopy_(mrelas, &ww[1], &c__0, &ww[1], &c__1);
    scalr = -one;
    if (ind[j] == 2) {
	scalr = one;
    }
    i__ = j - *nvars__;
    rz[j] = -scalr * duals[i__];
    ww[i__] = scalr;
    if (! (*stpedg)) {
	goto L20030;
    }
    trans = FALSE_;
    la05bd_(&basmat[1], &ibrc[ibrc_offset], lbm, mrelas, &ipr[1], &iwr[1], &
	    wr[1], gg, &ww[1], &trans);
    rg[j] = ddot_(mrelas, &ww[1], &c__1, &ww[1], &c__1) + one;
L20030:
L20008:

L20005:
    rcost = rz[j];
    if (ibb[j] % 2 == 0) {
	rcost = -rcost;
    }
    if (ind[j] == 4) {
	rcost = -abs(rcost);
    }
    cnorm = one;
    if (j <= *nvars__) {
	cnorm = colnrm[j];
    }
    if (rcost + *erdnrm * *dulnrm * cnorm < zero) {
	++nnegrc;
    }
    j = j % (*mrelas + *nvars__) + 1;
    if (! (nnegrc >= *npp || j == *jstrt)) {
	goto L20033;
    }
    goto L20003;
L20033:
    goto L20002;
L20003:
    *jstrt = j;
    return 0;
} /* dpincw_ */

