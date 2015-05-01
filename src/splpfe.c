/* splpfe.f -- translated by f2c (version 12.02.01).
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

/* DECK SPLPFE */
/* Subroutine */ int splpfe_(integer *mrelas, integer *nvars__, integer *lmx, 
	integer *lbm, integer *ienter, integer *ibasis, integer *imat, 
	integer *ibrc, integer *ipr, integer *iwr, integer *ind, integer *ibb,
	 real *erdnrm, real *eps, real *gg, real *dulnrm, real *dirnrm, real *
	amat, real *basmat, real *csc, real *wr, real *ww, real *bl, real *bu,
	 real *rz, real *rg, real *colnrm, real *duals, logical *found)
{
    /* System generated locals */
    integer ibrc_dim1, ibrc_offset, i__1;
    real r__1, r__2;

    /* Local variables */
    static integer i__, j, il1, iu1, ihi;
    static real one;
    static integer lpg, key, n20002, n20050;
    static real rmax;
    static integer ilow;
    static real zero;
    static integer ipage;
    extern /* Subroutine */ int la05bs_(real *, integer *, integer *, integer 
	    *, integer *, integer *, real *, real *, real *, logical *);
    extern integer iploc_(integer *, real *, integer *);
    static real cnorm, ratio;
    static logical trans;
    extern doublereal sasum_(integer *, real *, integer *);
    static real rcost;
    extern /* Subroutine */ int scopy_(integer *, real *, integer *, real *, 
	    integer *), prwpge_(integer *, integer *, integer *, real *, 
	    integer *);

/* ***BEGIN PROLOGUE  SPLPFE */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to SPLP */
/* ***LIBRARY   SLATEC */
/* ***TYPE      SINGLE PRECISION (SPLPFE-S, DPLPFE-D) */
/* ***AUTHOR  (UNKNOWN) */
/* ***DESCRIPTION */

/*     THE EDITING REQUIRED TO CONVERT THIS SUBROUTINE FROM SINGLE TO */
/*     DOUBLE PRECISION INVOLVES THE FOLLOWING CHARACTER STRING CHANGES. */

/*     USE AN EDITING COMMAND (CHANGE) /STRING-1/(TO)STRING-2/. */
/*     /REAL (12 BLANKS)/DOUBLE PRECISION/,/SASUM/DASUM/, */
/*     /SCOPY/DCOPY/. */

/*     THIS SUBPROGRAM IS PART OF THE SPLP( ) PACKAGE. */
/*     IT IMPLEMENTS THE PROCEDURE (FIND VARIABLE TO ENTER BASIS */
/*     AND GET SEARCH DIRECTION). */
/*     REVISED 811130-1100 */
/*     REVISED YYMMDD-HHMM */

/* ***SEE ALSO  SPLP */
/* ***ROUTINES CALLED  IPLOC, LA05BS, PRWPGE, SASUM, SCOPY */
/* ***REVISION HISTORY  (YYMMDD) */
/*   811215  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890605  Removed unreferenced labels.  (WRB) */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900328  Added TYPE section.  (WRB) */
/* ***END PROLOGUE  SPLPFE */
/* ***FIRST EXECUTABLE STATEMENT  SPLPFE */
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
    --bl;
    --bu;
    --rz;
    --rg;
    --colnrm;
    --duals;

    /* Function Body */
    lpg = *lmx - (*nvars__ + 4);
    zero = 0.f;
    one = 1.f;
    rmax = zero;
    *found = FALSE_;
    i__ = *mrelas + 1;
    n20002 = *mrelas + *nvars__;
    goto L20003;
L20002:
    ++i__;
L20003:
    if (n20002 - i__ < 0) {
	goto L20004;
    }
    j = ibasis[i__];

/*     IF J=IBASIS(I) .LT. 0 THEN THE VARIABLE LEFT AT A ZERO LEVEL */
/*     AND IS NOT CONSIDERED AS A CANDIDATE TO ENTER. */
    if (! (j > 0)) {
	goto L20006;
    }

/*     DO NOT CONSIDER VARIABLES CORRESPONDING TO UNBOUNDED STEP LENGTHS. */
    if (! (ibb[j] == 0)) {
	goto L20009;
    }
    goto L20002;
L20009:

/*     IF A VARIABLE CORRESPONDS TO AN EQUATION(IND=3 AND BL=BU), */
/*     THEN DO NOT CONSIDER IT AS A CANDIDATE TO ENTER. */
    if (! (ind[j] == 3)) {
	goto L20012;
    }
    if (! (bu[j] - bl[j] <= *eps * ((r__1 = bl[j], dabs(r__1)) + (r__2 = bu[j]
	    , dabs(r__2))))) {
	goto L20015;
    }
    goto L20002;
L20015:
L20012:
    rcost = rz[j];

/*     IF VARIABLE IS AT UPPER BOUND IT CAN ONLY DECREASE.  THIS */
/*     ACCOUNTS FOR THE POSSIBLE CHANGE OF SIGN. */
    if (ibb[j] % 2 == 0) {
	rcost = -rcost;
    }

/*     IF THE VARIABLE IS FREE, USE THE NEGATIVE MAGNITUDE OF THE */
/*     REDUCED COST FOR THAT VARIABLE. */
    if (ind[j] == 4) {
	rcost = -dabs(rcost);
    }
    cnorm = one;
    if (j <= *nvars__) {
	cnorm = colnrm[j];
    }

/*     TEST FOR NEGATIVITY OF REDUCED COSTS. */
    if (! (rcost + *erdnrm * *dulnrm * cnorm < zero)) {
	goto L20018;
    }
    *found = TRUE_;
/* Computing 2nd power */
    r__1 = rcost;
    ratio = r__1 * r__1 / rg[j];
    if (! (ratio > rmax)) {
	goto L20021;
    }
    rmax = ratio;
    *ienter = i__;
L20021:
L20018:
L20006:
    goto L20002;

/*     USE COL. CHOSEN TO COMPUTE SEARCH DIRECTION. */
L20004:
    if (! (*found)) {
	goto L20024;
    }
    j = ibasis[*ienter];
    ww[1] = zero;
    scopy_(mrelas, &ww[1], &c__0, &ww[1], &c__1);
    if (! (j <= *nvars__)) {
	goto L20027;
    }
    if (! (j == 1)) {
	goto L20030;
    }
    ilow = *nvars__ + 5;
    goto L20031;
L20030:
    ilow = imat[j + 3] + 1;
L20031:
    il1 = iploc_(&ilow, &amat[1], &imat[1]);
    if (! (il1 >= *lmx - 1)) {
	goto L20033;
    }
    ilow += 2;
    il1 = iploc_(&ilow, &amat[1], &imat[1]);
L20033:
    ipage = (i__1 = imat[*lmx - 1], abs(i__1));
    ihi = imat[j + 4] - (ilow - il1);
L20036:
/* Computing MIN */
    i__1 = *lmx - 2;
    iu1 = min(i__1,ihi);
    if (! (il1 > iu1)) {
	goto L20038;
    }
    goto L20037;
L20038:
    i__1 = iu1;
    for (i__ = il1; i__ <= i__1; ++i__) {
	ww[imat[i__]] = amat[i__] * csc[j];
/* L30: */
    }
    if (! (ihi <= *lmx - 2)) {
	goto L20041;
    }
    goto L20037;
L20041:
    ++ipage;
    key = 1;
    prwpge_(&key, &ipage, &lpg, &amat[1], &imat[1]);
    il1 = *nvars__ + 5;
    ihi -= lpg;
    goto L20036;
L20037:
    goto L20028;
L20027:
    if (! (ind[j] == 2)) {
	goto L20044;
    }
    ww[j - *nvars__] = one;
    goto L20045;
L20044:
    ww[j - *nvars__] = -one;
L20045:

/*     COMPUTE SEARCH DIRECTION. */
L20028:
    trans = FALSE_;
    la05bs_(&basmat[1], &ibrc[ibrc_offset], lbm, mrelas, &ipr[1], &iwr[1], &
	    wr[1], gg, &ww[1], &trans);

/*     THE SEARCH DIRECTION REQUIRES THE FOLLOWING SIGN CHANGE IF EITHER */
/*     VARIABLE ENTERING IS AT ITS UPPER BOUND OR IS FREE AND HAS */
/*     POSITIVE REDUCED COST. */
    if (! (ibb[j] % 2 == 0 || ind[j] == 4 && rz[j] > zero)) {
	goto L20047;
    }
    i__ = 1;
    n20050 = *mrelas;
    goto L20051;
L20050:
    ++i__;
L20051:
    if (n20050 - i__ < 0) {
	goto L20052;
    }
    ww[i__] = -ww[i__];
    goto L20050;
L20052:
L20047:
    *dirnrm = sasum_(mrelas, &ww[1], &c__1);

/*     COPY CONTENTS OF WR(*) TO DUALS(*) FOR USE IN */
/*     ADD-DROP (EXCHANGE) STEP, LA05CS( ). */
    scopy_(mrelas, &wr[1], &c__1, &duals[1], &c__1);
L20024:
    return 0;
} /* splpfe_ */

