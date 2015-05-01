/* splpmu.f -- translated by f2c (version 12.02.01).
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

/* DECK SPLPMU */
/* Subroutine */ int splpmu_(integer *mrelas, integer *nvars__, integer *lmx, 
	integer *lbm, integer *nredc, integer *info, integer *ienter, integer 
	*ileave, integer *iopt, integer *npp, integer *jstrt, integer *ibasis,
	 integer *imat, integer *ibrc, integer *ipr, integer *iwr, integer *
	ind, integer *ibb, real *anorm, real *eps, real *uu, real *gg, real *
	rprnrm, real *erdnrm, real *dulnrm, real *theta, real *costsc, real *
	xlamda, real *rhsnrm, real *amat, real *basmat, real *csc, real *wr, 
	real *rprim, real *ww, real *bu, real *bl, real *rhs, real *erd, real 
	*erp, real *rz, real *rg, real *colnrm, real *costs, real *primal, 
	real *duals, logical *singlr, logical *redbas, logical *zerolv, 
	logical *stpedg)
{
    /* Format strings */
    static char fmt_20009[] = "";
    static char fmt_20013[] = "";
    static char fmt_20017[] = "";
    static char fmt_20028[] = "";
    static char fmt_20042[] = "";
    static char fmt_20074[] = "";
    static char fmt_20109[] = "";

    /* System generated locals */
    integer ibrc_dim1, ibrc_offset, i__1;
    real r__1, r__2, r__3, r__4;

    /* Local variables */
    static integer i__, j, k;
    static real gq, wp;
    static integer il1, iu1;
    static real aij;
    static integer ihi;
    static real one;
    static integer lpg, key;
    static real rzj, two;
    static integer n20002, n20121, n20018, ibas, nerr;
    extern doublereal sdot_(integer *, real *, integer *, real *, integer *);
    static integer ilow;
    static real zero;
    static integer npr001, npr003;
    static real gamma, alpha;
    extern /* Subroutine */ int la05bs_(real *, integer *, integer *, integer 
	    *, integer *, integer *, real *, real *, real *, logical *), 
	    la05cs_(real *, integer *, integer *, integer *, integer *, 
	    integer *, real *, real *, real *, integer *);
    static integer ipage;
    static real scalr;
    extern integer iploc_(integer *, real *, integer *);
    static real cnorm;
    static logical trans;
    extern doublereal sasum_(integer *, real *, integer *);
    static real rcost;
    extern /* Subroutine */ int scopy_(integer *, real *, integer *, real *, 
	    integer *);
    static integer iplace;
    static logical pagepl;
    static integer nnegrc;
    extern /* Subroutine */ int splpdm_(integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, real *, 
	    real *, real *, real *, real *, real *, real *, real *, logical *,
	     logical *), prwpge_(integer *, integer *, integer *, real *, 
	    integer *), xermsg_(char *, char *, char *, integer *, integer *, 
	    ftnlen, ftnlen, ftnlen), pnnzrs_(integer *, real *, integer *, 
	    real *, integer *, integer *);

    /* Assigned format variables */
    static char *npr001_fmt, *npr003_fmt;

/* ***BEGIN PROLOGUE  SPLPMU */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to SPLP */
/* ***LIBRARY   SLATEC */
/* ***TYPE      SINGLE PRECISION (SPLPMU-S, DPLPMU-D) */
/* ***AUTHOR  (UNKNOWN) */
/* ***DESCRIPTION */

/*     THE EDITING REQUIRED TO CONVERT THIS SUBROUTINE FROM SINGLE TO */
/*     DOUBLE PRECISION INVOLVES THE FOLLOWING CHARACTER STRING CHANGES. */

/*     USE AN EDITING COMMAND (CHANGE) /STRING-1/(TO)STRING-2/. */
/*     /REAL (12 BLANKS)/DOUBLE PRECISION/, */
/*     /SASUM/DASUM/,/SCOPY/DCOPY/,/SDOT/DDOT/, */
/*     /.E0/.D0/ */

/*     THIS SUBPROGRAM IS FROM THE SPLP( ) PACKAGE.  IT PERFORMS THE */
/*     TASKS OF UPDATING THE PRIMAL SOLUTION, EDGE WEIGHTS, REDUCED */
/*     COSTS, AND MATRIX DECOMPOSITION. */
/*     IT IS THE MAIN PART OF THE PROCEDURE (MAKE MOVE AND UPDATE). */

/*     REVISED 821122-1100 */
/*     REVISED YYMMDD */

/* ***SEE ALSO  SPLP */
/* ***ROUTINES CALLED  IPLOC, LA05BS, LA05CS, PNNZRS, PRWPGE, SASUM, */
/*                    SCOPY, SDOT, SPLPDM, XERMSG */
/* ***REVISION HISTORY  (YYMMDD) */
/*   811215  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890605  Removed unreferenced labels.  (WRB) */
/*   890606  Removed unused COMMON block LA05DS.  (WRB) */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ) */
/*   900328  Added TYPE section.  (WRB) */
/* ***END PROLOGUE  SPLPMU */

/* ***FIRST EXECUTABLE STATEMENT  SPLPMU */
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
    --rprim;
    --ww;
    --bu;
    --bl;
    --rhs;
    --erd;
    --erp;
    --rz;
    --rg;
    --colnrm;
    --costs;
    --primal;
    --duals;

    /* Function Body */
    zero = 0.f;
    one = 1.f;
    two = 2.f;
    lpg = *lmx - (*nvars__ + 4);

/*     UPDATE THE PRIMAL SOLUTION WITH A MULTIPLE OF THE SEARCH */
/*     DIRECTION. */
    i__ = 1;
    n20002 = *mrelas;
    goto L20003;
L20002:
    ++i__;
L20003:
    if (n20002 - i__ < 0) {
	goto L20004;
    }
    rprim[i__] -= *theta * ww[i__];
    goto L20002;

/*     IF EJECTED VARIABLE IS LEAVING AT AN UPPER BOUND,  THEN */
/*     TRANSLATE RIGHT HAND SIDE. */
L20004:
    if (! (*ileave < 0)) {
	goto L20006;
    }
    ibas = ibasis[abs(*ileave)];
    scalr = rprim[abs(*ileave)];
    npr001 = 0;
    npr001_fmt = fmt_20009;
    goto L30001;
L20009:
    ibb[ibas] = (i__1 = ibb[ibas], abs(i__1)) + 1;

/*     IF ENTERING VARIABLE IS RESTRICTED TO ITS UPPER BOUND, TRANSLATE */
/*     RIGHT HAND SIDE.  IF THE VARIABLE DECREASED FROM ITS UPPER */
/*     BOUND, A SIGN CHANGE IS REQUIRED IN THE TRANSLATION. */
L20006:
    if (! (*ienter == *ileave)) {
	goto L20010;
    }
    ibas = ibasis[*ienter];
    scalr = *theta;
    if (ibb[ibas] % 2 == 0) {
	scalr = -scalr;
    }
    npr001 = 1;
    npr001_fmt = fmt_20013;
    goto L30001;
L20013:
    ++ibb[ibas];
    goto L20011;
L20010:
    ibas = ibasis[*ienter];

/*     IF ENTERING VARIABLE IS DECREASING FROM ITS UPPER BOUND, */
/*     COMPLEMENT ITS PRIMAL VALUE. */
    if (! (ind[ibas] == 3 && ibb[ibas] % 2 == 0)) {
	goto L20014;
    }
    scalr = -(bu[ibas] - bl[ibas]);
    if (ibas <= *nvars__) {
	scalr /= csc[ibas];
    }
    npr001 = 2;
    npr001_fmt = fmt_20017;
    goto L30001;
L20017:
    *theta = -scalr - *theta;
    ++ibb[ibas];
L20014:
    rprim[abs(*ileave)] = *theta;
    ibb[ibas] = -(i__1 = ibb[ibas], abs(i__1));
    i__ = ibasis[abs(*ileave)];
    ibb[i__] = (i__1 = ibb[i__], abs(i__1));
    if (primal[abs(*ileave) + *nvars__] > zero) {
	++ibb[i__];
    }

/*     INTERCHANGE COLUMN POINTERS TO NOTE EXCHANGE OF COLUMNS. */
L20011:
    ibas = ibasis[*ienter];
    ibasis[*ienter] = ibasis[abs(*ileave)];
    ibasis[abs(*ileave)] = ibas;

/*     IF VARIABLE WAS EXCHANGED AT A ZERO LEVEL, MARK IT SO THAT */
/*     IT CAN'T BE BROUGHT BACK IN.  THIS IS TO HELP PREVENT CYCLING. */
    if (*zerolv) {
	ibasis[*ienter] = -(i__1 = ibasis[*ienter], abs(i__1));
    }
/* Computing MAX */
    r__1 = *rprnrm, r__2 = sasum_(mrelas, &rprim[1], &c__1);
    *rprnrm = dmax(r__1,r__2);
    k = 1;
    n20018 = *mrelas;
    goto L20019;
L20018:
    ++k;
L20019:
    if (n20018 - k < 0) {
	goto L20020;
    }

/*     SEE IF VARIABLES THAT WERE CLASSIFIED AS INFEASIBLE HAVE NOW */
/*     BECOME FEASIBLE.  THIS MAY REQUIRED TRANSLATING UPPER BOUNDED */
/*     VARIABLES. */
    if (! (primal[k + *nvars__] != zero && (r__1 = rprim[k], dabs(r__1)) <= *
	    rprnrm * erp[k])) {
	goto L20022;
    }
    if (! (primal[k + *nvars__] > zero)) {
	goto L20025;
    }
    ibas = ibasis[k];
    scalr = -(bu[ibas] - bl[ibas]);
    if (ibas <= *nvars__) {
	scalr /= csc[ibas];
    }
    npr001 = 3;
    npr001_fmt = fmt_20028;
    goto L30001;
L20028:
    rprim[k] = -scalr;
    *rprnrm -= scalr;
L20025:
    primal[k + *nvars__] = zero;
L20022:
    goto L20018;

/*     UPDATE REDUCED COSTS, EDGE WEIGHTS, AND MATRIX DECOMPOSITION. */
L20020:
    if (! (*ienter != *ileave)) {
	goto L20029;
    }

/*     THE INCOMING VARIABLE IS ALWAYS CLASSIFIED AS FEASIBLE. */
    primal[abs(*ileave) + *nvars__] = zero;

    wp = ww[abs(*ileave)];
    gq = sdot_(mrelas, &ww[1], &c__1, &ww[1], &c__1) + one;

/*     COMPUTE INVERSE (TRANSPOSE) TIMES SEARCH DIRECTION. */
    trans = TRUE_;
    la05bs_(&basmat[1], &ibrc[ibrc_offset], lbm, mrelas, &ipr[1], &iwr[1], &
	    wr[1], gg, &ww[1], &trans);

/*     UPDATE THE MATRIX DECOMPOSITION.  COL. ABS(ILEAVE) IS LEAVING. */
/*     THE ARRAY DUALS(*) CONTAINS INTERMEDIATE RESULTS FOR THE */
/*     INCOMING COLUMN. */
    i__1 = abs(*ileave);
    la05cs_(&basmat[1], &ibrc[ibrc_offset], lbm, mrelas, &ipr[1], &iwr[1], &
	    duals[1], gg, uu, &i__1);
    *redbas = FALSE_;
    if (! (*gg < zero)) {
	goto L20032;
    }

/*     REDECOMPOSE BASIS MATRIX WHEN AN ERROR RETURN FROM */
/*     LA05CS( ) IS NOTED.  THIS WILL PROBABLY BE DUE TO */
/*     SPACE BEING EXHAUSTED, GG=-7. */
    splpdm_(mrelas, nvars__, lmx, lbm, nredc, info, iopt, &ibasis[1], &imat[1]
	    , &ibrc[ibrc_offset], &ipr[1], &iwr[1], &ind[1], &ibb[1], anorm, 
	    eps, uu, gg, &amat[1], &basmat[1], &csc[1], &wr[1], singlr, 
	    redbas);
    if (! (*singlr)) {
	goto L20035;
    }
    nerr = 26;
    xermsg_("SLATEC", "SPLPMU", "IN SPLP, MOVED TO A SINGULAR POINT.  THIS S"
	    "HOULD NOT HAPPEN.", &nerr, iopt, (ftnlen)6, (ftnlen)6, (ftnlen)60)
	    ;
    *info = -nerr;
    return 0;
L20035:
    goto L30002;
L20038:
L20032:

/*     IF STEEPEST EDGE PRICING IS USED, UPDATE REDUCED COSTS */
/*     AND EDGE WEIGHTS. */
    if (! (*stpedg)) {
	goto L20039;
    }

/*     COMPUTE COL. ABS(ILEAVE) OF THE NEW INVERSE (TRANSPOSE) MATRIX */
/*     HERE ABS(ILEAVE) POINTS TO THE EJECTED COLUMN. */
/*     USE ERD(*) FOR TEMP. STORAGE. */
    scopy_(mrelas, &zero, &c__0, &erd[1], &c__1);
    erd[abs(*ileave)] = one;
    trans = TRUE_;
    la05bs_(&basmat[1], &ibrc[ibrc_offset], lbm, mrelas, &ipr[1], &iwr[1], &
	    wr[1], gg, &erd[1], &trans);

/*     COMPUTE UPDATED DUAL VARIABLES IN DUALS(*). */
    npr003 = 0;
    npr003_fmt = fmt_20042;
    goto L30003;

/*     COMPUTE THE DOT PRODUCT OF COL. J OF THE NEW INVERSE (TRANSPOSE) */
/*     WITH EACH NON-BASIC COLUMN.  ALSO COMPUTE THE DOT PRODUCT OF THE */
/*     INVERSE (TRANSPOSE) OF NON-UPDATED MATRIX (TIMES) THE */
/*     SEARCH DIRECTION WITH EACH NON-BASIC COLUMN. */
/*     RECOMPUTE REDUCED COSTS. */
L20042:
    pagepl = TRUE_;
    i__1 = *nvars__ + *mrelas;
    scopy_(&i__1, &zero, &c__0, &rz[1], &c__1);
    nnegrc = 0;
    j = *jstrt;
L20043:
    if (! (ibb[j] <= 0)) {
	goto L20045;
    }
    pagepl = TRUE_;
    rg[j] = one;
    goto L20046;

/*     NONBASIC INDEPENDENT VARIABLES (COLUMN IN SPARSE MATRIX STORAGE) */
L20045:
    if (! (j <= *nvars__)) {
	goto L20048;
    }
    rzj = costs[j] * *costsc;
    alpha = zero;
    gamma = zero;

/*     COMPUTE THE DOT PRODUCT OF THE SPARSE MATRIX NONBASIC COLUMNS */
/*     WITH THREE VECTORS INVOLVED IN THE UPDATING STEP. */
    if (! (j == 1)) {
	goto L20051;
    }
    ilow = *nvars__ + 5;
    goto L20052;
L20051:
    ilow = imat[j + 3] + 1;
L20052:
    if (! pagepl) {
	goto L20054;
    }
    il1 = iploc_(&ilow, &amat[1], &imat[1]);
    if (! (il1 >= *lmx - 1)) {
	goto L20057;
    }
    ilow += 2;
    il1 = iploc_(&ilow, &amat[1], &imat[1]);
L20057:
    ipage = (i__1 = imat[*lmx - 1], abs(i__1));
    goto L20055;
L20054:
    il1 = ihi + 1;
L20055:
    ihi = imat[j + 4] - (ilow - il1);
L20060:
/* Computing MIN */
    i__1 = *lmx - 2;
    iu1 = min(i__1,ihi);
    if (! (il1 > iu1)) {
	goto L20062;
    }
    goto L20061;
L20062:
    i__1 = iu1;
    for (i__ = il1; i__ <= i__1; ++i__) {
	rzj -= amat[i__] * duals[imat[i__]];
	alpha += amat[i__] * erd[imat[i__]];
	gamma += amat[i__] * ww[imat[i__]];
/* L10: */
    }
    if (! (ihi <= *lmx - 2)) {
	goto L20065;
    }
    goto L20061;
L20065:
    ++ipage;
    key = 1;
    prwpge_(&key, &ipage, &lpg, &amat[1], &imat[1]);
    il1 = *nvars__ + 5;
    ihi -= lpg;
    goto L20060;
L20061:
    pagepl = ihi == *lmx - 2;
    rz[j] = rzj * csc[j];
    alpha *= csc[j];
    gamma *= csc[j];
/* Computing MAX */
/* Computing 2nd power */
    r__3 = alpha;
/* Computing 2nd power */
    r__4 = alpha;
    r__1 = rg[j] - two * alpha * gamma + r__3 * r__3 * gq, r__2 = one + r__4 *
	     r__4;
    rg[j] = dmax(r__1,r__2);

/*     NONBASIC DEPENDENT VARIABLES (COLUMNS DEFINED IMPLICITLY) */
    goto L20049;
L20048:
    pagepl = TRUE_;
    scalr = -one;
    if (ind[j] == 2) {
	scalr = one;
    }
    i__ = j - *nvars__;
    alpha = scalr * erd[i__];
    rz[j] = -scalr * duals[i__];
    gamma = scalr * ww[i__];
/* Computing MAX */
/* Computing 2nd power */
    r__3 = alpha;
/* Computing 2nd power */
    r__4 = alpha;
    r__1 = rg[j] - two * alpha * gamma + r__3 * r__3 * gq, r__2 = one + r__4 *
	     r__4;
    rg[j] = dmax(r__1,r__2);
L20049:
L20046:

    rcost = rz[j];
    if (ibb[j] % 2 == 0) {
	rcost = -rcost;
    }
    if (! (ind[j] == 3)) {
	goto L20068;
    }
    if (bu[j] == bl[j]) {
	rcost = zero;
    }
L20068:
    if (ind[j] == 4) {
	rcost = -dabs(rcost);
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
	goto L20071;
    }
    goto L20044;
L20071:
    goto L20043;
L20044:
    *jstrt = j;

/*     UPDATE THE EDGE WEIGHT FOR THE EJECTED VARIABLE. */
/* Computing 2nd power */
    r__1 = wp;
    rg[i__1 = ibasis[*ienter], abs(i__1)] = gq / (r__1 * r__1);

/*     IF MINIMUM REDUCED COST (DANTZIG) PRICING IS USED, */
/*     CALCULATE THE NEW REDUCED COSTS. */
    goto L20040;

/*     COMPUTE THE UPDATED DUALS IN DUALS(*). */
L20039:
    npr003 = 1;
    npr003_fmt = fmt_20074;
    goto L30003;
L20074:
    i__1 = *nvars__ + *mrelas;
    scopy_(&i__1, &zero, &c__0, &rz[1], &c__1);
    nnegrc = 0;
    j = *jstrt;
    pagepl = TRUE_;

L20075:
    if (! (ibb[j] <= 0)) {
	goto L20077;
    }
    pagepl = TRUE_;
    goto L20078;

/*     NONBASIC INDEPENDENT VARIABLE (COLUMN IN SPARSE MATRIX STORAGE) */
L20077:
    if (! (j <= *nvars__)) {
	goto L20080;
    }
    rz[j] = costs[j] * *costsc;
    if (! (j == 1)) {
	goto L20083;
    }
    ilow = *nvars__ + 5;
    goto L20084;
L20083:
    ilow = imat[j + 3] + 1;
L20084:
    if (! pagepl) {
	goto L20086;
    }
    il1 = iploc_(&ilow, &amat[1], &imat[1]);
    if (! (il1 >= *lmx - 1)) {
	goto L20089;
    }
    ilow += 2;
    il1 = iploc_(&ilow, &amat[1], &imat[1]);
L20089:
    ipage = (i__1 = imat[*lmx - 1], abs(i__1));
    goto L20087;
L20086:
    il1 = ihi + 1;
L20087:
    ihi = imat[j + 4] - (ilow - il1);
L20092:
/* Computing MIN */
    i__1 = *lmx - 2;
    iu1 = min(i__1,ihi);
    if (! (iu1 >= il1 && (iu1 - il1) % 2 == 0)) {
	goto L20094;
    }
    rz[j] -= amat[il1] * duals[imat[il1]];
    ++il1;
L20094:
    if (! (il1 > iu1)) {
	goto L20097;
    }
    goto L20093;
L20097:

/*     UNROLL THE DOT PRODUCT LOOP TO A DEPTH OF TWO.  (THIS IS DONE */
/*     FOR INCREASED EFFICIENCY). */
    i__1 = iu1;
    for (i__ = il1; i__ <= i__1; i__ += 2) {
	rz[j] = rz[j] - amat[i__] * duals[imat[i__]] - amat[i__ + 1] * duals[
		imat[i__ + 1]];
/* L40: */
    }
    if (! (ihi <= *lmx - 2)) {
	goto L20100;
    }
    goto L20093;
L20100:
    ++ipage;
    key = 1;
    prwpge_(&key, &ipage, &lpg, &amat[1], &imat[1]);
    il1 = *nvars__ + 5;
    ihi -= lpg;
    goto L20092;
L20093:
    pagepl = ihi == *lmx - 2;
    rz[j] *= csc[j];

/*     NONBASIC DEPENDENT VARIABLES (COLUMNS DEFINED IMPLICITLY) */
    goto L20081;
L20080:
    pagepl = TRUE_;
    scalr = -one;
    if (ind[j] == 2) {
	scalr = one;
    }
    i__ = j - *nvars__;
    rz[j] = -scalr * duals[i__];
L20081:
L20078:

    rcost = rz[j];
    if (ibb[j] % 2 == 0) {
	rcost = -rcost;
    }
    if (! (ind[j] == 3)) {
	goto L20103;
    }
    if (bu[j] == bl[j]) {
	rcost = zero;
    }
L20103:
    if (ind[j] == 4) {
	rcost = -dabs(rcost);
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
	goto L20106;
    }
    goto L20076;
L20106:
    goto L20075;
L20076:
    *jstrt = j;
L20040:
    goto L20030;

/*     THIS IS NECESSARY ONLY FOR PRINTING OF INTERMEDIATE RESULTS. */
L20029:
    npr003 = 2;
    npr003_fmt = fmt_20109;
    goto L30003;
L20109:
L20030:
    return 0;
/*     PROCEDURE (TRANSLATE RIGHT HAND SIDE) */

/*     PERFORM THE TRANSLATION ON THE RIGHT-HAND SIDE. */
L30001:
    if (! (ibas <= *nvars__)) {
	goto L20110;
    }
    i__ = 0;
L20113:
    pnnzrs_(&i__, &aij, &iplace, &amat[1], &imat[1], &ibas);
    if (! (i__ <= 0)) {
	goto L20115;
    }
    goto L20114;
L20115:
    rhs[i__] -= scalr * aij * csc[ibas];
    goto L20113;
L20114:
    goto L20111;
L20110:
    i__ = ibas - *nvars__;
    if (! (ind[ibas] == 2)) {
	goto L20118;
    }
    rhs[i__] -= scalr;
    goto L20119;
L20118:
    rhs[i__] += scalr;
L20119:
L20111:
/* Computing MAX */
    r__1 = *rhsnrm, r__2 = sasum_(mrelas, &rhs[1], &c__1);
    *rhsnrm = dmax(r__1,r__2);
    switch (npr001) {
	case 0: goto L20009;
	case 1: goto L20013;
	case 2: goto L20017;
	case 3: goto L20028;
    }
/*     PROCEDURE (COMPUTE NEW PRIMAL) */

/*     COPY RHS INTO WW(*), SOLVE SYSTEM. */
L30002:
    scopy_(mrelas, &rhs[1], &c__1, &ww[1], &c__1);
    trans = FALSE_;
    la05bs_(&basmat[1], &ibrc[ibrc_offset], lbm, mrelas, &ipr[1], &iwr[1], &
	    wr[1], gg, &ww[1], &trans);
    scopy_(mrelas, &ww[1], &c__1, &rprim[1], &c__1);
    *rprnrm = sasum_(mrelas, &rprim[1], &c__1);
    goto L20038;
/*     PROCEDURE (COMPUTE NEW DUALS) */

/*     SOLVE FOR DUAL VARIABLES. FIRST COPY COSTS INTO DUALS(*). */
L30003:
    i__ = 1;
    n20121 = *mrelas;
    goto L20122;
L20121:
    ++i__;
L20122:
    if (n20121 - i__ < 0) {
	goto L20123;
    }
    j = ibasis[i__];
    if (! (j <= *nvars__)) {
	goto L20125;
    }
    duals[i__] = *costsc * costs[j] * csc[j] + *xlamda * primal[i__ + *
	    nvars__];
    goto L20126;
L20125:
    duals[i__] = *xlamda * primal[i__ + *nvars__];
L20126:
    goto L20121;

L20123:
    trans = TRUE_;
    la05bs_(&basmat[1], &ibrc[ibrc_offset], lbm, mrelas, &ipr[1], &iwr[1], &
	    wr[1], gg, &duals[1], &trans);
    *dulnrm = sasum_(mrelas, &duals[1], &c__1);
    switch (npr003) {
	case 0: goto L20042;
	case 1: goto L20074;
	case 2: goto L20109;
    }
    abort();
    return 0;
} /* splpmu_ */

