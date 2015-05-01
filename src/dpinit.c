/* dpinit.f -- translated by f2c (version 12.02.01).
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

/* DECK DPINIT */
/* Subroutine */ int dpinit_(integer *mrelas, integer *nvars__, doublereal *
	costs, doublereal *bl, doublereal *bu, integer *ind, doublereal *
	primal, integer *info, doublereal *amat, doublereal *csc, doublereal *
	costsc, doublereal *colnrm, doublereal *xlamda, doublereal *anorm, 
	doublereal *rhs, doublereal *rhsnrm, integer *ibasis, integer *ibb, 
	integer *imat, logical *lopt)
{
    /* System generated locals */
    doublereal d__1, d__2, d__3;

    /* Local variables */
    static integer i__, j, ip;
    static doublereal aij, one;
    static integer n20041, n20007, n20070, n20019, n20028, n20056, n20066, 
	    n20074, n20078;
    static doublereal cmax, csum, zero, scalr;
    extern doublereal dasum_(integer *, doublereal *, integer *);
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *);
    static integer iplace;
    static logical colscp, minprb, contin, usrbas, cstscp;
    static doublereal testsc;
    extern /* Subroutine */ int dpnnzr_(integer *, doublereal *, integer *, 
	    doublereal *, integer *, integer *);

/* ***BEGIN PROLOGUE  DPINIT */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to DSPLP */
/* ***LIBRARY   SLATEC */
/* ***TYPE      DOUBLE PRECISION (SPINIT-S, DPINIT-D) */
/* ***AUTHOR  (UNKNOWN) */
/* ***DESCRIPTION */

/*     THE EDITING REQUIRED TO CONVERT THIS SUBROUTINE FROM SINGLE TO */
/*     DOUBLE PRECISION INVOLVES THE FOLLOWING CHARACTER STRING CHANGES. */

/*     USE AN EDITING COMMAND (CHANGE) /STRING-1/(TO)STRING-2/. */
/*     /REAL (12 BLANKS)/DOUBLE PRECISION/,/SCOPY/DCOPY/ */
/*     REVISED 810519-0900 */
/*     REVISED YYMMDD-HHMM */

/*     INITIALIZATION SUBROUTINE FOR DSPLP(*) PACKAGE. */

/* ***SEE ALSO  DSPLP */
/* ***ROUTINES CALLED  DASUM, DCOPY, DPNNZR */
/* ***REVISION HISTORY  (YYMMDD) */
/*   811215  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890605  Removed unreferenced labels.  (WRB) */
/*   891009  Removed unreferenced variable.  (WRB) */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900328  Added TYPE section.  (WRB) */
/* ***END PROLOGUE  DPINIT */

/* ***FIRST EXECUTABLE STATEMENT  DPINIT */
    /* Parameter adjustments */
    --lopt;
    --imat;
    --ibb;
    --ibasis;
    --rhs;
    --colnrm;
    --csc;
    --amat;
    --primal;
    --ind;
    --bu;
    --bl;
    --costs;

    /* Function Body */
    zero = 0.;
    one = 1.;
    contin = lopt[1];
    usrbas = lopt[2];
    colscp = lopt[5];
    cstscp = lopt[6];
    minprb = lopt[7];

/*     SCALE DATA. NORMALIZE BOUNDS. FORM COLUMN CHECK SUMS. */
    goto L30001;

/*     INITIALIZE ACTIVE BASIS MATRIX. */
L20002:
    goto L30002;
L20003:
    return 0;

/*     PROCEDURE (SCALE DATA. NORMALIZE BOUNDS. FORM COLUMN CHECK SUMS) */

/*     DO COLUMN SCALING IF NOT PROVIDED BY THE USER. */
L30001:
    if (colscp) {
	goto L20004;
    }
    j = 1;
    n20007 = *nvars__;
    goto L20008;
L20007:
    ++j;
L20008:
    if (n20007 - j < 0) {
	goto L20009;
    }
    cmax = zero;
    i__ = 0;
L20011:
    dpnnzr_(&i__, &aij, &iplace, &amat[1], &imat[1], &j);
    if (! (i__ == 0)) {
	goto L20013;
    }
    goto L20012;
L20013:
/* Computing MAX */
    d__1 = cmax, d__2 = abs(aij);
    cmax = max(d__1,d__2);
    goto L20011;
L20012:
    if (! (cmax == zero)) {
	goto L20016;
    }
    csc[j] = one;
    goto L20017;
L20016:
    csc[j] = one / cmax;
L20017:
    goto L20007;
L20009:

/*     FORM CHECK SUMS OF COLUMNS. COMPUTE MATRIX NORM OF SCALED MATRIX. */
L20004:
    *anorm = zero;
    j = 1;
    n20019 = *nvars__;
    goto L20020;
L20019:
    ++j;
L20020:
    if (n20019 - j < 0) {
	goto L20021;
    }
    primal[j] = zero;
    csum = zero;
    i__ = 0;
L20023:
    dpnnzr_(&i__, &aij, &iplace, &amat[1], &imat[1], &j);
    if (! (i__ <= 0)) {
	goto L20025;
    }
    goto L20024;
L20025:
    primal[j] += aij;
    csum += abs(aij);
    goto L20023;
L20024:
    if (ind[j] == 2) {
	csc[j] = -csc[j];
    }
    primal[j] *= csc[j];
    colnrm[j] = (d__1 = csc[j] * csum, abs(d__1));
/* Computing MAX */
    d__1 = *anorm, d__2 = colnrm[j];
    *anorm = max(d__1,d__2);
    goto L20019;

/*     IF THE USER HAS NOT PROVIDED COST VECTOR SCALING THEN SCALE IT */
/*     USING THE MAX. NORM OF THE TRANSFORMED COST VECTOR, IF NONZERO. */
L20021:
    testsc = zero;
    j = 1;
    n20028 = *nvars__;
    goto L20029;
L20028:
    ++j;
L20029:
    if (n20028 - j < 0) {
	goto L20030;
    }
/* Computing MAX */
    d__2 = testsc, d__3 = (d__1 = csc[j] * costs[j], abs(d__1));
    testsc = max(d__2,d__3);
    goto L20028;
L20030:
    if (cstscp) {
	goto L20032;
    }
    if (! (testsc > zero)) {
	goto L20035;
    }
    *costsc = one / testsc;
    goto L20036;
L20035:
    *costsc = one;
L20036:
L20032:
    *xlamda = (*costsc + *costsc) * testsc;
    if (*xlamda == zero) {
	*xlamda = one;
    }

/*     IF MAXIMIZATION PROBLEM, THEN CHANGE SIGN OF COSTSC AND LAMDA */
/*     =WEIGHT FOR PENALTY-FEASIBILITY METHOD. */
    if (minprb) {
	goto L20038;
    }
    *costsc = -(*costsc);
L20038:
    goto L20002;
/* :CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC */
/*     PROCEDURE (INITIALIZE RHS(*),IBASIS(*), AND IBB(*)) */

/*     INITIALLY SET RIGHT-HAND SIDE VECTOR TO ZERO. */
L30002:
    dcopy_(mrelas, &zero, &c__0, &rhs[1], &c__1);

/*     TRANSLATE RHS ACCORDING TO CLASSIFICATION OF INDEPENDENT VARIABLES */
    j = 1;
    n20041 = *nvars__;
    goto L20042;
L20041:
    ++j;
L20042:
    if (n20041 - j < 0) {
	goto L20043;
    }
    if (! (ind[j] == 1)) {
	goto L20045;
    }
    scalr = -bl[j];
    goto L20046;
L20045:
    if (! (ind[j] == 2)) {
	goto L10001;
    }
    scalr = -bu[j];
    goto L20046;
L10001:
    if (! (ind[j] == 3)) {
	goto L10002;
    }
    scalr = -bl[j];
    goto L20046;
L10002:
    if (! (ind[j] == 4)) {
	goto L10003;
    }
    scalr = zero;
L10003:
L20046:
    if (! (scalr != zero)) {
	goto L20048;
    }
    i__ = 0;
L20051:
    dpnnzr_(&i__, &aij, &iplace, &amat[1], &imat[1], &j);
    if (! (i__ <= 0)) {
	goto L20053;
    }
    goto L20052;
L20053:
    rhs[i__] = scalr * aij + rhs[i__];
    goto L20051;
L20052:
L20048:
    goto L20041;

/*     TRANSLATE RHS ACCORDING TO CLASSIFICATION OF DEPENDENT VARIABLES. */
L20043:
    i__ = *nvars__ + 1;
    n20056 = *nvars__ + *mrelas;
    goto L20057;
L20056:
    ++i__;
L20057:
    if (n20056 - i__ < 0) {
	goto L20058;
    }
    if (! (ind[i__] == 1)) {
	goto L20060;
    }
    scalr = bl[i__];
    goto L20061;
L20060:
    if (! (ind[i__] == 2)) {
	goto L10004;
    }
    scalr = bu[i__];
    goto L20061;
L10004:
    if (! (ind[i__] == 3)) {
	goto L10005;
    }
    scalr = bl[i__];
    goto L20061;
L10005:
    if (! (ind[i__] == 4)) {
	goto L10006;
    }
    scalr = zero;
L10006:
L20061:
    rhs[i__ - *nvars__] += scalr;
    goto L20056;
L20058:
    *rhsnrm = dasum_(mrelas, &rhs[1], &c__1);

/*     IF THIS IS NOT A CONTINUATION OR THE USER HAS NOT PROVIDED THE */
/*     INITIAL BASIS, THEN THE INITIAL BASIS IS COMPRISED OF THE */
/*     DEPENDENT VARIABLES. */
    if (contin || usrbas) {
	goto L20063;
    }
    j = 1;
    n20066 = *mrelas;
    goto L20067;
L20066:
    ++j;
L20067:
    if (n20066 - j < 0) {
	goto L20068;
    }
    ibasis[j] = *nvars__ + j;
    goto L20066;
L20068:

/*     DEFINE THE ARRAY IBB(*) */
L20063:
    j = 1;
    n20070 = *nvars__ + *mrelas;
    goto L20071;
L20070:
    ++j;
L20071:
    if (n20070 - j < 0) {
	goto L20072;
    }
    ibb[j] = 1;
    goto L20070;
L20072:
    j = 1;
    n20074 = *mrelas;
    goto L20075;
L20074:
    ++j;
L20075:
    if (n20074 - j < 0) {
	goto L20076;
    }
    ibb[ibasis[j]] = -1;
    goto L20074;

/*     DEFINE THE REST OF IBASIS(*) */
L20076:
    ip = *mrelas;
    j = 1;
    n20078 = *nvars__ + *mrelas;
    goto L20079;
L20078:
    ++j;
L20079:
    if (n20078 - j < 0) {
	goto L20080;
    }
    if (! (ibb[j] > 0)) {
	goto L20082;
    }
    ++ip;
    ibasis[ip] = j;
L20082:
    goto L20078;
L20080:
    goto L20003;
} /* dpinit_ */

