/* dplpfl.f -- translated by f2c (version 12.02.01).
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

/* DECK DPLPFL */
/* Subroutine */ int dplpfl_(integer *mrelas, integer *nvars__, integer *
	ienter, integer *ileave, integer *ibasis, integer *ind, integer *ibb, 
	doublereal *theta, doublereal *dirnrm, doublereal *rprnrm, doublereal 
	*csc, doublereal *ww, doublereal *bl, doublereal *bu, doublereal *erp,
	 doublereal *rprim, doublereal *primal, logical *finite, logical *
	zerolv)
{
    /* System generated locals */
    doublereal d__1;

    /* Local variables */
    static integer i__, j, n20005, n20036;
    static doublereal zero, bound, ratio;

/* ***BEGIN PROLOGUE  DPLPFL */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to DSPLP */
/* ***LIBRARY   SLATEC */
/* ***TYPE      DOUBLE PRECISION (SPLPFL-S, DPLPFL-D) */
/* ***AUTHOR  (UNKNOWN) */
/* ***DESCRIPTION */

/*     THE EDITING REQUIRED TO CONVERT THIS SUBROUTINE FROM SINGLE TO */
/*     DOUBLE PRECISION INVOLVES THE FOLLOWING CHARACTER STRING CHANGES. */

/*     USE AN EDITING COMMAND (CHANGE) /STRING-1/(TO)STRING-2/. */
/*     /REAL (12 BLANKS)/DOUBLE PRECISION/. */

/*     THIS SUBPROGRAM IS PART OF THE DSPLP( ) PACKAGE. */
/*     IT IMPLEMENTS THE PROCEDURE (CHOOSE VARIABLE TO LEAVE BASIS). */
/*     REVISED 811130-1045 */
/*     REVISED YYMMDD-HHMM */

/* ***SEE ALSO  DSPLP */
/* ***ROUTINES CALLED  (NONE) */
/* ***REVISION HISTORY  (YYMMDD) */
/*   811215  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890605  Removed unreferenced labels.  (WRB) */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900328  Added TYPE section.  (WRB) */
/* ***END PROLOGUE  DPLPFL */
/* ***FIRST EXECUTABLE STATEMENT  DPLPFL */
    /* Parameter adjustments */
    --primal;
    --rprim;
    --erp;
    --bu;
    --bl;
    --ww;
    --csc;
    --ibb;
    --ind;
    --ibasis;

    /* Function Body */
    zero = 0.;

/*     SEE IF THE ENTERING VARIABLE IS RESTRICTING THE STEP LENGTH */
/*     BECAUSE OF AN UPPER BOUND. */
    *finite = FALSE_;
    j = ibasis[*ienter];
    if (! (ind[j] == 3)) {
	goto L20002;
    }
    *theta = bu[j] - bl[j];
    if (j <= *nvars__) {
	*theta /= csc[j];
    }
    *finite = TRUE_;
    *ileave = *ienter;

/*     NOW USE THE BASIC VARIABLES TO POSSIBLY RESTRICT THE STEP */
/*     LENGTH EVEN FURTHER. */
L20002:
    i__ = 1;
    n20005 = *mrelas;
    goto L20006;
L20005:
    ++i__;
L20006:
    if (n20005 - i__ < 0) {
	goto L20007;
    }
    j = ibasis[i__];

/*     IF THIS IS A FREE VARIABLE, DO NOT USE IT TO */
/*     RESTRICT THE STEP LENGTH. */
    if (! (ind[j] == 4)) {
	goto L20009;
    }
    goto L20005;

/*     IF DIRECTION COMPONENT IS ABOUT ZERO, IGNORE IT FOR COMPUTING */
/*     THE STEP LENGTH. */
L20009:
    if (! ((d__1 = ww[i__], abs(d__1)) <= *dirnrm * erp[i__])) {
	goto L20012;
    }
    goto L20005;
L20012:
    if (! (ww[i__] > zero)) {
	goto L20015;
    }

/*     IF RPRIM(I) IS ESSENTIALLY ZERO, SET RATIO TO ZERO AND EXIT LOOP. */
    if (! ((d__1 = rprim[i__], abs(d__1)) <= *rprnrm * erp[i__])) {
	goto L20018;
    }
    *theta = zero;
    *ileave = i__;
    *finite = TRUE_;
    goto L20008;

/*     THE VALUE OF RPRIM(I) WILL DECREASE ONLY TO ITS LOWER BOUND OR */
/*     ONLY TO ITS UPPER BOUND.  IF IT DECREASES TO ITS */
/*     UPPER BOUND, THEN RPRIM(I) HAS ALREADY BEEN TRANSLATED */
/*     TO ITS UPPER BOUND AND NOTHING NEEDS TO BE DONE TO IBB(J). */
L20018:
    if (! (rprim[i__] > zero)) {
	goto L10001;
    }
    ratio = rprim[i__] / ww[i__];
    if (*finite) {
	goto L20021;
    }
    *ileave = i__;
    *theta = ratio;
    *finite = TRUE_;
    goto L20022;
L20021:
    if (! (ratio < *theta)) {
	goto L10002;
    }
    *ileave = i__;
    *theta = ratio;
L10002:
L20022:
    goto L20019;

/*     THE VALUE RPRIM(I).LT.ZERO WILL NOT RESTRICT THE STEP. */
L10001:

/*     THE DIRECTION COMPONENT IS NEGATIVE, THEREFORE THE VARIABLE WILL */
/*     INCREASE. */
L20019:
    goto L20016;

/*     IF THE VARIABLE IS LESS THAN ITS LOWER BOUND, IT CAN */
/*     INCREASE ONLY TO ITS LOWER BOUND. */
L20015:
    if (! (primal[i__ + *nvars__] < zero)) {
	goto L20024;
    }
    ratio = rprim[i__] / ww[i__];
    if (ratio < zero) {
	ratio = zero;
    }
    if (*finite) {
	goto L20027;
    }
    *ileave = i__;
    *theta = ratio;
    *finite = TRUE_;
    goto L20028;
L20027:
    if (! (ratio < *theta)) {
	goto L10003;
    }
    *ileave = i__;
    *theta = ratio;
L10003:
L20028:

/*     IF THE BASIC VARIABLE IS FEASIBLE AND IS NOT AT ITS UPPER BOUND, */
/*     THEN IT CAN INCREASE TO ITS UPPER BOUND. */
    goto L20025;
L20024:
    if (! (ind[j] == 3 && primal[i__ + *nvars__] == zero)) {
	goto L10004;
    }
    bound = bu[j] - bl[j];
    if (j <= *nvars__) {
	bound /= csc[j];
    }
    ratio = (bound - rprim[i__]) / (-ww[i__]);
    if (*finite) {
	goto L20030;
    }
    *ileave = -i__;
    *theta = ratio;
    *finite = TRUE_;
    goto L20031;
L20030:
    if (! (ratio < *theta)) {
	goto L10005;
    }
    *ileave = -i__;
    *theta = ratio;
L10005:
L20031:
L10004:
L20025:
L20016:
    goto L20005;
L20007:

/*     IF STEP LENGTH IS FINITE, SEE IF STEP LENGTH IS ABOUT ZERO. */
L20008:
    if (! (*finite)) {
	goto L20033;
    }
    *zerolv = TRUE_;
    i__ = 1;
    n20036 = *mrelas;
    goto L20037;
L20036:
    ++i__;
L20037:
    if (n20036 - i__ < 0) {
	goto L20038;
    }
    *zerolv = *zerolv && (d__1 = *theta * ww[i__], abs(d__1)) <= erp[i__] * *
	    rprnrm;
    if (*zerolv) {
	goto L20040;
    }
    goto L20039;
L20040:
    goto L20036;
L20038:
L20039:
L20033:
    return 0;
} /* dplpfl_ */

