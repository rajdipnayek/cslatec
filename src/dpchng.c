/* dpchng.f -- translated by f2c (version 12.02.01).
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

/* DECK DPCHNG */
/* Subroutine */ int dpchng_(integer *ii, doublereal *xval, integer *iplace, 
	doublereal *sx, integer *ix, integer *ircx)
{
    /* Initialized data */

    static doublereal zero = 0.;
    static doublereal one = 1.;

    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer i__, j, k, jj, il, ll, np, lpg, ipl, key, lmx, n20055, 
	    iend, nerr, iopt;
    extern integer idloc_(integer *, doublereal *, integer *);
    static integer ilast;
    static doublereal sxval;
    extern /* Subroutine */ int dprwpg_(integer *, integer *, integer *, 
	    doublereal *, integer *);
    static integer ixlast;
    extern /* Subroutine */ int xermsg_(char *, char *, char *, integer *, 
	    integer *, ftnlen, ftnlen, ftnlen);
    static integer istart, jstart;
    static doublereal sxlast;

/* ***BEGIN PROLOGUE  DPCHNG */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to DSPLP */
/* ***LIBRARY   SLATEC */
/* ***TYPE      DOUBLE PRECISION (PCHNGS-S, DPCHNG-D) */
/* ***AUTHOR  Hanson, R. J., (SNLA) */
/*           Wisniewski, J. A., (SNLA) */
/* ***DESCRIPTION */

/*     SUBROUTINE DPCHNG CHANGES ELEMENT II IN VECTOR +/- IRCX TO THE */
/*     VALUE XVAL. */
/*     DPCHNG LIMITS THE TYPE OF STORAGE TO A SEQUENTIAL SCHEME. */
/*     SPARSE MATRIX ELEMENT ALTERATION SUBROUTINE. */

/*            II THE ABSOLUTE VALUE OF THIS INTEGER IS THE SUBSCRIPT FOR */
/*               THE ELEMENT TO BE CHANGED. */
/*          XVAL NEW VALUE OF THE MATRIX ELEMENT BEING CHANGED. */
/*     IPLACE POINTER INFORMATION WHICH IS MAINTAINED BY THE PACKAGE. */
/*   SX(*),IX(*) THE WORK ARRAYS WHICH ARE USED TO STORE THE SPARSE */
/*               MATRIX. THESE ARRAYS ARE AUTOMATICALLY MAINTAINED BY THE */
/*               PACKAGE FOR THE USER. */
/*          IRCX POINTS TO THE VECTOR OF THE MATRIX BEING UPDATED. */
/*               A NEGATIVE VALUE OF IRCX INDICATES THAT ROW -IRCX IS */
/*               BEING UPDATED.  A POSITIVE VALUE OF IRCX INDICATES THAT */
/*               COLUMN IRCX IS BEING UPDATED.  A ZERO VALUE OF IRCX IS */
/*               AN ERROR. */

/*     SINCE DATA ITEMS ARE KEPT SORTED IN THE SEQUENTIAL DATA STRUCTURE, */
/*     CHANGING A MATRIX ELEMENT CAN REQUIRE THE MOVEMENT OF ALL THE DATA */
/*     ITEMS IN THE MATRIX. FOR THIS REASON, IT IS SUGGESTED THAT DATA */
/*     ITEMS BE ADDED A COL. AT A TIME, IN ASCENDING COL. SEQUENCE. */
/*     FURTHERMORE, SINCE DELETING ITEMS FROM THE DATA STRUCTURE MAY ALSO */
/*     REQUIRE MOVING LARGE AMOUNTS OF DATA, ZERO ELEMENTS ARE EXPLICITLY */
/*     STORED IN THE MATRIX. */

/*     THIS SUBROUTINE IS A MODIFICATION OF THE SUBROUTINE LCHNGS, */
/*     SANDIA LABS. REPT. SAND78-0785. */
/*     MODIFICATIONS BY K.L. HIEBERT AND R.J. HANSON */
/*     REVISED 811130-1000 */
/*     REVISED YYMMDD-HHMM */

/* ***SEE ALSO  DSPLP */
/* ***ROUTINES CALLED  DPRWPG, IDLOC, XERMSG */
/* ***REVISION HISTORY  (YYMMDD) */
/*   811215  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890606  Changed references from IPLOC to IDLOC.  (WRB) */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ) */
/*   900328  Added TYPE section.  (WRB) */
/*   910403  Updated AUTHOR and DESCRIPTION sections.  (WRB) */
/* ***END PROLOGUE  DPCHNG */
    /* Parameter adjustments */
    --ix;
    --sx;

    /* Function Body */
/* ***FIRST EXECUTABLE STATEMENT  DPCHNG */
    iopt = 1;

/*     DETERMINE NULL-CASES.. */
    if (*ii == 0) {
	return 0;
    }

/*     CHECK VALIDITY OF ROW/COL. INDEX. */

    if (! (*ircx == 0)) {
	goto L20002;
    }
    nerr = 55;
    xermsg_("SLATEC", "DPCHNG", "IRCX=0", &nerr, &iopt, (ftnlen)6, (ftnlen)6, 
	    (ftnlen)6);
L20002:
    lmx = ix[1];

/*     LMX IS THE LENGTH OF THE IN-MEMORY STORAGE AREA. */

    if (! (*ircx < 0)) {
	goto L20005;
    }

/*     CHECK SUBSCRIPTS OF THE ROW. THE ROW NUMBER MUST BE .LE. M AND */
/*     THE INDEX MUST BE .LE. N. */

    if (! (ix[2] < -(*ircx) || ix[3] < abs(*ii))) {
	goto L20008;
    }
    nerr = 55;
    xermsg_("SLATEC", "DPCHNG", "SUBSCRIPTS FOR ARRAY ELEMENT TO BE ACCESSED"
	    " WERE OUT OF BOUNDS", &nerr, &iopt, (ftnlen)6, (ftnlen)6, (ftnlen)
	    62);
L20008:
    goto L20006;

/*     CHECK SUBSCRIPTS OF THE COLUMN. THE COL. NUMBER MUST BE .LE. N AND */
/*     THE INDEX MUST BE .LE. M. */

L20005:
    if (! (ix[3] < *ircx || ix[2] < abs(*ii))) {
	goto L20011;
    }
    nerr = 55;
    xermsg_("SLATEC", "DPCHNG", "SUBSCRIPTS FOR ARRAY ELEMENT TO BE ACCESSED"
	    " WERE OUT OF BOUNDS", &nerr, &iopt, (ftnlen)6, (ftnlen)6, (ftnlen)
	    62);
L20011:

/*     SET I TO BE THE ELEMENT OF ROW/COLUMN J TO BE CHANGED. */

L20006:
    if (! (*ircx > 0)) {
	goto L20014;
    }
    i__ = abs(*ii);
    j = abs(*ircx);
    goto L20015;
L20014:
    i__ = abs(*ircx);
    j = abs(*ii);

/*     THE INTEGER LL POINTS TO THE START OF THE MATRIX ELEMENT DATA. */

L20015:
    ll = ix[3] + 4;
    *ii = abs(*ii);
    lpg = lmx - ll;

/*     SET IPLACE TO START OUR SCAN FOR THE ELEMENT AT THE BEGINNING */
/*     OF THE VECTOR. */

    if (! (j == 1)) {
	goto L20017;
    }
    *iplace = ll + 1;
    goto L20018;
L20017:
    *iplace = ix[j + 3] + 1;

/*     IEND POINTS TO THE LAST ELEMENT OF THE VECTOR TO BE SCANNED. */

L20018:
    iend = ix[j + 4];

/*     SCAN THROUGH SEVERAL PAGES, IF NECESSARY, TO FIND MATRIX ELEMENT. */

    ipl = idloc_(iplace, &sx[1], &ix[1]);
    np = (i__1 = ix[lmx - 1], abs(i__1));
    goto L20021;
L20020:
    if (ilast == iend) {
	goto L20022;
    }

/*     THE VIRTUAL END OF DATA FOR THIS PAGE IS ILAST. */

L20021:
/* Computing MIN */
    i__1 = iend, i__2 = np * lpg + ll - 2;
    ilast = min(i__1,i__2);

/*     THE RELATIVE END OF DATA FOR THIS PAGE IS IL. */
/*     SEARCH FOR A MATRIX VALUE WITH AN INDEX .GE. I ON THE PRESENT */
/*     PAGE. */

    il = idloc_(&ilast, &sx[1], &ix[1]);
/* Computing MIN */
    i__1 = il, i__2 = lmx - 2;
    il = min(i__1,i__2);
L20023:
    if (ipl >= il || ix[ipl] >= i__) {
	goto L20024;
    }
    ++ipl;
    goto L20023;

/*     SET IPLACE AND STORE DATA ITEM IF FOUND. */

L20024:
    if (! (ix[ipl] == i__ && ipl <= il)) {
	goto L20025;
    }
    sx[ipl] = *xval;
    sx[lmx] = one;
    return 0;

/*     EXIT FROM LOOP IF ITEM WAS FOUND. */

L20025:
    if (ix[ipl] > i__ && ipl <= il) {
	ilast = iend;
    }
    if (! (ilast != iend)) {
	goto L20028;
    }
    ipl = ll + 1;
    ++np;
L20028:
    goto L20020;

/*     INSERT NEW DATA ITEM INTO LOCATION AT IPLACE(IPL). */

L20022:
    if (! (ipl > il || ipl == il && i__ > ix[ipl])) {
	goto L20031;
    }
    ipl = il + 1;
    if (ipl == lmx - 1) {
	ipl += 2;
    }
L20031:
    *iplace = (np - 1) * lpg + ipl;

/*     GO TO A NEW PAGE, IF NECESSARY, TO INSERT THE ITEM. */

    if (! (ipl <= lmx || ix[lmx - 1] >= 0)) {
	goto L20034;
    }
    ipl = idloc_(iplace, &sx[1], &ix[1]);
L20034:
    iend = ix[ll];
    np = (i__1 = ix[lmx - 1], abs(i__1));
    sxval = *xval;

/*     LOOP THROUGH ALL SUBSEQUENT PAGES OF THE MATRIX MOVING DATA DOWN. */
/*     THIS IS NECESSARY TO MAKE ROOM FOR THE NEW MATRIX ELEMENT AND */
/*     KEEP THE ENTRIES SORTED. */

    goto L20038;
L20037:
    if (ix[lmx - 1] <= 0) {
	goto L20039;
    }
L20038:
/* Computing MIN */
    i__1 = iend, i__2 = np * lpg + ll - 2;
    ilast = min(i__1,i__2);
    il = idloc_(&ilast, &sx[1], &ix[1]);
/* Computing MIN */
    i__1 = il, i__2 = lmx - 2;
    il = min(i__1,i__2);
    sxlast = sx[il];
    ixlast = ix[il];
    istart = ipl + 1;
    if (! (istart <= il)) {
	goto L20040;
    }
    k = istart + il;
    i__1 = il;
    for (jj = istart; jj <= i__1; ++jj) {
	sx[k - jj] = sx[k - jj - 1];
	ix[k - jj] = ix[k - jj - 1];
/* L50: */
    }
    sx[lmx] = one;
L20040:
    if (! (ipl <= lmx)) {
	goto L20043;
    }
    sx[ipl] = sxval;
    ix[ipl] = i__;
    sxval = sxlast;
    i__ = ixlast;
    sx[lmx] = one;
    if (! (ix[lmx - 1] > 0)) {
	goto L20046;
    }
    ipl = ll + 1;
    ++np;
L20046:
L20043:
    goto L20037;
L20039:
    np = (i__1 = ix[lmx - 1], abs(i__1));

/*     DETERMINE IF A NEW PAGE IS TO BE CREATED FOR THE LAST ELEMENT */
/*     MOVED DOWN. */

    ++il;
    if (! (il == lmx - 1)) {
	goto L20049;
    }

/*     CREATE A NEW PAGE. */

    ix[lmx - 1] = np;

/*     WRITE THE OLD PAGE. */

    sx[lmx] = zero;
    key = 2;
    dprwpg_(&key, &np, &lpg, &sx[1], &ix[1]);
    sx[lmx] = one;

/*     STORE LAST ELEMENT MOVED DOWN IN A NEW PAGE. */

    ipl = ll + 1;
    ++np;
    ix[lmx - 1] = -np;
    sx[ipl] = sxval;
    ix[ipl] = i__;
    goto L20050;

/*     LAST ELEMENT MOVED REMAINED ON THE OLD PAGE. */

L20049:
    if (! (ipl != il)) {
	goto L20052;
    }
    sx[il] = sxval;
    ix[il] = i__;
    sx[lmx] = one;
L20052:

/*     INCREMENT POINTERS TO LAST ELEMENT IN VECTORS J,J+1,... . */

L20050:
    jstart = j + 4;
    jj = jstart;
    n20055 = ll;
    goto L20056;
L20055:
    ++jj;
L20056:
    if (n20055 - jj < 0) {
	goto L20057;
    }
    ++ix[jj];
    if ((ix[jj] - ll) % lpg == lpg - 1) {
	ix[jj] += 2;
    }
    goto L20055;

/*     IPLACE POINTS TO THE INSERTED DATA ITEM. */

L20057:
    ipl = idloc_(iplace, &sx[1], &ix[1]);
    return 0;
} /* dpchng_ */

