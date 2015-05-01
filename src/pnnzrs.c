/* pnnzrs.f -- translated by f2c (version 12.02.01).
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

/* DECK PNNZRS */
/* Subroutine */ int pnnzrs_(integer *i__, real *xval, integer *iplace, real *
	sx, integer *ix, integer *ircx)
{
    /* Initialized data */

    static real zero = 0.f;

    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer j, l, i1, ii, il, ll, np, lpg, ipl, lmx, n20046, iend, 
	    nerr, iopt, idiff;
    extern integer iploc_(integer *, real *, integer *);
    static integer ilast, ipploc;
    extern /* Subroutine */ int xermsg_(char *, char *, char *, integer *, 
	    integer *, ftnlen, ftnlen, ftnlen);
    static integer istart;

/* ***BEGIN PROLOGUE  PNNZRS */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to SPLP */
/* ***LIBRARY   SLATEC */
/* ***TYPE      SINGLE PRECISION (PNNZRS-S, DPNNZR-D) */
/* ***AUTHOR  Hanson, R. J., (SNLA) */
/*           Wisniewski, J. A., (SNLA) */
/* ***DESCRIPTION */

/*     PNNZRS LIMITS THE TYPE OF STORAGE TO A SEQUENTIAL SCHEME. */
/*     SPARSE MATRIX NON ZERO RETRIEVAL SUBROUTINE. */

/*     SUBROUTINE PNNZRS() GETS THE NEXT NONZERO VALUE IN ROW OR COLUMN */
/*     +/- IRCX WITH AN INDEX GREATER THAN THE VALUE OF I. */

/*             I ABSOLUTE VALUE OF THIS SUBSCRIPT IS TO BE EXCEEDED */
/*               IN THE SEARCH FOR THE NEXT NONZERO VALUE. A NEGATIVE */
/*               OR ZERO VALUE OF I CAUSES THE SEARCH TO START AT */
/*               THE BEGINNING OF THE VECTOR.  A POSITIVE VALUE */
/*               OF I CAUSES THE SEARCH TO CONTINUE FROM THE LAST PLACE */
/*               ACCESSED. ON OUTPUT, THE ARGUMENT I */
/*               CONTAINS THE VALUE OF THE SUBSCRIPT FOUND.  AN OUTPUT */
/*               VALUE OF I EQUAL TO ZERO INDICATES THAT ALL COMPONENTS */
/*               WITH AN INDEX GREATER THAN THE INPUT VALUE OF I ARE */
/*               ZERO. */
/*          XVAL VALUE OF THE NONZERO ELEMENT FOUND.  ON OUTPUT, */
/*               XVAL=0. WHENEVER I=0. */
/*     IPLACE POINTER INFORMATION WHICH IS MAINTAINED BY THE PACKAGE. */
/*   SX(*),IX(*) THE WORK ARRAYS WHICH ARE USED TO STORE THE SPARSE */
/*               MATRIX.  THESE ARRAY CONTENTS ARE AUTOMATICALLY */
/*               MAINTAINED BY THE PACKAGE FOR THE USER. */
/*          IRCX POINTS TO THE VECTOR OF THE MATRIX BEING SCANNED.  A */
/*               NEGATIVE VALUE OF IRCX INDICATES THAT ROW -IRCX IS TO BE */
/*               SCANNED.  A POSITIVE VALUE OF IRCX INDICATES THAT */
/*               COLUMN IRCX IS TO BE SCANNED.  A ZERO VALUE OF IRCX IS */
/*               AN ERROR. */

/*     THIS SUBROUTINE IS A MODIFICATION OF THE SUBROUTINE LNNZRS, */
/*     SANDIA LABS. REPT. SAND78-0785. */
/*     MODIFICATIONS BY K.L. HIEBERT AND R.J. HANSON */
/*     REVISED 811130-1000 */
/*     REVISED YYMMDD-HHMM */

/* ***SEE ALSO  SPLP */
/* ***ROUTINES CALLED  IPLOC, XERMSG */
/* ***REVISION HISTORY  (YYMMDD) */
/*   811215  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890605  Removed unreferenced labels.  (WRB) */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ) */
/*   900328  Added TYPE section.  (WRB) */
/*   910403  Updated AUTHOR and DESCRIPTION sections.  (WRB) */
/* ***END PROLOGUE  PNNZRS */
    /* Parameter adjustments */
    --ix;
    --sx;

    /* Function Body */
/* ***FIRST EXECUTABLE STATEMENT  PNNZRS */
    iopt = 1;

/*     CHECK VALIDITY OF ROW/COL. INDEX. */

    if (! (*ircx == 0)) {
	goto L20002;
    }
    nerr = 55;
    xermsg_("SLATEC", "PNNZRS", "IRCX=0.", &nerr, &iopt, (ftnlen)6, (ftnlen)6,
	     (ftnlen)7);

/*     LMX IS THE LENGTH OF THE IN-MEMORY STORAGE AREA. */

L20002:
    lmx = ix[1];
    if (! (*ircx < 0)) {
	goto L20005;
    }

/*     CHECK SUBSCRIPTS OF THE ROW. THE ROW NUMBER MUST BE .LE. M AND */
/*     THE INDEX MUST BE .LE. N. */

    if (! (ix[2] < -(*ircx) || ix[3] < abs(*i__))) {
	goto L20008;
    }
    nerr = 55;
    xermsg_("SLATEC", "PNNZRS", "SUBSCRIPTS FOR ARRAY ELEMENT TO BE ACCESSED"
	    " WERE OUT OF BOUNDS.", &nerr, &iopt, (ftnlen)6, (ftnlen)6, (
	    ftnlen)63);
L20008:
    l = ix[3];
    goto L20006;

/*     CHECK SUBSCRIPTS OF THE COLUMN. THE COL. NUMBER MUST BE .LE. N AND */
/*     THE INDEX MUST BE .LE. M. */

L20005:
    if (! (*ircx > ix[3] || abs(*i__) > ix[2])) {
	goto L20011;
    }
    nerr = 55;
    xermsg_("SLATEC", "PNNZRS", "SUBSCRIPTS FOR ARRAY ELEMENT TO BE ACCESSED"
	    " WERE OUT OF BOUNDS.", &nerr, &iopt, (ftnlen)6, (ftnlen)6, (
	    ftnlen)63);
L20011:
    l = ix[2];

/*     HERE L IS THE LARGEST POSSIBLE SUBSCRIPT WITHIN THE VECTOR. */

L20006:
    j = abs(*ircx);
    ll = ix[3] + 4;
    lpg = lmx - ll;
    if (! (*ircx > 0)) {
	goto L20014;
    }

/*     SEARCHING FOR THE NEXT NONZERO IN A COLUMN. */

/*     INITIALIZE STARTING LOCATIONS.. */
    if (! (*i__ <= 0)) {
	goto L20017;
    }
    if (! (j == 1)) {
	goto L20020;
    }
    *iplace = ll + 1;
    goto L20021;
L20020:
    *iplace = ix[j + 3] + 1;
L20021:

/*     THE CASE I.LE.0 SIGNALS THAT THE SCAN FOR THE ENTRY */
/*     IS TO BEGIN AT THE START OF THE VECTOR. */

L20017:
    *i__ = abs(*i__);
    if (! (j == 1)) {
	goto L20023;
    }
    istart = ll + 1;
    goto L20024;
L20023:
    istart = ix[j + 3] + 1;
L20024:
    iend = ix[j + 4];

/*     VALIDATE IPLACE. SET TO START OF VECTOR IF OUT OF RANGE. */

    if (! (istart > *iplace || *iplace > iend)) {
	goto L20026;
    }
    if (! (j == 1)) {
	goto L20029;
    }
    *iplace = ll + 1;
    goto L20030;
L20029:
    *iplace = ix[j + 3] + 1;
L20030:

/*     SCAN THROUGH SEVERAL PAGES, IF NECESSARY, TO FIND MATRIX ENTRY. */

L20026:
    ipl = iploc_(iplace, &sx[1], &ix[1]);

/*     FIX UP IPLACE AND IPL IF THEY POINT TO PAGING DATA. */
/*     THIS IS NECESSARY BECAUSE THERE IS CONTROL INFORMATION AT THE */
/*     END OF EACH PAGE. */

    idiff = lmx - ipl;
    if (! (idiff <= 1 && ix[lmx - 1] > 0)) {
	goto L20032;
    }

/*     UPDATE THE RELATIVE ADDRESS IN A NEW PAGE. */

    *iplace = *iplace + idiff + 1;
    ipl = iploc_(iplace, &sx[1], &ix[1]);
L20032:
    np = (i__1 = ix[lmx - 1], abs(i__1));
    goto L20036;
L20035:
    if (ilast == iend) {
	goto L20037;
    }
L20036:
/* Computing MIN */
    i__1 = iend, i__2 = np * lpg + ll - 2;
    ilast = min(i__1,i__2);

/*     THE VIRTUAL END OF THE DATA FOR THIS PAGE IS ILAST. */

    il = iploc_(&ilast, &sx[1], &ix[1]);
/* Computing MIN */
    i__1 = il, i__2 = lmx - 2;
    il = min(i__1,i__2);

/*     THE RELATIVE END OF DATA FOR THIS PAGE IS IL. */
/*     SEARCH FOR A NONZERO VALUE WITH AN INDEX .GT. I ON THE PRESENT */
/*     PAGE. */

L20038:
    if (ipl >= il || ix[ipl] > *i__ && sx[ipl] != zero) {
	goto L20039;
    }
    ++ipl;
    goto L20038;

/*     TEST IF WE HAVE FOUND THE NEXT NONZERO. */

L20039:
    if (! (ix[ipl] > *i__ && sx[ipl] != zero && ipl <= il)) {
	goto L20040;
    }
    *i__ = ix[ipl];
    *xval = sx[ipl];
    *iplace = (np - 1) * lpg + ipl;
    return 0;

/*     UPDATE TO SCAN THE NEXT PAGE. */
L20040:
    ipl = ll + 1;
    ++np;
    goto L20035;

/*     NO DATA WAS FOUND. END OF VECTOR ENCOUNTERED. */

L20037:
    *i__ = 0;
    *xval = zero;
    ++il;
    if (il == lmx - 1) {
	il += 2;
    }

/*     IF A NEW ITEM WOULD BE INSERTED, IPLACE POINTS TO THE PLACE */
/*     TO PUT IT. */

    *iplace = (np - 1) * lpg + il;
    return 0;

/*     SEARCH A ROW FOR THE NEXT NONZERO. */
/*     FIND ELEMENT J=ABS(IRCX) IN ROWS ABS(I)+1,...,L. */

L20014:
    *i__ = abs(*i__);

/*     CHECK FOR END OF VECTOR. */

    if (! (*i__ == l)) {
	goto L20043;
    }
    *i__ = 0;
    *xval = zero;
    return 0;
L20043:
    i1 = *i__ + 1;
    ii = i1;
    n20046 = l;
    goto L20047;
L20046:
    ++ii;
L20047:
    if (n20046 - ii < 0) {
	goto L20048;
    }

/*     INITIALIZE IPPLOC FOR ORTHOGONAL SCAN. */
/*     LOOK FOR J AS A SUBSCRIPT IN ROWS II, II=I+1,...,L. */

    if (! (ii == 1)) {
	goto L20050;
    }
    ipploc = ll + 1;
    goto L20051;
L20050:
    ipploc = ix[ii + 3] + 1;
L20051:
    iend = ix[ii + 4];

/*     SCAN THROUGH SEVERAL PAGES, IF NECESSARY, TO FIND MATRIX ENTRY. */

    ipl = iploc_(&ipploc, &sx[1], &ix[1]);

/*     FIX UP IPPLOC AND IPL TO POINT TO MATRIX DATA. */

    idiff = lmx - ipl;
    if (! (idiff <= 1 && ix[lmx - 1] > 0)) {
	goto L20053;
    }
    ipploc = ipploc + idiff + 1;
    ipl = iploc_(&ipploc, &sx[1], &ix[1]);
L20053:
    np = (i__1 = ix[lmx - 1], abs(i__1));
    goto L20057;
L20056:
    if (ilast == iend) {
	goto L20058;
    }
L20057:
/* Computing MIN */
    i__1 = iend, i__2 = np * lpg + ll - 2;
    ilast = min(i__1,i__2);
    il = iploc_(&ilast, &sx[1], &ix[1]);
/* Computing MIN */
    i__1 = il, i__2 = lmx - 2;
    il = min(i__1,i__2);
L20059:
    if (ipl >= il || ix[ipl] >= j) {
	goto L20060;
    }
    ++ipl;
    goto L20059;

/*     TEST IF WE HAVE FOUND THE NEXT NONZERO. */

L20060:
    if (! (ix[ipl] == j && sx[ipl] != zero && ipl <= il)) {
	goto L20061;
    }
    *i__ = ii;
    *xval = sx[ipl];
    return 0;
L20061:
    if (ix[ipl] >= j) {
	ilast = iend;
    }
    ipl = ll + 1;
    ++np;
    goto L20056;
L20058:
    goto L20046;

/*     ORTHOGONAL SCAN FAILED. THE VALUE J WAS NOT A SUBSCRIPT */
/*     IN ANY ROW. */

L20048:
    *i__ = 0;
    *xval = zero;
    return 0;
} /* pnnzrs_ */

