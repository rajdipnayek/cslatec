/* dpopt.f -- translated by f2c (version 12.02.01).
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

static integer c__4 = 4;

/* DECK DPOPT */
/* Subroutine */ int dpopt_(doublereal *prgopt, integer *mrelas, integer *
	nvars__, integer *info, doublereal *csc, integer *ibasis, doublereal *
	ropt, integer *intopt, logical *lopt)
{
    /* System generated locals */
    doublereal d__1, d__2;

    /* Local variables */
    static integer i__, j, idg;
    static doublereal one;
    static integer lds;
    static doublereal eps;
    static integer key, npp, n20043, n20053, n20096;
    static doublereal abig;
    static integer last, lprg, nerr;
    static doublereal tune;
    static integer iopt;
    static doublereal zero;
    static integer next, itbrc, isave, itest;
    static doublereal tolls;
    extern doublereal d1mach_(integer *);
    static integer iadbig, ipagef;
    static doublereal asmall, tolabs;
    static logical colscp, savedt, minprb, stpedg;
    static integer ictmax;
    static logical contin;
    static doublereal costsc;
    static logical usrbas, cstscp;
    static integer ictopt;
    extern /* Subroutine */ int xermsg_(char *, char *, char *, integer *, 
	    integer *, ftnlen, ftnlen, ftnlen);
    static integer kprint, mxitlp;
    static logical sizeup;

/* ***BEGIN PROLOGUE  DPOPT */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to DSPLP */
/* ***LIBRARY   SLATEC */
/* ***TYPE      DOUBLE PRECISION (SPOPT-S, DPOPT-D) */
/* ***AUTHOR  (UNKNOWN) */
/* ***DESCRIPTION */

/*     THE EDITING REQUIRED TO CONVERT THIS SUBROUTINE FROM SINGLE TO */
/*     DOUBLE PRECISION INVOLVES THE FOLLOWING CHARACTER STRING CHANGES. */

/*     USE AN EDITING COMMAND (CHANGE) /STRING-1/(TO)STRING-2/. */
/*     /REAL (12 BLANKS)/DOUBLE PRECISION/,/R1MACH/D1MACH/,/E0/D0/ */

/*     REVISED 821122-1045 */
/*     REVISED YYMMDD-HHMM */

/*     THIS SUBROUTINE PROCESSES THE OPTION VECTOR, PRGOPT(*), */
/*     AND VALIDATES ANY MODIFIED DATA. */

/* ***SEE ALSO  DSPLP */
/* ***ROUTINES CALLED  D1MACH, XERMSG */
/* ***REVISION HISTORY  (YYMMDD) */
/*   811215  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890605  Removed unreferenced labels.  (WRB) */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ) */
/*   900328  Added TYPE section.  (WRB) */
/*   900510  Fixed an error message.  (RWC) */
/* ***END PROLOGUE  DPOPT */

/* ***FIRST EXECUTABLE STATEMENT  DPOPT */
    /* Parameter adjustments */
    --lopt;
    --intopt;
    --ropt;
    --ibasis;
    --csc;
    --prgopt;

    /* Function Body */
    iopt = 1;
    zero = 0.;
    one = 1.;
    goto L30001;
L20002:
    goto L30002;

L20003:
    lopt[1] = contin;
    lopt[2] = usrbas;
    lopt[3] = sizeup;
    lopt[4] = savedt;
    lopt[5] = colscp;
    lopt[6] = cstscp;
    lopt[7] = minprb;
    lopt[8] = stpedg;

    intopt[1] = idg;
    intopt[2] = ipagef;
    intopt[3] = isave;
    intopt[4] = mxitlp;
    intopt[5] = kprint;
    intopt[6] = itbrc;
    intopt[7] = npp;
    intopt[8] = lprg;

    ropt[1] = eps;
    ropt[2] = asmall;
    ropt[3] = abig;
    ropt[4] = costsc;
    ropt[5] = tolls;
    ropt[6] = tune;
    ropt[7] = tolabs;
    return 0;


/*     PROCEDURE (INITIALIZE PARAMETERS AND PROCESS USER OPTIONS) */
L30001:
    contin = FALSE_;
    usrbas = FALSE_;
    sizeup = FALSE_;
    savedt = FALSE_;
    colscp = FALSE_;
    cstscp = FALSE_;
    minprb = TRUE_;
    stpedg = TRUE_;

/*     GET THE MACHINE REL. FLOATING POINT ACCURACY VALUE FROM THE */
/*     LIBRARY SUBPROGRAM, D1MACH( ). */
    eps = d1mach_(&c__4);
    tolls = d1mach_(&c__4);
    tune = one;
    tolabs = zero;

/*     DEFINE NOMINAL FILE NUMBERS FOR MATRIX PAGES AND DATA SAVING. */
    ipagef = 1;
    isave = 2;
    itbrc = 10;
    mxitlp = (*nvars__ + *mrelas) * 3;
    kprint = 0;
    idg = -4;
    npp = *nvars__;
    lprg = 0;

    last = 1;
    iadbig = 10000;
    ictmax = 1000;
    ictopt = 0;
L20004:
    next = (integer) prgopt[last];
    if (! (next <= 0 || next > iadbig)) {
	goto L20006;
    }

/*     THE CHECKS FOR SMALL OR LARGE VALUES OF NEXT ARE TO PREVENT */
/*     WORKING WITH UNDEFINED DATA. */
    nerr = 14;
    xermsg_("SLATEC", "DPOPT", "IN DSPLP, THE USER OPTION ARRAY HAS UNDEFINE"
	    "D DATA.", &nerr, &iopt, (ftnlen)6, (ftnlen)5, (ftnlen)51);
    *info = -nerr;
    return 0;
L20006:
    if (! (next == 1)) {
	goto L10001;
    }
    goto L20005;
L10001:
    if (! (ictopt > ictmax)) {
	goto L10002;
    }
    nerr = 15;
    xermsg_("SLATEC", "DPOPT", "IN DSPLP, OPTION ARRAY PROCESSING IS CYCLING."
	    , &nerr, &iopt, (ftnlen)6, (ftnlen)5, (ftnlen)45);
    *info = -nerr;
    return 0;
L10002:
    key = (integer) prgopt[last + 1];

/*     IF KEY = 50, THIS IS TO BE A MAXIMIZATION PROBLEM */
/*     INSTEAD OF A MINIMIZATION PROBLEM. */
    if (! (key == 50)) {
	goto L20010;
    }
    minprb = prgopt[last + 2] == zero;
    lds = 3;
    goto L20009;
L20010:

/*     IF KEY = 51, THE LEVEL OF OUTPUT IS BEING MODIFIED. */
/*     KPRINT = 0, NO OUTPUT */
/*            = 1, SUMMARY OUTPUT */
/*            = 2, LOTS OF OUTPUT */
/*            = 3, EVEN MORE OUTPUT */
    if (! (key == 51)) {
	goto L20013;
    }
    kprint = (integer) prgopt[last + 2];
    lds = 3;
    goto L20009;
L20013:

/*     IF KEY = 52, REDEFINE THE FORMAT AND PRECISION USED */
/*     IN THE OUTPUT. */
    if (! (key == 52)) {
	goto L20016;
    }
    if (prgopt[last + 2] != zero) {
	idg = (integer) prgopt[last + 3];
    }
    lds = 4;
    goto L20009;
L20016:

/*     IF KEY = 53, THE ALLOTTED SPACE FOR THE SPARSE MATRIX */
/*     STORAGE AND/OR SPARSE EQUATION SOLVING HAS BEEN CHANGED. */
/*    (PROCESSED IN DSPLP(). THIS IS TO COMPUTE THE LENGTH OF PRGOPT(*).) */
    if (! (key == 53)) {
	goto L20019;
    }
    lds = 5;
    goto L20009;
L20019:

/*     IF KEY = 54, REDEFINE THE FILE NUMBER WHERE THE PAGES */
/*     FOR THE SPARSE MATRIX ARE STORED. */
    if (! (key == 54)) {
	goto L20022;
    }
    if (prgopt[last + 2] != zero) {
	ipagef = (integer) prgopt[last + 3];
    }
    lds = 4;
    goto L20009;
L20022:

/*     IF KEY = 55,  A CONTINUATION FOR A PROBLEM MAY BE REQUESTED. */
    if (! (key == 55)) {
	goto L20025;
    }
    contin = prgopt[last + 2] != zero;
    lds = 3;
    goto L20009;
L20025:

/*     IF KEY = 56, REDEFINE THE FILE NUMBER WHERE THE SAVED DATA */
/*     WILL BE STORED. */
    if (! (key == 56)) {
	goto L20028;
    }
    if (prgopt[last + 2] != zero) {
	isave = (integer) prgopt[last + 3];
    }
    lds = 4;
    goto L20009;
L20028:

/*     IF KEY = 57, SAVE DATA (ON EXTERNAL FILE)  AT MXITLP ITERATIONS OR */
/*     THE OPTIMUM, WHICHEVER COMES FIRST. */
    if (! (key == 57)) {
	goto L20031;
    }
    savedt = prgopt[last + 2] != zero;
    lds = 3;
    goto L20009;
L20031:

/*     IF KEY = 58,  SEE IF PROBLEM IS TO RUN ONLY A GIVEN */
/*     NUMBER OF ITERATIONS. */
    if (! (key == 58)) {
	goto L20034;
    }
    if (prgopt[last + 2] != zero) {
	mxitlp = (integer) prgopt[last + 3];
    }
    lds = 4;
    goto L20009;
L20034:

/*     IF KEY = 59,  SEE IF USER PROVIDES THE BASIS INDICES. */
    if (! (key == 59)) {
	goto L20037;
    }
    usrbas = prgopt[last + 2] != zero;
    if (! usrbas) {
	goto L20040;
    }
    i__ = 1;
    n20043 = *mrelas;
    goto L20044;
L20043:
    ++i__;
L20044:
    if (n20043 - i__ < 0) {
	goto L20045;
    }
    ibasis[i__] = (integer) prgopt[last + 2 + i__];
    goto L20043;
L20045:
L20040:
    lds = *mrelas + 3;
    goto L20009;
L20037:

/*     IF KEY = 60,  SEE IF USER HAS PROVIDED SCALING OF COLUMNS. */
    if (! (key == 60)) {
	goto L20047;
    }
    colscp = prgopt[last + 2] != zero;
    if (! colscp) {
	goto L20050;
    }
    j = 1;
    n20053 = *nvars__;
    goto L20054;
L20053:
    ++j;
L20054:
    if (n20053 - j < 0) {
	goto L20055;
    }
    csc[j] = (d__1 = prgopt[last + 2 + j], abs(d__1));
    goto L20053;
L20055:
L20050:
    lds = *nvars__ + 3;
    goto L20009;
L20047:

/*     IF KEY = 61,  SEE IF USER HAS PROVIDED SCALING OF COSTS. */
    if (! (key == 61)) {
	goto L20057;
    }
    cstscp = prgopt[last + 2] != zero;
    if (cstscp) {
	costsc = prgopt[last + 3];
    }
    lds = 4;
    goto L20009;
L20057:

/*     IF KEY = 62,  SEE IF SIZE PARAMETERS ARE PROVIDED WITH THE DATA. */
/*     THESE WILL BE CHECKED AGAINST THE MATRIX ELEMENT SIZES LATER. */
    if (! (key == 62)) {
	goto L20060;
    }
    sizeup = prgopt[last + 2] != zero;
    if (! sizeup) {
	goto L20063;
    }
    asmall = prgopt[last + 3];
    abig = prgopt[last + 4];
L20063:
    lds = 5;
    goto L20009;
L20060:

/*     IF KEY = 63, SEE IF TOLERANCE FOR LINEAR SYSTEM RESIDUAL ERROR IS */
/*     PROVIDED. */
    if (! (key == 63)) {
	goto L20066;
    }
    if (prgopt[last + 2] != zero) {
/* Computing MAX */
	d__1 = eps, d__2 = prgopt[last + 3];
	tolls = max(d__1,d__2);
    }
    lds = 4;
    goto L20009;
L20066:

/*     IF KEY = 64,  SEE IF MINIMUM REDUCED COST OR STEEPEST EDGE */
/*     DESCENT IS TO BE USED FOR SELECTING VARIABLES TO ENTER BASIS. */
    if (! (key == 64)) {
	goto L20069;
    }
    stpedg = prgopt[last + 2] == zero;
    lds = 3;
    goto L20009;
L20069:

/*     IF KEY = 65, SET THE NUMBER OF ITERATIONS BETWEEN RECALCULATING */
/*     THE ERROR IN THE PRIMAL SOLUTION. */
    if (! (key == 65)) {
	goto L20072;
    }
    if (prgopt[last + 2] != zero) {
/* Computing MAX */
	d__1 = one, d__2 = prgopt[last + 3];
	itbrc = (integer) max(d__1,d__2);
    }
    lds = 4;
    goto L20009;
L20072:

/*     IF KEY = 66, SET THE NUMBER OF NEGATIVE REDUCED COSTS TO BE FOUND */
/*     IN THE PARTIAL PRICING STRATEGY. */
    if (! (key == 66)) {
	goto L20075;
    }
    if (! (prgopt[last + 2] != zero)) {
	goto L20078;
    }
/* Computing MAX */
    d__1 = prgopt[last + 3];
    npp = (integer) max(d__1,one);
    npp = min(npp,*nvars__);
L20078:
    lds = 4;
    goto L20009;
L20075:
/*     IF KEY = 67, CHANGE THE TUNING PARAMETER TO APPLY TO THE ERROR */
/*     ESTIMATES FOR THE PRIMAL AND DUAL SYSTEMS. */
    if (! (key == 67)) {
	goto L20081;
    }
    if (! (prgopt[last + 2] != zero)) {
	goto L20084;
    }
    tune = (d__1 = prgopt[last + 3], abs(d__1));
L20084:
    lds = 4;
    goto L20009;
L20081:
    if (! (key == 68)) {
	goto L20087;
    }
    lds = 6;
    goto L20009;
L20087:

/*     RESET THE ABSOLUTE TOLERANCE TO BE USED ON THE FEASIBILITY */
/*     DECISION PROVIDED THE RELATIVE ERROR TEST FAILED. */
    if (! (key == 69)) {
	goto L20090;
    }
    if (prgopt[last + 2] != zero) {
	tolabs = prgopt[last + 3];
    }
    lds = 4;
    goto L20009;
L20090:

L20009:
    ++ictopt;
    last = next;
    lprg += lds;
    goto L20004;
L20005:
    goto L20002;

/*     PROCEDURE (VALIDATE OPTIONALLY MODIFIED DATA) */

/*     IF USER HAS DEFINED THE BASIS, CHECK FOR VALIDITY OF INDICES. */
L30002:
    if (! usrbas) {
	goto L20093;
    }
    i__ = 1;
    n20096 = *mrelas;
    goto L20097;
L20096:
    ++i__;
L20097:
    if (n20096 - i__ < 0) {
	goto L20098;
    }
    itest = ibasis[i__];
    if (! (itest <= 0 || itest > *nvars__ + *mrelas)) {
	goto L20100;
    }
    nerr = 16;
    xermsg_("SLATEC", "DPOPT", "IN DSPLP, AN INDEX OF USER-SUPPLIED BASIS IS"
	    " OUT OF RANGE.", &nerr, &iopt, (ftnlen)6, (ftnlen)5, (ftnlen)58);
    *info = -nerr;
    return 0;
L20100:
    goto L20096;
L20098:
L20093:

/*     IF USER HAS PROVIDED SIZE PARAMETERS, MAKE SURE THEY ARE ORDERED */
/*     AND POSITIVE. */
    if (! sizeup) {
	goto L20103;
    }
    if (! (asmall <= zero || abig < asmall)) {
	goto L20106;
    }
    nerr = 17;
    xermsg_("SLATEC", "DPOPT", "IN DSPLP, SIZE PARAMETERS FOR MATRIX MUST BE"
	    " SMALLEST AND LARGEST MAGNITUDES OF NONZERO ENTRIES.", &nerr, &
	    iopt, (ftnlen)6, (ftnlen)5, (ftnlen)96);
    *info = -nerr;
    return 0;
L20106:
L20103:

/*     THE NUMBER OF ITERATIONS OF REV. SIMPLEX STEPS MUST BE POSITIVE. */
    if (! (mxitlp <= 0)) {
	goto L20109;
    }
    nerr = 18;
    xermsg_("SLATEC", "DPOPT", "IN DSPLP, THE NUMBER OF REVISED SIMPLEX STEP"
	    "S BETWEEN CHECK-POINTS MUST BE POSITIVE.", &nerr, &iopt, (ftnlen)
	    6, (ftnlen)5, (ftnlen)84);
    *info = -nerr;
    return 0;
L20109:

/*     CHECK THAT SAVE AND PAGE FILE NUMBERS ARE DEFINED AND NOT EQUAL. */
    if (! (isave <= 0 || ipagef <= 0 || isave == ipagef)) {
	goto L20112;
    }
    nerr = 19;
    xermsg_("SLATEC", "DPOPT", "IN DSPLP, FILE NUMBERS FOR SAVED DATA AND MA"
	    "TRIX PAGES MUST BE POSITIVE AND NOT EQUAL.", &nerr, &iopt, (
	    ftnlen)6, (ftnlen)5, (ftnlen)86);
    *info = -nerr;
    return 0;
L20112:
    goto L20003;
} /* dpopt_ */

