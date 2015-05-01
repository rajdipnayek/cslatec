/* dlsod.f -- translated by f2c (version 12.02.01).
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
    doublereal told, rowns[210], el0, h__, hmin, hmxi, hu, x, u;
    integer iquit, init, lyh, lewt, lacor, lsavf, lwm, ksteps, ibegin, itol, 
	    iinteg, itstop, ijac, iband, iowns[6], ier, jstart, kflag, ldum, 
	    meth, miter, maxord, n, nq, nst, nfe, nje, nqu;
} ddebd1_;

#define ddebd1_1 ddebd1_

/* Table of constant values */

static integer c__4 = 4;
static integer c__1 = 1;
static integer c__6 = 6;
static integer c__2 = 2;
static integer c__7 = 7;
static integer c__8 = 8;
static doublereal c_b41 = 1.;
static integer c__14 = 14;
static integer c__9 = 9;
static integer c__3 = 3;
static integer c__10 = 10;
static integer c__5 = 5;
static integer c__11 = 11;
static integer c__12 = 12;
static integer c__0 = 0;

/* DECK DLSOD */
/* Subroutine */ int dlsod_(S_fp df, integer *neq, doublereal *t, doublereal *
	y, doublereal *tout, doublereal *rtol, doublereal *atol, integer *
	idid, doublereal *ypout, doublereal *yh, doublereal *yh1, doublereal *
	ewt, doublereal *savf, doublereal *acor, doublereal *wm, integer *iwm,
	 U_fp djac, logical *intout, doublereal *tstop, doublereal *tolfac, 
	doublereal *delsgn, doublereal *rpar, integer *ipar)
{
    /* Initialized data */

    static integer maxnum = 500;

    /* System generated locals */
    address a__1[2], a__2[7], a__3[6], a__4[8], a__5[3], a__6[5];
    integer yh_dim1, yh_offset, i__1[2], i__2, i__3[7], i__4[6], i__5[8], 
	    i__6[3], i__7[5];
    doublereal d__1, d__2, d__3, d__4;
    char ch__1[108], ch__2[216], ch__3[208], ch__4[112], ch__5[128], ch__6[
	    159];

    /* Local variables */
    static integer k, l;
    static doublereal ha, dt, big, del, tol;
    static integer ltol;
    static char xern1[8], xern3[16], xern4[16];
    extern /* Subroutine */ int dstod_(integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, doublereal *, doublereal *, doublereal *,
	     doublereal *, integer *, S_fp, U_fp, doublereal *, integer *);
    extern doublereal d1mach_(integer *);
    static doublereal absdel;
    static integer intflg;
    extern /* Subroutine */ int dintyd_(doublereal *, integer *, doublereal *,
	     integer *, doublereal *, integer *);
    static integer natolp;
    extern /* Subroutine */ int xermsg_(char *, char *, char *, integer *, 
	    integer *, ftnlen, ftnlen, ftnlen), dhstrt_(S_fp, integer *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, integer *, doublereal *, doublereal *, doublereal *,
	     doublereal *, doublereal *, doublereal *, doublereal *, integer *
	    , doublereal *);
    extern doublereal dvnrms_(integer *, doublereal *, doublereal *);
    static integer nrtolp;

    /* Fortran I/O blocks */
    static icilist io___3 = { 0, xern1, 0, "(I8)", 8, 1 };
    static icilist io___7 = { 0, xern1, 0, "(I8)", 8, 1 };
    static icilist io___9 = { 0, xern3, 0, "(1PE15.6)", 16, 1 };
    static icilist io___10 = { 0, xern1, 0, "(I8)", 8, 1 };
    static icilist io___11 = { 0, xern3, 0, "(1PE15.6)", 16, 1 };
    static icilist io___12 = { 0, xern3, 0, "(1PE15.6)", 16, 1 };
    static icilist io___14 = { 0, xern4, 0, "(1PE15.6)", 16, 1 };
    static icilist io___15 = { 0, xern3, 0, "(1PE15.6)", 16, 1 };
    static icilist io___16 = { 0, xern3, 0, "(1PE15.6)", 16, 1 };
    static icilist io___17 = { 0, xern4, 0, "(1PE15.6)", 16, 1 };
    static icilist io___18 = { 0, xern3, 0, "(1PE15.6)", 16, 1 };


/* ***BEGIN PROLOGUE  DLSOD */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to DDEBDF */
/* ***LIBRARY   SLATEC */
/* ***TYPE      DOUBLE PRECISION (LSOD-S, DLSOD-D) */
/* ***AUTHOR  (UNKNOWN) */
/* ***DESCRIPTION */

/*   DDEBDF  merely allocates storage for  DLSOD  to relieve the user of */
/*   the inconvenience of a long call list.  Consequently  DLSOD  is used */
/*   as described in the comments for  DDEBDF . */

/* ***SEE ALSO  DDEBDF */
/* ***ROUTINES CALLED  D1MACH, DHSTRT, DINTYD, DSTOD, DVNRMS, XERMSG */
/* ***COMMON BLOCKS    DDEBD1 */
/* ***REVISION HISTORY  (YYMMDD) */
/*   820301  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890831  Modified array declarations.  (WRB) */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900328  Added TYPE section.  (WRB) */
/*   900510  Convert XERRWV calls to XERMSG calls.  (RWC) */
/* ***END PROLOGUE  DLSOD */






/*     .................................................................. */

/*       THE EXPENSE OF SOLVING THE PROBLEM IS MONITORED BY COUNTING THE */
/*       NUMBER OF  STEPS ATTEMPTED. WHEN THIS EXCEEDS  MAXNUM, THE */
/*       COUNTER IS RESET TO ZERO AND THE USER IS INFORMED ABOUT POSSIBLE */
/*       EXCESSIVE WORK. */

    /* Parameter adjustments */
    yh_dim1 = *neq;
    yh_offset = 1 + yh_dim1;
    yh -= yh_offset;
    --y;
    --rtol;
    --atol;
    --ypout;
    --yh1;
    --ewt;
    --savf;
    --acor;
    --wm;
    --iwm;
    --rpar;
    --ipar;

    /* Function Body */

/*     .................................................................. */

/* ***FIRST EXECUTABLE STATEMENT  DLSOD */
    if (ddebd1_1.ibegin == 0) {

/*        ON THE FIRST CALL , PERFORM INITIALIZATION -- */
/*        DEFINE THE MACHINE UNIT ROUNDOFF QUANTITY  U  BY CALLING THE */
/*        FUNCTION ROUTINE D1MACH. THE USER MUST MAKE SURE THAT THE */
/*        VALUES SET IN D1MACH ARE RELEVANT TO THE COMPUTER BEING USED. */

	ddebd1_1.u = d1mach_(&c__4);
/*                          -- SET ASSOCIATED MACHINE DEPENDENT PARAMETER */
	wm[1] = sqrt(ddebd1_1.u);
/*                          -- SET TERMINATION FLAG */
	ddebd1_1.iquit = 0;
/*                          -- SET INITIALIZATION INDICATOR */
	ddebd1_1.init = 0;
/*                          -- SET COUNTER FOR ATTEMPTED STEPS */
	ddebd1_1.ksteps = 0;
/*                          -- SET INDICATOR FOR INTERMEDIATE-OUTPUT */
	*intout = FALSE_;
/*                          -- SET START INDICATOR FOR DSTOD CODE */
	ddebd1_1.jstart = 0;
/*                          -- SET BDF METHOD INDICATOR */
	ddebd1_1.meth = 2;
/*                          -- SET MAXIMUM ORDER FOR BDF METHOD */
	ddebd1_1.maxord = 5;
/*                          -- SET ITERATION MATRIX INDICATOR */

	if (ddebd1_1.ijac == 0 && ddebd1_1.iband == 0) {
	    ddebd1_1.miter = 2;
	}
	if (ddebd1_1.ijac == 1 && ddebd1_1.iband == 0) {
	    ddebd1_1.miter = 1;
	}
	if (ddebd1_1.ijac == 0 && ddebd1_1.iband == 1) {
	    ddebd1_1.miter = 5;
	}
	if (ddebd1_1.ijac == 1 && ddebd1_1.iband == 1) {
	    ddebd1_1.miter = 4;
	}

/*                          -- SET OTHER NECESSARY ITEMS IN COMMON BLOCK */
	ddebd1_1.n = *neq;
	ddebd1_1.nst = 0;
	ddebd1_1.nje = 0;
	ddebd1_1.hmxi = 0.;
	ddebd1_1.nq = 1;
	ddebd1_1.h__ = 1.;
/*                          -- RESET IBEGIN FOR SUBSEQUENT CALLS */
	ddebd1_1.ibegin = 1;
    }

/*     .................................................................. */

/*      CHECK VALIDITY OF INPUT PARAMETERS ON EACH ENTRY */

    if (*neq < 1) {
	s_wsfi(&io___3);
	do_fio(&c__1, (char *)&(*neq), (ftnlen)sizeof(integer));
	e_wsfi();
/* Writing concatenation */
	i__1[0] = 100, a__1[0] = "IN DDEBDF, THE NUMBER OF EQUATIONS MUST BE"
		" A POSITIVE INTEGER.$$YOU HAVE CALLED THE CODE WITH NEQ = ";
	i__1[1] = 8, a__1[1] = xern1;
	s_cat(ch__1, a__1, i__1, &c__2, (ftnlen)108);
	xermsg_("SLATEC", "DLSOD", ch__1, &c__6, &c__1, (ftnlen)6, (ftnlen)5, 
		(ftnlen)108);
	*idid = -33;
    }

    nrtolp = 0;
    natolp = 0;
    i__2 = *neq;
    for (k = 1; k <= i__2; ++k) {
	if (nrtolp <= 0) {
	    if (rtol[k] < 0.f) {
		s_wsfi(&io___7);
		do_fio(&c__1, (char *)&k, (ftnlen)sizeof(integer));
		e_wsfi();
		s_wsfi(&io___9);
		do_fio(&c__1, (char *)&rtol[k], (ftnlen)sizeof(doublereal));
		e_wsfi();
/* Writing concatenation */
		i__3[0] = 99, a__2[0] = "IN DDEBDF, THE RELATIVE ERROR TOLER"
			"ANCES MUST BE NON-NEGATIVE.$$YOU HAVE CALLED THE COD"
			"E WITH RTOL(";
		i__3[1] = 8, a__2[1] = xern1;
		i__3[2] = 4, a__2[2] = ") = ";
		i__3[3] = 16, a__2[3] = xern3;
		i__3[4] = 9, a__2[4] = "$$IN THE ";
		i__3[5] = 44, a__2[5] = "CASE OF VECTOR ERROR TOLERANCES, NO"
			" FURTHER ";
		i__3[6] = 36, a__2[6] = "CHECKING OF RTOL COMPONENTS IS DONE."
			;
		s_cat(ch__2, a__2, i__3, &c__7, (ftnlen)216);
		xermsg_("SLATEC", "DLSOD", ch__2, &c__7, &c__1, (ftnlen)6, (
			ftnlen)5, (ftnlen)216);
		*idid = -33;
		if (natolp > 0) {
		    goto L70;
		}
		nrtolp = 1;
	    } else if (natolp > 0) {
		goto L50;
	    }
	}

	if (atol[k] < 0.f) {
	    s_wsfi(&io___10);
	    do_fio(&c__1, (char *)&k, (ftnlen)sizeof(integer));
	    e_wsfi();
	    s_wsfi(&io___11);
	    do_fio(&c__1, (char *)&atol[k], (ftnlen)sizeof(doublereal));
	    e_wsfi();
/* Writing concatenation */
	    i__4[0] = 99, a__3[0] = "IN DDEBDF, THE ABSOLUTE ERROR TOLERANCE"
		    "S MUST BE NON-NEGATIVE.$$YOU HAVE CALLED THE CODE WITH A"
		    "TOL(";
	    i__4[1] = 8, a__3[1] = xern1;
	    i__4[2] = 4, a__3[2] = ") = ";
	    i__4[3] = 16, a__3[3] = xern3;
	    i__4[4] = 53, a__3[4] = "$$IN THE CASE OF VECTOR ERROR TOLERANCE"
		    "S, NO FURTHER ";
	    i__4[5] = 36, a__3[5] = "CHECKING OF ATOL COMPONENTS IS DONE.";
	    s_cat(ch__2, a__3, i__4, &c__6, (ftnlen)216);
	    xermsg_("SLATEC", "DLSOD", ch__2, &c__8, &c__1, (ftnlen)6, (
		    ftnlen)5, (ftnlen)216);
	    *idid = -33;
	    if (nrtolp > 0) {
		goto L70;
	    }
	    natolp = 1;
	}
L50:
	if (ddebd1_1.itol == 0) {
	    goto L70;
	}
/* L60: */
    }

L70:
    if (ddebd1_1.itstop == 1) {
	d__3 = *tout - *t;
	d__4 = *tstop - *t;
	if (d_sign(&c_b41, &d__3) != d_sign(&c_b41, &d__4) || (d__1 = *tout - 
		*t, abs(d__1)) > (d__2 = *tstop - *t, abs(d__2))) {
	    s_wsfi(&io___12);
	    do_fio(&c__1, (char *)&(*tout), (ftnlen)sizeof(doublereal));
	    e_wsfi();
	    s_wsfi(&io___14);
	    do_fio(&c__1, (char *)&(*tstop), (ftnlen)sizeof(doublereal));
	    e_wsfi();
/* Writing concatenation */
	    i__5[0] = 48, a__4[0] = "IN DDEBDF, YOU HAVE CALLED THE CODE WIT"
		    "H TOUT = ";
	    i__5[1] = 16, a__4[1] = xern3;
	    i__5[2] = 15, a__4[2] = "$$BUT YOU HAVE ";
	    i__5[3] = 51, a__4[3] = "ALSO TOLD THE CODE NOT TO INTEGRATE PAS"
		    "T THE POINT ";
	    i__5[4] = 8, a__4[4] = "TSTOP = ";
	    i__5[5] = 16, a__4[5] = xern4;
	    i__5[6] = 26, a__4[6] = " BY SETTING INFO(4) = 1.$$";
	    i__5[7] = 28, a__4[7] = "THESE INSTRUCTIONS CONFLICT.";
	    s_cat(ch__3, a__4, i__5, &c__8, (ftnlen)208);
	    xermsg_("SLATEC", "DLSOD", ch__3, &c__14, &c__1, (ftnlen)6, (
		    ftnlen)5, (ftnlen)208);
	    *idid = -33;
	}
    }

/*        CHECK SOME CONTINUATION POSSIBILITIES */

    if (ddebd1_1.init != 0) {
	if (*t == *tout) {
	    s_wsfi(&io___15);
	    do_fio(&c__1, (char *)&(*t), (ftnlen)sizeof(doublereal));
	    e_wsfi();
/* Writing concatenation */
	    i__6[0] = 52, a__5[0] = "IN DDEBDF, YOU HAVE CALLED THE CODE WIT"
		    "H T = TOUT = ";
	    i__6[1] = 16, a__5[1] = xern3;
	    i__6[2] = 44, a__5[2] = "$$THIS IS NOT ALLOWED ON CONTINUATION C"
		    "ALLS.";
	    s_cat(ch__4, a__5, i__6, &c__3, (ftnlen)112);
	    xermsg_("SLATEC", "DLSOD", ch__4, &c__9, &c__1, (ftnlen)6, (
		    ftnlen)5, (ftnlen)112);
	    *idid = -33;
	}

	if (*t != ddebd1_1.told) {
	    s_wsfi(&io___16);
	    do_fio(&c__1, (char *)&ddebd1_1.told, (ftnlen)sizeof(doublereal));
	    e_wsfi();
	    s_wsfi(&io___17);
	    do_fio(&c__1, (char *)&(*t), (ftnlen)sizeof(doublereal));
	    e_wsfi();
/* Writing concatenation */
	    i__7[0] = 48, a__6[0] = "IN DDEBDF, YOU HAVE CHANGED THE VALUE O"
		    "F T FROM ";
	    i__7[1] = 16, a__6[1] = xern3;
	    i__7[2] = 4, a__6[2] = " TO ";
	    i__7[3] = 16, a__6[3] = xern4;
	    i__7[4] = 44, a__6[4] = "  THIS IS NOT ALLOWED ON CONTINUATION C"
		    "ALLS.";
	    s_cat(ch__5, a__6, i__7, &c__5, (ftnlen)128);
	    xermsg_("SLATEC", "DLSOD", ch__5, &c__10, &c__1, (ftnlen)6, (
		    ftnlen)5, (ftnlen)128);
	    *idid = -33;
	}

	if (ddebd1_1.init != 1) {
	    if (*delsgn * (*tout - *t) < 0.) {
		s_wsfi(&io___18);
		do_fio(&c__1, (char *)&(*tout), (ftnlen)sizeof(doublereal));
		e_wsfi();
/* Writing concatenation */
		i__7[0] = 43, a__6[0] = "IN DDEBDF, BY CALLING THE CODE WITH"
			" TOUT = ";
		i__7[1] = 16, a__6[1] = xern3;
		i__7[2] = 34, a__6[2] = " YOU ARE ATTEMPTING TO CHANGE THE ";
		i__7[3] = 47, a__6[3] = "DIRECTION OF INTEGRATION.$$THIS IS "
			"NOT ALLOWED ";
		i__7[4] = 19, a__6[4] = "WITHOUT RESTARTING.";
		s_cat(ch__6, a__6, i__7, &c__5, (ftnlen)159);
		xermsg_("SLATEC", "DLSOD", ch__6, &c__11, &c__1, (ftnlen)6, (
			ftnlen)5, (ftnlen)159);
		*idid = -33;
	    }
	}
    }

    if (*idid == -33) {
	if (ddebd1_1.iquit != -33) {
/*                       INVALID INPUT DETECTED */
	    ddebd1_1.iquit = -33;
	    ddebd1_1.ibegin = -1;
	} else {
	    xermsg_("SLATEC", "DLSOD", "IN DDEBDF, INVALID INPUT WAS DETECTE"
		    "D ON SUCCESSIVE ENTRIES.  IT IS IMPOSSIBLE TO PROCEED BE"
		    "CAUSE YOU HAVE NOT CORRECTED THE PROBLEM, SO EXECUTION I"
		    "S BEING TERMINATED.", &c__12, &c__2, (ftnlen)6, (ftnlen)5,
		     (ftnlen)167);
	}
	return 0;
    }

/*        ............................................................... */

/*             RTOL = ATOL = 0. IS ALLOWED AS VALID INPUT AND INTERPRETED */
/*             AS ASKING FOR THE MOST ACCURATE SOLUTION POSSIBLE. IN THIS */
/*             CASE, THE RELATIVE ERROR TOLERANCE RTOL IS RESET TO THE */
/*             SMALLEST VALUE 100*U WHICH IS LIKELY TO BE REASONABLE FOR */
/*             THIS METHOD AND MACHINE */

    i__2 = *neq;
    for (k = 1; k <= i__2; ++k) {
	if (rtol[k] + atol[k] > 0.) {
	    goto L170;
	}
	rtol[k] = ddebd1_1.u * 100.;
	*idid = -2;
L170:
/*     ...EXIT */
	if (ddebd1_1.itol == 0) {
	    goto L190;
	}
/* L180: */
    }
L190:

    if (*idid != -2) {
	goto L200;
    }
/*        RTOL=ATOL=0 ON INPUT, SO RTOL IS CHANGED TO A */
/*                                 SMALL POSITIVE VALUE */
    ddebd1_1.ibegin = -1;
    goto L460;
L200:
/*        BEGIN BLOCK PERMITTING ...EXITS TO 450 */
/*           BEGIN BLOCK PERMITTING ...EXITS TO 430 */
/*              BEGIN BLOCK PERMITTING ...EXITS TO 260 */
/*                 BEGIN BLOCK PERMITTING ...EXITS TO 230 */

/*                    BRANCH ON STATUS OF INITIALIZATION INDICATOR */
/*                           INIT=0 MEANS INITIAL DERIVATIVES AND */
/*                           NOMINAL STEP SIZE */
/*                                  AND DIRECTION NOT YET SET */
/*                           INIT=1 MEANS NOMINAL STEP SIZE AND */
/*                           DIRECTION NOT YET SET INIT=2 MEANS NO */
/*                           FURTHER INITIALIZATION REQUIRED */

    if (ddebd1_1.init == 0) {
	goto L210;
    }
/*                 ......EXIT */
    if (ddebd1_1.init == 1) {
	goto L230;
    }
/*              .........EXIT */
    goto L260;
L210:

/*                    ................................................ */

/*                         MORE INITIALIZATION -- */
/*                                             -- EVALUATE INITIAL */
/*                                             DERIVATIVES */

    ddebd1_1.init = 1;
    (*df)(t, &y[1], &yh[(yh_dim1 << 1) + 1], &rpar[1], &ipar[1]);
    ddebd1_1.nfe = 1;
/*                 ...EXIT */
    if (*t != *tout) {
	goto L230;
    }
    *idid = 2;
    i__2 = *neq;
    for (l = 1; l <= i__2; ++l) {
	ypout[l] = yh[l + (yh_dim1 << 1)];
/* L220: */
    }
    ddebd1_1.told = *t;
/*        ............EXIT */
    goto L450;
L230:

/*                 -- COMPUTE INITIAL STEP SIZE */
/*                 -- SAVE SIGN OF INTEGRATION DIRECTION */
/*                 -- SET INDEPENDENT AND DEPENDENT VARIABLES */
/*                                      X AND YH(*) FOR DSTOD */

    ltol = 1;
    i__2 = *neq;
    for (l = 1; l <= i__2; ++l) {
	if (ddebd1_1.itol == 1) {
	    ltol = l;
	}
	tol = rtol[ltol] * (d__1 = y[l], abs(d__1)) + atol[ltol];
	if (tol == 0.) {
	    goto L390;
	}
	ewt[l] = tol;
/* L240: */
    }

    big = sqrt(d1mach_(&c__2));
    dhstrt_((S_fp)df, neq, t, tout, &y[1], &yh[(yh_dim1 << 1) + 1], &ewt[1], &
	    c__1, &ddebd1_1.u, &big, &yh[yh_dim1 * 3 + 1], &yh[(yh_dim1 << 2) 
	    + 1], &yh[yh_dim1 * 5 + 1], &yh[yh_dim1 * 6 + 1], &rpar[1], &ipar[
	    1], &ddebd1_1.h__);

    d__1 = *tout - *t;
    *delsgn = d_sign(&c_b41, &d__1);
    ddebd1_1.x = *t;
    i__2 = *neq;
    for (l = 1; l <= i__2; ++l) {
	yh[l + yh_dim1] = y[l];
	yh[l + (yh_dim1 << 1)] = ddebd1_1.h__ * yh[l + (yh_dim1 << 1)];
/* L250: */
    }
    ddebd1_1.init = 2;
L260:

/*              ...................................................... */

/*                 ON EACH CALL SET INFORMATION WHICH DETERMINES THE */
/*                 ALLOWED INTERVAL OF INTEGRATION BEFORE RETURNING */
/*                 WITH AN ANSWER AT TOUT */

    del = *tout - *t;
    absdel = abs(del);

/*              ...................................................... */

/*                 IF ALREADY PAST OUTPUT POINT, INTERPOLATE AND */
/*                 RETURN */

L270:
/*                 BEGIN BLOCK PERMITTING ...EXITS TO 400 */
/*                    BEGIN BLOCK PERMITTING ...EXITS TO 380 */
    if ((d__1 = ddebd1_1.x - *t, abs(d__1)) < absdel) {
	goto L290;
    }
    dintyd_(tout, &c__0, &yh[yh_offset], neq, &y[1], &intflg);
    dintyd_(tout, &c__1, &yh[yh_offset], neq, &ypout[1], &intflg);
    *idid = 3;
    if (ddebd1_1.x != *tout) {
	goto L280;
    }
    *idid = 2;
    *intout = FALSE_;
L280:
    *t = *tout;
    ddebd1_1.told = *t;
/*        ..................EXIT */
    goto L450;
L290:

/*                       IF CANNOT GO PAST TSTOP AND SUFFICIENTLY */
/*                       CLOSE, EXTRAPOLATE AND RETURN */

    if (ddebd1_1.itstop != 1) {
	goto L310;
    }
    if ((d__1 = *tstop - ddebd1_1.x, abs(d__1)) >= ddebd1_1.u * 100. * abs(
	    ddebd1_1.x)) {
	goto L310;
    }
    dt = *tout - ddebd1_1.x;
    i__2 = *neq;
    for (l = 1; l <= i__2; ++l) {
	y[l] = yh[l + yh_dim1] + dt / ddebd1_1.h__ * yh[l + (yh_dim1 << 1)];
/* L300: */
    }
    (*df)(tout, &y[1], &ypout[1], &rpar[1], &ipar[1]);
    ++ddebd1_1.nfe;
    *idid = 3;
    *t = *tout;
    ddebd1_1.told = *t;
/*        ..................EXIT */
    goto L450;
L310:

    if (ddebd1_1.iinteg == 0 || ! (*intout)) {
	goto L320;
    }

/*                          INTERMEDIATE-OUTPUT MODE */

    *idid = 1;
    goto L370;
L320:

/*                       ............................................. */

/*                            MONITOR NUMBER OF STEPS ATTEMPTED */

    if (ddebd1_1.ksteps <= maxnum) {
	goto L330;
    }

/*                          A SIGNIFICANT AMOUNT OF WORK HAS BEEN */
/*                          EXPENDED */
    *idid = -1;
    ddebd1_1.ksteps = 0;
    ddebd1_1.ibegin = -1;
    goto L370;
L330:

/*                          .......................................... */

/*                             LIMIT STEP SIZE AND SET WEIGHT VECTOR */

    ddebd1_1.hmin = ddebd1_1.u * 100. * abs(ddebd1_1.x);
/* Computing MAX */
    d__1 = abs(ddebd1_1.h__);
    ha = max(d__1,ddebd1_1.hmin);
    if (ddebd1_1.itstop == 1) {
/* Computing MIN */
	d__2 = ha, d__3 = (d__1 = *tstop - ddebd1_1.x, abs(d__1));
	ha = min(d__2,d__3);
    }
    ddebd1_1.h__ = d_sign(&ha, &ddebd1_1.h__);
    ltol = 1;
    i__2 = *neq;
    for (l = 1; l <= i__2; ++l) {
	if (ddebd1_1.itol == 1) {
	    ltol = l;
	}
	ewt[l] = rtol[ltol] * (d__1 = yh[l + yh_dim1], abs(d__1)) + atol[ltol]
		;
/*                    .........EXIT */
	if (ewt[l] <= 0.) {
	    goto L380;
	}
/* L340: */
    }
    *tolfac = ddebd1_1.u * dvnrms_(neq, &yh[yh_offset], &ewt[1]);
/*                 .........EXIT */
    if (*tolfac <= 1.) {
	goto L400;
    }

/*                          TOLERANCES TOO SMALL */
    *idid = -2;
    *tolfac *= 2.;
    rtol[1] = *tolfac * rtol[1];
    atol[1] = *tolfac * atol[1];
    if (ddebd1_1.itol == 0) {
	goto L360;
    }
    i__2 = *neq;
    for (l = 2; l <= i__2; ++l) {
	rtol[l] = *tolfac * rtol[l];
	atol[l] = *tolfac * atol[l];
/* L350: */
    }
L360:
    ddebd1_1.ibegin = -1;
L370:
/*           ............EXIT */
    goto L430;
L380:

/*                    RELATIVE ERROR CRITERION INAPPROPRIATE */
L390:
    *idid = -3;
    ddebd1_1.ibegin = -1;
/*           .........EXIT */
    goto L430;
L400:

/*                 ................................................... */

/*                      TAKE A STEP */

    dstod_(neq, &y[1], &yh[yh_offset], neq, &yh1[1], &ewt[1], &savf[1], &acor[
	    1], &wm[1], &iwm[1], (S_fp)df, (U_fp)djac, &rpar[1], &ipar[1]);

    ddebd1_1.jstart = -2;
    *intout = TRUE_;
    if (ddebd1_1.kflag == 0) {
	goto L270;
    }

/*              ...................................................... */

    if (ddebd1_1.kflag == -1) {
	goto L410;
    }

/*                 REPEATED CORRECTOR CONVERGENCE FAILURES */
    *idid = -6;
    ddebd1_1.ibegin = -1;
    goto L420;
L410:

/*                 REPEATED ERROR TEST FAILURES */
    *idid = -7;
    ddebd1_1.ibegin = -1;
L420:
L430:

/*           ......................................................... */

/*                                  STORE VALUES BEFORE RETURNING TO */
/*                                  DDEBDF */
    i__2 = *neq;
    for (l = 1; l <= i__2; ++l) {
	y[l] = yh[l + yh_dim1];
	ypout[l] = yh[l + (yh_dim1 << 1)] / ddebd1_1.h__;
/* L440: */
    }
    *t = ddebd1_1.x;
    ddebd1_1.told = *t;
    *intout = FALSE_;
L450:
L460:
    return 0;
} /* dlsod_ */

