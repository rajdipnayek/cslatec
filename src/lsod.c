/* lsod.f -- translated by f2c (version 12.02.01).
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
    real told, rowns[210], el0, h__, hmin, hmxi, hu, x, u;
    integer iquit, init, lyh, lewt, lacor, lsavf, lwm, ksteps, ibegin, itol, 
	    iinteg, itstop, ijac, iband, iowns[6], ier, jstart, kflag, ldum, 
	    meth, miter, maxord, n, nq, nst, nfe, nje, nqu;
} debdf1_;

#define debdf1_1 debdf1_

/* Table of constant values */

static integer c__4 = 4;
static integer c__1 = 1;
static integer c__6 = 6;
static integer c__2 = 2;
static integer c__7 = 7;
static integer c__8 = 8;
static real c_b41 = 1.f;
static integer c__14 = 14;
static integer c__9 = 9;
static integer c__3 = 3;
static integer c__10 = 10;
static integer c__5 = 5;
static integer c__11 = 11;
static integer c__12 = 12;
static integer c__0 = 0;

/* DECK LSOD */
/* Subroutine */ int lsod_(S_fp f, integer *neq, real *t, real *y, real *tout,
	 real *rtol, real *atol, integer *idid, real *ypout, real *yh, real *
	yh1, real *ewt, real *savf, real *acor, real *wm, integer *iwm, U_fp 
	jac, logical *intout, real *tstop, real *tolfac, real *delsgn, real *
	rpar, integer *ipar)
{
    /* Initialized data */

    static integer maxnum = 500;

    /* System generated locals */
    address a__1[2], a__2[7], a__3[6], a__4[8], a__5[3], a__6[5];
    integer yh_dim1, yh_offset, i__1[2], i__2, i__3[7], i__4[6], i__5[8], 
	    i__6[3], i__7[5];
    real r__1, r__2, r__3, r__4;
    char ch__1[107], ch__2[215], ch__3[207], ch__4[111], ch__5[127], ch__6[
	    158];

    /* Local variables */
    static integer k, l;
    static real ha, dt, big, del, tol;
    extern /* Subroutine */ int stod_(integer *, real *, real *, integer *, 
	    real *, real *, real *, real *, real *, integer *, S_fp, U_fp, 
	    real *, integer *);
    static integer ltol;
    static char xern1[8], xern3[16], xern4[16];
    extern /* Subroutine */ int intyd_(real *, integer *, real *, integer *, 
	    real *, integer *);
    extern doublereal r1mach_(integer *);
    static real absdel;
    static integer intflg, natolp;
    extern /* Subroutine */ int xermsg_(char *, char *, char *, integer *, 
	    integer *, ftnlen, ftnlen, ftnlen), hstart_(S_fp, integer *, real 
	    *, real *, real *, real *, real *, integer *, real *, real *, 
	    real *, real *, real *, real *, real *, integer *, real *);
    static integer nrtolp;
    extern doublereal vnwrms_(integer *, real *, real *);

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


/* ***BEGIN PROLOGUE  LSOD */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to DEBDF */
/* ***LIBRARY   SLATEC */
/* ***TYPE      SINGLE PRECISION (LSOD-S, DLSOD-D) */
/* ***AUTHOR  (UNKNOWN) */
/* ***DESCRIPTION */

/*   DEBDF  merely allocates storage for  LSOD  to relieve the user of */
/*   the inconvenience of a long call list.  Consequently  LSOD  is used */
/*   as described in the comments for  DEBDF . */

/* ***SEE ALSO  DEBDF */
/* ***ROUTINES CALLED  HSTART, INTYD, R1MACH, STOD, VNWRMS, XERMSG */
/* ***COMMON BLOCKS    DEBDF1 */
/* ***REVISION HISTORY  (YYMMDD) */
/*   800901  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890831  Modified array declarations.  (WRB) */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900328  Added TYPE section.  (WRB) */
/*   900510  Convert XERRWV calls to XERMSG calls.  (RWC) */
/* ***END PROLOGUE  LSOD */






/* ....................................................................... */

/*  THE EXPENSE OF SOLVING THE PROBLEM IS MONITORED BY COUNTING THE */
/*  NUMBER OF  STEPS ATTEMPTED. WHEN THIS EXCEEDS  MAXNUM, THE COUNTER */
/*  IS RESET TO ZERO AND THE USER IS INFORMED ABOUT POSSIBLE EXCESSIVE */
/*  WORK. */

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

/* ....................................................................... */

/* ***FIRST EXECUTABLE STATEMENT  LSOD */
    if (debdf1_1.ibegin == 0) {

/*        ON THE FIRST CALL , PERFORM INITIALIZATION -- */
/*        DEFINE THE MACHINE UNIT ROUNDOFF QUANTITY  U  BY CALLING THE */
/*        FUNCTION ROUTINE R1MACH. THE USER MUST MAKE SURE THAT THE */
/*        VALUES SET IN R1MACH ARE RELEVANT TO THE COMPUTER BEING USED. */

	debdf1_1.u = r1mach_(&c__4);
/*                          -- SET ASSOCIATED MACHINE DEPENDENT PARAMETER */
	wm[1] = sqrt(debdf1_1.u);
/*                          -- SET TERMINATION FLAG */
	debdf1_1.iquit = 0;
/*                          -- SET INITIALIZATION INDICATOR */
	debdf1_1.init = 0;
/*                          -- SET COUNTER FOR ATTEMPTED STEPS */
	debdf1_1.ksteps = 0;
/*                          -- SET INDICATOR FOR INTERMEDIATE-OUTPUT */
	*intout = FALSE_;
/*                          -- SET START INDICATOR FOR STOD CODE */
	debdf1_1.jstart = 0;
/*                          -- SET BDF METHOD INDICATOR */
	debdf1_1.meth = 2;
/*                          -- SET MAXIMUM ORDER FOR BDF METHOD */
	debdf1_1.maxord = 5;
/*                          -- SET ITERATION MATRIX INDICATOR */

	if (debdf1_1.ijac == 0 && debdf1_1.iband == 0) {
	    debdf1_1.miter = 2;
	}
	if (debdf1_1.ijac == 1 && debdf1_1.iband == 0) {
	    debdf1_1.miter = 1;
	}
	if (debdf1_1.ijac == 0 && debdf1_1.iband == 1) {
	    debdf1_1.miter = 5;
	}
	if (debdf1_1.ijac == 1 && debdf1_1.iband == 1) {
	    debdf1_1.miter = 4;
	}

/*                          -- SET OTHER NECESSARY ITEMS IN COMMON BLOCK */
	debdf1_1.n = *neq;
	debdf1_1.nst = 0;
	debdf1_1.nje = 0;
	debdf1_1.hmxi = 0.f;
	debdf1_1.nq = 1;
	debdf1_1.h__ = 1.f;
/*                          -- RESET IBEGIN FOR SUBSEQUENT CALLS */
	debdf1_1.ibegin = 1;
    }

/* ....................................................................... */

/*      CHECK VALIDITY OF INPUT PARAMETERS ON EACH ENTRY */

    if (*neq < 1) {
	s_wsfi(&io___3);
	do_fio(&c__1, (char *)&(*neq), (ftnlen)sizeof(integer));
	e_wsfi();
/* Writing concatenation */
	i__1[0] = 99, a__1[0] = "IN DEBDF, THE NUMBER OF EQUATIONS MUST BE A"
		" POSITIVE INTEGER.$$YOU HAVE CALLED THE CODE WITH NEQ = ";
	i__1[1] = 8, a__1[1] = xern1;
	s_cat(ch__1, a__1, i__1, &c__2, (ftnlen)107);
	xermsg_("SLATEC", "LSOD", ch__1, &c__6, &c__1, (ftnlen)6, (ftnlen)4, (
		ftnlen)107);
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
		do_fio(&c__1, (char *)&rtol[k], (ftnlen)sizeof(real));
		e_wsfi();
/* Writing concatenation */
		i__3[0] = 98, a__2[0] = "IN DEBDF, THE RELATIVE ERROR TOLERA"
			"NCES MUST BE NON-NEGATIVE.$$YOU HAVE CALLED THE CODE"
			" WITH RTOL(";
		i__3[1] = 8, a__2[1] = xern1;
		i__3[2] = 4, a__2[2] = ") = ";
		i__3[3] = 16, a__2[3] = xern3;
		i__3[4] = 9, a__2[4] = "$$IN THE ";
		i__3[5] = 44, a__2[5] = "CASE OF VECTOR ERROR TOLERANCES, NO"
			" FURTHER ";
		i__3[6] = 36, a__2[6] = "CHECKING OF RTOL COMPONENTS IS DONE."
			;
		s_cat(ch__2, a__2, i__3, &c__7, (ftnlen)215);
		xermsg_("SLATEC", "LSOD", ch__2, &c__7, &c__1, (ftnlen)6, (
			ftnlen)4, (ftnlen)215);
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
	    do_fio(&c__1, (char *)&atol[k], (ftnlen)sizeof(real));
	    e_wsfi();
/* Writing concatenation */
	    i__4[0] = 98, a__3[0] = "IN DEBDF, THE ABSOLUTE ERROR TOLERANCES"
		    " MUST BE NON-NEGATIVE.$$YOU HAVE CALLED THE CODE WITH AT"
		    "OL(";
	    i__4[1] = 8, a__3[1] = xern1;
	    i__4[2] = 4, a__3[2] = ") = ";
	    i__4[3] = 16, a__3[3] = xern3;
	    i__4[4] = 53, a__3[4] = "$$IN THE CASE OF VECTOR ERROR TOLERANCE"
		    "S, NO FURTHER ";
	    i__4[5] = 36, a__3[5] = "CHECKING OF ATOL COMPONENTS IS DONE.";
	    s_cat(ch__2, a__3, i__4, &c__6, (ftnlen)215);
	    xermsg_("SLATEC", "LSOD", ch__2, &c__8, &c__1, (ftnlen)6, (ftnlen)
		    4, (ftnlen)215);
	    *idid = -33;
	    if (nrtolp > 0) {
		goto L70;
	    }
	    natolp = 1;
	}
L50:
	if (debdf1_1.itol == 0) {
	    goto L70;
	}
/* L60: */
    }

L70:
    if (debdf1_1.itstop == 1) {
	r__3 = *tout - *t;
	r__4 = *tstop - *t;
	if (r_sign(&c_b41, &r__3) != r_sign(&c_b41, &r__4) || (r__1 = *tout - 
		*t, dabs(r__1)) > (r__2 = *tstop - *t, dabs(r__2))) {
	    s_wsfi(&io___12);
	    do_fio(&c__1, (char *)&(*tout), (ftnlen)sizeof(real));
	    e_wsfi();
	    s_wsfi(&io___14);
	    do_fio(&c__1, (char *)&(*tstop), (ftnlen)sizeof(real));
	    e_wsfi();
/* Writing concatenation */
	    i__5[0] = 47, a__4[0] = "IN DEBDF, YOU HAVE CALLED THE CODE WITH"
		    " TOUT = ";
	    i__5[1] = 16, a__4[1] = xern3;
	    i__5[2] = 15, a__4[2] = "$$BUT YOU HAVE ";
	    i__5[3] = 51, a__4[3] = "ALSO TOLD THE CODE NOT TO INTEGRATE PAS"
		    "T THE POINT ";
	    i__5[4] = 8, a__4[4] = "TSTOP = ";
	    i__5[5] = 16, a__4[5] = xern4;
	    i__5[6] = 26, a__4[6] = " BY SETTING INFO(4) = 1.  ";
	    i__5[7] = 28, a__4[7] = "THESE INSTRUCTIONS CONFLICT.";
	    s_cat(ch__3, a__4, i__5, &c__8, (ftnlen)207);
	    xermsg_("SLATEC", "LSOD", ch__3, &c__14, &c__1, (ftnlen)6, (
		    ftnlen)4, (ftnlen)207);
	    *idid = -33;
	}
    }

/*        CHECK SOME CONTINUATION POSSIBILITIES */

    if (debdf1_1.init != 0) {
	if (*t == *tout) {
	    s_wsfi(&io___15);
	    do_fio(&c__1, (char *)&(*t), (ftnlen)sizeof(real));
	    e_wsfi();
/* Writing concatenation */
	    i__6[0] = 51, a__5[0] = "IN DEBDF, YOU HAVE CALLED THE CODE WITH"
		    " T = TOUT = ";
	    i__6[1] = 16, a__5[1] = xern3;
	    i__6[2] = 44, a__5[2] = "  THIS IS NOT ALLOWED ON CONTINUATION C"
		    "ALLS.";
	    s_cat(ch__4, a__5, i__6, &c__3, (ftnlen)111);
	    xermsg_("SLATEC", "LSOD", ch__4, &c__9, &c__1, (ftnlen)6, (ftnlen)
		    4, (ftnlen)111);
	    *idid = -33;
	}

	if (*t != debdf1_1.told) {
	    s_wsfi(&io___16);
	    do_fio(&c__1, (char *)&debdf1_1.told, (ftnlen)sizeof(real));
	    e_wsfi();
	    s_wsfi(&io___17);
	    do_fio(&c__1, (char *)&(*t), (ftnlen)sizeof(real));
	    e_wsfi();
/* Writing concatenation */
	    i__7[0] = 47, a__6[0] = "IN DEBDF, YOU HAVE CHANGED THE VALUE OF"
		    " T FROM ";
	    i__7[1] = 16, a__6[1] = xern3;
	    i__7[2] = 4, a__6[2] = " TO ";
	    i__7[3] = 16, a__6[3] = xern4;
	    i__7[4] = 44, a__6[4] = "  THIS IS NOT ALLOWED ON CONTINUATION C"
		    "ALLS.";
	    s_cat(ch__5, a__6, i__7, &c__5, (ftnlen)127);
	    xermsg_("SLATEC", "LSOD", ch__5, &c__10, &c__1, (ftnlen)6, (
		    ftnlen)4, (ftnlen)127);
	    *idid = -33;
	}

	if (debdf1_1.init != 1) {
	    if (*delsgn * (*tout - *t) < 0.f) {
		s_wsfi(&io___18);
		do_fio(&c__1, (char *)&(*tout), (ftnlen)sizeof(real));
		e_wsfi();
/* Writing concatenation */
		i__7[0] = 42, a__6[0] = "IN DEBDF, BY CALLING THE CODE WITH "
			"TOUT = ";
		i__7[1] = 16, a__6[1] = xern3;
		i__7[2] = 34, a__6[2] = " YOU ARE ATTEMPTING TO CHANGE THE ";
		i__7[3] = 27, a__6[3] = "DIRECTION OF INTEGRATION.$$";
		i__7[4] = 39, a__6[4] = "THIS IS NOT ALLOWED WITHOUT RESTART"
			"ING.";
		s_cat(ch__6, a__6, i__7, &c__5, (ftnlen)158);
		xermsg_("SLATEC", "LSOD", ch__6, &c__11, &c__1, (ftnlen)6, (
			ftnlen)4, (ftnlen)158);
		*idid = -33;
	    }
	}
    }

    if (*idid == -33) {
	if (debdf1_1.iquit != -33) {
/*                       INVALID INPUT DETECTED */
	    debdf1_1.iquit = -33;
	    debdf1_1.ibegin = -1;
	} else {
	    xermsg_("SLATEC", "LSOD", "IN DEBDF, INVALID INPUT WAS DETECTED "
		    "ON SUCCESSIVE ENTRIES.  IT IS IMPOSSIBLE TO PROCEED BECA"
		    "USE YOU HAVE NOT CORRECTED THE PROBLEM, SO EXECUTION IS "
		    "BEING TERMINATED.", &c__12, &c__2, (ftnlen)6, (ftnlen)4, (
		    ftnlen)166);
	}
	return 0;
    }

/* ....................................................................... */

/*     RTOL = ATOL = 0. IS ALLOWED AS VALID INPUT AND INTERPRETED AS */
/*     ASKING FOR THE MOST ACCURATE SOLUTION POSSIBLE. IN THIS CASE, */
/*     THE RELATIVE ERROR TOLERANCE RTOL IS RESET TO THE SMALLEST VALUE */
/*     100*U WHICH IS LIKELY TO BE REASONABLE FOR THIS METHOD AND MACHINE */

    i__2 = *neq;
    for (k = 1; k <= i__2; ++k) {
	if (rtol[k] + atol[k] > 0.f) {
	    goto L160;
	}
	rtol[k] = debdf1_1.u * 100.f;
	*idid = -2;
L160:
	if (debdf1_1.itol == 0) {
	    goto L180;
	}
/* L170: */
    }

L180:
    if (*idid != -2) {
	goto L190;
    }
/*                       RTOL=ATOL=0 ON INPUT, SO RTOL IS CHANGED TO A */
/*                                                SMALL POSITIVE VALUE */
    debdf1_1.ibegin = -1;
    return 0;

/*     BRANCH ON STATUS OF INITIALIZATION INDICATOR */
/*            INIT=0 MEANS INITIAL DERIVATIVES AND NOMINAL STEP SIZE */
/*                   AND DIRECTION NOT YET SET */
/*            INIT=1 MEANS NOMINAL STEP SIZE AND DIRECTION NOT YET SET */
/*            INIT=2 MEANS NO FURTHER INITIALIZATION REQUIRED */

L190:
    if (debdf1_1.init == 0) {
	goto L200;
    }
    if (debdf1_1.init == 1) {
	goto L220;
    }
    goto L240;

/* ....................................................................... */

/*     MORE INITIALIZATION -- */
/*                         -- EVALUATE INITIAL DERIVATIVES */

L200:
    debdf1_1.init = 1;
    (*f)(t, &y[1], &yh[(yh_dim1 << 1) + 1], &rpar[1], &ipar[1]);
    debdf1_1.nfe = 1;
    if (*t != *tout) {
	goto L220;
    }
    *idid = 2;
    i__2 = *neq;
    for (l = 1; l <= i__2; ++l) {
/* L210: */
	ypout[l] = yh[l + (yh_dim1 << 1)];
    }
    debdf1_1.told = *t;
    return 0;

/*                         -- COMPUTE INITIAL STEP SIZE */
/*                         -- SAVE SIGN OF INTEGRATION DIRECTION */
/*                         -- SET INDEPENDENT AND DEPENDENT VARIABLES */
/*                                              X AND YH(*) FOR STOD */

L220:
    ltol = 1;
    i__2 = *neq;
    for (l = 1; l <= i__2; ++l) {
	if (debdf1_1.itol == 1) {
	    ltol = l;
	}
	tol = rtol[ltol] * (r__1 = y[l], dabs(r__1)) + atol[ltol];
	if (tol == 0.f) {
	    goto L380;
	}
/* L225: */
	ewt[l] = tol;
    }

    big = sqrt(r1mach_(&c__2));
    hstart_((S_fp)f, neq, t, tout, &y[1], &yh[(yh_dim1 << 1) + 1], &ewt[1], &
	    c__1, &debdf1_1.u, &big, &yh[yh_dim1 * 3 + 1], &yh[(yh_dim1 << 2) 
	    + 1], &yh[yh_dim1 * 5 + 1], &yh[yh_dim1 * 6 + 1], &rpar[1], &ipar[
	    1], &debdf1_1.h__);

    r__1 = *tout - *t;
    *delsgn = r_sign(&c_b41, &r__1);
    debdf1_1.x = *t;
    i__2 = *neq;
    for (l = 1; l <= i__2; ++l) {
	yh[l + yh_dim1] = y[l];
/* L230: */
	yh[l + (yh_dim1 << 1)] = debdf1_1.h__ * yh[l + (yh_dim1 << 1)];
    }
    debdf1_1.init = 2;

/* ....................................................................... */

/*   ON EACH CALL SET INFORMATION WHICH DETERMINES THE ALLOWED INTERVAL */
/*   OF INTEGRATION BEFORE RETURNING WITH AN ANSWER AT TOUT */

L240:
    del = *tout - *t;
    absdel = dabs(del);

/* ....................................................................... */

/*   IF ALREADY PAST OUTPUT POINT, INTERPOLATE AND RETURN */

L250:
    if ((r__1 = debdf1_1.x - *t, dabs(r__1)) < absdel) {
	goto L270;
    }
    intyd_(tout, &c__0, &yh[yh_offset], neq, &y[1], &intflg);
    intyd_(tout, &c__1, &yh[yh_offset], neq, &ypout[1], &intflg);
    *idid = 3;
    if (debdf1_1.x != *tout) {
	goto L260;
    }
    *idid = 2;
    *intout = FALSE_;
L260:
    *t = *tout;
    debdf1_1.told = *t;
    return 0;

/*   IF CANNOT GO PAST TSTOP AND SUFFICIENTLY CLOSE, */
/*   EXTRAPOLATE AND RETURN */

L270:
    if (debdf1_1.itstop != 1) {
	goto L290;
    }
    if ((r__1 = *tstop - debdf1_1.x, dabs(r__1)) >= debdf1_1.u * 100.f * dabs(
	    debdf1_1.x)) {
	goto L290;
    }
    dt = *tout - debdf1_1.x;
    i__2 = *neq;
    for (l = 1; l <= i__2; ++l) {
/* L280: */
	y[l] = yh[l + yh_dim1] + dt / debdf1_1.h__ * yh[l + (yh_dim1 << 1)];
    }
    (*f)(tout, &y[1], &ypout[1], &rpar[1], &ipar[1]);
    ++debdf1_1.nfe;
    *idid = 3;
    *t = *tout;
    debdf1_1.told = *t;
    return 0;

L290:
    if (debdf1_1.iinteg == 0 || ! (*intout)) {
	goto L300;
    }

/*   INTERMEDIATE-OUTPUT MODE */

    *idid = 1;
    goto L500;

/* ....................................................................... */

/*     MONITOR NUMBER OF STEPS ATTEMPTED */

L300:
    if (debdf1_1.ksteps <= maxnum) {
	goto L330;
    }

/*                       A SIGNIFICANT AMOUNT OF WORK HAS BEEN EXPENDED */
    *idid = -1;
    debdf1_1.ksteps = 0;
    debdf1_1.ibegin = -1;
    goto L500;

/* ....................................................................... */

/*   LIMIT STEP SIZE AND SET WEIGHT VECTOR */

L330:
    debdf1_1.hmin = debdf1_1.u * 100.f * dabs(debdf1_1.x);
/* Computing MAX */
    r__1 = dabs(debdf1_1.h__);
    ha = dmax(r__1,debdf1_1.hmin);
    if (debdf1_1.itstop != 1) {
	goto L340;
    }
/* Computing MIN */
    r__2 = ha, r__3 = (r__1 = *tstop - debdf1_1.x, dabs(r__1));
    ha = dmin(r__2,r__3);
L340:
    debdf1_1.h__ = r_sign(&ha, &debdf1_1.h__);
    ltol = 1;
    i__2 = *neq;
    for (l = 1; l <= i__2; ++l) {
	if (debdf1_1.itol == 1) {
	    ltol = l;
	}
	ewt[l] = rtol[ltol] * (r__1 = yh[l + yh_dim1], dabs(r__1)) + atol[
		ltol];
	if (ewt[l] <= 0.f) {
	    goto L380;
	}
/* L350: */
    }
    *tolfac = debdf1_1.u * vnwrms_(neq, &yh[yh_offset], &ewt[1]);
    if (*tolfac <= 1.f) {
	goto L400;
    }

/*                       TOLERANCES TOO SMALL */
    *idid = -2;
    *tolfac *= 2.f;
    rtol[1] = *tolfac * rtol[1];
    atol[1] = *tolfac * atol[1];
    if (debdf1_1.itol == 0) {
	goto L370;
    }
    i__2 = *neq;
    for (l = 2; l <= i__2; ++l) {
	rtol[l] = *tolfac * rtol[l];
/* L360: */
	atol[l] = *tolfac * atol[l];
    }
L370:
    debdf1_1.ibegin = -1;
    goto L500;

/*                       RELATIVE ERROR CRITERION INAPPROPRIATE */
L380:
    *idid = -3;
    debdf1_1.ibegin = -1;
    goto L500;

/* ....................................................................... */

/*     TAKE A STEP */

L400:
    stod_(neq, &y[1], &yh[yh_offset], neq, &yh1[1], &ewt[1], &savf[1], &acor[
	    1], &wm[1], &iwm[1], (S_fp)f, (U_fp)jac, &rpar[1], &ipar[1]);

    debdf1_1.jstart = -2;
    *intout = TRUE_;
    if (debdf1_1.kflag == 0) {
	goto L250;
    }

/* ....................................................................... */

    if (debdf1_1.kflag == -1) {
	goto L450;
    }

/*                       REPEATED CORRECTOR CONVERGENCE FAILURES */
    *idid = -6;
    debdf1_1.ibegin = -1;
    goto L500;

/*                       REPEATED ERROR TEST FAILURES */
L450:
    *idid = -7;
    debdf1_1.ibegin = -1;

/* ....................................................................... */

/*                       STORE VALUES BEFORE RETURNING TO DEBDF */
L500:
    i__2 = *neq;
    for (l = 1; l <= i__2; ++l) {
	y[l] = yh[l + yh_dim1];
/* L555: */
	ypout[l] = yh[l + (yh_dim1 << 1)] / debdf1_1.h__;
    }
    *t = debdf1_1.x;
    debdf1_1.told = *t;
    *intout = FALSE_;
    return 0;
} /* lsod_ */

