/* ddes.f -- translated by f2c (version 12.02.01).
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
static integer c__1 = 1;
static integer c__3 = 3;
static integer c__2 = 2;
static integer c__5 = 5;
static integer c__14 = 14;
static integer c__6 = 6;
static integer c__7 = 7;
static integer c__8 = 8;
static doublereal c_b71 = 1.;
static integer c__9 = 9;
static integer c__10 = 10;
static integer c__11 = 11;
static integer c__12 = 12;

/* DECK DDES */
/* Subroutine */ int ddes_(S_fp df, integer *neq, doublereal *t, doublereal *
	y, doublereal *tout, integer *info, doublereal *rtol, doublereal *
	atol, integer *idid, doublereal *ypout, doublereal *yp, doublereal *
	yy, doublereal *wt, doublereal *p, doublereal *phi, doublereal *alpha,
	 doublereal *beta, doublereal *psi, doublereal *v, doublereal *w, 
	doublereal *sig, doublereal *g, doublereal *gi, doublereal *h__, 
	doublereal *eps, doublereal *x, doublereal *xold, doublereal *hold, 
	doublereal *told, doublereal *delsgn, doublereal *tstop, doublereal *
	twou, doublereal *fouru, logical *start, logical *phase1, logical *
	nornd, logical *stiff, logical *intout, integer *ns, integer *kord, 
	integer *kold, integer *init, integer *ksteps, integer *kle4, integer 
	*iquit, integer *kprev, integer *ivc, integer *iv, integer *kgi, 
	doublereal *rpar, integer *ipar)
{
    /* Initialized data */

    static integer maxnum = 500;

    /* System generated locals */
    address a__1[2], a__2[6], a__3[7], a__4[3], a__5[5];
    integer phi_dim1, phi_offset, i__1[2], i__2, i__3[6], i__4[7], i__5[3], 
	    i__6[5];
    doublereal d__1, d__2, d__3, d__4, d__5;
    char ch__1[221], ch__2[144], ch__3[166], ch__4[172], ch__5[114], ch__6[
	    223], ch__7[197], ch__8[113], ch__9[128], ch__10[159];

    /* Local variables */
    static doublereal a;
    static integer k, l;
    static doublereal u, ha, dt, del;
    static integer ltol;
    static char xern1[8], xern3[16], xern4[16];
    static logical crash;
    extern /* Subroutine */ int dintp_(doublereal *, doublereal *, doublereal 
	    *, doublereal *, doublereal *, integer *, integer *, doublereal *,
	     integer *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *);
    extern doublereal d1mach_(integer *);
    static doublereal absdel;
    static integer natolp;
    extern /* Subroutine */ int dsteps_(S_fp, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, logical *,
	     doublereal *, integer *, integer *, logical *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, logical *, integer *, logical *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, integer *, 
	    integer *, integer *, doublereal *, doublereal *, integer *), 
	    xermsg_(char *, char *, char *, integer *, integer *, ftnlen, 
	    ftnlen, ftnlen);
    static integer nrtolp;

    /* Fortran I/O blocks */
    static icilist io___4 = { 0, xern1, 0, "(I8)", 8, 1 };
    static icilist io___5 = { 0, xern1, 0, "(I8)", 8, 1 };
    static icilist io___6 = { 0, xern1, 0, "(I8)", 8, 1 };
    static icilist io___7 = { 0, xern1, 0, "(I8)", 8, 1 };
    static icilist io___8 = { 0, xern1, 0, "(I8)", 8, 1 };
    static icilist io___12 = { 0, xern1, 0, "(I8)", 8, 1 };
    static icilist io___14 = { 0, xern3, 0, "(1PE15.6)", 16, 1 };
    static icilist io___15 = { 0, xern1, 0, "(I8)", 8, 1 };
    static icilist io___16 = { 0, xern3, 0, "(1PE15.6)", 16, 1 };
    static icilist io___17 = { 0, xern3, 0, "(1PE15.6)", 16, 1 };
    static icilist io___19 = { 0, xern4, 0, "(1PE15.6)", 16, 1 };
    static icilist io___20 = { 0, xern3, 0, "(1PE15.6)", 16, 1 };
    static icilist io___21 = { 0, xern3, 0, "(1PE15.6)", 16, 1 };
    static icilist io___22 = { 0, xern4, 0, "(1PE15.6)", 16, 1 };
    static icilist io___23 = { 0, xern3, 0, "(1PE15.6)", 16, 1 };


/* ***BEGIN PROLOGUE  DDES */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to DDEABM */
/* ***LIBRARY   SLATEC */
/* ***TYPE      DOUBLE PRECISION (DES-S, DDES-D) */
/* ***AUTHOR  Watts, H. A., (SNLA) */
/* ***DESCRIPTION */

/*   DDEABM merely allocates storage for DDES to relieve the user of the */
/*   inconvenience of a long call list.  Consequently  DDES  is used as */
/*   described in the comments for  DDEABM . */

/* ***SEE ALSO  DDEABM */
/* ***ROUTINES CALLED  D1MACH, DINTP, DSTEPS, XERMSG */
/* ***REVISION HISTORY  (YYMMDD) */
/*   820301  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890831  Modified array declarations.  (WRB) */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900328  Added TYPE section.  (WRB) */
/*   900510  Convert XERRWV calls to XERMSG calls, cvt GOTOs to */
/*           IF-THEN-ELSE.  (RWC) */
/*   910722  Updated AUTHOR section.  (ALS) */
/* ***END PROLOGUE  DDES */




/* ....................................................................... */

/*  THE EXPENSE OF SOLVING THE PROBLEM IS MONITORED BY COUNTING THE */
/*  NUMBER OF  STEPS ATTEMPTED. WHEN THIS EXCEEDS  MAXNUM, THE COUNTER */
/*  IS RESET TO ZERO AND THE USER IS INFORMED ABOUT POSSIBLE EXCESSIVE */
/*  WORK. */

    /* Parameter adjustments */
    phi_dim1 = *neq;
    phi_offset = 1 + phi_dim1;
    phi -= phi_offset;
    --y;
    --info;
    --rtol;
    --atol;
    --ypout;
    --yp;
    --yy;
    --wt;
    --p;
    --alpha;
    --beta;
    --psi;
    --v;
    --w;
    --sig;
    --g;
    --gi;
    --iv;
    --rpar;
    --ipar;

    /* Function Body */

/* ....................................................................... */

/* ***FIRST EXECUTABLE STATEMENT  DDES */
    if (info[1] == 0) {

/* ON THE FIRST CALL , PERFORM INITIALIZATION -- */
/*        DEFINE THE MACHINE UNIT ROUNDOFF QUANTITY  U  BY CALLING THE */
/*        FUNCTION ROUTINE  D1MACH. THE USER MUST MAKE SURE THAT THE */
/*        VALUES SET IN D1MACH ARE RELEVANT TO THE COMPUTER BEING USED. */

	u = d1mach_(&c__4);
/*                       -- SET ASSOCIATED MACHINE DEPENDENT PARAMETERS */
	*twou = u * 2.;
	*fouru = u * 4.;
/*                       -- SET TERMINATION FLAG */
	*iquit = 0;
/*                       -- SET INITIALIZATION INDICATOR */
	*init = 0;
/*                       -- SET COUNTER FOR ATTEMPTED STEPS */
	*ksteps = 0;
/*                       -- SET INDICATOR FOR INTERMEDIATE-OUTPUT */
	*intout = FALSE_;
/*                       -- SET INDICATOR FOR STIFFNESS DETECTION */
	*stiff = FALSE_;
/*                       -- SET STEP COUNTER FOR STIFFNESS DETECTION */
	*kle4 = 0;
/*                       -- SET INDICATORS FOR STEPS CODE */
	*start = TRUE_;
	*phase1 = TRUE_;
	*nornd = TRUE_;
/*                       -- RESET INFO(1) FOR SUBSEQUENT CALLS */
	info[1] = 1;
    }

/* ....................................................................... */

/*      CHECK VALIDITY OF INPUT PARAMETERS ON EACH ENTRY */

    if (info[1] != 0 && info[1] != 1) {
	s_wsfi(&io___4);
	do_fio(&c__1, (char *)&info[1], (ftnlen)sizeof(integer));
	e_wsfi();
/* Writing concatenation */
	i__1[0] = 213, a__1[0] = "IN DDEABM, INFO(1) MUST BE SET TO 0 FOR TH"
		"E START OF A NEW PROBLEM, AND MUST BE SET TO 1 FOLLOWING AN "
		"INTERRUPTED TASK.  YOU ARE ATTEMPTING TO CONTINUE THE INTEGR"
		"ATION ILLEGALLY BY CALLING THE CODE WITH INFO(1) = ";
	i__1[1] = 8, a__1[1] = xern1;
	s_cat(ch__1, a__1, i__1, &c__2, (ftnlen)221);
	xermsg_("SLATEC", "DDES", ch__1, &c__3, &c__1, (ftnlen)6, (ftnlen)4, (
		ftnlen)221);
	*idid = -33;
    }

    if (info[2] != 0 && info[2] != 1) {
	s_wsfi(&io___5);
	do_fio(&c__1, (char *)&info[2], (ftnlen)sizeof(integer));
	e_wsfi();
/* Writing concatenation */
	i__1[0] = 136, a__1[0] = "IN DDEABM, INFO(2) MUST BE 0 OR 1 INDICATI"
		"NG SCALAR AND VECTOR ERROR TOLERANCES, RESPECTIVELY.  YOU HA"
		"VE CALLED THE CODE WITH INFO(2) = ";
	i__1[1] = 8, a__1[1] = xern1;
	s_cat(ch__2, a__1, i__1, &c__2, (ftnlen)144);
	xermsg_("SLATEC", "DDES", ch__2, &c__4, &c__1, (ftnlen)6, (ftnlen)4, (
		ftnlen)144);
	*idid = -33;
    }

    if (info[3] != 0 && info[3] != 1) {
	s_wsfi(&io___6);
	do_fio(&c__1, (char *)&info[3], (ftnlen)sizeof(integer));
	e_wsfi();
/* Writing concatenation */
	i__1[0] = 158, a__1[0] = "IN DDEABM, INFO(3) MUST BE 0 OR 1 INDICATI"
		"NG THE INTERVAL OR INTERMEDIATE-OUTPUT MODE OF INTEGRATION, "
		"RESPECTIVELY.  YOU HAVE CALLED THE CODE WITH  INFO(3) = ";
	i__1[1] = 8, a__1[1] = xern1;
	s_cat(ch__3, a__1, i__1, &c__2, (ftnlen)166);
	xermsg_("SLATEC", "DDES", ch__3, &c__5, &c__1, (ftnlen)6, (ftnlen)4, (
		ftnlen)166);
	*idid = -33;
    }

    if (info[4] != 0 && info[4] != 1) {
	s_wsfi(&io___7);
	do_fio(&c__1, (char *)&info[4], (ftnlen)sizeof(integer));
	e_wsfi();
/* Writing concatenation */
	i__1[0] = 164, a__1[0] = "IN DDEABM, INFO(4) MUST BE 0 OR 1 INDICATI"
		"NG WHETHER OR NOT THE INTEGRATION INTERVAL IS TO BE RESTRICT"
		"ED BY A POINT TSTOP.  YOU HAVE CALLED THE CODE WITH INFO(4) "
		"= ";
	i__1[1] = 8, a__1[1] = xern1;
	s_cat(ch__4, a__1, i__1, &c__2, (ftnlen)172);
	xermsg_("SLATEC", "DDES", ch__4, &c__14, &c__1, (ftnlen)6, (ftnlen)4, 
		(ftnlen)172);
	*idid = -33;
    }

    if (*neq < 1) {
	s_wsfi(&io___8);
	do_fio(&c__1, (char *)&(*neq), (ftnlen)sizeof(integer));
	e_wsfi();
/* Writing concatenation */
	i__1[0] = 106, a__1[0] = "IN DDEABM,  THE NUMBER OF EQUATIONS NEQ MU"
		"ST BE A POSITIVE INTEGER.  YOU HAVE CALLED THE CODE WITH  NE"
		"Q = ";
	i__1[1] = 8, a__1[1] = xern1;
	s_cat(ch__5, a__1, i__1, &c__2, (ftnlen)114);
	xermsg_("SLATEC", "DDES", ch__5, &c__6, &c__1, (ftnlen)6, (ftnlen)4, (
		ftnlen)114);
	*idid = -33;
    }

    nrtolp = 0;
    natolp = 0;
    i__2 = *neq;
    for (k = 1; k <= i__2; ++k) {
	if (nrtolp == 0 && rtol[k] < 0.) {
	    s_wsfi(&io___12);
	    do_fio(&c__1, (char *)&k, (ftnlen)sizeof(integer));
	    e_wsfi();
	    s_wsfi(&io___14);
	    do_fio(&c__1, (char *)&rtol[k], (ftnlen)sizeof(doublereal));
	    e_wsfi();
/* Writing concatenation */
	    i__3[0] = 105, a__2[0] = "IN DDEABM, THE RELATIVE ERROR TOLERANC"
		    "ES RTOL MUST BE NON-NEGATIVE.  YOU HAVE CALLED THE CODE "
		    "WITH  RTOL(";
	    i__3[1] = 8, a__2[1] = xern1;
	    i__3[2] = 4, a__2[2] = ") = ";
	    i__3[3] = 16, a__2[3] = xern3;
	    i__3[4] = 43, a__2[4] = ".  IN THE CASE OF VECTOR ERROR TOLERANC"
		    "ES, ";
	    i__3[5] = 47, a__2[5] = "NO FURTHER CHECKING OF RTOL COMPONENTS "
		    "IS DONE.";
	    s_cat(ch__6, a__2, i__3, &c__6, (ftnlen)223);
	    xermsg_("SLATEC", "DDES", ch__6, &c__7, &c__1, (ftnlen)6, (ftnlen)
		    4, (ftnlen)223);
	    *idid = -33;
	    nrtolp = 1;
	}

	if (natolp == 0 && atol[k] < 0.) {
	    s_wsfi(&io___15);
	    do_fio(&c__1, (char *)&k, (ftnlen)sizeof(integer));
	    e_wsfi();
	    s_wsfi(&io___16);
	    do_fio(&c__1, (char *)&atol[k], (ftnlen)sizeof(doublereal));
	    e_wsfi();
/* Writing concatenation */
	    i__3[0] = 105, a__2[0] = "IN DDEABM, THE ABSOLUTE ERROR TOLERANC"
		    "ES ATOL MUST BE NON-NEGATIVE.  YOU HAVE CALLED THE CODE "
		    "WITH  ATOL(";
	    i__3[1] = 8, a__2[1] = xern1;
	    i__3[2] = 4, a__2[2] = ") = ";
	    i__3[3] = 16, a__2[3] = xern3;
	    i__3[4] = 43, a__2[4] = ".  IN THE CASE OF VECTOR ERROR TOLERANC"
		    "ES, ";
	    i__3[5] = 47, a__2[5] = "NO FURTHER CHECKING OF ATOL COMPONENTS "
		    "IS DONE.";
	    s_cat(ch__6, a__2, i__3, &c__6, (ftnlen)223);
	    xermsg_("SLATEC", "DDES", ch__6, &c__8, &c__1, (ftnlen)6, (ftnlen)
		    4, (ftnlen)223);
	    *idid = -33;
	    natolp = 1;
	}

	if (info[2] == 0) {
	    goto L100;
	}
	if (natolp > 0 && nrtolp > 0) {
	    goto L100;
	}
/* L90: */
    }

L100:
    if (info[4] == 1) {
	d__3 = *tout - *t;
	d__4 = *tstop - *t;
	if (d_sign(&c_b71, &d__3) != d_sign(&c_b71, &d__4) || (d__1 = *tout - 
		*t, abs(d__1)) > (d__2 = *tstop - *t, abs(d__2))) {
	    s_wsfi(&io___17);
	    do_fio(&c__1, (char *)&(*tout), (ftnlen)sizeof(doublereal));
	    e_wsfi();
	    s_wsfi(&io___19);
	    do_fio(&c__1, (char *)&(*tstop), (ftnlen)sizeof(doublereal));
	    e_wsfi();
/* Writing concatenation */
	    i__4[0] = 49, a__3[0] = "IN DDEABM, YOU HAVE CALLED THE CODE WIT"
		    "H  TOUT = ";
	    i__4[1] = 16, a__3[1] = xern3;
	    i__4[2] = 5, a__3[2] = " BUT ";
	    i__4[3] = 49, a__3[3] = "YOU HAVE ALSO TOLD THE CODE (INFO(4) = "
		    "1) NOT TO ";
	    i__4[4] = 33, a__3[4] = "INTEGRATE PAST THE POINT TSTOP = ";
	    i__4[5] = 16, a__3[5] = xern4;
	    i__4[6] = 29, a__3[6] = " THESE INSTRUCTIONS CONFLICT.";
	    s_cat(ch__7, a__3, i__4, &c__7, (ftnlen)197);
	    xermsg_("SLATEC", "DDES", ch__7, &c__14, &c__1, (ftnlen)6, (
		    ftnlen)4, (ftnlen)197);
	    *idid = -33;
	}
    }

/*     CHECK SOME CONTINUATION POSSIBILITIES */

    if (*init != 0) {
	if (*t == *tout) {
	    s_wsfi(&io___20);
	    do_fio(&c__1, (char *)&(*t), (ftnlen)sizeof(doublereal));
	    e_wsfi();
/* Writing concatenation */
	    i__5[0] = 53, a__4[0] = "IN DDEABM, YOU HAVE CALLED THE CODE WIT"
		    "H  T = TOUT = ";
	    i__5[1] = 16, a__4[1] = xern3;
	    i__5[2] = 44, a__4[2] = "$$THIS IS NOT ALLOWED ON CONTINUATION C"
		    "ALLS.";
	    s_cat(ch__8, a__4, i__5, &c__3, (ftnlen)113);
	    xermsg_("SLATEC", "DDES", ch__8, &c__9, &c__1, (ftnlen)6, (ftnlen)
		    4, (ftnlen)113);
	    *idid = -33;
	}

	if (*t != *told) {
	    s_wsfi(&io___21);
	    do_fio(&c__1, (char *)&(*told), (ftnlen)sizeof(doublereal));
	    e_wsfi();
	    s_wsfi(&io___22);
	    do_fio(&c__1, (char *)&(*t), (ftnlen)sizeof(doublereal));
	    e_wsfi();
/* Writing concatenation */
	    i__6[0] = 48, a__5[0] = "IN DDEABM, YOU HAVE CHANGED THE VALUE O"
		    "F T FROM ";
	    i__6[1] = 16, a__5[1] = xern3;
	    i__6[2] = 4, a__5[2] = " TO ";
	    i__6[3] = 16, a__5[3] = xern4;
	    i__6[4] = 44, a__5[4] = "  THIS IS NOT ALLOWED ON CONTINUATION C"
		    "ALLS.";
	    s_cat(ch__9, a__5, i__6, &c__5, (ftnlen)128);
	    xermsg_("SLATEC", "DDES", ch__9, &c__10, &c__1, (ftnlen)6, (
		    ftnlen)4, (ftnlen)128);
	    *idid = -33;
	}

	if (*init != 1) {
	    if (*delsgn * (*tout - *t) < 0.) {
		s_wsfi(&io___23);
		do_fio(&c__1, (char *)&(*tout), (ftnlen)sizeof(doublereal));
		e_wsfi();
/* Writing concatenation */
		i__6[0] = 43, a__5[0] = "IN DDEABM, BY CALLING THE CODE WITH"
			" TOUT = ";
		i__6[1] = 16, a__5[1] = xern3;
		i__6[2] = 47, a__5[2] = " YOU ARE ATTEMPTING TO CHANGE THE D"
			"IRECTION OF ";
		i__6[3] = 42, a__5[3] = "INTEGRATION.$$THIS IS NOT ALLOWED W"
			"ITHOUT ";
		i__6[4] = 11, a__5[4] = "RESTARTING.";
		s_cat(ch__10, a__5, i__6, &c__5, (ftnlen)159);
		xermsg_("SLATEC", "DDES", ch__10, &c__11, &c__1, (ftnlen)6, (
			ftnlen)4, (ftnlen)159);
		*idid = -33;
	    }
	}
    }

/*     INVALID INPUT DETECTED */

    if (*idid == -33) {
	if (*iquit != -33) {
	    *iquit = -33;
	    info[1] = -1;
	} else {
	    xermsg_("SLATEC", "DDES", "IN DDEABM, INVALID INPUT WAS DETECTED"
		    " ON SUCCESSIVE ENTRIES.  IT IS IMPOSSIBLE TO PROCEED BEC"
		    "AUSE YOU HAVE NOT CORRECTED THE PROBLEM, SO EXECUTION IS"
		    " BEING TERMINATED.", &c__12, &c__2, (ftnlen)6, (ftnlen)4, 
		    (ftnlen)167);
	}
	return 0;
    }

/* ....................................................................... */

/*     RTOL = ATOL = 0. IS ALLOWED AS VALID INPUT AND INTERPRETED AS */
/*     ASKING FOR THE MOST ACCURATE SOLUTION POSSIBLE. IN THIS CASE, */
/*     THE RELATIVE ERROR TOLERANCE RTOL IS RESET TO THE SMALLEST VALUE */
/*     FOURU WHICH IS LIKELY TO BE REASONABLE FOR THIS METHOD AND MACHINE */

    i__2 = *neq;
    for (k = 1; k <= i__2; ++k) {
	if (rtol[k] + atol[k] > 0.) {
	    goto L170;
	}
	rtol[k] = *fouru;
	*idid = -2;
L170:
	if (info[2] == 0) {
	    goto L190;
	}
/* L180: */
    }

L190:
    if (*idid != -2) {
	goto L200;
    }
/*                       RTOL=ATOL=0 ON INPUT, SO RTOL IS CHANGED TO A */
/*                                                SMALL POSITIVE VALUE */
    info[1] = -1;
    return 0;

/*     BRANCH ON STATUS OF INITIALIZATION INDICATOR */
/*            INIT=0 MEANS INITIAL DERIVATIVES AND NOMINAL STEP SIZE */
/*                   AND DIRECTION NOT YET SET */
/*            INIT=1 MEANS NOMINAL STEP SIZE AND DIRECTION NOT YET SET */
/*            INIT=2 MEANS NO FURTHER INITIALIZATION REQUIRED */

L200:
    if (*init == 0) {
	goto L210;
    }
    if (*init == 1) {
	goto L220;
    }
    goto L240;

/* ....................................................................... */

/*     MORE INITIALIZATION -- */
/*                         -- EVALUATE INITIAL DERIVATIVES */

L210:
    *init = 1;
    a = *t;
    (*df)(&a, &y[1], &yp[1], &rpar[1], &ipar[1]);
    if (*t != *tout) {
	goto L220;
    }
    *idid = 2;
    i__2 = *neq;
    for (l = 1; l <= i__2; ++l) {
/* L215: */
	ypout[l] = yp[l];
    }
    *told = *t;
    return 0;

/*                         -- SET INDEPENDENT AND DEPENDENT VARIABLES */
/*                                              X AND YY(*) FOR STEPS */
/*                         -- SET SIGN OF INTEGRATION DIRECTION */
/*                         -- INITIALIZE THE STEP SIZE */

L220:
    *init = 2;
    *x = *t;
    i__2 = *neq;
    for (l = 1; l <= i__2; ++l) {
/* L230: */
	yy[l] = y[l];
    }
    d__1 = *tout - *t;
    *delsgn = d_sign(&c_b71, &d__1);
/* Computing MAX */
    d__3 = *fouru * abs(*x), d__4 = (d__1 = *tout - *x, abs(d__1));
    d__2 = max(d__3,d__4);
    d__5 = *tout - *x;
    *h__ = d_sign(&d__2, &d__5);

/* ....................................................................... */

/*   ON EACH CALL SET INFORMATION WHICH DETERMINES THE ALLOWED INTERVAL */
/*   OF INTEGRATION BEFORE RETURNING WITH AN ANSWER AT TOUT */

L240:
    del = *tout - *t;
    absdel = abs(del);

/* ....................................................................... */

/*   IF ALREADY PAST OUTPUT POINT, INTERPOLATE AND RETURN */

L250:
    if ((d__1 = *x - *t, abs(d__1)) < absdel) {
	goto L260;
    }
    dintp_(x, &yy[1], tout, &y[1], &ypout[1], neq, kold, &phi[phi_offset], 
	    ivc, &iv[1], kgi, &gi[1], &alpha[1], &g[1], &w[1], xold, &p[1]);
    *idid = 3;
    if (*x != *tout) {
	goto L255;
    }
    *idid = 2;
    *intout = FALSE_;
L255:
    *t = *tout;
    *told = *t;
    return 0;

/*   IF CANNOT GO PAST TSTOP AND SUFFICIENTLY CLOSE, */
/*   EXTRAPOLATE AND RETURN */

L260:
    if (info[4] != 1) {
	goto L280;
    }
    if ((d__1 = *tstop - *x, abs(d__1)) >= *fouru * abs(*x)) {
	goto L280;
    }
    dt = *tout - *x;
    i__2 = *neq;
    for (l = 1; l <= i__2; ++l) {
/* L270: */
	y[l] = yy[l] + dt * yp[l];
    }
    (*df)(tout, &y[1], &ypout[1], &rpar[1], &ipar[1]);
    *idid = 3;
    *t = *tout;
    *told = *t;
    return 0;

L280:
    if (info[3] == 0 || ! (*intout)) {
	goto L300;
    }

/*   INTERMEDIATE-OUTPUT MODE */

    *idid = 1;
    i__2 = *neq;
    for (l = 1; l <= i__2; ++l) {
	y[l] = yy[l];
/* L290: */
	ypout[l] = yp[l];
    }
    *t = *x;
    *told = *t;
    *intout = FALSE_;
    return 0;

/* ....................................................................... */

/*     MONITOR NUMBER OF STEPS ATTEMPTED */

L300:
    if (*ksteps <= maxnum) {
	goto L330;
    }

/*                       A SIGNIFICANT AMOUNT OF WORK HAS BEEN EXPENDED */
    *idid = -1;
    *ksteps = 0;
    if (! (*stiff)) {
	goto L310;
    }

/*                       PROBLEM APPEARS TO BE STIFF */
    *idid = -4;
    *stiff = FALSE_;
    *kle4 = 0;

L310:
    i__2 = *neq;
    for (l = 1; l <= i__2; ++l) {
	y[l] = yy[l];
/* L320: */
	ypout[l] = yp[l];
    }
    *t = *x;
    *told = *t;
    info[1] = -1;
    *intout = FALSE_;
    return 0;

/* ....................................................................... */

/*   LIMIT STEP SIZE, SET WEIGHT VECTOR AND TAKE A STEP */

L330:
    ha = abs(*h__);
    if (info[4] != 1) {
	goto L340;
    }
/* Computing MIN */
    d__2 = ha, d__3 = (d__1 = *tstop - *x, abs(d__1));
    ha = min(d__2,d__3);
L340:
    *h__ = d_sign(&ha, h__);
    *eps = 1.;
    ltol = 1;
    i__2 = *neq;
    for (l = 1; l <= i__2; ++l) {
	if (info[2] == 1) {
	    ltol = l;
	}
	wt[l] = rtol[ltol] * (d__1 = yy[l], abs(d__1)) + atol[ltol];
	if (wt[l] <= 0.) {
	    goto L360;
	}
/* L350: */
    }
    goto L380;

/*                       RELATIVE ERROR CRITERION INAPPROPRIATE */
L360:
    *idid = -3;
    i__2 = *neq;
    for (l = 1; l <= i__2; ++l) {
	y[l] = yy[l];
/* L370: */
	ypout[l] = yp[l];
    }
    *t = *x;
    *told = *t;
    info[1] = -1;
    *intout = FALSE_;
    return 0;

L380:
    dsteps_((S_fp)df, neq, &yy[1], x, h__, eps, &wt[1], start, hold, kord, 
	    kold, &crash, &phi[phi_offset], &p[1], &yp[1], &psi[1], &alpha[1],
	     &beta[1], &sig[1], &v[1], &w[1], &g[1], phase1, ns, nornd, 
	    ksteps, twou, fouru, xold, kprev, ivc, &iv[1], kgi, &gi[1], &rpar[
	    1], &ipar[1]);

/* ....................................................................... */

    if (! crash) {
	goto L420;
    }

/*                       TOLERANCES TOO SMALL */
    *idid = -2;
    rtol[1] = *eps * rtol[1];
    atol[1] = *eps * atol[1];
    if (info[2] == 0) {
	goto L400;
    }
    i__2 = *neq;
    for (l = 2; l <= i__2; ++l) {
	rtol[l] = *eps * rtol[l];
/* L390: */
	atol[l] = *eps * atol[l];
    }
L400:
    i__2 = *neq;
    for (l = 1; l <= i__2; ++l) {
	y[l] = yy[l];
/* L410: */
	ypout[l] = yp[l];
    }
    *t = *x;
    *told = *t;
    info[1] = -1;
    *intout = FALSE_;
    return 0;

/*   (STIFFNESS TEST) COUNT NUMBER OF CONSECUTIVE STEPS TAKEN WITH THE */
/*   ORDER OF THE METHOD BEING LESS OR EQUAL TO FOUR */

L420:
    ++(*kle4);
    if (*kold > 4) {
	*kle4 = 0;
    }
    if (*kle4 >= 50) {
	*stiff = TRUE_;
    }
    *intout = TRUE_;
    goto L250;
} /* ddes_ */

