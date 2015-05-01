/* des.f -- translated by f2c (version 12.02.01).
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
static real c_b71 = 1.f;
static integer c__9 = 9;
static integer c__10 = 10;
static integer c__11 = 11;
static integer c__12 = 12;

/* DECK DES */
/* Subroutine */ int des_(S_fp f, integer *neq, real *t, real *y, real *tout, 
	integer *info, real *rtol, real *atol, integer *idid, real *ypout, 
	real *yp, real *yy, real *wt, real *p, real *phi, real *alpha, real *
	beta, real *psi, real *v, real *w, real *sig, real *g, real *gi, real 
	*h__, real *eps, real *x, real *xold, real *hold, real *told, real *
	delsgn, real *tstop, real *twou, real *fouru, logical *start, logical 
	*phase1, logical *nornd, logical *stiff, logical *intout, integer *ns,
	 integer *kord, integer *kold, integer *init, integer *ksteps, 
	integer *kle4, integer *iquit, integer *kprev, integer *ivc, integer *
	iv, integer *kgi, real *rpar, integer *ipar)
{
    /* Initialized data */

    static integer maxnum = 500;

    /* System generated locals */
    address a__1[2], a__2[6], a__3[8], a__4[4], a__5[5];
    integer phi_dim1, phi_offset, i__1[2], i__2, i__3[6], i__4[8], i__5[4], 
	    i__6[5];
    real r__1, r__2, r__3, r__4, r__5;
    char ch__1[220], ch__2[143], ch__3[165], ch__4[171], ch__5[112], ch__6[
	    222], ch__7[196], ch__8[127], ch__9[158];

    /* Local variables */
    static real a;
    static integer k, l;
    static real u, ha, dt, del;
    static integer ltol;
    static char xern1[8], xern3[16], xern4[16];
    static logical crash;
    extern /* Subroutine */ int steps_(S_fp, integer *, real *, real *, real *
	    , real *, real *, logical *, real *, integer *, integer *, 
	    logical *, real *, real *, real *, real *, real *, real *, real *,
	     real *, real *, real *, logical *, integer *, logical *, integer 
	    *, real *, real *, real *, integer *, integer *, integer *, 
	    integer *, real *, real *, integer *);
    extern doublereal r1mach_(integer *);
    static real absdel;
    static integer natolp;
    extern /* Subroutine */ int xermsg_(char *, char *, char *, integer *, 
	    integer *, ftnlen, ftnlen, ftnlen);
    static integer nrtolp;
    extern /* Subroutine */ int sintrp_(real *, real *, real *, real *, real *
	    , integer *, integer *, real *, integer *, integer *, integer *, 
	    real *, real *, real *, real *, real *, real *);

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


/* ***BEGIN PROLOGUE  DES */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to DEABM */
/* ***LIBRARY   SLATEC */
/* ***TYPE      SINGLE PRECISION (DES-S, DDES-D) */
/* ***AUTHOR  Watts, H. A., (SNLA) */
/* ***DESCRIPTION */

/*   DEABM  merely allocates storage for  DES  to relieve the user of the */
/*   inconvenience of a long call list.  Consequently  DES  is used as */
/*   described in the comments for  DEABM . */

/* ***SEE ALSO  DEABM */
/* ***ROUTINES CALLED  R1MACH, SINTRP, STEPS, XERMSG */
/* ***REVISION HISTORY  (YYMMDD) */
/*   800501  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890831  Modified array declarations.  (WRB) */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900328  Added TYPE section.  (WRB) */
/*   900510  Convert XERRWV calls to XERMSG calls, replace GOTOs with */
/*           IF-THEN-ELSEs.  (RWC) */
/*   910722  Updated AUTHOR section.  (ALS) */
/* ***END PROLOGUE  DES */




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

/* ***FIRST EXECUTABLE STATEMENT  DES */
    if (info[1] == 0) {

/* ON THE FIRST CALL , PERFORM INITIALIZATION -- */
/*        DEFINE THE MACHINE UNIT ROUNDOFF QUANTITY  U  BY CALLING THE */
/*        FUNCTION ROUTINE  R1MACH. THE USER MUST MAKE SURE THAT THE */
/*        VALUES SET IN R1MACH ARE RELEVANT TO THE COMPUTER BEING USED. */

	u = r1mach_(&c__4);
/*                       -- SET ASSOCIATED MACHINE DEPENDENT PARAMETERS */
	*twou = u * 2.f;
	*fouru = u * 4.f;
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
	i__1[0] = 212, a__1[0] = "IN DEABM, INFO(1) MUST BE SET TO 0 FOR THE"
		" START OF A NEW PROBLEM, AND MUST BE SET TO 1 FOLLOWING AN I"
		"NTERRUPTED TASK.  YOU ARE ATTEMPTING TO CONTINUE THE INTEGRA"
		"TION ILLEGALLY BY CALLING THE CODE WITH INFO(1) = ";
	i__1[1] = 8, a__1[1] = xern1;
	s_cat(ch__1, a__1, i__1, &c__2, (ftnlen)220);
	xermsg_("SLATEC", "DES", ch__1, &c__3, &c__1, (ftnlen)6, (ftnlen)3, (
		ftnlen)220);
	*idid = -33;
    }

    if (info[2] != 0 && info[2] != 1) {
	s_wsfi(&io___5);
	do_fio(&c__1, (char *)&info[2], (ftnlen)sizeof(integer));
	e_wsfi();
/* Writing concatenation */
	i__1[0] = 135, a__1[0] = "IN DEABM, INFO(2) MUST BE 0 OR 1 INDICATIN"
		"G SCALAR AND VECTOR ERROR TOLERANCES, RESPECTIVELY.  YOU HAV"
		"E CALLED THE CODE WITH INFO(2) = ";
	i__1[1] = 8, a__1[1] = xern1;
	s_cat(ch__2, a__1, i__1, &c__2, (ftnlen)143);
	xermsg_("SLATEC", "DES", ch__2, &c__4, &c__1, (ftnlen)6, (ftnlen)3, (
		ftnlen)143);
	*idid = -33;
    }

    if (info[3] != 0 && info[3] != 1) {
	s_wsfi(&io___6);
	do_fio(&c__1, (char *)&info[3], (ftnlen)sizeof(integer));
	e_wsfi();
/* Writing concatenation */
	i__1[0] = 157, a__1[0] = "IN DEABM, INFO(3) MUST BE 0 OR 1 INDICATIN"
		"G THE INTERVAL OR INTERMEDIATE-OUTPUT MODE OF INTEGRATION, R"
		"ESPECTIVELY.  YOU HAVE CALLED THE CODE WITH  INFO(3) = ";
	i__1[1] = 8, a__1[1] = xern1;
	s_cat(ch__3, a__1, i__1, &c__2, (ftnlen)165);
	xermsg_("SLATEC", "DES", ch__3, &c__5, &c__1, (ftnlen)6, (ftnlen)3, (
		ftnlen)165);
	*idid = -33;
    }

    if (info[4] != 0 && info[4] != 1) {
	s_wsfi(&io___7);
	do_fio(&c__1, (char *)&info[4], (ftnlen)sizeof(integer));
	e_wsfi();
/* Writing concatenation */
	i__1[0] = 163, a__1[0] = "IN DEABM, INFO(4) MUST BE 0 OR 1 INDICATIN"
		"G WHETHER OR NOT THE INTEGRATION INTERVAL IS TO BE RESTRICTE"
		"D BY A POINT TSTOP.  YOU HAVE CALLED THE CODE WITH INFO(4) = "
		;
	i__1[1] = 8, a__1[1] = xern1;
	s_cat(ch__4, a__1, i__1, &c__2, (ftnlen)171);
	xermsg_("SLATEC", "DES", ch__4, &c__14, &c__1, (ftnlen)6, (ftnlen)3, (
		ftnlen)171);
	*idid = -33;
    }

    if (*neq < 1) {
	s_wsfi(&io___8);
	do_fio(&c__1, (char *)&(*neq), (ftnlen)sizeof(integer));
	e_wsfi();
/* Writing concatenation */
	i__1[0] = 104, a__1[0] = "IN DEABM, THE NUMBER OF EQUATIONS NEQ MUST"
		" BE A POSITIVE INTEGER.  YOU HAVE CALLED THE CODE WITH  NEQ "
		"= ";
	i__1[1] = 8, a__1[1] = xern1;
	s_cat(ch__5, a__1, i__1, &c__2, (ftnlen)112);
	xermsg_("SLATEC", "DES", ch__5, &c__6, &c__1, (ftnlen)6, (ftnlen)3, (
		ftnlen)112);
	*idid = -33;
    }

    nrtolp = 0;
    natolp = 0;
    i__2 = *neq;
    for (k = 1; k <= i__2; ++k) {
	if (nrtolp == 0 && rtol[k] < 0.f) {
	    s_wsfi(&io___12);
	    do_fio(&c__1, (char *)&k, (ftnlen)sizeof(integer));
	    e_wsfi();
	    s_wsfi(&io___14);
	    do_fio(&c__1, (char *)&rtol[k], (ftnlen)sizeof(real));
	    e_wsfi();
/* Writing concatenation */
	    i__3[0] = 104, a__2[0] = "IN DEABM, THE RELATIVE ERROR TOLERANCE"
		    "S RTOL MUST BE NON-NEGATIVE.  YOU HAVE CALLED THE CODE W"
		    "ITH  RTOL(";
	    i__3[1] = 8, a__2[1] = xern1;
	    i__3[2] = 4, a__2[2] = ") = ";
	    i__3[3] = 16, a__2[3] = xern3;
	    i__3[4] = 43, a__2[4] = ".  IN THE CASE OF VECTOR ERROR TOLERANC"
		    "ES, ";
	    i__3[5] = 47, a__2[5] = "NO FURTHER CHECKING OF RTOL COMPONENTS "
		    "IS DONE.";
	    s_cat(ch__6, a__2, i__3, &c__6, (ftnlen)222);
	    xermsg_("SLATEC", "DES", ch__6, &c__7, &c__1, (ftnlen)6, (ftnlen)
		    3, (ftnlen)222);
	    *idid = -33;
	    nrtolp = 1;
	}

	if (natolp == 0 && atol[k] < 0.f) {
	    s_wsfi(&io___15);
	    do_fio(&c__1, (char *)&k, (ftnlen)sizeof(integer));
	    e_wsfi();
	    s_wsfi(&io___16);
	    do_fio(&c__1, (char *)&atol[k], (ftnlen)sizeof(real));
	    e_wsfi();
/* Writing concatenation */
	    i__3[0] = 104, a__2[0] = "IN DEABM, THE ABSOLUTE ERROR TOLERANCE"
		    "S ATOL MUST BE NON-NEGATIVE.  YOU HAVE CALLED THE CODE W"
		    "ITH  ATOL(";
	    i__3[1] = 8, a__2[1] = xern1;
	    i__3[2] = 4, a__2[2] = ") = ";
	    i__3[3] = 16, a__2[3] = xern3;
	    i__3[4] = 43, a__2[4] = ".  IN THE CASE OF VECTOR ERROR TOLERANC"
		    "ES, ";
	    i__3[5] = 47, a__2[5] = "NO FURTHER CHECKING OF ATOL COMPONENTS "
		    "IS DONE.";
	    s_cat(ch__6, a__2, i__3, &c__6, (ftnlen)222);
	    xermsg_("SLATEC", "DES", ch__6, &c__8, &c__1, (ftnlen)6, (ftnlen)
		    3, (ftnlen)222);
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
	r__3 = *tout - *t;
	r__4 = *tstop - *t;
	if (r_sign(&c_b71, &r__3) != r_sign(&c_b71, &r__4) || (r__1 = *tout - 
		*t, dabs(r__1)) > (r__2 = *tstop - *t, dabs(r__2))) {
	    s_wsfi(&io___17);
	    do_fio(&c__1, (char *)&(*tout), (ftnlen)sizeof(real));
	    e_wsfi();
	    s_wsfi(&io___19);
	    do_fio(&c__1, (char *)&(*tstop), (ftnlen)sizeof(real));
	    e_wsfi();
/* Writing concatenation */
	    i__4[0] = 48, a__3[0] = "IN DEABM, YOU HAVE CALLED THE CODE WITH"
		    "  TOUT = ";
	    i__4[1] = 16, a__3[1] = xern3;
	    i__4[2] = 14, a__3[2] = " BUT YOU HAVE ";
	    i__4[3] = 50, a__3[3] = "ALSO TOLD THE CODE (INFO(4) = 1) NOT TO"
		    " INTEGRATE ";
	    i__4[4] = 23, a__3[4] = "PAST THE POINT TSTOP = ";
	    i__4[5] = 16, a__3[5] = xern4;
	    i__4[6] = 7, a__3[6] = " THESE ";
	    i__4[7] = 22, a__3[7] = "INSTRUCTIONS CONFLICT.";
	    s_cat(ch__7, a__3, i__4, &c__8, (ftnlen)196);
	    xermsg_("SLATEC", "DES", ch__7, &c__14, &c__1, (ftnlen)6, (ftnlen)
		    3, (ftnlen)196);
	    *idid = -33;
	}
    }

/*     CHECK SOME CONTINUATION POSSIBILITIES */

    if (*init != 0) {
	if (*t == *tout) {
	    s_wsfi(&io___20);
	    do_fio(&c__1, (char *)&(*t), (ftnlen)sizeof(real));
	    e_wsfi();
/* Writing concatenation */
	    i__5[0] = 52, a__4[0] = "IN DEABM, YOU HAVE CALLED THE CODE WITH"
		    "  T = TOUT = ";
	    i__5[1] = 16, a__4[1] = xern3;
	    i__5[2] = 14, a__4[2] = "$$THIS IS NOT ";
	    i__5[3] = 30, a__4[3] = "ALLOWED ON CONTINUATION CALLS.";
	    s_cat(ch__5, a__4, i__5, &c__4, (ftnlen)112);
	    xermsg_("SLATEC", "DES", ch__5, &c__9, &c__1, (ftnlen)6, (ftnlen)
		    3, (ftnlen)112);
	    *idid = -33;
	}

	if (*t != *told) {
	    s_wsfi(&io___21);
	    do_fio(&c__1, (char *)&(*told), (ftnlen)sizeof(real));
	    e_wsfi();
	    s_wsfi(&io___22);
	    do_fio(&c__1, (char *)&(*t), (ftnlen)sizeof(real));
	    e_wsfi();
/* Writing concatenation */
	    i__6[0] = 47, a__5[0] = "IN DEABM, YOU HAVE CHANGED THE VALUE OF"
		    " T FROM ";
	    i__6[1] = 16, a__5[1] = xern3;
	    i__6[2] = 4, a__5[2] = " TO ";
	    i__6[3] = 16, a__5[3] = xern4;
	    i__6[4] = 44, a__5[4] = "  THIS IS NOT ALLOWED ON CONTINUATION C"
		    "ALLS.";
	    s_cat(ch__8, a__5, i__6, &c__5, (ftnlen)127);
	    xermsg_("SLATEC", "DES", ch__8, &c__10, &c__1, (ftnlen)6, (ftnlen)
		    3, (ftnlen)127);
	    *idid = -33;
	}

	if (*init != 1) {
	    if (*delsgn * (*tout - *t) < 0.f) {
		s_wsfi(&io___23);
		do_fio(&c__1, (char *)&(*tout), (ftnlen)sizeof(real));
		e_wsfi();
/* Writing concatenation */
		i__3[0] = 42, a__2[0] = "IN DEABM, BY CALLING THE CODE WITH "
			"TOUT = ";
		i__3[1] = 16, a__2[1] = xern3;
		i__3[2] = 9, a__2[2] = " YOU ARE ";
		i__3[3] = 38, a__2[3] = "ATTEMPTING TO CHANGE THE DIRECTION "
			"OF ";
		i__3[4] = 42, a__2[4] = "INTEGRATION.$$THIS IS NOT ALLOWED W"
			"ITHOUT ";
		i__3[5] = 11, a__2[5] = "RESTARTING.";
		s_cat(ch__9, a__2, i__3, &c__6, (ftnlen)158);
		xermsg_("SLATEC", "DES", ch__9, &c__11, &c__1, (ftnlen)6, (
			ftnlen)3, (ftnlen)158);
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
	    xermsg_("SLATEC", "DES", "IN DEABM, INVALID INPUT WAS DETECTED O"
		    "N SUCCESSIVE ENTRIES.  IT IS IMPOSSIBLE TO PROCEED BECAU"
		    "SE YOU HAVE NOT CORRECTED THE PROBLEM, SO EXECUTION IS B"
		    "EING TERMINATED.", &c__12, &c__2, (ftnlen)6, (ftnlen)3, (
		    ftnlen)166);
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
	if (rtol[k] + atol[k] > 0.f) {
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
    (*f)(&a, &y[1], &yp[1], &rpar[1], &ipar[1]);
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
    r__1 = *tout - *t;
    *delsgn = r_sign(&c_b71, &r__1);
/* Computing MAX */
    r__3 = *fouru * dabs(*x), r__4 = (r__1 = *tout - *x, dabs(r__1));
    r__2 = dmax(r__3,r__4);
    r__5 = *tout - *x;
    *h__ = r_sign(&r__2, &r__5);

/* ....................................................................... */

/*   ON EACH CALL SET INFORMATION WHICH DETERMINES THE ALLOWED INTERVAL */
/*   OF INTEGRATION BEFORE RETURNING WITH AN ANSWER AT TOUT */

L240:
    del = *tout - *t;
    absdel = dabs(del);

/* ....................................................................... */

/*   IF ALREADY PAST OUTPUT POINT, INTERPOLATE AND RETURN */

L250:
    if ((r__1 = *x - *t, dabs(r__1)) < absdel) {
	goto L260;
    }
    sintrp_(x, &yy[1], tout, &y[1], &ypout[1], neq, kold, &phi[phi_offset], 
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
    if ((r__1 = *tstop - *x, dabs(r__1)) >= *fouru * dabs(*x)) {
	goto L280;
    }
    dt = *tout - *x;
    i__2 = *neq;
    for (l = 1; l <= i__2; ++l) {
/* L270: */
	y[l] = yy[l] + dt * yp[l];
    }
    (*f)(tout, &y[1], &ypout[1], &rpar[1], &ipar[1]);
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
    ha = dabs(*h__);
    if (info[4] != 1) {
	goto L340;
    }
/* Computing MIN */
    r__2 = ha, r__3 = (r__1 = *tstop - *x, dabs(r__1));
    ha = dmin(r__2,r__3);
L340:
    *h__ = r_sign(&ha, h__);
    *eps = 1.f;
    ltol = 1;
    i__2 = *neq;
    for (l = 1; l <= i__2; ++l) {
	if (info[2] == 1) {
	    ltol = l;
	}
	wt[l] = rtol[ltol] * (r__1 = yy[l], dabs(r__1)) + atol[ltol];
	if (wt[l] <= 0.f) {
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
    steps_((S_fp)f, neq, &yy[1], x, h__, eps, &wt[1], start, hold, kord, kold,
	     &crash, &phi[phi_offset], &p[1], &yp[1], &psi[1], &alpha[1], &
	    beta[1], &sig[1], &v[1], &w[1], &g[1], phase1, ns, nornd, ksteps, 
	    twou, fouru, xold, kprev, ivc, &iv[1], kgi, &gi[1], &rpar[1], &
	    ipar[1]);

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
} /* des_ */

