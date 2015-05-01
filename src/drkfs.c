/* drkfs.f -- translated by f2c (version 12.02.01).
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
static integer c__6 = 6;
static integer c__7 = 7;
static integer c__8 = 8;
static integer c__9 = 9;
static integer c__10 = 10;
static integer c__11 = 11;
static integer c__12 = 12;
static doublereal c_b112 = 1.;
static doublereal c_b115 = .375;
static doublereal c_b139 = .2;

/* DECK DRKFS */
/* Subroutine */ int drkfs_(S_fp df, integer *neq, doublereal *t, doublereal *
	y, doublereal *tout, integer *info, doublereal *rtol, doublereal *
	atol, integer *idid, doublereal *h__, doublereal *tolfac, doublereal *
	yp, doublereal *f1, doublereal *f2, doublereal *f3, doublereal *f4, 
	doublereal *f5, doublereal *ys, doublereal *told, doublereal *dtsign, 
	doublereal *u26, doublereal *rer, integer *init, integer *ksteps, 
	integer *kop, integer *iquit, logical *stiff, logical *nonstf, 
	integer *ntstep, integer *nstifs, doublereal *rpar, integer *ipar)
{
    /* Initialized data */

    static doublereal remin = 1e-12;
    static integer mxstep = 500;
    static integer mxkop = 100;

    /* System generated locals */
    address a__1[2], a__2[6], a__3[4], a__4[5];
    integer i__1[2], i__2, i__3[6], i__4[4], i__5[5];
    doublereal d__1, d__2, d__3;
    char ch__1[222], ch__2[144], ch__3[166], ch__4[112], ch__5[223], ch__6[
	    113], ch__7[128], ch__8[159];

    /* Local variables */
    static doublereal a;
    static integer k;
    static doublereal s, u, ee, dt, es, et, dy, big, ute, tol, hmin, yavg;
    static integer ktol;
    static char xern1[8], xern3[16], xern4[16];
    extern /* Subroutine */ int dfehl_(S_fp, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, integer *);
    static doublereal eeoet;
    extern doublereal d1mach_(integer *);
    static logical hfaild;
    static doublereal estiff;
    static integer natolp;
    extern doublereal dhvnrm_(doublereal *, integer *);
    extern /* Subroutine */ int xermsg_(char *, char *, char *, integer *, 
	    integer *, ftnlen, ftnlen, ftnlen), dhstrt_(S_fp, integer *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, integer *, doublereal *, doublereal *, doublereal *,
	     doublereal *, doublereal *, doublereal *, doublereal *, integer *
	    , doublereal *);
    static doublereal esttol;
    static integer nrtolp;
    static logical output;

    /* Fortran I/O blocks */
    static icilist io___6 = { 0, xern1, 0, "(I8)", 8, 1 };
    static icilist io___7 = { 0, xern1, 0, "(I8)", 8, 1 };
    static icilist io___8 = { 0, xern1, 0, "(I8)", 8, 1 };
    static icilist io___9 = { 0, xern1, 0, "(I8)", 8, 1 };
    static icilist io___13 = { 0, xern1, 0, "(I8)", 8, 1 };
    static icilist io___15 = { 0, xern3, 0, "(1PE15.6)", 16, 1 };
    static icilist io___16 = { 0, xern1, 0, "(I8)", 8, 1 };
    static icilist io___17 = { 0, xern3, 0, "(1PE15.6)", 16, 1 };
    static icilist io___18 = { 0, xern3, 0, "(1PE15.6)", 16, 1 };
    static icilist io___19 = { 0, xern3, 0, "(1PE15.6)", 16, 1 };
    static icilist io___21 = { 0, xern4, 0, "(1PE15.6)", 16, 1 };
    static icilist io___22 = { 0, xern3, 0, "(1PE15.6)", 16, 1 };


/* ***BEGIN PROLOGUE  DRKFS */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to DDERKF */
/* ***LIBRARY   SLATEC */
/* ***TYPE      DOUBLE PRECISION (DERKFS-S, DRKFS-D) */
/* ***AUTHOR  Watts, H. A., (SNLA) */
/* ***DESCRIPTION */

/*     Fehlberg Fourth-Fifth Order Runge-Kutta Method */
/* ********************************************************************** */

/*     DRKFS integrates a system of first order ordinary differential */
/*     equations as described in the comments for DDERKF . */

/*     The arrays YP,F1,F2,F3,F4,F5,and YS  (of length at least NEQ) */
/*     appear in the call list for variable dimensioning purposes. */

/*     The variables H,TOLFAC,TOLD,DTSIGN,U26,RER,INIT,KSTEPS,KOP,IQUIT, */
/*     STIFF,NONSTF,NTSTEP, and NSTIFS are used internally by the code */
/*     and appear in the call list to eliminate local retention of */
/*     variables between calls. Accordingly, these variables and the */
/*     array YP should not be altered. */
/*     Items of possible interest are */
/*         H  - An appropriate step size to be used for the next step */
/*         TOLFAC - Factor of change in the tolerances */
/*         YP - Derivative of solution vector at T */
/*         KSTEPS - Counter on the number of steps attempted */

/* ********************************************************************** */

/* ***SEE ALSO  DDERKF */
/* ***ROUTINES CALLED  D1MACH, DFEHL, DHSTRT, DHVNRM, XERMSG */
/* ***REVISION HISTORY  (YYMMDD) */
/*   820301  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890831  Modified array declarations.  (WRB) */
/*   891024  Changed references from DVNORM to DHVNRM.  (WRB) */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900328  Added TYPE section.  (WRB) */
/*   900510  Convert XERRWV calls to XERMSG calls, change GOTOs to */
/*           IF-THEN-ELSEs.  (RWC) */
/*   910722  Updated AUTHOR section.  (ALS) */
/* ***END PROLOGUE  DRKFS */




/*     .................................................................. */

/*       A FIFTH ORDER METHOD WILL GENERALLY NOT BE CAPABLE OF DELIVERING */
/*       ACCURACIES NEAR LIMITING PRECISION ON COMPUTERS WITH LONG */
/*       WORDLENGTHS. TO PROTECT AGAINST LIMITING PRECISION DIFFICULTIES */
/*       ARISING FROM UNREASONABLE ACCURACY REQUESTS, AN APPROPRIATE */
/*       TOLERANCE THRESHOLD REMIN IS ASSIGNED FOR THIS METHOD. THIS */
/*       VALUE SHOULD NOT BE CHANGED ACROSS DIFFERENT MACHINES. */

    /* Parameter adjustments */
    --ipar;
    --rpar;
    --ys;
    --f5;
    --f4;
    --f3;
    --f2;
    --f1;
    --yp;
    --atol;
    --rtol;
    --info;
    --y;

    /* Function Body */

/*     .................................................................. */

/*       THE EXPENSE OF SOLVING THE PROBLEM IS MONITORED BY COUNTING THE */
/*       NUMBER OF  STEPS ATTEMPTED. WHEN THIS EXCEEDS  MXSTEP, THE */
/*       COUNTER IS RESET TO ZERO AND THE USER IS INFORMED ABOUT POSSIBLE */
/*       EXCESSIVE WORK. */


/*     .................................................................. */

/*       INEFFICIENCY CAUSED BY TOO FREQUENT OUTPUT IS MONITORED BY */
/*       COUNTING THE NUMBER OF STEP SIZES WHICH ARE SEVERELY SHORTENED */
/*       DUE SOLELY TO THE CHOICE OF OUTPUT POINTS. WHEN THE NUMBER OF */
/*       ABUSES EXCEED MXKOP, THE COUNTER IS RESET TO ZERO AND THE USER */
/*       IS INFORMED ABOUT POSSIBLE MISUSE OF THE CODE. */


/*     .................................................................. */

/* ***FIRST EXECUTABLE STATEMENT  DRKFS */
    if (info[1] == 0) {

/* ON THE FIRST CALL , PERFORM INITIALIZATION -- */
/*        DEFINE THE MACHINE UNIT ROUNDOFF QUANTITY  U  BY CALLING THE */
/*        FUNCTION ROUTINE  D1MACH. THE USER MUST MAKE SURE THAT THE */
/*        VALUES SET IN D1MACH ARE RELEVANT TO THE COMPUTER BEING USED. */

	u = d1mach_(&c__4);
/*                       -- SET ASSOCIATED MACHINE DEPENDENT PARAMETERS */
	*u26 = u * 26.;
	*rer = u * 2. + remin;
/*                       -- SET TERMINATION FLAG */
	*iquit = 0;
/*                       -- SET INITIALIZATION INDICATOR */
	*init = 0;
/*                       -- SET COUNTER FOR IMPACT OF OUTPUT POINTS */
	*kop = 0;
/*                       -- SET COUNTER FOR ATTEMPTED STEPS */
	*ksteps = 0;
/*                       -- SET INDICATORS FOR STIFFNESS DETECTION */
	*stiff = FALSE_;
	*nonstf = FALSE_;
/*                       -- SET STEP COUNTERS FOR STIFFNESS DETECTION */
	*ntstep = 0;
	*nstifs = 0;
/*                       -- RESET INFO(1) FOR SUBSEQUENT CALLS */
	info[1] = 1;
    }

/* ....................................................................... */

/*        CHECK VALIDITY OF INPUT PARAMETERS ON EACH ENTRY */

    if (info[1] != 0 && info[1] != 1) {
	s_wsfi(&io___6);
	do_fio(&c__1, (char *)&info[1], (ftnlen)sizeof(integer));
	e_wsfi();
/* Writing concatenation */
	i__1[0] = 214, a__1[0] = "IN DDERKF, INFO(1) MUST BE SET TO 0 FOR TH"
		"E START OF A NEW PROBLEM, AND MUST BE SET TO 1 FOLLOWING AN "
		"INTERRUPTED TASK.  YOU ARE ATTEMPTING TO CONTINUE THE INTEGR"
		"ATION ILLEGALLY BY CALLING THE CODE WITH  INFO(1) = ";
	i__1[1] = 8, a__1[1] = xern1;
	s_cat(ch__1, a__1, i__1, &c__2, (ftnlen)222);
	xermsg_("SLATEC", "DRKFS", ch__1, &c__3, &c__1, (ftnlen)6, (ftnlen)5, 
		(ftnlen)222);
	*idid = -33;
    }

    if (info[2] != 0 && info[2] != 1) {
	s_wsfi(&io___7);
	do_fio(&c__1, (char *)&info[2], (ftnlen)sizeof(integer));
	e_wsfi();
/* Writing concatenation */
	i__1[0] = 136, a__1[0] = "IN DDERKF, INFO(2) MUST BE 0 OR 1 INDICATI"
		"NG SCALAR AND VECTOR ERROR TOLERANCES, RESPECTIVELY.  YOU HA"
		"VE CALLED THE CODE WITH INFO(2) = ";
	i__1[1] = 8, a__1[1] = xern1;
	s_cat(ch__2, a__1, i__1, &c__2, (ftnlen)144);
	xermsg_("SLATEC", "DRKFS", ch__2, &c__4, &c__1, (ftnlen)6, (ftnlen)5, 
		(ftnlen)144);
	*idid = -33;
    }

    if (info[3] != 0 && info[3] != 1) {
	s_wsfi(&io___8);
	do_fio(&c__1, (char *)&info[3], (ftnlen)sizeof(integer));
	e_wsfi();
/* Writing concatenation */
	i__1[0] = 158, a__1[0] = "IN DDERKF, INFO(3) MUST BE 0 OR 1 INDICATI"
		"NG THE INTERVAL OR INTERMEDIATE-OUTPUT MODE OF INTEGRATION, "
		"RESPECTIVELY.  YOU HAVE CALLED THE CODE WITH  INFO(3) = ";
	i__1[1] = 8, a__1[1] = xern1;
	s_cat(ch__3, a__1, i__1, &c__2, (ftnlen)166);
	xermsg_("SLATEC", "DRKFS", ch__3, &c__5, &c__1, (ftnlen)6, (ftnlen)5, 
		(ftnlen)166);
	*idid = -33;
    }

    if (*neq < 1) {
	s_wsfi(&io___9);
	do_fio(&c__1, (char *)&(*neq), (ftnlen)sizeof(integer));
	e_wsfi();
/* Writing concatenation */
	i__1[0] = 104, a__1[0] = "IN DDERKF, THE NUMBER OF EQUATIONS NEQ MUS"
		"T BE A POSITIVE INTEGER.  YOU HAVE CALLED THE CODE WITH NEQ "
		"= ";
	i__1[1] = 8, a__1[1] = xern1;
	s_cat(ch__4, a__1, i__1, &c__2, (ftnlen)112);
	xermsg_("SLATEC", "DRKFS", ch__4, &c__6, &c__1, (ftnlen)6, (ftnlen)5, 
		(ftnlen)112);
	*idid = -33;
    }

    nrtolp = 0;
    natolp = 0;
    i__2 = *neq;
    for (k = 1; k <= i__2; ++k) {
	if (nrtolp == 0 && rtol[k] < 0.) {
	    s_wsfi(&io___13);
	    do_fio(&c__1, (char *)&k, (ftnlen)sizeof(integer));
	    e_wsfi();
	    s_wsfi(&io___15);
	    do_fio(&c__1, (char *)&rtol[k], (ftnlen)sizeof(doublereal));
	    e_wsfi();
/* Writing concatenation */
	    i__3[0] = 105, a__2[0] = "IN DDERKF, THE RELATIVE ERROR TOLERANC"
		    "ES RTOL MUST BE NON-NEGATIVE.  YOU HAVE CALLED THE CODE "
		    "WITH  RTOL(";
	    i__3[1] = 8, a__2[1] = xern1;
	    i__3[2] = 4, a__2[2] = ") = ";
	    i__3[3] = 16, a__2[3] = xern3;
	    i__3[4] = 43, a__2[4] = ".  IN THE CASE OF VECTOR ERROR TOLERANC"
		    "ES, ";
	    i__3[5] = 47, a__2[5] = "NO FURTHER CHECKING OF RTOL COMPONENTS "
		    "IS DONE.";
	    s_cat(ch__5, a__2, i__3, &c__6, (ftnlen)223);
	    xermsg_("SLATEC", "DRKFS", ch__5, &c__7, &c__1, (ftnlen)6, (
		    ftnlen)5, (ftnlen)223);
	    *idid = -33;
	    nrtolp = 1;
	}

	if (natolp == 0 && atol[k] < 0.) {
	    s_wsfi(&io___16);
	    do_fio(&c__1, (char *)&k, (ftnlen)sizeof(integer));
	    e_wsfi();
	    s_wsfi(&io___17);
	    do_fio(&c__1, (char *)&atol[k], (ftnlen)sizeof(doublereal));
	    e_wsfi();
/* Writing concatenation */
	    i__3[0] = 105, a__2[0] = "IN DDERKF, THE ABSOLUTE ERROR TOLERANC"
		    "ES ATOL MUST BE NON-NEGATIVE.  YOU HAVE CALLED THE CODE "
		    "WITH  ATOL(";
	    i__3[1] = 8, a__2[1] = xern1;
	    i__3[2] = 4, a__2[2] = ") = ";
	    i__3[3] = 16, a__2[3] = xern3;
	    i__3[4] = 43, a__2[4] = ".  IN THE CASE OF VECTOR ERROR TOLERANC"
		    "ES, ";
	    i__3[5] = 47, a__2[5] = "NO FURTHER CHECKING OF ATOL COMPONENTS "
		    "IS DONE.";
	    s_cat(ch__5, a__2, i__3, &c__6, (ftnlen)223);
	    xermsg_("SLATEC", "DRKFS", ch__5, &c__8, &c__1, (ftnlen)6, (
		    ftnlen)5, (ftnlen)223);
	    *idid = -33;
	    natolp = 1;
	}

	if (info[2] == 0) {
	    goto L20;
	}
	if (natolp > 0 && nrtolp > 0) {
	    goto L20;
	}
/* L10: */
    }


/*     CHECK SOME CONTINUATION POSSIBILITIES */

L20:
    if (*init != 0) {
	if (*t == *tout) {
	    s_wsfi(&io___18);
	    do_fio(&c__1, (char *)&(*t), (ftnlen)sizeof(doublereal));
	    e_wsfi();
/* Writing concatenation */
	    i__4[0] = 53, a__3[0] = "IN DDERKF, YOU HAVE CALLED THE CODE WIT"
		    "H  T = TOUT = ";
	    i__4[1] = 16, a__3[1] = xern3;
	    i__4[2] = 14, a__3[2] = "$$THIS IS NOT ";
	    i__4[3] = 30, a__3[3] = "ALLOWED ON CONTINUATION CALLS.";
	    s_cat(ch__6, a__3, i__4, &c__4, (ftnlen)113);
	    xermsg_("SLATEC", "DRKFS", ch__6, &c__9, &c__1, (ftnlen)6, (
		    ftnlen)5, (ftnlen)113);
	    *idid = -33;
	}

	if (*t != *told) {
	    s_wsfi(&io___19);
	    do_fio(&c__1, (char *)&(*told), (ftnlen)sizeof(doublereal));
	    e_wsfi();
	    s_wsfi(&io___21);
	    do_fio(&c__1, (char *)&(*t), (ftnlen)sizeof(doublereal));
	    e_wsfi();
/* Writing concatenation */
	    i__5[0] = 48, a__4[0] = "IN DDERKF, YOU HAVE CHANGED THE VALUE O"
		    "F T FROM ";
	    i__5[1] = 16, a__4[1] = xern3;
	    i__5[2] = 4, a__4[2] = " TO ";
	    i__5[3] = 16, a__4[3] = xern4;
	    i__5[4] = 44, a__4[4] = "$$THIS IS NOT ALLOWED ON CONTINUATION C"
		    "ALLS.";
	    s_cat(ch__7, a__4, i__5, &c__5, (ftnlen)128);
	    xermsg_("SLATEC", "DRKFS", ch__7, &c__10, &c__1, (ftnlen)6, (
		    ftnlen)5, (ftnlen)128);
	    *idid = -33;
	}

	if (*init != 1) {
	    if (*dtsign * (*tout - *t) < 0.) {
		s_wsfi(&io___22);
		do_fio(&c__1, (char *)&(*tout), (ftnlen)sizeof(doublereal));
		e_wsfi();
/* Writing concatenation */
		i__5[0] = 43, a__4[0] = "IN DDERKF, BY CALLING THE CODE WITH"
			" TOUT = ";
		i__5[1] = 16, a__4[1] = xern3;
		i__5[2] = 34, a__4[2] = " YOU ARE ATTEMPTING TO CHANGE THE ";
		i__5[3] = 47, a__4[3] = "DIRECTION OF INTEGRATION.$$THIS IS "
			"NOT ALLOWED ";
		i__5[4] = 19, a__4[4] = "WITHOUT RESTARTING.";
		s_cat(ch__8, a__4, i__5, &c__5, (ftnlen)159);
		xermsg_("SLATEC", "DRKFS", ch__8, &c__11, &c__1, (ftnlen)6, (
			ftnlen)5, (ftnlen)159);
		*idid = -33;
	    }
	}
    }

/*     INVALID INPUT DETECTED */

    if (*idid == -33) {
	if (*iquit != -33) {
	    *iquit = -33;
	    goto L540;
	} else {
	    xermsg_("SLATEC", "DRKFS", "IN DDERKF, INVALID INPUT WAS DETECTE"
		    "D ON SUCCESSIVE ENTRIES.  IT IS IMPOSSIBLE TO PROCEED BE"
		    "CAUSE YOU HAVE NOT CORRECTED THE PROBLEM, SO EXECUTION I"
		    "S BEING TERMINATED.", &c__12, &c__2, (ftnlen)6, (ftnlen)5,
		     (ftnlen)167);
	    return 0;
	}
    }

/*           ............................................................ */

/*                RTOL = ATOL = 0. IS ALLOWED AS VALID INPUT AND */
/*                INTERPRETED AS ASKING FOR THE MOST ACCURATE SOLUTION */
/*                POSSIBLE. IN THIS CASE, THE RELATIVE ERROR TOLERANCE */
/*                RTOL IS RESET TO THE SMALLEST VALUE RER WHICH IS LIKELY */
/*                TO BE REASONABLE FOR THIS METHOD AND MACHINE. */

    i__2 = *neq;
    for (k = 1; k <= i__2; ++k) {
	if (rtol[k] + atol[k] > 0.) {
	    goto L180;
	}
	rtol[k] = *rer;
	*idid = -2;
L180:
/*           ...EXIT */
	if (info[2] == 0) {
	    goto L200;
	}
/* L190: */
    }
L200:

    if (*idid != -2) {
	goto L210;
    }

/*              RTOL=ATOL=0 ON INPUT, SO RTOL WAS CHANGED TO A */
/*                                       SMALL POSITIVE VALUE */
    *tolfac = 1.;
    goto L530;
L210:

/*                       BRANCH ON STATUS OF INITIALIZATION INDICATOR */
/*                              INIT=0 MEANS INITIAL DERIVATIVES AND */
/*                              STARTING STEP SIZE */
/*                                     NOT YET COMPUTED */
/*                              INIT=1 MEANS STARTING STEP SIZE NOT YET */
/*                              COMPUTED INIT=2 MEANS NO FURTHER */
/*                              INITIALIZATION REQUIRED */

    if (*init == 0) {
	goto L220;
    }
/*                    ......EXIT */
    if (*init == 1) {
	goto L240;
    }
/*                 .........EXIT */
    goto L260;
L220:

/*                       ................................................ */

/*                            MORE INITIALIZATION -- */
/*                                                -- EVALUATE INITIAL */
/*                                                DERIVATIVES */

    *init = 1;
    a = *t;
    (*df)(&a, &y[1], &yp[1], &rpar[1], &ipar[1]);
    if (*t != *tout) {
	goto L230;
    }

/*                          INTERVAL MODE */
    *idid = 2;
    *t = *tout;
    *told = *t;
/*     .....................EXIT */
    goto L560;
L230:
L240:

/*                    -- SET SIGN OF INTEGRATION DIRECTION  AND */
/*                    -- ESTIMATE STARTING STEP SIZE */

    *init = 2;
    d__1 = *tout - *t;
    *dtsign = d_sign(&c_b112, &d__1);
    u = d1mach_(&c__4);
    big = sqrt(d1mach_(&c__2));
    ute = pow_dd(&u, &c_b115);
    dy = ute * dhvnrm_(&y[1], neq);
    if (dy == 0.) {
	dy = ute;
    }
    ktol = 1;
    i__2 = *neq;
    for (k = 1; k <= i__2; ++k) {
	if (info[2] == 1) {
	    ktol = k;
	}
	tol = rtol[ktol] * (d__1 = y[k], abs(d__1)) + atol[ktol];
	if (tol == 0.) {
	    tol = dy * rtol[ktol];
	}
	f1[k] = tol;
/* L250: */
    }

    dhstrt_((S_fp)df, neq, t, tout, &y[1], &yp[1], &f1[1], &c__4, &u, &big, &
	    f2[1], &f3[1], &f4[1], &f5[1], &rpar[1], &ipar[1], h__);
L260:

/*                 ...................................................... */

/*                      SET STEP SIZE FOR INTEGRATION IN THE DIRECTION */
/*                      FROM T TO TOUT AND SET OUTPUT POINT INDICATOR */

    dt = *tout - *t;
    *h__ = d_sign(h__, &dt);
    output = FALSE_;

/*                 TEST TO SEE IF DDERKF IS BEING SEVERELY IMPACTED BY */
/*                 TOO MANY OUTPUT POINTS */

    if (abs(*h__) >= abs(dt) * 2.) {
	++(*kop);
    }
    if (*kop <= mxkop) {
	goto L270;
    }

/*                    UNNECESSARY FREQUENCY OF OUTPUT IS RESTRICTING */
/*                                              THE STEP SIZE CHOICE */
    *idid = -5;
    *kop = 0;
    goto L510;
L270:

    if (abs(dt) > *u26 * abs(*t)) {
	goto L290;
    }

/*                       IF TOO CLOSE TO OUTPUT POINT,EXTRAPOLATE AND */
/*                       RETURN */

    i__2 = *neq;
    for (k = 1; k <= i__2; ++k) {
	y[k] += dt * yp[k];
/* L280: */
    }
    a = *tout;
    (*df)(&a, &y[1], &yp[1], &rpar[1], &ipar[1]);
    ++(*ksteps);
    goto L500;
L290:
/*                       BEGIN BLOCK PERMITTING ...EXITS TO 490 */

/*                          ********************************************* */
/*                          ********************************************* */
/*                               STEP BY STEP INTEGRATION */

L300:
/*                             BEGIN BLOCK PERMITTING ...EXITS TO 480 */
    hfaild = FALSE_;

/*                                TO PROTECT AGAINST IMPOSSIBLE ACCURACY */
/*                                REQUESTS, COMPUTE A TOLERANCE FACTOR */
/*                                BASED ON THE REQUESTED ERROR TOLERANCE */
/*                                AND A LEVEL OF ACCURACY ACHIEVABLE AT */
/*                                LIMITING PRECISION */

    *tolfac = 0.;
    ktol = 1;
    i__2 = *neq;
    for (k = 1; k <= i__2; ++k) {
	if (info[2] == 1) {
	    ktol = k;
	}
	et = rtol[ktol] * (d__1 = y[k], abs(d__1)) + atol[ktol];
	if (et > 0.) {
	    goto L310;
	}
/* Computing MAX */
	d__1 = *tolfac, d__2 = *rer / rtol[ktol];
	*tolfac = max(d__1,d__2);
	goto L320;
L310:
/* Computing MAX */
	d__2 = *tolfac, d__3 = (d__1 = y[k], abs(d__1)) * (*rer / et);
	*tolfac = max(d__2,d__3);
L320:
/* L330: */
	;
    }
    if (*tolfac <= 1.) {
	goto L340;
    }

/*                          REQUESTED ERROR UNATTAINABLE DUE TO LIMITED */
/*                                                  PRECISION AVAILABLE */
    *tolfac *= 2.;
    *idid = -2;
/*              .....................EXIT */
    goto L520;
L340:

/*                                SET SMALLEST ALLOWABLE STEP SIZE */

    hmin = *u26 * abs(*t);

/*                                ADJUST STEP SIZE IF NECESSARY TO HIT */
/*                                THE OUTPUT POINT -- LOOK AHEAD TWO */
/*                                STEPS TO AVOID DRASTIC CHANGES IN THE */
/*                                STEP SIZE AND THUS LESSEN THE IMPACT OF */
/*                                OUTPUT POINTS ON THE CODE.  STRETCH THE */
/*                                STEP SIZE BY, AT MOST, AN AMOUNT EQUAL */
/*                                TO THE SAFETY FACTOR OF 9/10. */

    dt = *tout - *t;
    if (abs(dt) >= abs(*h__) * 2.) {
	goto L370;
    }
    if (abs(dt) > abs(*h__) / .9) {
	goto L350;
    }

/*                                      THE NEXT STEP, IF SUCCESSFUL, */
/*                                      WILL COMPLETE THE INTEGRATION TO */
/*                                      THE OUTPUT POINT */

    output = TRUE_;
    *h__ = dt;
    goto L360;
L350:

    *h__ = dt * .5;
L360:
L370:


/*                                *************************************** */
/*                                     CORE INTEGRATOR FOR TAKING A */
/*                                     SINGLE STEP */
/*                                *************************************** */
/*                                     TO AVOID PROBLEMS WITH ZERO */
/*                                     CROSSINGS, RELATIVE ERROR IS */
/*                                     MEASURED USING THE AVERAGE OF THE */
/*                                     MAGNITUDES OF THE SOLUTION AT THE */
/*                                     BEGINNING AND END OF A STEP. */
/*                                     THE ERROR ESTIMATE FORMULA HAS */
/*                                     BEEN GROUPED TO CONTROL LOSS OF */
/*                                     SIGNIFICANCE. */
/*                                     LOCAL ERROR ESTIMATES FOR A FIRST */
/*                                     ORDER METHOD USING THE SAME */
/*                                     STEP SIZE AS THE FEHLBERG METHOD */
/*                                     ARE CALCULATED AS PART OF THE */
/*                                     TEST FOR STIFFNESS. */
/*                                     TO DISTINGUISH THE VARIOUS */
/*                                     ARGUMENTS, H IS NOT PERMITTED */
/*                                     TO BECOME SMALLER THAN 26 UNITS OF */
/*                                     ROUNDOFF IN T.  PRACTICAL LIMITS */
/*                                     ON THE CHANGE IN THE STEP SIZE ARE */
/*                                     ENFORCED TO SMOOTH THE STEP SIZE */
/*                                     SELECTION PROCESS AND TO AVOID */
/*                                     EXCESSIVE CHATTERING ON PROBLEMS */
/*                                     HAVING DISCONTINUITIES.  TO */
/*                                     PREVENT UNNECESSARY FAILURES, THE */
/*                                     CODE USES 9/10 THE STEP SIZE */
/*                                     IT ESTIMATES WILL SUCCEED. */
/*                                     AFTER A STEP FAILURE, THE STEP */
/*                                     SIZE IS NOT ALLOWED TO INCREASE */
/*                                     FOR THE NEXT ATTEMPTED STEP. THIS */
/*                                     MAKES THE CODE MORE EFFICIENT ON */
/*                                     PROBLEMS HAVING DISCONTINUITIES */
/*                                     AND MORE EFFECTIVE IN GENERAL */
/*                                     SINCE LOCAL EXTRAPOLATION IS BEING */
/*                                     USED AND EXTRA CAUTION SEEMS */
/*                                     WARRANTED. */
/*                                ....................................... */

/*                                     MONITOR NUMBER OF STEPS ATTEMPTED */

L380:
    if (*ksteps <= mxstep) {
	goto L390;
    }

/*                                      A SIGNIFICANT AMOUNT OF WORK HAS */
/*                                      BEEN EXPENDED */
    *idid = -1;
    *ksteps = 0;
/*              ........................EXIT */
    if (! (*stiff)) {
	goto L520;
    }

/*                                      PROBLEM APPEARS TO BE STIFF */
    *idid = -4;
    *stiff = FALSE_;
    *nonstf = FALSE_;
    *ntstep = 0;
    *nstifs = 0;
/*              ........................EXIT */
    goto L520;
L390:

/*                                   ADVANCE AN APPROXIMATE SOLUTION OVER */
/*                                   ONE STEP OF LENGTH H */

    dfehl_((S_fp)df, neq, t, &y[1], h__, &yp[1], &f1[1], &f2[1], &f3[1], &f4[
	    1], &f5[1], &ys[1], &rpar[1], &ipar[1]);
    ++(*ksteps);

/*                                   .................................... */

/*                                        COMPUTE AND TEST ALLOWABLE */
/*                                        TOLERANCES VERSUS LOCAL ERROR */
/*                                        ESTIMATES.  NOTE THAT RELATIVE */
/*                                        ERROR IS MEASURED WITH RESPECT */
/*                                        TO THE AVERAGE OF THE */
/*                                        MAGNITUDES OF THE SOLUTION AT */
/*                                        THE BEGINNING AND END OF THE */
/*                                        STEP.  LOCAL ERROR ESTIMATES */
/*                                        FOR A SPECIAL FIRST ORDER */
/*                                        METHOD ARE CALCULATED ONLY WHEN */
/*                                        THE STIFFNESS DETECTION IS */
/*                                        TURNED ON. */

    eeoet = 0.;
    estiff = 0.;
    ktol = 1;
    i__2 = *neq;
    for (k = 1; k <= i__2; ++k) {
	yavg = ((d__1 = y[k], abs(d__1)) + (d__2 = ys[k], abs(d__2))) * .5;
	if (info[2] == 1) {
	    ktol = k;
	}
	et = rtol[ktol] * yavg + atol[ktol];
	if (et > 0.) {
	    goto L400;
	}

/*           PURE RELATIVE ERROR INAPPROPRIATE WHEN SOLUTION */
/*                                                  VANISHES */
	*idid = -3;
/*              ...........................EXIT */
	goto L520;
L400:

	ee = (d__1 = yp[k] * -2090. + (f3[k] * 21970. - f4[k] * 15048.) + (f2[
		k] * 22528. - f5[k] * 27360.), abs(d__1));
	if (*stiff || *nonstf) {
	    goto L410;
	}
	es = (d__1 = *h__ * (yp[k] * .055455 - f1[k] * .035493 - f2[k] * 
		.036571 + f3[k] * .023107 - f4[k] * .009515 + f5[k] * .003017)
		, abs(d__1));
/* Computing MAX */
	d__1 = estiff, d__2 = es / et;
	estiff = max(d__1,d__2);
L410:
/* Computing MAX */
	d__1 = eeoet, d__2 = ee / et;
	eeoet = max(d__1,d__2);
/* L420: */
    }

    esttol = abs(*h__) * eeoet / 752400.;

/*                                ...EXIT */
    if (esttol <= 1.) {
	goto L440;
    }

/*                                   .................................... */

/*                                        UNSUCCESSFUL STEP */

    if (abs(*h__) > hmin) {
	goto L430;
    }

/*                             REQUESTED ERROR UNATTAINABLE AT SMALLEST */
/*                                                  ALLOWABLE STEP SIZE */
    *tolfac = esttol * 1.69;
    *idid = -2;
/*              ........................EXIT */
    goto L520;
L430:

/*                                   REDUCE THE STEP SIZE , TRY AGAIN */
/*                                   THE DECREASE IS LIMITED TO A FACTOR */
/*                                   OF 1/10 */

    hfaild = TRUE_;
    output = FALSE_;
    s = .1;
    if (esttol < 59049.) {
	s = .9 / pow_dd(&esttol, &c_b139);
    }
/* Computing MAX */
    d__2 = s * abs(*h__);
    d__1 = max(d__2,hmin);
    *h__ = d_sign(&d__1, h__);
    goto L380;
L440:

/*                                ....................................... */

/*                                SUCCESSFUL STEP */
/*                                                  STORE SOLUTION AT T+H */
/*                                                  AND EVALUATE */
/*                                                  DERIVATIVES THERE */

    *t += *h__;
    i__2 = *neq;
    for (k = 1; k <= i__2; ++k) {
	y[k] = ys[k];
/* L450: */
    }
    a = *t;
    (*df)(&a, &y[1], &yp[1], &rpar[1], &ipar[1]);

/*                                CHOOSE NEXT STEP SIZE */
/*                                THE INCREASE IS LIMITED TO A FACTOR OF */
/*                                5 IF STEP FAILURE HAS JUST OCCURRED, */
/*                                NEXT */
/*                                   STEP SIZE IS NOT ALLOWED TO INCREASE */

    s = 5.;
    if (esttol > 1.889568e-4) {
	s = .9 / pow_dd(&esttol, &c_b139);
    }
    if (hfaild) {
	s = min(s,1.);
    }
/* Computing MAX */
    d__2 = s * abs(*h__);
    d__1 = max(d__2,hmin);
    *h__ = d_sign(&d__1, h__);

/*                                ....................................... */

/*                                     CHECK FOR STIFFNESS (IF NOT */
/*                                     ALREADY DETECTED) */

/*                                     IN A SEQUENCE OF 50 SUCCESSFUL */
/*                                     STEPS BY THE FEHLBERG METHOD, 25 */
/*                                     SUCCESSFUL STEPS BY THE FIRST */
/*                                     ORDER METHOD INDICATES STIFFNESS */
/*                                     AND TURNS THE TEST OFF. IF 26 */
/*                                     FAILURES BY THE FIRST ORDER METHOD */
/*                                     OCCUR, THE TEST IS TURNED OFF */
/*                                     UNTIL THIS SEQUENCE OF 50 STEPS BY */
/*                                     THE FEHLBERG METHOD IS COMPLETED. */

/*                             ...EXIT */
    if (*stiff) {
	goto L480;
    }
    *ntstep = (*ntstep + 1) % 50;
    if (*ntstep == 1) {
	*nonstf = FALSE_;
    }
/*                             ...EXIT */
    if (*nonstf) {
	goto L480;
    }
    if (estiff > 1.) {
	goto L460;
    }

/*                                   SUCCESSFUL STEP WITH FIRST ORDER */
/*                                   METHOD */
    ++(*nstifs);
/*                                   TURN TEST OFF AFTER 25 INDICATIONS */
/*                                   OF STIFFNESS */
    if (*nstifs == 25) {
	*stiff = TRUE_;
    }
    goto L470;
L460:

/*                                UNSUCCESSFUL STEP WITH FIRST ORDER */
/*                                METHOD */
    if (*ntstep - *nstifs <= 25) {
	goto L470;
    }
/*               TURN STIFFNESS DETECTION OFF FOR THIS BLOCK OF */
/*                                                  FIFTY STEPS */
    *nonstf = TRUE_;
/*                                   RESET STIFF STEP COUNTER */
    *nstifs = 0;
L470:
L480:

/*                             ****************************************** */
/*                                  END OF CORE INTEGRATOR */
/*                             ****************************************** */


/*                                  SHOULD WE TAKE ANOTHER STEP */

/*                       ......EXIT */
    if (output) {
	goto L490;
    }
    if (info[3] == 0) {
	goto L300;
    }

/*                          ********************************************* */
/*                          ********************************************* */

/*                               INTEGRATION SUCCESSFULLY COMPLETED */

/*                                           ONE-STEP MODE */
    *idid = 1;
    *told = *t;
/*     .....................EXIT */
    goto L560;
L490:
L500:

/*                    INTERVAL MODE */
    *idid = 2;
    *t = *tout;
    *told = *t;
/*     ...............EXIT */
    goto L560;
L510:
L520:
L530:
L540:

/*        INTEGRATION TASK INTERRUPTED */

    info[1] = -1;
    *told = *t;
/*     ...EXIT */
    if (*idid != -2) {
	goto L560;
    }

/*        THE ERROR TOLERANCES ARE INCREASED TO VALUES */
/*                WHICH ARE APPROPRIATE FOR CONTINUING */
    rtol[1] = *tolfac * rtol[1];
    atol[1] = *tolfac * atol[1];
/*     ...EXIT */
    if (info[2] == 0) {
	goto L560;
    }
    i__2 = *neq;
    for (k = 2; k <= i__2; ++k) {
	rtol[k] = *tolfac * rtol[k];
	atol[k] = *tolfac * atol[k];
/* L550: */
    }
L560:
    return 0;
} /* drkfs_ */

