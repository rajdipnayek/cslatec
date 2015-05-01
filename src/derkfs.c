/* derkfs.f -- translated by f2c (version 12.02.01).
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
static real c_b110 = 1.f;
static doublereal c_b113 = .375;
static doublereal c_b129 = .2;

/* DECK DERKFS */
/* Subroutine */ int derkfs_(S_fp f, integer *neq, real *t, real *y, real *
	tout, integer *info, real *rtol, real *atol, integer *idid, real *h__,
	 real *tolfac, real *yp, real *f1, real *f2, real *f3, real *f4, real 
	*f5, real *ys, real *told, real *dtsign, real *u26, real *rer, 
	integer *init, integer *ksteps, integer *kop, integer *iquit, logical 
	*stiff, logical *nonstf, integer *ntstep, integer *nstifs, real *rpar,
	 integer *ipar)
{
    /* Initialized data */

    static real remin = 1e-12f;
    static integer mxstep = 500;
    static integer mxkop = 100;

    /* System generated locals */
    address a__1[2], a__2[6], a__3[4], a__4[5];
    integer i__1[2], i__2, i__3[6], i__4[4], i__5[5];
    real r__1, r__2, r__3;
    doublereal d__1;
    char ch__1[221], ch__2[143], ch__3[156], ch__4[111], ch__5[222], ch__6[
	    112], ch__7[127], ch__8[158];

    /* Local variables */
    static real a;
    static integer k;
    static real s, u, ee, dt, et, es, dy, big, ute, tol, hmin, yavg;
    static integer ktol;
    static char xern1[8], xern3[16], xern4[16];
    static real eeoet;
    extern doublereal hvnrm_(real *, integer *), r1mach_(integer *);
    static logical hfaild;
    extern /* Subroutine */ int defehl_(S_fp, integer *, real *, real *, real 
	    *, real *, real *, real *, real *, real *, real *, real *, real *,
	     integer *);
    static real estiff;
    static integer natolp;
    extern /* Subroutine */ int xermsg_(char *, char *, char *, integer *, 
	    integer *, ftnlen, ftnlen, ftnlen), hstart_(S_fp, integer *, real 
	    *, real *, real *, real *, real *, integer *, real *, real *, 
	    real *, real *, real *, real *, real *, integer *, real *);
    static real esttol;
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


/* ***BEGIN PROLOGUE  DERKFS */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to DERKF */
/* ***LIBRARY   SLATEC */
/* ***TYPE      SINGLE PRECISION (DERKFS-S, DRKFS-D) */
/* ***AUTHOR  Watts, H. A., (SNLA) */
/* ***DESCRIPTION */

/*     Fehlberg Fourth-Fifth order Runge-Kutta Method */
/* ********************************************************************** */

/*     DERKFS integrates a system of first order ordinary differential */
/*     equations as described in the comments for DERKF . */

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

/* ***SEE ALSO  DERKF */
/* ***ROUTINES CALLED  DEFEHL, HSTART, HVNRM, R1MACH, XERMSG */
/* ***REVISION HISTORY  (YYMMDD) */
/*   800501  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890831  Modified array declarations.  (WRB) */
/*   891024  Changed references from VNORM to HVNRM.  (WRB) */
/*   891024  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900328  Added TYPE section.  (WRB) */
/*   900510  Convert XERRWV calls to XERMSG calls, replace GOTOs with */
/*           IF-THEN-ELSEs.  (RWC) */
/*   910722  Updated AUTHOR section.  (ALS) */
/* ***END PROLOGUE  DERKFS */




/* ....................................................................... */

/*  A FIFTH ORDER METHOD WILL GENERALLY NOT BE CAPABLE OF DELIVERING */
/*  ACCURACIES NEAR LIMITING PRECISION ON COMPUTERS WITH LONG */
/*  WORDLENGTHS. TO PROTECT AGAINST LIMITING PRECISION DIFFICULTIES */
/*  ARISING FROM UNREASONABLE ACCURACY REQUESTS, AN APPROPRIATE */
/*  TOLERANCE THRESHOLD REMIN IS ASSIGNED FOR THIS METHOD. THIS VALUE */
/*  SHOULD NOT BE CHANGED ACROSS DIFFERENT MACHINES. */

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

/* ....................................................................... */

/*  THE EXPENSE OF SOLVING THE PROBLEM IS MONITORED BY COUNTING THE */
/*  NUMBER OF  STEPS ATTEMPTED. WHEN THIS EXCEEDS  MXSTEP, THE COUNTER */
/*  IS RESET TO ZERO AND THE USER IS INFORMED ABOUT POSSIBLE EXCESSIVE */
/*  WORK. */


/* ....................................................................... */

/*  INEFFICIENCY CAUSED BY TOO FREQUENT OUTPUT IS MONITORED BY COUNTING */
/*  THE NUMBER OF STEP SIZES WHICH ARE SEVERELY SHORTENED DUE SOLELY TO */
/*  THE CHOICE OF OUTPUT POINTS. WHEN THE NUMBER OF ABUSES EXCEED MXKOP, */
/*  THE COUNTER IS RESET TO ZERO AND THE USER IS INFORMED ABOUT POSSIBLE */
/*  MISUSE OF THE CODE. */


/* ....................................................................... */

/* ***FIRST EXECUTABLE STATEMENT  DERKFS */
    if (info[1] == 0) {

/* ON THE FIRST CALL , PERFORM INITIALIZATION -- */
/*        DEFINE THE MACHINE UNIT ROUNDOFF QUANTITY  U  BY CALLING THE */
/*        FUNCTION ROUTINE  R1MACH. THE USER MUST MAKE SURE THAT THE */
/*        VALUES SET IN R1MACH ARE RELEVANT TO THE COMPUTER BEING USED. */

	u = r1mach_(&c__4);
/*                       -- SET ASSOCIATED MACHINE DEPENDENT PARAMETERS */
	*u26 = u * 26.f;
	*rer = u * 2.f + remin;
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
	i__1[0] = 213, a__1[0] = "IN DERKF, INFO(1) MUST BE SET TO 0 FOR THE"
		" START OF A NEW PROBLEM, AND MUST BE SET TO 1 FOLLOWING AN I"
		"NTERRUPTED TASK.  YOU ARE ATTEMPTING TO CONTINUE THE INTEGRA"
		"TION ILLEGALLY BY CALLING THE CODE WITH  INFO(1) = ";
	i__1[1] = 8, a__1[1] = xern1;
	s_cat(ch__1, a__1, i__1, &c__2, (ftnlen)221);
	xermsg_("SLATEC", "DERKFS", ch__1, &c__3, &c__1, (ftnlen)6, (ftnlen)6,
		 (ftnlen)221);
	*idid = -33;
    }

    if (info[2] != 0 && info[2] != 1) {
	s_wsfi(&io___7);
	do_fio(&c__1, (char *)&info[2], (ftnlen)sizeof(integer));
	e_wsfi();
/* Writing concatenation */
	i__1[0] = 135, a__1[0] = "IN DERKF, INFO(2) MUST BE 0 OR 1 INDICATIN"
		"G SCALAR AND VECTOR ERROR TOLERANCES, RESPECTIVELY.  YOU HAV"
		"E CALLED THE CODE WITH INFO(2) = ";
	i__1[1] = 8, a__1[1] = xern1;
	s_cat(ch__2, a__1, i__1, &c__2, (ftnlen)143);
	xermsg_("SLATEC", "DERKFS", ch__2, &c__4, &c__1, (ftnlen)6, (ftnlen)6,
		 (ftnlen)143);
	*idid = -33;
    }

    if (info[3] != 0 && info[3] != 1) {
	s_wsfi(&io___8);
	do_fio(&c__1, (char *)&info[3], (ftnlen)sizeof(integer));
	e_wsfi();
/* Writing concatenation */
	i__1[0] = 148, a__1[0] = "IN DERKF, INFO(3) MUST BE 0 OR 1 INDICATIN"
		"G THE OR INTERMEDIATE-OUTPUT MODE OF INTEGRATION, RESPECTIVE"
		"LY.  YOU HAVE CALLED THE CODE WITH  INFO(3) = ";
	i__1[1] = 8, a__1[1] = xern1;
	s_cat(ch__3, a__1, i__1, &c__2, (ftnlen)156);
	xermsg_("SLATEC", "DERKFS", ch__3, &c__5, &c__1, (ftnlen)6, (ftnlen)6,
		 (ftnlen)156);
	*idid = -33;
    }

    if (*neq < 1) {
	s_wsfi(&io___9);
	do_fio(&c__1, (char *)&(*neq), (ftnlen)sizeof(integer));
	e_wsfi();
/* Writing concatenation */
	i__1[0] = 103, a__1[0] = "IN DERKF, THE NUMBER OF EQUATIONS NEQ MUST"
		" BE A POSITIVE INTEGER.  YOU HAVE CALLED THE CODE WITH NEQ = "
		;
	i__1[1] = 8, a__1[1] = xern1;
	s_cat(ch__4, a__1, i__1, &c__2, (ftnlen)111);
	xermsg_("SLATEC", "DERKFS", ch__4, &c__6, &c__1, (ftnlen)6, (ftnlen)6,
		 (ftnlen)111);
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
	    do_fio(&c__1, (char *)&rtol[k], (ftnlen)sizeof(real));
	    e_wsfi();
/* Writing concatenation */
	    i__3[0] = 104, a__2[0] = "IN DERKF, THE RELATIVE ERROR TOLERANCE"
		    "S RTOL MUST BE NON-NEGATIVE.  YOU HAVE CALLED THE CODE W"
		    "ITH  RTOL(";
	    i__3[1] = 8, a__2[1] = xern1;
	    i__3[2] = 4, a__2[2] = ") = ";
	    i__3[3] = 16, a__2[3] = xern3;
	    i__3[4] = 43, a__2[4] = ".  IN THE CASE OF VECTOR ERROR TOLERANC"
		    "ES, ";
	    i__3[5] = 47, a__2[5] = "NO FURTHER CHECKING OF RTOL COMPONENTS "
		    "IS DONE.";
	    s_cat(ch__5, a__2, i__3, &c__6, (ftnlen)222);
	    xermsg_("SLATEC", "DERKFS", ch__5, &c__7, &c__1, (ftnlen)6, (
		    ftnlen)6, (ftnlen)222);
	    *idid = -33;
	    nrtolp = 1;
	}

	if (natolp == 0 && atol[k] < 0.) {
	    s_wsfi(&io___16);
	    do_fio(&c__1, (char *)&k, (ftnlen)sizeof(integer));
	    e_wsfi();
	    s_wsfi(&io___17);
	    do_fio(&c__1, (char *)&atol[k], (ftnlen)sizeof(real));
	    e_wsfi();
/* Writing concatenation */
	    i__3[0] = 104, a__2[0] = "IN DERKF, THE ABSOLUTE ERROR TOLERANCE"
		    "S ATOL MUST BE NON-NEGATIVE.  YOU HAVE CALLED THE CODE W"
		    "ITH  ATOL(";
	    i__3[1] = 8, a__2[1] = xern1;
	    i__3[2] = 4, a__2[2] = ") = ";
	    i__3[3] = 16, a__2[3] = xern3;
	    i__3[4] = 43, a__2[4] = ".  IN THE CASE OF VECTOR ERROR TOLERANC"
		    "ES, ";
	    i__3[5] = 47, a__2[5] = "NO FURTHER CHECKING OF ATOL COMPONENTS "
		    "IS DONE.";
	    s_cat(ch__5, a__2, i__3, &c__6, (ftnlen)222);
	    xermsg_("SLATEC", "DERKFS", ch__5, &c__8, &c__1, (ftnlen)6, (
		    ftnlen)6, (ftnlen)222);
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
	    do_fio(&c__1, (char *)&(*t), (ftnlen)sizeof(real));
	    e_wsfi();
/* Writing concatenation */
	    i__4[0] = 52, a__3[0] = "IN DERKF, YOU HAVE CALLED THE CODE WITH"
		    "  T = TOUT = ";
	    i__4[1] = 16, a__3[1] = xern3;
	    i__4[2] = 14, a__3[2] = "$$THIS IS NOT ";
	    i__4[3] = 30, a__3[3] = "ALLOWED ON CONTINUATION CALLS.";
	    s_cat(ch__6, a__3, i__4, &c__4, (ftnlen)112);
	    xermsg_("SLATEC", "DERKFS", ch__6, &c__9, &c__1, (ftnlen)6, (
		    ftnlen)6, (ftnlen)112);
	    *idid = -33;
	}

	if (*t != *told) {
	    s_wsfi(&io___19);
	    do_fio(&c__1, (char *)&(*told), (ftnlen)sizeof(real));
	    e_wsfi();
	    s_wsfi(&io___21);
	    do_fio(&c__1, (char *)&(*t), (ftnlen)sizeof(real));
	    e_wsfi();
/* Writing concatenation */
	    i__5[0] = 47, a__4[0] = "IN DERKF, YOU HAVE CHANGED THE VALUE OF"
		    " T FROM ";
	    i__5[1] = 16, a__4[1] = xern3;
	    i__5[2] = 4, a__4[2] = " TO ";
	    i__5[3] = 16, a__4[3] = xern4;
	    i__5[4] = 44, a__4[4] = "$$THIS IS NOT ALLOWED ON CONTINUATION C"
		    "ALLS.";
	    s_cat(ch__7, a__4, i__5, &c__5, (ftnlen)127);
	    xermsg_("SLATEC", "DERKFS", ch__7, &c__10, &c__1, (ftnlen)6, (
		    ftnlen)6, (ftnlen)127);
	    *idid = -33;
	}

	if (*init != 1) {
	    if (*dtsign * (*tout - *t) < 0.) {
		s_wsfi(&io___22);
		do_fio(&c__1, (char *)&(*tout), (ftnlen)sizeof(real));
		e_wsfi();
/* Writing concatenation */
		i__5[0] = 42, a__4[0] = "IN DERKF, BY CALLING THE CODE WITH "
			"TOUT = ";
		i__5[1] = 16, a__4[1] = xern3;
		i__5[2] = 20, a__4[2] = " YOU ARE ATTEMPTING ";
		i__5[3] = 49, a__4[3] = "TO CHANGE THE DIRECTION OF INTEGRAT"
			"ION.$$THIS IS ";
		i__5[4] = 31, a__4[4] = "NOT ALLOWED WITHOUT RESTARTING.";
		s_cat(ch__8, a__4, i__5, &c__5, (ftnlen)158);
		xermsg_("SLATEC", "DERKFS", ch__8, &c__11, &c__1, (ftnlen)6, (
			ftnlen)6, (ftnlen)158);
		*idid = -33;
	    }
	}
    }

/*     INVALID INPUT DETECTED */

    if (*idid == -33) {
	if (*iquit != -33) {
	    *iquit = -33;
	    goto L909;
	} else {
	    xermsg_("SLATEC", "DERKFS", "IN DERKF, INVALID INPUT WAS DETECTE"
		    "D ON SUCCESSIVE ENTRIES.  IT IS IMPOSSIBLE TO PROCEED BE"
		    "CAUSE YOU HAVE NOT CORRECTED THE PROBLEM, SO EXECUTION I"
		    "S BEING TERMINATED.", &c__12, &c__2, (ftnlen)6, (ftnlen)6,
		     (ftnlen)166);
	    return 0;
	}
    }

/* ....................................................................... */

/*     RTOL = ATOL = 0. IS ALLOWED AS VALID INPUT AND INTERPRETED AS */
/*     ASKING FOR THE MOST ACCURATE SOLUTION POSSIBLE. IN THIS CASE, */
/*     THE RELATIVE ERROR TOLERANCE RTOL IS RESET TO THE SMALLEST VALUE */
/*     RER WHICH IS LIKELY TO BE REASONABLE FOR THIS METHOD AND MACHINE. */

    i__2 = *neq;
    for (k = 1; k <= i__2; ++k) {
	if (rtol[k] + atol[k] > 0.f) {
	    goto L45;
	}
	rtol[k] = *rer;
	*idid = -2;
L45:
	if (info[2] == 0) {
	    goto L55;
	}
/* L50: */
    }

L55:
    if (*idid != -2) {
	goto L60;
    }

/*                       RTOL=ATOL=0 ON INPUT, SO RTOL WAS CHANGED TO A */
/*                                                SMALL POSITIVE VALUE */
    *tolfac = 1.f;
    goto L909;

/*     BRANCH ON STATUS OF INITIALIZATION INDICATOR */
/*            INIT=0 MEANS INITIAL DERIVATIVES AND STARTING STEP SIZE */
/*                   NOT YET COMPUTED */
/*            INIT=1 MEANS STARTING STEP SIZE NOT YET COMPUTED */
/*            INIT=2 MEANS NO FURTHER INITIALIZATION REQUIRED */

L60:
    if (*init == 0) {
	goto L65;
    }
    if (*init == 1) {
	goto L70;
    }
    goto L80;

/* ....................................................................... */

/*     MORE INITIALIZATION -- */
/*                         -- EVALUATE INITIAL DERIVATIVES */

L65:
    *init = 1;
    a = *t;
    (*f)(&a, &y[1], &yp[1], &rpar[1], &ipar[1]);
    if (*t == *tout) {
	goto L666;
    }

/*                         -- SET SIGN OF INTEGRATION DIRECTION  AND */
/*                         -- ESTIMATE STARTING STEP SIZE */

L70:
    *init = 2;
    r__1 = *tout - *t;
    *dtsign = r_sign(&c_b110, &r__1);
    u = r1mach_(&c__4);
    big = sqrt(r1mach_(&c__2));
    d__1 = (doublereal) u;
    ute = pow_dd(&d__1, &c_b113);
    dy = ute * hvnrm_(&y[1], neq);
    if (dy == 0.f) {
	dy = ute;
    }
    ktol = 1;
    i__2 = *neq;
    for (k = 1; k <= i__2; ++k) {
	if (info[2] == 1) {
	    ktol = k;
	}
	tol = rtol[ktol] * (r__1 = y[k], dabs(r__1)) + atol[ktol];
	if (tol == 0.f) {
	    tol = dy * rtol[ktol];
	}
/* L75: */
	f1[k] = tol;
    }

    hstart_((S_fp)f, neq, t, tout, &y[1], &yp[1], &f1[1], &c__4, &u, &big, &
	    f2[1], &f3[1], &f4[1], &f5[1], &rpar[1], &ipar[1], h__);

/* ....................................................................... */

/*     SET STEP SIZE FOR INTEGRATION IN THE DIRECTION FROM T TO TOUT */
/*     AND SET OUTPUT POINT INDICATOR */

L80:
    dt = *tout - *t;
    *h__ = r_sign(h__, &dt);
    output = FALSE_;

/*     TEST TO SEE IF DERKF IS BEING SEVERELY IMPACTED BY TOO MANY */
/*     OUTPUT POINTS */

    if (dabs(*h__) >= dabs(dt) * 2.f) {
	++(*kop);
    }
    if (*kop <= mxkop) {
	goto L85;
    }

/*                       UNNECESSARY FREQUENCY OF OUTPUT IS RESTRICTING */
/*                                                 THE STEP SIZE CHOICE */
    *idid = -5;
    *kop = 0;
    goto L909;

L85:
    if (dabs(dt) > *u26 * dabs(*t)) {
	goto L100;
    }

/*     IF TOO CLOSE TO OUTPUT POINT,EXTRAPOLATE AND RETURN */

    i__2 = *neq;
    for (k = 1; k <= i__2; ++k) {
/* L90: */
	y[k] += dt * yp[k];
    }
    a = *tout;
    (*f)(&a, &y[1], &yp[1], &rpar[1], &ipar[1]);
    ++(*ksteps);
    goto L666;

/* ********************************************************************** */
/* ********************************************************************** */
/*     STEP BY STEP INTEGRATION */

L100:
    hfaild = FALSE_;

/*     TO PROTECT AGAINST IMPOSSIBLE ACCURACY REQUESTS, COMPUTE A */
/*     TOLERANCE FACTOR BASED ON THE REQUESTED ERROR TOLERANCE AND A */
/*     LEVEL OF ACCURACY ACHIEVABLE AT LIMITING PRECISION */

    *tolfac = 0.f;
    ktol = 1;
    i__2 = *neq;
    for (k = 1; k <= i__2; ++k) {
	if (info[2] == 1) {
	    ktol = k;
	}
	et = rtol[ktol] * (r__1 = y[k], dabs(r__1)) + atol[ktol];
	if (et > 0.f) {
	    goto L120;
	}
/* Computing MAX */
	r__1 = *tolfac, r__2 = *rer / rtol[ktol];
	*tolfac = dmax(r__1,r__2);
	goto L125;
L120:
/* Computing MAX */
	r__2 = *tolfac, r__3 = (r__1 = y[k], dabs(r__1)) * (*rer / et);
	*tolfac = dmax(r__2,r__3);
L125:
	;
    }
    if (*tolfac <= 1.f) {
	goto L150;
    }

/*                       REQUESTED ERROR UNATTAINABLE DUE TO LIMITED */
/*                                               PRECISION AVAILABLE */
    *tolfac *= 2.f;
    *idid = -2;
    goto L909;

/*     SET SMALLEST ALLOWABLE STEP SIZE */

L150:
    hmin = *u26 * dabs(*t);

/*     ADJUST STEP SIZE IF NECESSARY TO HIT THE OUTPUT POINT -- */
/*     LOOK AHEAD TWO STEPS TO AVOID DRASTIC CHANGES IN THE STEP SIZE AND */
/*     THUS LESSEN THE IMPACT OF OUTPUT POINTS ON THE CODE. */
/*     STRETCH THE STEP SIZE BY, AT MOST, AN AMOUNT EQUAL TO THE */
/*     SAFETY FACTOR OF 9/10. */

    dt = *tout - *t;
    if (dabs(dt) >= dabs(*h__) * 2.f) {
	goto L200;
    }
    if (dabs(dt) > dabs(*h__) / .9f) {
	goto L175;
    }

/*     THE NEXT STEP, IF SUCCESSFUL, WILL COMPLETE THE INTEGRATION TO */
/*     THE OUTPUT POINT */

    output = TRUE_;
    *h__ = dt;
    goto L200;

L175:
    *h__ = dt * .5f;


/* ********************************************************************** */
/*     CORE INTEGRATOR FOR TAKING A SINGLE STEP */
/* ********************************************************************** */
/*     TO AVOID PROBLEMS WITH ZERO CROSSINGS, RELATIVE ERROR IS MEASURED */
/*     USING THE AVERAGE OF THE MAGNITUDES OF THE SOLUTION AT THE */
/*     BEGINNING AND END OF A STEP. */
/*     THE ERROR ESTIMATE FORMULA HAS BEEN GROUPED TO CONTROL LOSS OF */
/*     SIGNIFICANCE. */
/*     LOCAL ERROR ESTIMATES FOR A FIRST ORDER METHOD USING THE SAME */
/*     STEP SIZE AS THE FEHLBERG METHOD ARE CALCULATED AS PART OF THE */
/*     TEST FOR STIFFNESS. */
/*     TO DISTINGUISH THE VARIOUS ARGUMENTS, H IS NOT PERMITTED */
/*     TO BECOME SMALLER THAN 26 UNITS OF ROUNDOFF IN T. */
/*     PRACTICAL LIMITS ON THE CHANGE IN THE STEP SIZE ARE ENFORCED TO */
/*     SMOOTH THE STEP SIZE SELECTION PROCESS AND TO AVOID EXCESSIVE */
/*     CHATTERING ON PROBLEMS HAVING DISCONTINUITIES. */
/*     TO PREVENT UNNECESSARY FAILURES, THE CODE USES 9/10 THE STEP SIZE */
/*     IT ESTIMATES WILL SUCCEED. */
/*     AFTER A STEP FAILURE, THE STEP SIZE IS NOT ALLOWED TO INCREASE FOR */
/*     THE NEXT ATTEMPTED STEP. THIS MAKES THE CODE MORE EFFICIENT ON */
/*     PROBLEMS HAVING DISCONTINUITIES AND MORE EFFECTIVE IN GENERAL */
/*     SINCE LOCAL EXTRAPOLATION IS BEING USED AND EXTRA CAUTION SEEMS */
/*     WARRANTED. */
/* ....................................................................... */

/*     MONITOR NUMBER OF STEPS ATTEMPTED */

L200:
    if (*ksteps <= mxstep) {
	goto L222;
    }

/*                       A SIGNIFICANT AMOUNT OF WORK HAS BEEN EXPENDED */
    *idid = -1;
    *ksteps = 0;
    if (! (*stiff)) {
	goto L909;
    }

/*                       PROBLEM APPEARS TO BE STIFF */
    *idid = -4;
    *stiff = FALSE_;
    *nonstf = FALSE_;
    *ntstep = 0;
    *nstifs = 0;
    goto L909;

/*     ADVANCE AN APPROXIMATE SOLUTION OVER ONE STEP OF LENGTH H */

L222:
    defehl_((S_fp)f, neq, t, &y[1], h__, &yp[1], &f1[1], &f2[1], &f3[1], &f4[
	    1], &f5[1], &ys[1], &rpar[1], &ipar[1]);
    ++(*ksteps);

/* ....................................................................... */

/*     COMPUTE AND TEST ALLOWABLE TOLERANCES VERSUS LOCAL ERROR */
/*     ESTIMATES.  NOTE THAT RELATIVE ERROR IS MEASURED WITH RESPECT TO */
/*     THE AVERAGE OF THE MAGNITUDES OF THE SOLUTION AT THE BEGINNING */
/*     AND END OF THE STEP. */
/*     LOCAL ERROR ESTIMATES FOR A SPECIAL FIRST ORDER METHOD ARE */
/*     CALCULATED ONLY WHEN THE STIFFNESS DETECTION IS TURNED ON. */

    eeoet = 0.f;
    estiff = 0.f;
    ktol = 1;
    i__2 = *neq;
    for (k = 1; k <= i__2; ++k) {
	yavg = ((r__1 = y[k], dabs(r__1)) + (r__2 = ys[k], dabs(r__2))) * .5f;
	if (info[2] == 1) {
	    ktol = k;
	}
	et = rtol[ktol] * yavg + atol[ktol];
	if (et > 0.f) {
	    goto L325;
	}

/*                       PURE RELATIVE ERROR INAPPROPRIATE WHEN SOLUTION */
/*                                                              VANISHES */
	*idid = -3;
	goto L909;

L325:
	ee = (r__1 = yp[k] * -2090.f + (f3[k] * 21970.f - f4[k] * 15048.f) + (
		f2[k] * 22528.f - f5[k] * 27360.f), dabs(r__1));
	if (*stiff || *nonstf) {
	    goto L350;
	}
	es = (r__1 = *h__ * (yp[k] * .055455f - f1[k] * .035493f - f2[k] * 
		.036571f + f3[k] * .023107f - f4[k] * .009515f + f5[k] * 
		.003017f), dabs(r__1));
/* Computing MAX */
	r__1 = estiff, r__2 = es / et;
	estiff = dmax(r__1,r__2);
L350:
/* Computing MAX */
	r__1 = eeoet, r__2 = ee / et;
	eeoet = dmax(r__1,r__2);
    }

    esttol = dabs(*h__) * eeoet / 752400.f;

    if (esttol <= 1.f) {
	goto L500;
    }

/* ....................................................................... */

/*     UNSUCCESSFUL STEP */

    if (dabs(*h__) > hmin) {
	goto L400;
    }

/*                       REQUESTED ERROR UNATTAINABLE AT SMALLEST */
/*                                            ALLOWABLE STEP SIZE */
    *tolfac = esttol * 1.69f;
    *idid = -2;
    goto L909;

/*                       REDUCE THE STEP SIZE , TRY AGAIN */
/*                       THE DECREASE IS LIMITED TO A FACTOR OF 1/10 */

L400:
    hfaild = TRUE_;
    output = FALSE_;
    s = .1f;
    if (esttol < 59049.f) {
	d__1 = (doublereal) esttol;
	s = .9f / pow_dd(&d__1, &c_b129);
    }
/* Computing MAX */
    r__2 = s * dabs(*h__);
    r__1 = dmax(r__2,hmin);
    *h__ = r_sign(&r__1, h__);
    goto L200;

/* ....................................................................... */

/*     SUCCESSFUL STEP */
/*                       STORE SOLUTION AT T+H */
/*                       AND EVALUATE DERIVATIVES THERE */

L500:
    *t += *h__;
    i__2 = *neq;
    for (k = 1; k <= i__2; ++k) {
/* L525: */
	y[k] = ys[k];
    }
    a = *t;
    (*f)(&a, &y[1], &yp[1], &rpar[1], &ipar[1]);

/*                       CHOOSE NEXT STEP SIZE */
/*                       THE INCREASE IS LIMITED TO A FACTOR OF 5 */
/*                       IF STEP FAILURE HAS JUST OCCURRED, NEXT */
/*                          STEP SIZE IS NOT ALLOWED TO INCREASE */

    s = 5.f;
    if (esttol > 1.889568e-4f) {
	d__1 = (doublereal) esttol;
	s = .9f / pow_dd(&d__1, &c_b129);
    }
    if (hfaild) {
	s = dmin(s,1.f);
    }
/* Computing MAX */
    r__2 = s * dabs(*h__);
    r__1 = dmax(r__2,hmin);
    *h__ = r_sign(&r__1, h__);

/* ....................................................................... */

/*     CHECK FOR STIFFNESS (IF NOT ALREADY DETECTED) */

/*     IN A SEQUENCE OF 50 SUCCESSFUL STEPS BY THE FEHLBERG METHOD, 25 */
/*     SUCCESSFUL STEPS BY THE FIRST ORDER METHOD INDICATES STIFFNESS */
/*     AND TURNS THE TEST OFF. IF 26 FAILURES BY THE FIRST ORDER METHOD */
/*     OCCUR, THE TEST IS TURNED OFF UNTIL THIS SEQUENCE OF 50 STEPS */
/*     BY THE FEHLBERG METHOD IS COMPLETED. */

    if (*stiff) {
	goto L600;
    }
    *ntstep = (*ntstep + 1) % 50;
    if (*ntstep == 1) {
	*nonstf = FALSE_;
    }
    if (*nonstf) {
	goto L600;
    }
    if (estiff > 1.f) {
	goto L550;
    }

/*                       SUCCESSFUL STEP WITH FIRST ORDER METHOD */
    ++(*nstifs);
/*                       TURN TEST OFF AFTER 25 INDICATIONS OF STIFFNESS */
    if (*nstifs == 25) {
	*stiff = TRUE_;
    }
    goto L600;

/*                       UNSUCCESSFUL STEP WITH FIRST ORDER METHOD */
L550:
    if (*ntstep - *nstifs <= 25) {
	goto L600;
    }
/*                       TURN STIFFNESS DETECTION OFF FOR THIS BLOCK OF */
/*                                                          FIFTY STEPS */
    *nonstf = TRUE_;
/*                       RESET STIFF STEP COUNTER */
    *nstifs = 0;

/* ********************************************************************** */
/*     END OF CORE INTEGRATOR */
/* ********************************************************************** */


/*     SHOULD WE TAKE ANOTHER STEP */

L600:
    if (output) {
	goto L666;
    }
    if (info[3] == 0) {
	goto L100;
    }

/* ********************************************************************** */
/* ********************************************************************** */

/*     INTEGRATION SUCCESSFULLY COMPLETED */

/*                 ONE-STEP MODE */
    *idid = 1;
    *told = *t;
    return 0;

/*                 INTERVAL MODE */
L666:
    *idid = 2;
    *t = *tout;
    *told = *t;
    return 0;

/*     INTEGRATION TASK INTERRUPTED */

L909:
    info[1] = -1;
    *told = *t;
    if (*idid != -2) {
	return 0;
    }

/*                       THE ERROR TOLERANCES ARE INCREASED TO VALUES */
/*                               WHICH ARE APPROPRIATE FOR CONTINUING */
    rtol[1] = *tolfac * rtol[1];
    atol[1] = *tolfac * atol[1];
    if (info[2] == 0) {
	return 0;
    }
    i__2 = *neq;
    for (k = 2; k <= i__2; ++k) {
	rtol[k] = *tolfac * rtol[k];
/* L939: */
	atol[k] = *tolfac * atol[k];
    }
    return 0;
} /* derkfs_ */

