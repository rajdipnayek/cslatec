/* sdajac.f -- translated by f2c (version 12.02.01).
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

/* DECK SDAJAC */
/* Subroutine */ int sdajac_(integer *neq, real *x, real *y, real *yprime, 
	real *delta, real *cj, real *h__, integer *ier, real *wt, real *e, 
	real *wm, integer *iwm, S_fp res, integer *ires, real *uround, S_fp 
	jac, real *rpar, integer *ipar, integer *ntemp)
{
    /* System generated locals */
    integer i__1, i__2, i__3, i__4, i__5;
    real r__1, r__2, r__3, r__4, r__5;

    /* Local variables */
    static integer i__, j, k, l, n, i1, i2, ii, mba;
    static real del;
    static integer meb1, nrow;
    static real squr;
    static integer npdm1, mband;
    extern /* Subroutine */ int sgbfa_(real *, integer *, integer *, integer *
	    , integer *, integer *, integer *), sgefa_(real *, integer *, 
	    integer *, integer *, integer *);
    static integer lenpd, isave, msave;
    static real ysave;
    static integer mtype, meband;
    static real delinv;
    static integer ipsave;
    static real ypsave;

/* ***BEGIN PROLOGUE  SDAJAC */
/* ***SUBSIDIARY */
/* ***PURPOSE  Compute the iteration matrix for SDASSL and form the */
/*            LU-decomposition. */
/* ***LIBRARY   SLATEC (DASSL) */
/* ***TYPE      SINGLE PRECISION (SDAJAC-S, DDAJAC-D) */
/* ***AUTHOR  Petzold, Linda R., (LLNL) */
/* ***DESCRIPTION */
/* ----------------------------------------------------------------------- */
/*     THIS ROUTINE COMPUTES THE ITERATION MATRIX */
/*     PD=DG/DY+CJ*DG/DYPRIME (WHERE G(X,Y,YPRIME)=0). */
/*     HERE PD IS COMPUTED BY THE USER-SUPPLIED */
/*     ROUTINE JAC IF IWM(MTYPE) IS 1 OR 4, AND */
/*     IT IS COMPUTED BY NUMERICAL FINITE DIFFERENCING */
/*     IF IWM(MTYPE)IS 2 OR 5 */
/*     THE PARAMETERS HAVE THE FOLLOWING MEANINGS. */
/*     Y        = ARRAY CONTAINING PREDICTED VALUES */
/*     YPRIME   = ARRAY CONTAINING PREDICTED DERIVATIVES */
/*     DELTA    = RESIDUAL EVALUATED AT (X,Y,YPRIME) */
/*                (USED ONLY IF IWM(MTYPE)=2 OR 5) */
/*     CJ       = SCALAR PARAMETER DEFINING ITERATION MATRIX */
/*     H        = CURRENT STEPSIZE IN INTEGRATION */
/*     IER      = VARIABLE WHICH IS .NE. 0 */
/*                IF ITERATION MATRIX IS SINGULAR, */
/*                AND 0 OTHERWISE. */
/*     WT       = VECTOR OF WEIGHTS FOR COMPUTING NORMS */
/*     E        = WORK SPACE (TEMPORARY) OF LENGTH NEQ */
/*     WM       = REAL WORK SPACE FOR MATRICES. ON */
/*                OUTPUT IT CONTAINS THE LU DECOMPOSITION */
/*                OF THE ITERATION MATRIX. */
/*     IWM      = INTEGER WORK SPACE CONTAINING */
/*                MATRIX INFORMATION */
/*     RES      = NAME OF THE EXTERNAL USER-SUPPLIED ROUTINE */
/*                TO EVALUATE THE RESIDUAL FUNCTION G(X,Y,YPRIME) */
/*     IRES     = FLAG WHICH IS EQUAL TO ZERO IF NO ILLEGAL VALUES */
/*                IN RES, AND LESS THAN ZERO OTHERWISE.  (IF IRES */
/*                IS LESS THAN ZERO, THE MATRIX WAS NOT COMPLETED) */
/*                IN THIS CASE (IF IRES .LT. 0), THEN IER = 0. */
/*     UROUND   = THE UNIT ROUNDOFF ERROR OF THE MACHINE BEING USED. */
/*     JAC      = NAME OF THE EXTERNAL USER-SUPPLIED ROUTINE */
/*                TO EVALUATE THE ITERATION MATRIX (THIS ROUTINE */
/*                IS ONLY USED IF IWM(MTYPE) IS 1 OR 4) */
/* ----------------------------------------------------------------------- */
/* ***ROUTINES CALLED  SGBFA, SGEFA */
/* ***REVISION HISTORY  (YYMMDD) */
/*   830315  DATE WRITTEN */
/*   901009  Finished conversion to SLATEC 4.0 format (F.N.Fritsch) */
/*   901010  Modified three MAX calls to be all on one line.  (FNF) */
/*   901019  Merged changes made by C. Ulrich with SLATEC 4.0 format. */
/*   901026  Added explicit declarations for all variables and minor */
/*           cosmetic changes to prologue.  (FNF) */
/*   901101  Corrected PURPOSE.  (FNF) */
/* ***END PROLOGUE  SDAJAC */





/* ***FIRST EXECUTABLE STATEMENT  SDAJAC */
    /* Parameter adjustments */
    --ipar;
    --rpar;
    --iwm;
    --wm;
    --e;
    --wt;
    --delta;
    --yprime;
    --y;

    /* Function Body */
    *ier = 0;
    npdm1 = 0;
    mtype = iwm[4];
    switch (mtype) {
	case 1:  goto L100;
	case 2:  goto L200;
	case 3:  goto L300;
	case 4:  goto L400;
	case 5:  goto L500;
    }


/*     DENSE USER-SUPPLIED MATRIX */
L100:
    lenpd = *neq * *neq;
    i__1 = lenpd;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L110: */
	wm[npdm1 + i__] = 0.f;
    }
    (*jac)(x, &y[1], &yprime[1], &wm[1], cj, &rpar[1], &ipar[1]);
    goto L230;


/*     DENSE FINITE-DIFFERENCE-GENERATED MATRIX */
L200:
    *ires = 0;
    nrow = npdm1;
    squr = sqrt(*uround);
    i__1 = *neq;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing MAX */
	r__4 = (r__1 = y[i__], dabs(r__1)), r__5 = (r__2 = *h__ * yprime[i__],
		 dabs(r__2)), r__4 = max(r__4,r__5), r__5 = (r__3 = wt[i__], 
		dabs(r__3));
	del = squr * dmax(r__4,r__5);
	r__1 = *h__ * yprime[i__];
	del = r_sign(&del, &r__1);
	del = y[i__] + del - y[i__];
	ysave = y[i__];
	ypsave = yprime[i__];
	y[i__] += del;
	yprime[i__] += *cj * del;
	(*res)(x, &y[1], &yprime[1], &e[1], ires, &rpar[1], &ipar[1]);
	if (*ires < 0) {
	    return 0;
	}
	delinv = 1.f / del;
	i__2 = *neq;
	for (l = 1; l <= i__2; ++l) {
/* L220: */
	    wm[nrow + l] = (e[l] - delta[l]) * delinv;
	}
	nrow += *neq;
	y[i__] = ysave;
	yprime[i__] = ypsave;
/* L210: */
    }


/*     DO DENSE-MATRIX LU DECOMPOSITION ON PD */
L230:
    sgefa_(&wm[1], neq, neq, &iwm[21], ier);
    return 0;


/*     DUMMY SECTION FOR IWM(MTYPE)=3 */
L300:
    return 0;


/*     BANDED USER-SUPPLIED MATRIX */
L400:
    lenpd = ((iwm[1] << 1) + iwm[2] + 1) * *neq;
    i__1 = lenpd;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L410: */
	wm[npdm1 + i__] = 0.f;
    }
    (*jac)(x, &y[1], &yprime[1], &wm[1], cj, &rpar[1], &ipar[1]);
    meband = (iwm[1] << 1) + iwm[2] + 1;
    goto L550;


/*     BANDED FINITE-DIFFERENCE-GENERATED MATRIX */
L500:
    mband = iwm[1] + iwm[2] + 1;
    mba = min(mband,*neq);
    meband = mband + iwm[1];
    meb1 = meband - 1;
    msave = *neq / mband + 1;
    isave = *ntemp - 1;
    ipsave = isave + msave;
    *ires = 0;
    squr = sqrt(*uround);
    i__1 = mba;
    for (j = 1; j <= i__1; ++j) {
	i__2 = *neq;
	i__3 = mband;
	for (n = j; i__3 < 0 ? n >= i__2 : n <= i__2; n += i__3) {
	    k = (n - j) / mband + 1;
	    wm[isave + k] = y[n];
	    wm[ipsave + k] = yprime[n];
/* Computing MAX */
	    r__4 = (r__1 = y[n], dabs(r__1)), r__5 = (r__2 = *h__ * yprime[n],
		     dabs(r__2)), r__4 = max(r__4,r__5), r__5 = (r__3 = wt[n],
		     dabs(r__3));
	    del = squr * dmax(r__4,r__5);
	    r__1 = *h__ * yprime[n];
	    del = r_sign(&del, &r__1);
	    del = y[n] + del - y[n];
	    y[n] += del;
/* L510: */
	    yprime[n] += *cj * del;
	}
	(*res)(x, &y[1], &yprime[1], &e[1], ires, &rpar[1], &ipar[1]);
	if (*ires < 0) {
	    return 0;
	}
	i__3 = *neq;
	i__2 = mband;
	for (n = j; i__2 < 0 ? n >= i__3 : n <= i__3; n += i__2) {
	    k = (n - j) / mband + 1;
	    y[n] = wm[isave + k];
	    yprime[n] = wm[ipsave + k];
/* Computing MAX */
	    r__4 = (r__1 = y[n], dabs(r__1)), r__5 = (r__2 = *h__ * yprime[n],
		     dabs(r__2)), r__4 = max(r__4,r__5), r__5 = (r__3 = wt[n],
		     dabs(r__3));
	    del = squr * dmax(r__4,r__5);
	    r__1 = *h__ * yprime[n];
	    del = r_sign(&del, &r__1);
	    del = y[n] + del - y[n];
	    delinv = 1.f / del;
/* Computing MAX */
	    i__4 = 1, i__5 = n - iwm[2];
	    i1 = max(i__4,i__5);
/* Computing MIN */
	    i__4 = *neq, i__5 = n + iwm[1];
	    i2 = min(i__4,i__5);
	    ii = n * meb1 - iwm[1] + npdm1;
	    i__4 = i2;
	    for (i__ = i1; i__ <= i__4; ++i__) {
/* L520: */
		wm[ii + i__] = (e[i__] - delta[i__]) * delinv;
	    }
/* L530: */
	}
/* L540: */
    }


/*     DO LU DECOMPOSITION OF BANDED PD */
L550:
    sgbfa_(&wm[1], &meband, neq, &iwm[1], &iwm[2], &iwm[21], ier);
    return 0;
/* ------END OF SUBROUTINE SDAJAC------ */
} /* sdajac_ */

