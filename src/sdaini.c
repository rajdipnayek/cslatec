/* sdaini.f -- translated by f2c (version 12.02.01).
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

/* DECK SDAINI */
/* Subroutine */ int sdaini_(real *x, real *y, real *yprime, integer *neq, 
	S_fp res, U_fp jac, real *h__, real *wt, integer *idid, real *rpar, 
	integer *ipar, real *phi, real *delta, real *e, real *wm, integer *
	iwm, real *hmin, real *uround, integer *nonneg, integer *ntemp)
{
    /* Initialized data */

    static integer maxit = 10;
    static integer mjac = 5;
    static real damp = .75f;

    /* System generated locals */
    integer phi_dim1, phi_offset, i__1;
    real r__1, r__2;
    doublereal d__1, d__2;

    /* Local variables */
    static integer i__, m;
    static real r__, s, cj;
    static integer ncf, nef, ier, nsf;
    static real err, rate;
    static integer ires;
    static real xold;
    static integer jcalc;
    static real ynorm;
    extern /* Subroutine */ int sdajac_(integer *, real *, real *, real *, 
	    real *, real *, real *, integer *, real *, real *, real *, 
	    integer *, S_fp, integer *, real *, U_fp, real *, integer *, 
	    integer *);
    static logical convgd;
    static real delnrm;
    extern doublereal sdanrm_(integer *, real *, real *, real *, integer *);
    static real oldnrm;
    extern /* Subroutine */ int sdaslv_(integer *, real *, real *, integer *);

/* ***BEGIN PROLOGUE  SDAINI */
/* ***SUBSIDIARY */
/* ***PURPOSE  Initialization routine for SDASSL. */
/* ***LIBRARY   SLATEC (DASSL) */
/* ***TYPE      SINGLE PRECISION (SDAINI-S, DDAINI-D) */
/* ***AUTHOR  Petzold, Linda R., (LLNL) */
/* ***DESCRIPTION */
/* ----------------------------------------------------------------- */
/*     SDAINI TAKES ONE STEP OF SIZE H OR SMALLER */
/*     WITH THE BACKWARD EULER METHOD, TO */
/*     FIND YPRIME.  X AND Y ARE UPDATED TO BE CONSISTENT WITH THE */
/*     NEW STEP.  A MODIFIED DAMPED NEWTON ITERATION IS USED TO */
/*     SOLVE THE CORRECTOR ITERATION. */

/*     THE INITIAL GUESS FOR YPRIME IS USED IN THE */
/*     PREDICTION, AND IN FORMING THE ITERATION */
/*     MATRIX, BUT IS NOT INVOLVED IN THE */
/*     ERROR TEST. THIS MAY HAVE TROUBLE */
/*     CONVERGING IF THE INITIAL GUESS IS NO */
/*     GOOD, OR IF G(X,Y,YPRIME) DEPENDS */
/*     NONLINEARLY ON YPRIME. */

/*     THE PARAMETERS REPRESENT: */
/*     X --         INDEPENDENT VARIABLE */
/*     Y --         SOLUTION VECTOR AT X */
/*     YPRIME --    DERIVATIVE OF SOLUTION VECTOR */
/*     NEQ --       NUMBER OF EQUATIONS */
/*     H --         STEPSIZE. IMDER MAY USE A STEPSIZE */
/*                  SMALLER THAN H. */
/*     WT --        VECTOR OF WEIGHTS FOR ERROR */
/*                  CRITERION */
/*     IDID --      COMPLETION CODE WITH THE FOLLOWING MEANINGS */
/*                  IDID= 1 -- YPRIME WAS FOUND SUCCESSFULLY */
/*                  IDID=-12 -- SDAINI FAILED TO FIND YPRIME */
/*     RPAR,IPAR -- REAL AND INTEGER PARAMETER ARRAYS */
/*                  THAT ARE NOT ALTERED BY SDAINI */
/*     PHI --       WORK SPACE FOR SDAINI */
/*     DELTA,E --   WORK SPACE FOR SDAINI */
/*     WM,IWM --    REAL AND INTEGER ARRAYS STORING */
/*                  MATRIX INFORMATION */

/* ----------------------------------------------------------------- */
/* ***ROUTINES CALLED  SDAJAC, SDANRM, SDASLV */
/* ***REVISION HISTORY  (YYMMDD) */
/*   830315  DATE WRITTEN */
/*   901009  Finished conversion to SLATEC 4.0 format (F.N.Fritsch) */
/*   901019  Merged changes made by C. Ulrich with SLATEC 4.0 format. */
/*   901026  Added explicit declarations for all variables and minor */
/*           cosmetic changes to prologue.  (FNF) */
/*   901030  Minor corrections to declarations.  (FNF) */
/* ***END PROLOGUE  SDAINI */





    /* Parameter adjustments */
    --y;
    --yprime;
    phi_dim1 = *neq;
    phi_offset = 1 + phi_dim1;
    phi -= phi_offset;
    --wt;
    --rpar;
    --ipar;
    --delta;
    --e;
    --wm;
    --iwm;

    /* Function Body */


/* --------------------------------------------------- */
/*     BLOCK 1. */
/*     INITIALIZATIONS. */
/* --------------------------------------------------- */

/* ***FIRST EXECUTABLE STATEMENT  SDAINI */
    *idid = 1;
    nef = 0;
    ncf = 0;
    nsf = 0;
    xold = *x;
    ynorm = sdanrm_(neq, &y[1], &wt[1], &rpar[1], &ipar[1]);

/*     SAVE Y AND YPRIME IN PHI */
    i__1 = *neq;
    for (i__ = 1; i__ <= i__1; ++i__) {
	phi[i__ + phi_dim1] = y[i__];
/* L100: */
	phi[i__ + (phi_dim1 << 1)] = yprime[i__];
    }


/* ---------------------------------------------------- */
/*     BLOCK 2. */
/*     DO ONE BACKWARD EULER STEP. */
/* ---------------------------------------------------- */

/*     SET UP FOR START OF CORRECTOR ITERATION */
L200:
    cj = 1.f / *h__;
    *x += *h__;

/*     PREDICT SOLUTION AND DERIVATIVE */
    i__1 = *neq;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L250: */
	y[i__] += *h__ * yprime[i__];
    }

    jcalc = -1;
    m = 0;
    convgd = TRUE_;


/*     CORRECTOR LOOP. */
L300:
    ++iwm[12];
    ires = 0;

    (*res)(x, &y[1], &yprime[1], &delta[1], &ires, &rpar[1], &ipar[1]);
    if (ires < 0) {
	goto L430;
    }


/*     EVALUATE THE ITERATION MATRIX */
    if (jcalc != -1) {
	goto L310;
    }
    ++iwm[13];
    jcalc = 0;
    sdajac_(neq, x, &y[1], &yprime[1], &delta[1], &cj, h__, &ier, &wt[1], &e[
	    1], &wm[1], &iwm[1], (S_fp)res, &ires, uround, (U_fp)jac, &rpar[1]
	    , &ipar[1], ntemp);

    s = 1e6f;
    if (ires < 0) {
	goto L430;
    }
    if (ier != 0) {
	goto L430;
    }
    nsf = 0;



/*     MULTIPLY RESIDUAL BY DAMPING FACTOR */
L310:
    i__1 = *neq;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L320: */
	delta[i__] *= damp;
    }

/*     COMPUTE A NEW ITERATE (BACK SUBSTITUTION) */
/*     STORE THE CORRECTION IN DELTA */

    sdaslv_(neq, &delta[1], &wm[1], &iwm[1]);

/*     UPDATE Y AND YPRIME */
    i__1 = *neq;
    for (i__ = 1; i__ <= i__1; ++i__) {
	y[i__] -= delta[i__];
/* L330: */
	yprime[i__] -= cj * delta[i__];
    }

/*     TEST FOR CONVERGENCE OF THE ITERATION. */

    delnrm = sdanrm_(neq, &delta[1], &wt[1], &rpar[1], &ipar[1]);
    if (delnrm <= *uround * 100.f * ynorm) {
	goto L400;
    }

    if (m > 0) {
	goto L340;
    }
    oldnrm = delnrm;
    goto L350;

L340:
    d__1 = (doublereal) (delnrm / oldnrm);
    d__2 = (doublereal) (1.f / m);
    rate = pow_dd(&d__1, &d__2);
    if (rate > .9f) {
	goto L430;
    }
    s = rate / (1.f - rate);

L350:
    if (s * delnrm <= .33f) {
	goto L400;
    }


/*     THE CORRECTOR HAS NOT YET CONVERGED. UPDATE */
/*     M AND AND TEST WHETHER THE MAXIMUM */
/*     NUMBER OF ITERATIONS HAVE BEEN TRIED. */
/*     EVERY MJAC ITERATIONS, GET A NEW */
/*     ITERATION MATRIX. */

    ++m;
    if (m >= maxit) {
	goto L430;
    }

    if (m / mjac * mjac == m) {
	jcalc = -1;
    }
    goto L300;


/*     THE ITERATION HAS CONVERGED. */
/*     CHECK NONNEGATIVITY CONSTRAINTS */
L400:
    if (*nonneg == 0) {
	goto L450;
    }
    i__1 = *neq;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L410: */
/* Computing MIN */
	r__1 = y[i__];
	delta[i__] = dmin(r__1,0.f);
    }

    delnrm = sdanrm_(neq, &delta[1], &wt[1], &rpar[1], &ipar[1]);
    if (delnrm > .33f) {
	goto L430;
    }

    i__1 = *neq;
    for (i__ = 1; i__ <= i__1; ++i__) {
	y[i__] -= delta[i__];
/* L420: */
	yprime[i__] -= cj * delta[i__];
    }
    goto L450;


/*     EXITS FROM CORRECTOR LOOP. */
L430:
    convgd = FALSE_;
L450:
    if (! convgd) {
	goto L600;
    }



/* ----------------------------------------------------- */
/*     BLOCK 3. */
/*     THE CORRECTOR ITERATION CONVERGED. */
/*     DO ERROR TEST. */
/* ----------------------------------------------------- */

    i__1 = *neq;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L510: */
	e[i__] = y[i__] - phi[i__ + phi_dim1];
    }
    err = sdanrm_(neq, &e[1], &wt[1], &rpar[1], &ipar[1]);

    if (err <= 1.f) {
	return 0;
    }



/* -------------------------------------------------------- */
/*     BLOCK 4. */
/*     THE BACKWARD EULER STEP FAILED. RESTORE X, Y */
/*     AND YPRIME TO THEIR ORIGINAL VALUES. */
/*     REDUCE STEPSIZE AND TRY AGAIN, IF */
/*     POSSIBLE. */
/* --------------------------------------------------------- */

L600:
    *x = xold;
    i__1 = *neq;
    for (i__ = 1; i__ <= i__1; ++i__) {
	y[i__] = phi[i__ + phi_dim1];
/* L610: */
	yprime[i__] = phi[i__ + (phi_dim1 << 1)];
    }

    if (convgd) {
	goto L640;
    }
    if (ier == 0) {
	goto L620;
    }
    ++nsf;
    *h__ *= .25f;
    if (nsf < 3 && dabs(*h__) >= *hmin) {
	goto L690;
    }
    *idid = -12;
    return 0;
L620:
    if (ires > -2) {
	goto L630;
    }
    *idid = -12;
    return 0;
L630:
    ++ncf;
    *h__ *= .25f;
    if (ncf < 10 && dabs(*h__) >= *hmin) {
	goto L690;
    }
    *idid = -12;
    return 0;

L640:
    ++nef;
    r__ = .9f / (err * 2.f + 1e-4f);
/* Computing MAX */
    r__1 = .1f, r__2 = dmin(.5f,r__);
    r__ = dmax(r__1,r__2);
    *h__ *= r__;
    if (dabs(*h__) >= *hmin && nef < 10) {
	goto L690;
    }
    *idid = -12;
    return 0;
L690:
    goto L200;

/* -------------END OF SUBROUTINE SDAINI---------------------- */
} /* sdaini_ */

