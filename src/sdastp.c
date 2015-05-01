/* sdastp.f -- translated by f2c (version 12.02.01).
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

/* DECK SDASTP */
/* Subroutine */ int sdastp_(real *x, real *y, real *yprime, integer *neq, 
	S_fp res, U_fp jac, real *h__, real *wt, integer *jstart, integer *
	idid, real *rpar, integer *ipar, real *phi, real *delta, real *e, 
	real *wm, integer *iwm, real *alpha, real *beta, real *gamma, real *
	psi, real *sigma, real *cj, real *cjold, real *hold, real *s, real *
	hmin, real *uround, integer *iphase, integer *jcalc, integer *k, 
	integer *kold, integer *ns, integer *nonneg, integer *ntemp)
{
    /* Initialized data */

    static integer maxit = 4;
    static real xrate = .25f;

    /* System generated locals */
    integer phi_dim1, phi_offset, i__1, i__2;
    real r__1, r__2;
    doublereal d__1, d__2;

    /* Local variables */
    static integer i__, j, m;
    static real r__;
    static integer j1;
    static real ck;
    static integer km1, kp1, kp2, ncf, nef, ier;
    static real erk;
    static integer nsf;
    static real err, est;
    static integer nsp1;
    static real rate, hnew;
    static integer ires, knew;
    static real terk, xold, erkm1, erkm2, erkp1, temp1, temp2;
    static integer kdiff;
    static real enorm, pnorm, alpha0, terkm1, terkm2, terkp1;
    extern /* Subroutine */ int sdajac_(integer *, real *, real *, real *, 
	    real *, real *, real *, integer *, real *, real *, real *, 
	    integer *, S_fp, integer *, real *, U_fp, real *, integer *, 
	    integer *);
    static real alphas, cjlast, delnrm;
    static logical convgd;
    extern doublereal sdanrm_(integer *, real *, real *, real *, integer *);
    static real oldnrm;
    extern /* Subroutine */ int sdaslv_(integer *, real *, real *, integer *),
	     sdatrp_(real *, real *, real *, real *, integer *, integer *, 
	    real *, real *);

/* ***BEGIN PROLOGUE  SDASTP */
/* ***SUBSIDIARY */
/* ***PURPOSE  Perform one step of the SDASSL integration. */
/* ***LIBRARY   SLATEC (DASSL) */
/* ***TYPE      SINGLE PRECISION (SDASTP-S, DDASTP-D) */
/* ***AUTHOR  Petzold, Linda R., (LLNL) */
/* ***DESCRIPTION */
/* ----------------------------------------------------------------------- */
/*     SDASTP SOLVES A SYSTEM OF DIFFERENTIAL/ */
/*     ALGEBRAIC EQUATIONS OF THE FORM */
/*     G(X,Y,YPRIME) = 0,  FOR ONE STEP (NORMALLY */
/*     FROM X TO X+H). */

/*     THE METHODS USED ARE MODIFIED DIVIDED */
/*     DIFFERENCE,FIXED LEADING COEFFICIENT */
/*     FORMS OF BACKWARD DIFFERENTIATION */
/*     FORMULAS. THE CODE ADJUSTS THE STEPSIZE */
/*     AND ORDER TO CONTROL THE LOCAL ERROR PER */
/*     STEP. */


/*     THE PARAMETERS REPRESENT */
/*     X  --        INDEPENDENT VARIABLE */
/*     Y  --        SOLUTION VECTOR AT X */
/*     YPRIME --    DERIVATIVE OF SOLUTION VECTOR */
/*                  AFTER SUCCESSFUL STEP */
/*     NEQ --       NUMBER OF EQUATIONS TO BE INTEGRATED */
/*     RES --       EXTERNAL USER-SUPPLIED SUBROUTINE */
/*                  TO EVALUATE THE RESIDUAL.  THE CALL IS */
/*                  CALL RES(X,Y,YPRIME,DELTA,IRES,RPAR,IPAR) */
/*                  X,Y,YPRIME ARE INPUT.  DELTA IS OUTPUT. */
/*                  ON INPUT, IRES=0.  RES SHOULD ALTER IRES ONLY */
/*                  IF IT ENCOUNTERS AN ILLEGAL VALUE OF Y OR A */
/*                  STOP CONDITION.  SET IRES=-1 IF AN INPUT VALUE */
/*                  OF Y IS ILLEGAL, AND SDASTP WILL TRY TO SOLVE */
/*                  THE PROBLEM WITHOUT GETTING IRES = -1.  IF */
/*                  IRES=-2, SDASTP RETURNS CONTROL TO THE CALLING */
/*                  PROGRAM WITH IDID = -11. */
/*     JAC --       EXTERNAL USER-SUPPLIED ROUTINE TO EVALUATE */
/*                  THE ITERATION MATRIX (THIS IS OPTIONAL) */
/*                  THE CALL IS OF THE FORM */
/*                  CALL JAC(X,Y,YPRIME,PD,CJ,RPAR,IPAR) */
/*                  PD IS THE MATRIX OF PARTIAL DERIVATIVES, */
/*                  PD=DG/DY+CJ*DG/DYPRIME */
/*     H --         APPROPRIATE STEP SIZE FOR NEXT STEP. */
/*                  NORMALLY DETERMINED BY THE CODE */
/*     WT --        VECTOR OF WEIGHTS FOR ERROR CRITERION. */
/*     JSTART --    INTEGER VARIABLE SET 0 FOR */
/*                  FIRST STEP, 1 OTHERWISE. */
/*     IDID --      COMPLETION CODE WITH THE FOLLOWING MEANINGS: */
/*                  IDID= 1 -- THE STEP WAS COMPLETED SUCCESSFULLY */
/*                  IDID=-6 -- THE ERROR TEST FAILED REPEATEDLY */
/*                  IDID=-7 -- THE CORRECTOR COULD NOT CONVERGE */
/*                  IDID=-8 -- THE ITERATION MATRIX IS SINGULAR */
/*                  IDID=-9 -- THE CORRECTOR COULD NOT CONVERGE. */
/*                             THERE WERE REPEATED ERROR TEST */
/*                             FAILURES ON THIS STEP. */
/*                  IDID=-10-- THE CORRECTOR COULD NOT CONVERGE */
/*                             BECAUSE IRES WAS EQUAL TO MINUS ONE */
/*                  IDID=-11-- IRES EQUAL TO -2 WAS ENCOUNTERED, */
/*                             AND CONTROL IS BEING RETURNED TO */
/*                             THE CALLING PROGRAM */
/*     RPAR,IPAR -- REAL AND INTEGER PARAMETER ARRAYS THAT */
/*                  ARE USED FOR COMMUNICATION BETWEEN THE */
/*                  CALLING PROGRAM AND EXTERNAL USER ROUTINES */
/*                  THEY ARE NOT ALTERED BY SDASTP */
/*     PHI --       ARRAY OF DIVIDED DIFFERENCES USED BY */
/*                  SDASTP. THE LENGTH IS NEQ*(K+1),WHERE */
/*                  K IS THE MAXIMUM ORDER */
/*     DELTA,E --   WORK VECTORS FOR SDASTP OF LENGTH NEQ */
/*     WM,IWM --    REAL AND INTEGER ARRAYS STORING */
/*                  MATRIX INFORMATION SUCH AS THE MATRIX */
/*                  OF PARTIAL DERIVATIVES,PERMUTATION */
/*                  VECTOR, AND VARIOUS OTHER INFORMATION. */

/*     THE OTHER PARAMETERS ARE INFORMATION */
/*     WHICH IS NEEDED INTERNALLY BY SDASTP TO */
/*     CONTINUE FROM STEP TO STEP. */

/* ----------------------------------------------------------------------- */
/* ***ROUTINES CALLED  SDAJAC, SDANRM, SDASLV, SDATRP */
/* ***REVISION HISTORY  (YYMMDD) */
/*   830315  DATE WRITTEN */
/*   901009  Finished conversion to SLATEC 4.0 format (F.N.Fritsch) */
/*   901019  Merged changes made by C. Ulrich with SLATEC 4.0 format. */
/*   901026  Added explicit declarations for all variables and minor */
/*           cosmetic changes to prologue.  (FNF) */
/* ***END PROLOGUE  SDASTP */





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
    --alpha;
    --beta;
    --gamma;
    --psi;
    --sigma;

    /* Function Body */





/* ----------------------------------------------------------------------- */
/*     BLOCK 1. */
/*     INITIALIZE. ON THE FIRST CALL,SET */
/*     THE ORDER TO 1 AND INITIALIZE */
/*     OTHER VARIABLES. */
/* ----------------------------------------------------------------------- */

/*     INITIALIZATIONS FOR ALL CALLS */
/* ***FIRST EXECUTABLE STATEMENT  SDASTP */
    *idid = 1;
    xold = *x;
    ncf = 0;
    nsf = 0;
    nef = 0;
    if (*jstart != 0) {
	goto L120;
    }

/*     IF THIS IS THE FIRST STEP,PERFORM */
/*     OTHER INITIALIZATIONS */
    iwm[14] = 0;
    iwm[15] = 0;
    *k = 1;
    *kold = 0;
    *hold = 0.f;
    *jstart = 1;
    psi[1] = *h__;
    *cjold = 1.f / *h__;
    *cj = *cjold;
    *s = 100.f;
    *jcalc = -1;
    delnrm = 1.f;
    *iphase = 0;
    *ns = 0;
L120:





/* ----------------------------------------------------------------------- */
/*     BLOCK 2 */
/*     COMPUTE COEFFICIENTS OF FORMULAS FOR */
/*     THIS STEP. */
/* ----------------------------------------------------------------------- */
L200:
    kp1 = *k + 1;
    kp2 = *k + 2;
    km1 = *k - 1;
    xold = *x;
    if (*h__ != *hold || *k != *kold) {
	*ns = 0;
    }
/* Computing MIN */
    i__1 = *ns + 1, i__2 = *kold + 2;
    *ns = min(i__1,i__2);
    nsp1 = *ns + 1;
    if (kp1 < *ns) {
	goto L230;
    }

    beta[1] = 1.f;
    alpha[1] = 1.f;
    temp1 = *h__;
    gamma[1] = 0.f;
    sigma[1] = 1.f;
    i__1 = kp1;
    for (i__ = 2; i__ <= i__1; ++i__) {
	temp2 = psi[i__ - 1];
	psi[i__ - 1] = temp1;
	beta[i__] = beta[i__ - 1] * psi[i__ - 1] / temp2;
	temp1 = temp2 + *h__;
	alpha[i__] = *h__ / temp1;
	sigma[i__] = (i__ - 1) * sigma[i__ - 1] * alpha[i__];
	gamma[i__] = gamma[i__ - 1] + alpha[i__ - 1] / *h__;
/* L210: */
    }
    psi[kp1] = temp1;
L230:

/*     COMPUTE ALPHAS, ALPHA0 */
    alphas = 0.f;
    alpha0 = 0.f;
    i__1 = *k;
    for (i__ = 1; i__ <= i__1; ++i__) {
	alphas -= 1.f / i__;
	alpha0 -= alpha[i__];
/* L240: */
    }

/*     COMPUTE LEADING COEFFICIENT CJ */
    cjlast = *cj;
    *cj = -alphas / *h__;

/*     COMPUTE VARIABLE STEPSIZE ERROR COEFFICIENT CK */
    ck = (r__1 = alpha[kp1] + alphas - alpha0, dabs(r__1));
/* Computing MAX */
    r__1 = ck, r__2 = alpha[kp1];
    ck = dmax(r__1,r__2);

/*     DECIDE WHETHER NEW JACOBIAN IS NEEDED */
    temp1 = (1.f - xrate) / (xrate + 1.f);
    temp2 = 1.f / temp1;
    if (*cj / *cjold < temp1 || *cj / *cjold > temp2) {
	*jcalc = -1;
    }
    if (*cj != cjlast) {
	*s = 100.f;
    }

/*     CHANGE PHI TO PHI STAR */
    if (kp1 < nsp1) {
	goto L280;
    }
    i__1 = kp1;
    for (j = nsp1; j <= i__1; ++j) {
	i__2 = *neq;
	for (i__ = 1; i__ <= i__2; ++i__) {
/* L260: */
	    phi[i__ + j * phi_dim1] = beta[j] * phi[i__ + j * phi_dim1];
	}
/* L270: */
    }
L280:

/*     UPDATE TIME */
    *x += *h__;





/* ----------------------------------------------------------------------- */
/*     BLOCK 3 */
/*     PREDICT THE SOLUTION AND DERIVATIVE, */
/*     AND SOLVE THE CORRECTOR EQUATION */
/* ----------------------------------------------------------------------- */

/*     FIRST,PREDICT THE SOLUTION AND DERIVATIVE */
L300:
    i__1 = *neq;
    for (i__ = 1; i__ <= i__1; ++i__) {
	y[i__] = phi[i__ + phi_dim1];
/* L310: */
	yprime[i__] = 0.f;
    }
    i__1 = kp1;
    for (j = 2; j <= i__1; ++j) {
	i__2 = *neq;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    y[i__] += phi[i__ + j * phi_dim1];
/* L320: */
	    yprime[i__] += gamma[j] * phi[i__ + j * phi_dim1];
	}
/* L330: */
    }
    pnorm = sdanrm_(neq, &y[1], &wt[1], &rpar[1], &ipar[1]);



/*     SOLVE THE CORRECTOR EQUATION USING A */
/*     MODIFIED NEWTON SCHEME. */
    convgd = TRUE_;
    m = 0;
    ++iwm[12];
    ires = 0;
    (*res)(x, &y[1], &yprime[1], &delta[1], &ires, &rpar[1], &ipar[1]);
    if (ires < 0) {
	goto L380;
    }


/*     IF INDICATED,REEVALUATE THE */
/*     ITERATION MATRIX PD = DG/DY + CJ*DG/DYPRIME */
/*     (WHERE G(X,Y,YPRIME)=0). SET */
/*     JCALC TO 0 AS AN INDICATOR THAT */
/*     THIS HAS BEEN DONE. */
    if (*jcalc != -1) {
	goto L340;
    }
    ++iwm[13];
    *jcalc = 0;
    sdajac_(neq, x, &y[1], &yprime[1], &delta[1], cj, h__, &ier, &wt[1], &e[1]
	    , &wm[1], &iwm[1], (S_fp)res, &ires, uround, (U_fp)jac, &rpar[1], 
	    &ipar[1], ntemp);
    *cjold = *cj;
    *s = 100.f;
    if (ires < 0) {
	goto L380;
    }
    if (ier != 0) {
	goto L380;
    }
    nsf = 0;


/*     INITIALIZE THE ERROR ACCUMULATION VECTOR E. */
L340:
    i__1 = *neq;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L345: */
	e[i__] = 0.f;
    }


/*     CORRECTOR LOOP. */
L350:

/*     MULTIPLY RESIDUAL BY TEMP1 TO ACCELERATE CONVERGENCE */
    temp1 = 2.f / (*cj / *cjold + 1.f);
    i__1 = *neq;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L355: */
	delta[i__] *= temp1;
    }

/*     COMPUTE A NEW ITERATE (BACK-SUBSTITUTION). */
/*     STORE THE CORRECTION IN DELTA. */
    sdaslv_(neq, &delta[1], &wm[1], &iwm[1]);

/*     UPDATE Y, E, AND YPRIME */
    i__1 = *neq;
    for (i__ = 1; i__ <= i__1; ++i__) {
	y[i__] -= delta[i__];
	e[i__] -= delta[i__];
/* L360: */
	yprime[i__] -= *cj * delta[i__];
    }

/*     TEST FOR CONVERGENCE OF THE ITERATION */
    delnrm = sdanrm_(neq, &delta[1], &wt[1], &rpar[1], &ipar[1]);
    if (delnrm <= *uround * 100.f * pnorm) {
	goto L375;
    }
    if (m > 0) {
	goto L365;
    }
    oldnrm = delnrm;
    goto L367;
L365:
    d__1 = (doublereal) (delnrm / oldnrm);
    d__2 = (doublereal) (1.f / m);
    rate = pow_dd(&d__1, &d__2);
    if (rate > .9f) {
	goto L370;
    }
    *s = rate / (1.f - rate);
L367:
    if (*s * delnrm <= .33f) {
	goto L375;
    }

/*     THE CORRECTOR HAS NOT YET CONVERGED. */
/*     UPDATE M AND TEST WHETHER THE */
/*     MAXIMUM NUMBER OF ITERATIONS HAVE */
/*     BEEN TRIED. */
    ++m;
    if (m >= maxit) {
	goto L370;
    }

/*     EVALUATE THE RESIDUAL */
/*     AND GO BACK TO DO ANOTHER ITERATION */
    ++iwm[12];
    ires = 0;
    (*res)(x, &y[1], &yprime[1], &delta[1], &ires, &rpar[1], &ipar[1]);
    if (ires < 0) {
	goto L380;
    }
    goto L350;


/*     THE CORRECTOR FAILED TO CONVERGE IN MAXIT */
/*     ITERATIONS. IF THE ITERATION MATRIX */
/*     IS NOT CURRENT,RE-DO THE STEP WITH */
/*     A NEW ITERATION MATRIX. */
L370:
    if (*jcalc == 0) {
	goto L380;
    }
    *jcalc = -1;
    goto L300;


/*     THE ITERATION HAS CONVERGED.  IF NONNEGATIVITY OF SOLUTION IS */
/*     REQUIRED, SET THE SOLUTION NONNEGATIVE, IF THE PERTURBATION */
/*     TO DO IT IS SMALL ENOUGH.  IF THE CHANGE IS TOO LARGE, THEN */
/*     CONSIDER THE CORRECTOR ITERATION TO HAVE FAILED. */
L375:
    if (*nonneg == 0) {
	goto L390;
    }
    i__1 = *neq;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L377: */
/* Computing MIN */
	r__1 = y[i__];
	delta[i__] = dmin(r__1,0.f);
    }
    delnrm = sdanrm_(neq, &delta[1], &wt[1], &rpar[1], &ipar[1]);
    if (delnrm > .33f) {
	goto L380;
    }
    i__1 = *neq;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L378: */
	e[i__] -= delta[i__];
    }
    goto L390;


/*     EXITS FROM BLOCK 3 */
/*     NO CONVERGENCE WITH CURRENT ITERATION */
/*     MATRIX,OR SINGULAR ITERATION MATRIX */
L380:
    convgd = FALSE_;
L390:
    *jcalc = 1;
    if (! convgd) {
	goto L600;
    }





/* ----------------------------------------------------------------------- */
/*     BLOCK 4 */
/*     ESTIMATE THE ERRORS AT ORDERS K,K-1,K-2 */
/*     AS IF CONSTANT STEPSIZE WAS USED. ESTIMATE */
/*     THE LOCAL ERROR AT ORDER K AND TEST */
/*     WHETHER THE CURRENT STEP IS SUCCESSFUL. */
/* ----------------------------------------------------------------------- */

/*     ESTIMATE ERRORS AT ORDERS K,K-1,K-2 */
    enorm = sdanrm_(neq, &e[1], &wt[1], &rpar[1], &ipar[1]);
    erk = sigma[*k + 1] * enorm;
    terk = (*k + 1) * erk;
    est = erk;
    knew = *k;
    if (*k == 1) {
	goto L430;
    }
    i__1 = *neq;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L405: */
	delta[i__] = phi[i__ + kp1 * phi_dim1] + e[i__];
    }
    erkm1 = sigma[*k] * sdanrm_(neq, &delta[1], &wt[1], &rpar[1], &ipar[1]);
    terkm1 = *k * erkm1;
    if (*k > 2) {
	goto L410;
    }
    if (terkm1 <= terk * .5f) {
	goto L420;
    }
    goto L430;
L410:
    i__1 = *neq;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L415: */
	delta[i__] = phi[i__ + *k * phi_dim1] + delta[i__];
    }
    erkm2 = sigma[*k - 1] * sdanrm_(neq, &delta[1], &wt[1], &rpar[1], &ipar[1]
	    );
    terkm2 = (*k - 1) * erkm2;
    if (dmax(terkm1,terkm2) > terk) {
	goto L430;
    }
/*     LOWER THE ORDER */
L420:
    knew = *k - 1;
    est = erkm1;


/*     CALCULATE THE LOCAL ERROR FOR THE CURRENT STEP */
/*     TO SEE IF THE STEP WAS SUCCESSFUL */
L430:
    err = ck * enorm;
    if (err > 1.f) {
	goto L600;
    }





/* ----------------------------------------------------------------------- */
/*     BLOCK 5 */
/*     THE STEP IS SUCCESSFUL. DETERMINE */
/*     THE BEST ORDER AND STEPSIZE FOR */
/*     THE NEXT STEP. UPDATE THE DIFFERENCES */
/*     FOR THE NEXT STEP. */
/* ----------------------------------------------------------------------- */
    *idid = 1;
    ++iwm[11];
    kdiff = *k - *kold;
    *kold = *k;
    *hold = *h__;


/*     ESTIMATE THE ERROR AT ORDER K+1 UNLESS: */
/*        ALREADY DECIDED TO LOWER ORDER, OR */
/*        ALREADY USING MAXIMUM ORDER, OR */
/*        STEPSIZE NOT CONSTANT, OR */
/*        ORDER RAISED IN PREVIOUS STEP */
    if (knew == km1 || *k == iwm[3]) {
	*iphase = 1;
    }
    if (*iphase == 0) {
	goto L545;
    }
    if (knew == km1) {
	goto L540;
    }
    if (*k == iwm[3]) {
	goto L550;
    }
    if (kp1 >= *ns || kdiff == 1) {
	goto L550;
    }
    i__1 = *neq;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L510: */
	delta[i__] = e[i__] - phi[i__ + kp2 * phi_dim1];
    }
    erkp1 = 1.f / (*k + 2) * sdanrm_(neq, &delta[1], &wt[1], &rpar[1], &ipar[
	    1]);
    terkp1 = (*k + 2) * erkp1;
    if (*k > 1) {
	goto L520;
    }
    if (terkp1 >= terk * .5f) {
	goto L550;
    }
    goto L530;
L520:
    if (terkm1 <= dmin(terk,terkp1)) {
	goto L540;
    }
    if (terkp1 >= terk || *k == iwm[3]) {
	goto L550;
    }

/*     RAISE ORDER */
L530:
    *k = kp1;
    est = erkp1;
    goto L550;

/*     LOWER ORDER */
L540:
    *k = km1;
    est = erkm1;
    goto L550;

/*     IF IPHASE = 0, INCREASE ORDER BY ONE AND MULTIPLY STEPSIZE BY */
/*     FACTOR TWO */
L545:
    *k = kp1;
    hnew = *h__ * 2.f;
    *h__ = hnew;
    goto L575;


/*     DETERMINE THE APPROPRIATE STEPSIZE FOR */
/*     THE NEXT STEP. */
L550:
    hnew = *h__;
    temp2 = (real) (*k + 1);
    d__1 = (doublereal) (est * 2.f + 1e-4f);
    d__2 = (doublereal) (-1.f / temp2);
    r__ = pow_dd(&d__1, &d__2);
    if (r__ < 2.f) {
	goto L555;
    }
    hnew = *h__ * 2.f;
    goto L560;
L555:
    if (r__ > 1.f) {
	goto L560;
    }
/* Computing MAX */
    r__1 = .5f, r__2 = dmin(.9f,r__);
    r__ = dmax(r__1,r__2);
    hnew = *h__ * r__;
L560:
    *h__ = hnew;


/*     UPDATE DIFFERENCES FOR NEXT STEP */
L575:
    if (*kold == iwm[3]) {
	goto L585;
    }
    i__1 = *neq;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L580: */
	phi[i__ + kp2 * phi_dim1] = e[i__];
    }
L585:
    i__1 = *neq;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L590: */
	phi[i__ + kp1 * phi_dim1] += e[i__];
    }
    i__1 = kp1;
    for (j1 = 2; j1 <= i__1; ++j1) {
	j = kp1 - j1 + 1;
	i__2 = *neq;
	for (i__ = 1; i__ <= i__2; ++i__) {
/* L595: */
	    phi[i__ + j * phi_dim1] += phi[i__ + (j + 1) * phi_dim1];
	}
    }
    return 0;





/* ----------------------------------------------------------------------- */
/*     BLOCK 6 */
/*     THE STEP IS UNSUCCESSFUL. RESTORE X,PSI,PHI */
/*     DETERMINE APPROPRIATE STEPSIZE FOR */
/*     CONTINUING THE INTEGRATION, OR EXIT WITH */
/*     AN ERROR FLAG IF THERE HAVE BEEN MANY */
/*     FAILURES. */
/* ----------------------------------------------------------------------- */
L600:
    *iphase = 1;

/*     RESTORE X,PHI,PSI */
    *x = xold;
    if (kp1 < nsp1) {
	goto L630;
    }
    i__2 = kp1;
    for (j = nsp1; j <= i__2; ++j) {
	temp1 = 1.f / beta[j];
	i__1 = *neq;
	for (i__ = 1; i__ <= i__1; ++i__) {
/* L610: */
	    phi[i__ + j * phi_dim1] = temp1 * phi[i__ + j * phi_dim1];
	}
/* L620: */
    }
L630:
    i__2 = kp1;
    for (i__ = 2; i__ <= i__2; ++i__) {
/* L640: */
	psi[i__ - 1] = psi[i__] - *h__;
    }


/*     TEST WHETHER FAILURE IS DUE TO CORRECTOR ITERATION */
/*     OR ERROR TEST */
    if (convgd) {
	goto L660;
    }
    ++iwm[15];


/*     THE NEWTON ITERATION FAILED TO CONVERGE WITH */
/*     A CURRENT ITERATION MATRIX.  DETERMINE THE CAUSE */
/*     OF THE FAILURE AND TAKE APPROPRIATE ACTION. */
    if (ier == 0) {
	goto L650;
    }

/*     THE ITERATION MATRIX IS SINGULAR. REDUCE */
/*     THE STEPSIZE BY A FACTOR OF 4. IF */
/*     THIS HAPPENS THREE TIMES IN A ROW ON */
/*     THE SAME STEP, RETURN WITH AN ERROR FLAG */
    ++nsf;
    r__ = .25f;
    *h__ *= r__;
    if (nsf < 3 && dabs(*h__) >= *hmin) {
	goto L690;
    }
    *idid = -8;
    goto L675;


/*     THE NEWTON ITERATION FAILED TO CONVERGE FOR A REASON */
/*     OTHER THAN A SINGULAR ITERATION MATRIX.  IF IRES = -2, THEN */
/*     RETURN.  OTHERWISE, REDUCE THE STEPSIZE AND TRY AGAIN, UNLESS */
/*     TOO MANY FAILURES HAVE OCCURRED. */
L650:
    if (ires > -2) {
	goto L655;
    }
    *idid = -11;
    goto L675;
L655:
    ++ncf;
    r__ = .25f;
    *h__ *= r__;
    if (ncf < 10 && dabs(*h__) >= *hmin) {
	goto L690;
    }
    *idid = -7;
    if (ires < 0) {
	*idid = -10;
    }
    if (nef >= 3) {
	*idid = -9;
    }
    goto L675;


/*     THE NEWTON SCHEME CONVERGED, AND THE CAUSE */
/*     OF THE FAILURE WAS THE ERROR ESTIMATE */
/*     EXCEEDING THE TOLERANCE. */
L660:
    ++nef;
    ++iwm[14];
    if (nef > 1) {
	goto L665;
    }

/*     ON FIRST ERROR TEST FAILURE, KEEP CURRENT ORDER OR LOWER */
/*     ORDER BY ONE.  COMPUTE NEW STEPSIZE BASED ON DIFFERENCES */
/*     OF THE SOLUTION. */
    *k = knew;
    temp2 = (real) (*k + 1);
    d__1 = (doublereal) (est * 2.f + 1e-4f);
    d__2 = (doublereal) (-1.f / temp2);
    r__ = pow_dd(&d__1, &d__2) * .9f;
/* Computing MAX */
    r__1 = .25f, r__2 = dmin(.9f,r__);
    r__ = dmax(r__1,r__2);
    *h__ *= r__;
    if (dabs(*h__) >= *hmin) {
	goto L690;
    }
    *idid = -6;
    goto L675;

/*     ON SECOND ERROR TEST FAILURE, USE THE CURRENT ORDER OR */
/*     DECREASE ORDER BY ONE.  REDUCE THE STEPSIZE BY A FACTOR OF */
/*     FOUR. */
L665:
    if (nef > 2) {
	goto L670;
    }
    *k = knew;
    *h__ *= .25f;
    if (dabs(*h__) >= *hmin) {
	goto L690;
    }
    *idid = -6;
    goto L675;

/*     ON THIRD AND SUBSEQUENT ERROR TEST FAILURES, SET THE ORDER TO */
/*     ONE AND REDUCE THE STEPSIZE BY A FACTOR OF FOUR. */
L670:
    *k = 1;
    *h__ *= .25f;
    if (dabs(*h__) >= *hmin) {
	goto L690;
    }
    *idid = -6;
    goto L675;




/*     FOR ALL CRASHES, RESTORE Y TO ITS LAST VALUE, */
/*     INTERPOLATE TO FIND YPRIME AT LAST X, AND RETURN */
L675:
    sdatrp_(x, x, &y[1], &yprime[1], neq, k, &phi[phi_offset], &psi[1]);
    return 0;


/*     GO BACK AND TRY THIS STEP AGAIN */
L690:
    goto L200;

/* ------END OF SUBROUTINE SDASTP------ */
} /* sdastp_ */

