/* stod.f -- translated by f2c (version 12.02.01).
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
    real rownd, conit, crate, el[13], elco[156]	/* was [13][12] */, hold, rc, 
	    rmax, tesco[36]	/* was [3][12] */, el0, h__, hmin, hmxi, hu, 
	    tn, uround;
    integer iownd[7], ksteps, iod[6], ialth, ipup, lmax, meo, nqnyh, nstepj, 
	    ier, jstart, kflag, l, meth, miter, maxord, n, nq, nst, nfe, nje, 
	    nqu;
} debdf1_;

#define debdf1_1 debdf1_

/* DECK STOD */
/* Subroutine */ int stod_(integer *neq, real *y, real *yh, integer *nyh, 
	real *yh1, real *ewt, real *savf, real *acor, real *wm, integer *iwm, 
	S_fp f, U_fp jac, real *rpar, integer *ipar)
{
    /* System generated locals */
    integer yh_dim1, yh_offset, i__1, i__2;
    real r__1, r__2, r__3;
    doublereal d__1, d__2;

    /* Local variables */
    static integer i__, j, m;
    static real r__;
    static integer i1, jb;
    static real rh, del, ddn;
    static integer ncf;
    static real dsm, dup;
    extern /* Subroutine */ int cfod_(integer *, real *, real *), pjac_(
	    integer *, real *, real *, integer *, real *, real *, real *, 
	    real *, integer *, S_fp, U_fp, real *, integer *);
    static real dcon, delp, exdn, rhdn, told;
    static integer iret, newq;
    static real rhsm, exsm, rhup, exup;
    extern /* Subroutine */ int slvs_(real *, integer *, real *, real *);
    static integer iredo;
    extern doublereal vnwrms_(integer *, real *, real *);

/* ***BEGIN PROLOGUE  STOD */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to DEBDF */
/* ***LIBRARY   SLATEC */
/* ***TYPE      SINGLE PRECISION (STOD-S, DSTOD-D) */
/* ***AUTHOR  Watts, H. A., (SNLA) */
/* ***DESCRIPTION */

/*   STOD integrates a system of first order odes over one step in the */
/*   integrator package DEBDF. */
/* ---------------------------------------------------------------------- */
/* STOD  performs one step of the integration of an initial value */
/* problem for a system of ordinary differential equations. */
/* Note.. STOD  is independent of the value of the iteration method */
/* indicator MITER, when this is .NE. 0, and hence is independent */
/* of the type of chord method used, or the Jacobian structure. */
/* Communication with STOD  is done with the following variables.. */

/* Y      = An array of length .GE. n used as the Y argument in */
/*          all calls to F and JAC. */
/* NEQ    = Integer array containing problem size in NEQ(1), and */
/*          passed as the NEQ argument in all calls to F and JAC. */
/* YH     = An NYH by LMAX array containing the dependent variables */
/*          and their approximate scaled derivatives, where */
/*          LMAX = MAXORD + 1.  YH(I,J+1) contains the approximate */
/*          J-th derivative of Y(I), scaled by H**J/Factorial(j) */
/*          (J = 0,1,...,NQ).  On entry for the first step, the first */
/*          two columns of YH must be set from the initial values. */
/* NYH    = A constant integer .GE. N, the first dimension of YH. */
/* YH1    = A one-dimensional array occupying the same space as YH. */
/* EWT    = An array of N elements with which the estimated local */
/*          errors in YH are compared. */
/* SAVF   = An array of working storage, of length N. */
/* ACOR   = A work array of length N, used for the accumulated */
/*          corrections.  On a successful return, ACOR(I) contains */
/*          the estimated one-step local error in Y(I). */
/* WM,IWM = Real and integer work arrays associated with matrix */
/*          operations in chord iteration (MITER .NE. 0). */
/* PJAC   = Name of routine to evaluate and preprocess Jacobian matrix */
/*          if a chord method is being used. */
/* SLVS   = Name of routine to solve linear system in chord iteration. */
/* H      = The step size to be attempted on the next step. */
/*          H is altered by the error control algorithm during the */
/*          problem.  H can be either positive or negative, but its */
/*          sign must remain constant throughout the problem. */
/* HMIN   = The minimum absolute value of the step size H to be used. */
/* HMXI   = Inverse of the maximum absolute value of H to be used. */
/*          HMXI = 0.0 is allowed and corresponds to an infinite HMAX. */
/*          HMIN and HMXI may be changed at any time, but will not */
/*          take effect until the next change of H is considered. */
/* TN     = The independent variable. TN is updated on each step taken. */
/* JSTART = An integer used for input only, with the following */
/*          values and meanings.. */
/*               0  Perform the first step. */
/*           .GT.0  Take a new step continuing from the last. */
/*              -1  Take the next step with a new value of H, MAXORD, */
/*                    N, METH, MITER, and/or matrix parameters. */
/*              -2  Take the next step with a new value of H, */
/*                    but with other inputs unchanged. */
/*          On return, JSTART is set to 1 to facilitate continuation. */
/* KFLAG  = a completion code with the following meanings.. */
/*               0  The step was successful. */
/*              -1  The requested error could not be achieved. */
/*              -2  Corrector convergence could not be achieved. */
/*          A return with KFLAG = -1 or -2 means either */
/*          ABS(H) = HMIN or 10 consecutive failures occurred. */
/*          On a return with KFLAG negative, the values of TN and */
/*          the YH array are as of the beginning of the last */
/*          step, and H is the last step size attempted. */
/* MAXORD = The maximum order of integration method to be allowed. */
/* METH/MITER = The method flags.  See description in driver. */
/* N      = The number of first-order differential equations. */
/* ---------------------------------------------------------------------- */

/* ***SEE ALSO  DEBDF */
/* ***ROUTINES CALLED  CFOD, PJAC, SLVS, VNWRMS */
/* ***COMMON BLOCKS    DEBDF1 */
/* ***REVISION HISTORY  (YYMMDD) */
/*   800901  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900328  Added TYPE section.  (WRB) */
/*   910722  Updated AUTHOR section.  (ALS) */
/*   920422  Changed DIMENSION statement.  (WRB) */
/* ***END PROLOGUE  STOD */

/* LLL. OPTIMIZE */


/* ***FIRST EXECUTABLE STATEMENT  STOD */
    /* Parameter adjustments */
    --y;
    yh_dim1 = *nyh;
    yh_offset = 1 + yh_dim1;
    yh -= yh_offset;
    --yh1;
    --ewt;
    --savf;
    --acor;
    --wm;
    --iwm;
    --rpar;
    --ipar;

    /* Function Body */
    debdf1_1.kflag = 0;
    told = debdf1_1.tn;
    ncf = 0;
    if (debdf1_1.jstart > 0) {
	goto L200;
    }
    if (debdf1_1.jstart == -1) {
	goto L100;
    }
    if (debdf1_1.jstart == -2) {
	goto L160;
    }
/* ----------------------------------------------------------------------- */
/* ON THE FIRST CALL, THE ORDER IS SET TO 1, AND OTHER VARIABLES ARE */
/* INITIALIZED.  RMAX IS THE MAXIMUM RATIO BY WHICH H CAN BE INCREASED */
/* IN A SINGLE STEP.  IT IS INITIALLY 1.E4 TO COMPENSATE FOR THE SMALL */
/* INITIAL H, BUT THEN IS NORMALLY EQUAL TO 10.  IF A FAILURE */
/* OCCURS (IN CORRECTOR CONVERGENCE OR ERROR TEST), RMAX IS SET AT 2 */
/* FOR THE NEXT INCREASE. */
/* ----------------------------------------------------------------------- */
    debdf1_1.lmax = debdf1_1.maxord + 1;
    debdf1_1.nq = 1;
    debdf1_1.l = 2;
    debdf1_1.ialth = 2;
    debdf1_1.rmax = 1e4f;
    debdf1_1.rc = 0.f;
    debdf1_1.el0 = 1.f;
    debdf1_1.crate = .7f;
    delp = 0.f;
    debdf1_1.hold = debdf1_1.h__;
    debdf1_1.meo = debdf1_1.meth;
    debdf1_1.nstepj = 0;
    iret = 3;
    goto L140;
/* ----------------------------------------------------------------------- */
/* THE FOLLOWING BLOCK HANDLES PRELIMINARIES NEEDED WHEN JSTART = -1. */
/* IPUP IS SET TO MITER TO FORCE A MATRIX UPDATE. */
/* IF AN ORDER INCREASE IS ABOUT TO BE CONSIDERED (IALTH = 1), */
/* IALTH IS RESET TO 2 TO POSTPONE CONSIDERATION ONE MORE STEP. */
/* IF THE CALLER HAS CHANGED METH, CFOD  IS CALLED TO RESET */
/* THE COEFFICIENTS OF THE METHOD. */
/* IF THE CALLER HAS CHANGED MAXORD TO A VALUE LESS THAN THE CURRENT */
/* ORDER NQ, NQ IS REDUCED TO MAXORD, AND A NEW H CHOSEN ACCORDINGLY. */
/* IF H IS TO BE CHANGED, YH MUST BE RESCALED. */
/* IF H OR METH IS BEING CHANGED, IALTH IS RESET TO L = NQ + 1 */
/* TO PREVENT FURTHER CHANGES IN H FOR THAT MANY STEPS. */
/* ----------------------------------------------------------------------- */
L100:
    debdf1_1.ipup = debdf1_1.miter;
    debdf1_1.lmax = debdf1_1.maxord + 1;
    if (debdf1_1.ialth == 1) {
	debdf1_1.ialth = 2;
    }
    if (debdf1_1.meth == debdf1_1.meo) {
	goto L110;
    }
    cfod_(&debdf1_1.meth, debdf1_1.elco, debdf1_1.tesco);
    debdf1_1.meo = debdf1_1.meth;
    if (debdf1_1.nq > debdf1_1.maxord) {
	goto L120;
    }
    debdf1_1.ialth = debdf1_1.l;
    iret = 1;
    goto L150;
L110:
    if (debdf1_1.nq <= debdf1_1.maxord) {
	goto L160;
    }
L120:
    debdf1_1.nq = debdf1_1.maxord;
    debdf1_1.l = debdf1_1.lmax;
    i__1 = debdf1_1.l;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L125: */
	debdf1_1.el[i__ - 1] = debdf1_1.elco[i__ + debdf1_1.nq * 13 - 14];
    }
    debdf1_1.nqnyh = debdf1_1.nq * *nyh;
    debdf1_1.rc = debdf1_1.rc * debdf1_1.el[0] / debdf1_1.el0;
    debdf1_1.el0 = debdf1_1.el[0];
    debdf1_1.conit = .5f / (debdf1_1.nq + 2);
    ddn = vnwrms_(&debdf1_1.n, &savf[1], &ewt[1]) / debdf1_1.tesco[debdf1_1.l 
	    * 3 - 3];
    exdn = 1.f / debdf1_1.l;
    d__1 = (doublereal) ddn;
    d__2 = (doublereal) exdn;
    rhdn = 1.f / (pow_dd(&d__1, &d__2) * 1.3f + 1.3e-6f);
    rh = dmin(rhdn,1.f);
    iredo = 3;
    if (debdf1_1.h__ == debdf1_1.hold) {
	goto L170;
    }
/* Computing MIN */
    r__2 = rh, r__3 = (r__1 = debdf1_1.h__ / debdf1_1.hold, dabs(r__1));
    rh = dmin(r__2,r__3);
    debdf1_1.h__ = debdf1_1.hold;
    goto L175;
/* ----------------------------------------------------------------------- */
/* CFOD  IS CALLED TO GET ALL THE INTEGRATION COEFFICIENTS FOR THE */
/* CURRENT METH.  THEN THE EL VECTOR AND RELATED CONSTANTS ARE RESET */
/* WHENEVER THE ORDER NQ IS CHANGED, OR AT THE START OF THE PROBLEM. */
/* ----------------------------------------------------------------------- */
L140:
    cfod_(&debdf1_1.meth, debdf1_1.elco, debdf1_1.tesco);
L150:
    i__1 = debdf1_1.l;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L155: */
	debdf1_1.el[i__ - 1] = debdf1_1.elco[i__ + debdf1_1.nq * 13 - 14];
    }
    debdf1_1.nqnyh = debdf1_1.nq * *nyh;
    debdf1_1.rc = debdf1_1.rc * debdf1_1.el[0] / debdf1_1.el0;
    debdf1_1.el0 = debdf1_1.el[0];
    debdf1_1.conit = .5f / (debdf1_1.nq + 2);
    switch (iret) {
	case 1:  goto L160;
	case 2:  goto L170;
	case 3:  goto L200;
    }
/* ----------------------------------------------------------------------- */
/* IF H IS BEING CHANGED, THE H RATIO RH IS CHECKED AGAINST */
/* RMAX, HMIN, AND HMXI, AND THE YH ARRAY RESCALED.  IALTH IS SET TO */
/* L = NQ + 1 TO PREVENT A CHANGE OF H FOR THAT MANY STEPS, UNLESS */
/* FORCED BY A CONVERGENCE OR ERROR TEST FAILURE. */
/* ----------------------------------------------------------------------- */
L160:
    if (debdf1_1.h__ == debdf1_1.hold) {
	goto L200;
    }
    rh = debdf1_1.h__ / debdf1_1.hold;
    debdf1_1.h__ = debdf1_1.hold;
    iredo = 3;
    goto L175;
L170:
/* Computing MAX */
    r__1 = rh, r__2 = debdf1_1.hmin / dabs(debdf1_1.h__);
    rh = dmax(r__1,r__2);
L175:
    rh = dmin(rh,debdf1_1.rmax);
/* Computing MAX */
    r__1 = 1.f, r__2 = dabs(debdf1_1.h__) * debdf1_1.hmxi * rh;
    rh /= dmax(r__1,r__2);
    r__ = 1.f;
    i__1 = debdf1_1.l;
    for (j = 2; j <= i__1; ++j) {
	r__ *= rh;
	i__2 = debdf1_1.n;
	for (i__ = 1; i__ <= i__2; ++i__) {
/* L180: */
	    yh[i__ + j * yh_dim1] *= r__;
	}
    }
    debdf1_1.h__ *= rh;
    debdf1_1.rc *= rh;
    debdf1_1.ialth = debdf1_1.l;
    if (iredo == 0) {
	goto L680;
    }
/* ----------------------------------------------------------------------- */
/* THIS SECTION COMPUTES THE PREDICTED VALUES BY EFFECTIVELY */
/* MULTIPLYING THE YH ARRAY BY THE PASCAL TRIANGLE MATRIX. */
/* RC IS THE RATIO OF NEW TO OLD VALUES OF THE COEFFICIENT  H*EL(1). */
/* WHEN RC DIFFERS FROM 1 BY MORE THAN 30 PERCENT, IPUP IS SET TO MITER */
/* TO FORCE PJAC TO BE CALLED, IF A JACOBIAN IS INVOLVED. */
/* IN ANY CASE, PJAC IS CALLED AT LEAST EVERY 20-TH STEP. */
/* ----------------------------------------------------------------------- */
L200:
    if ((r__1 = debdf1_1.rc - 1.f, dabs(r__1)) > .3f) {
	debdf1_1.ipup = debdf1_1.miter;
    }
    if (debdf1_1.nst >= debdf1_1.nstepj + 20) {
	debdf1_1.ipup = debdf1_1.miter;
    }
    debdf1_1.tn += debdf1_1.h__;
    i1 = debdf1_1.nqnyh + 1;
    i__2 = debdf1_1.nq;
    for (jb = 1; jb <= i__2; ++jb) {
	i1 -= *nyh;
	i__1 = debdf1_1.nqnyh;
	for (i__ = i1; i__ <= i__1; ++i__) {
/* L210: */
	    yh1[i__] += yh1[i__ + *nyh];
	}
/* L215: */
    }
    ++debdf1_1.ksteps;
/* ----------------------------------------------------------------------- */
/* UP TO 3 CORRECTOR ITERATIONS ARE TAKEN.  A CONVERGENCE TEST IS */
/* MADE ON THE R.M.S. NORM OF EACH CORRECTION, WEIGHTED BY THE ERROR */
/* WEIGHT VECTOR EWT.  THE SUM OF THE CORRECTIONS IS ACCUMULATED IN THE */
/* VECTOR ACOR(I).  THE YH ARRAY IS NOT ALTERED IN THE CORRECTOR LOOP. */
/* ----------------------------------------------------------------------- */
L220:
    m = 0;
    i__2 = debdf1_1.n;
    for (i__ = 1; i__ <= i__2; ++i__) {
/* L230: */
	y[i__] = yh[i__ + yh_dim1];
    }
    (*f)(&debdf1_1.tn, &y[1], &savf[1], &rpar[1], &ipar[1]);
    ++debdf1_1.nfe;
    if (debdf1_1.ipup <= 0) {
	goto L250;
    }
/* ----------------------------------------------------------------------- */
/* IF INDICATED, THE MATRIX P = I - H*EL(1)*J IS REEVALUATED AND */
/* PREPROCESSED BEFORE STARTING THE CORRECTOR ITERATION.  IPUP IS SET */
/* TO 0 AS AN INDICATOR THAT THIS HAS BEEN DONE. */
/* ----------------------------------------------------------------------- */
    debdf1_1.ipup = 0;
    debdf1_1.rc = 1.f;
    debdf1_1.nstepj = debdf1_1.nst;
    debdf1_1.crate = .7f;
    pjac_(neq, &y[1], &yh[yh_offset], nyh, &ewt[1], &acor[1], &savf[1], &wm[1]
	    , &iwm[1], (S_fp)f, (U_fp)jac, &rpar[1], &ipar[1]);
    if (debdf1_1.ier != 0) {
	goto L430;
    }
L250:
    i__2 = debdf1_1.n;
    for (i__ = 1; i__ <= i__2; ++i__) {
/* L260: */
	acor[i__] = 0.f;
    }
L270:
    if (debdf1_1.miter != 0) {
	goto L350;
    }
/* ----------------------------------------------------------------------- */
/* IN THE CASE OF FUNCTIONAL ITERATION, UPDATE Y DIRECTLY FROM */
/* THE RESULT OF THE LAST FUNCTION EVALUATION. */
/* ----------------------------------------------------------------------- */
    i__2 = debdf1_1.n;
    for (i__ = 1; i__ <= i__2; ++i__) {
	savf[i__] = debdf1_1.h__ * savf[i__] - yh[i__ + (yh_dim1 << 1)];
/* L290: */
	y[i__] = savf[i__] - acor[i__];
    }
    del = vnwrms_(&debdf1_1.n, &y[1], &ewt[1]);
    i__2 = debdf1_1.n;
    for (i__ = 1; i__ <= i__2; ++i__) {
	y[i__] = yh[i__ + yh_dim1] + debdf1_1.el[0] * savf[i__];
/* L300: */
	acor[i__] = savf[i__];
    }
    goto L400;
/* ----------------------------------------------------------------------- */
/* IN THE CASE OF THE CHORD METHOD, COMPUTE THE CORRECTOR ERROR, */
/* AND SOLVE THE LINEAR SYSTEM WITH THAT AS RIGHT-HAND SIDE AND */
/* P AS COEFFICIENT MATRIX. */
/* ----------------------------------------------------------------------- */
L350:
    i__2 = debdf1_1.n;
    for (i__ = 1; i__ <= i__2; ++i__) {
/* L360: */
	y[i__] = debdf1_1.h__ * savf[i__] - (yh[i__ + (yh_dim1 << 1)] + acor[
		i__]);
    }
    slvs_(&wm[1], &iwm[1], &y[1], &savf[1]);
    if (debdf1_1.ier != 0) {
	goto L410;
    }
    del = vnwrms_(&debdf1_1.n, &y[1], &ewt[1]);
    i__2 = debdf1_1.n;
    for (i__ = 1; i__ <= i__2; ++i__) {
	acor[i__] += y[i__];
/* L380: */
	y[i__] = yh[i__ + yh_dim1] + debdf1_1.el[0] * acor[i__];
    }
/* ----------------------------------------------------------------------- */
/* TEST FOR CONVERGENCE.  IF M.GT.0, AN ESTIMATE OF THE CONVERGENCE */
/* RATE CONSTANT IS STORED IN CRATE, AND THIS IS USED IN THE TEST. */
/* ----------------------------------------------------------------------- */
L400:
    if (m != 0) {
/* Computing MAX */
	r__1 = debdf1_1.crate * .2f, r__2 = del / delp;
	debdf1_1.crate = dmax(r__1,r__2);
    }
/* Computing MIN */
    r__1 = 1.f, r__2 = debdf1_1.crate * 1.5f;
    dcon = del * dmin(r__1,r__2) / (debdf1_1.tesco[debdf1_1.nq * 3 - 2] * 
	    debdf1_1.conit);
    if (dcon <= 1.f) {
	goto L450;
    }
    ++m;
    if (m == 3) {
	goto L410;
    }
    if (m >= 2 && del > delp * 2.f) {
	goto L410;
    }
    delp = del;
    (*f)(&debdf1_1.tn, &y[1], &savf[1], &rpar[1], &ipar[1]);
    ++debdf1_1.nfe;
    goto L270;
/* ----------------------------------------------------------------------- */
/* THE CORRECTOR ITERATION FAILED TO CONVERGE IN 3 TRIES. */
/* IF MITER .NE. 0 AND THE JACOBIAN IS OUT OF DATE, PJAC IS CALLED FOR */
/* THE NEXT TRY.  OTHERWISE THE YH ARRAY IS RETRACTED TO ITS VALUES */
/* BEFORE PREDICTION, AND H IS REDUCED, IF POSSIBLE.  IF H CANNOT BE */
/* REDUCED OR 10 FAILURES HAVE OCCURRED, EXIT WITH KFLAG = -2. */
/* ----------------------------------------------------------------------- */
L410:
    if (debdf1_1.ipup == 0) {
	goto L430;
    }
    debdf1_1.ipup = debdf1_1.miter;
    goto L220;
L430:
    debdf1_1.tn = told;
    ++ncf;
    debdf1_1.rmax = 2.f;
    i1 = debdf1_1.nqnyh + 1;
    i__2 = debdf1_1.nq;
    for (jb = 1; jb <= i__2; ++jb) {
	i1 -= *nyh;
	i__1 = debdf1_1.nqnyh;
	for (i__ = i1; i__ <= i__1; ++i__) {
/* L440: */
	    yh1[i__] -= yh1[i__ + *nyh];
	}
/* L445: */
    }
    if (dabs(debdf1_1.h__) <= debdf1_1.hmin * 1.00001f) {
	goto L670;
    }
    if (ncf == 10) {
	goto L670;
    }
    rh = .25f;
    debdf1_1.ipup = debdf1_1.miter;
    iredo = 1;
    goto L170;
/* ----------------------------------------------------------------------- */
/* THE CORRECTOR HAS CONVERGED.  IPUP IS SET TO -1 IF MITER .NE. 0, */
/* TO SIGNAL THAT THE JACOBIAN INVOLVED MAY NEED UPDATING LATER. */
/* THE LOCAL ERROR TEST IS MADE AND CONTROL PASSES TO STATEMENT 500 */
/* IF IT FAILS. */
/* ----------------------------------------------------------------------- */
L450:
    if (debdf1_1.miter != 0) {
	debdf1_1.ipup = -1;
    }
    if (m == 0) {
	dsm = del / debdf1_1.tesco[debdf1_1.nq * 3 - 2];
    }
    if (m > 0) {
	dsm = vnwrms_(&debdf1_1.n, &acor[1], &ewt[1]) / debdf1_1.tesco[
		debdf1_1.nq * 3 - 2];
    }
    if (dsm > 1.f) {
	goto L500;
    }
/* ----------------------------------------------------------------------- */
/* AFTER A SUCCESSFUL STEP, UPDATE THE YH ARRAY. */
/* CONSIDER CHANGING H IF IALTH = 1.  OTHERWISE DECREASE IALTH BY 1. */
/* IF IALTH IS THEN 1 AND NQ .LT. MAXORD, THEN ACOR IS SAVED FOR */
/* USE IN A POSSIBLE ORDER INCREASE ON THE NEXT STEP. */
/* IF A CHANGE IN H IS CONSIDERED, AN INCREASE OR DECREASE IN ORDER */
/* BY ONE IS CONSIDERED ALSO.  A CHANGE IN H IS MADE ONLY IF IT IS BY A */
/* FACTOR OF AT LEAST 1.1.  IF NOT, IALTH IS SET TO 3 TO PREVENT */
/* TESTING FOR THAT MANY STEPS. */
/* ----------------------------------------------------------------------- */
    debdf1_1.kflag = 0;
    iredo = 0;
    ++debdf1_1.nst;
    debdf1_1.hu = debdf1_1.h__;
    debdf1_1.nqu = debdf1_1.nq;
    i__2 = debdf1_1.l;
    for (j = 1; j <= i__2; ++j) {
	i__1 = debdf1_1.n;
	for (i__ = 1; i__ <= i__1; ++i__) {
/* L470: */
	    yh[i__ + j * yh_dim1] += debdf1_1.el[j - 1] * acor[i__];
	}
    }
    --debdf1_1.ialth;
    if (debdf1_1.ialth == 0) {
	goto L520;
    }
    if (debdf1_1.ialth > 1) {
	goto L690;
    }
    if (debdf1_1.l == debdf1_1.lmax) {
	goto L690;
    }
    i__1 = debdf1_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L490: */
	yh[i__ + debdf1_1.lmax * yh_dim1] = acor[i__];
    }
    goto L690;
/* ----------------------------------------------------------------------- */
/* THE ERROR TEST FAILED.  KFLAG KEEPS TRACK OF MULTIPLE FAILURES. */
/* RESTORE TN AND THE YH ARRAY TO THEIR PREVIOUS VALUES, AND PREPARE */
/* TO TRY THE STEP AGAIN.  COMPUTE THE OPTIMUM STEP SIZE FOR THIS OR */
/* ONE LOWER ORDER.  AFTER 2 OR MORE FAILURES, H IS FORCED TO DECREASE */
/* BY A FACTOR OF 0.2 OR LESS. */
/* ----------------------------------------------------------------------- */
L500:
    --debdf1_1.kflag;
    debdf1_1.tn = told;
    i1 = debdf1_1.nqnyh + 1;
    i__1 = debdf1_1.nq;
    for (jb = 1; jb <= i__1; ++jb) {
	i1 -= *nyh;
	i__2 = debdf1_1.nqnyh;
	for (i__ = i1; i__ <= i__2; ++i__) {
/* L510: */
	    yh1[i__] -= yh1[i__ + *nyh];
	}
/* L515: */
    }
    debdf1_1.rmax = 2.f;
    if (dabs(debdf1_1.h__) <= debdf1_1.hmin * 1.00001f) {
	goto L660;
    }
    if (debdf1_1.kflag <= -3) {
	goto L640;
    }
    iredo = 2;
    rhup = 0.f;
    goto L540;
/* ----------------------------------------------------------------------- */
/* REGARDLESS OF THE SUCCESS OR FAILURE OF THE STEP, FACTORS */
/* RHDN, RHSM, AND RHUP ARE COMPUTED, BY WHICH H COULD BE MULTIPLIED */
/* AT ORDER NQ - 1, ORDER NQ, OR ORDER NQ + 1, RESPECTIVELY. */
/* IN THE CASE OF FAILURE, RHUP = 0.0 TO AVOID AN ORDER INCREASE. */
/* THE LARGEST OF THESE IS DETERMINED AND THE NEW ORDER CHOSEN */
/* ACCORDINGLY.  IF THE ORDER IS TO BE INCREASED, WE COMPUTE ONE */
/* ADDITIONAL SCALED DERIVATIVE. */
/* ----------------------------------------------------------------------- */
L520:
    rhup = 0.f;
    if (debdf1_1.l == debdf1_1.lmax) {
	goto L540;
    }
    i__1 = debdf1_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L530: */
	savf[i__] = acor[i__] - yh[i__ + debdf1_1.lmax * yh_dim1];
    }
    dup = vnwrms_(&debdf1_1.n, &savf[1], &ewt[1]) / debdf1_1.tesco[
	    debdf1_1.nq * 3 - 1];
    exup = 1.f / (debdf1_1.l + 1);
    d__1 = (doublereal) dup;
    d__2 = (doublereal) exup;
    rhup = 1.f / (pow_dd(&d__1, &d__2) * 1.4f + 1.4e-6f);
L540:
    exsm = 1.f / debdf1_1.l;
    d__1 = (doublereal) dsm;
    d__2 = (doublereal) exsm;
    rhsm = 1.f / (pow_dd(&d__1, &d__2) * 1.2f + 1.2e-6f);
    rhdn = 0.f;
    if (debdf1_1.nq == 1) {
	goto L560;
    }
    ddn = vnwrms_(&debdf1_1.n, &yh[debdf1_1.l * yh_dim1 + 1], &ewt[1]) / 
	    debdf1_1.tesco[debdf1_1.nq * 3 - 3];
    exdn = 1.f / debdf1_1.nq;
    d__1 = (doublereal) ddn;
    d__2 = (doublereal) exdn;
    rhdn = 1.f / (pow_dd(&d__1, &d__2) * 1.3f + 1.3e-6f);
L560:
    if (rhsm >= rhup) {
	goto L570;
    }
    if (rhup > rhdn) {
	goto L590;
    }
    goto L580;
L570:
    if (rhsm < rhdn) {
	goto L580;
    }
    newq = debdf1_1.nq;
    rh = rhsm;
    goto L620;
L580:
    newq = debdf1_1.nq - 1;
    rh = rhdn;
    if (debdf1_1.kflag < 0 && rh > 1.f) {
	rh = 1.f;
    }
    goto L620;
L590:
    newq = debdf1_1.l;
    rh = rhup;
    if (rh < 1.1f) {
	goto L610;
    }
    r__ = debdf1_1.el[debdf1_1.l - 1] / debdf1_1.l;
    i__1 = debdf1_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L600: */
	yh[i__ + (newq + 1) * yh_dim1] = acor[i__] * r__;
    }
    goto L630;
L610:
    debdf1_1.ialth = 3;
    goto L690;
L620:
    if (debdf1_1.kflag == 0 && rh < 1.1f) {
	goto L610;
    }
    if (debdf1_1.kflag <= -2) {
	rh = dmin(rh,.2f);
    }
/* ----------------------------------------------------------------------- */
/* IF THERE IS A CHANGE OF ORDER, RESET NQ, L, AND THE COEFFICIENTS. */
/* IN ANY CASE H IS RESET ACCORDING TO RH AND THE YH ARRAY IS RESCALED. */
/* THEN EXIT FROM 680 IF THE STEP WAS OK, OR REDO THE STEP OTHERWISE. */
/* ----------------------------------------------------------------------- */
    if (newq == debdf1_1.nq) {
	goto L170;
    }
L630:
    debdf1_1.nq = newq;
    debdf1_1.l = debdf1_1.nq + 1;
    iret = 2;
    goto L150;
/* ----------------------------------------------------------------------- */
/* CONTROL REACHES THIS SECTION IF 3 OR MORE FAILURES HAVE OCCURRED. */
/* IF 10 FAILURES HAVE OCCURRED, EXIT WITH KFLAG = -1. */
/* IT IS ASSUMED THAT THE DERIVATIVES THAT HAVE ACCUMULATED IN THE */
/* YH ARRAY HAVE ERRORS OF THE WRONG ORDER.  HENCE THE FIRST */
/* DERIVATIVE IS RECOMPUTED, AND THE ORDER IS SET TO 1.  THEN */
/* H IS REDUCED BY A FACTOR OF 10, AND THE STEP IS RETRIED, */
/* UNTIL IT SUCCEEDS OR H REACHES HMIN. */
/* ----------------------------------------------------------------------- */
L640:
    if (debdf1_1.kflag == -10) {
	goto L660;
    }
    rh = .1f;
/* Computing MAX */
    r__1 = debdf1_1.hmin / dabs(debdf1_1.h__);
    rh = dmax(r__1,rh);
    debdf1_1.h__ *= rh;
    i__1 = debdf1_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L645: */
	y[i__] = yh[i__ + yh_dim1];
    }
    (*f)(&debdf1_1.tn, &y[1], &savf[1], &rpar[1], &ipar[1]);
    ++debdf1_1.nfe;
    i__1 = debdf1_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L650: */
	yh[i__ + (yh_dim1 << 1)] = debdf1_1.h__ * savf[i__];
    }
    debdf1_1.ipup = debdf1_1.miter;
    debdf1_1.ialth = 5;
    if (debdf1_1.nq == 1) {
	goto L200;
    }
    debdf1_1.nq = 1;
    debdf1_1.l = 2;
    iret = 3;
    goto L150;
/* ----------------------------------------------------------------------- */
/* ALL RETURNS ARE MADE THROUGH THIS SECTION.  H IS SAVED IN HOLD */
/* TO ALLOW THE CALLER TO CHANGE H ON THE NEXT STEP. */
/* ----------------------------------------------------------------------- */
L660:
    debdf1_1.kflag = -1;
    goto L700;
L670:
    debdf1_1.kflag = -2;
    goto L700;
L680:
    debdf1_1.rmax = 10.f;
L690:
    r__ = 1.f / debdf1_1.tesco[debdf1_1.nqu * 3 - 2];
    i__1 = debdf1_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L695: */
	acor[i__] *= r__;
    }
L700:
    debdf1_1.hold = debdf1_1.h__;
    debdf1_1.jstart = 1;
    return 0;
/* ----------------------- END OF SUBROUTINE STOD  ----------------------- */
} /* stod_ */

