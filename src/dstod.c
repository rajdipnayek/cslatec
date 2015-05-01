/* dstod.f -- translated by f2c (version 12.02.01).
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
    doublereal rownd, conit, crate, el[13], elco[156]	/* was [13][12] */, 
	    hold, rc, rmax, tesco[36]	/* was [3][12] */, el0, h__, hmin, 
	    hmxi, hu, tn, uround;
    integer iownd[7], ksteps, iod[6], ialth, ipup, lmax, meo, nqnyh, nstepj, 
	    ier, jstart, kflag, l, meth, miter, maxord, n, nq, nst, nfe, nje, 
	    nqu;
} ddebd1_;

#define ddebd1_1 ddebd1_

/* DECK DSTOD */
/* Subroutine */ int dstod_(integer *neq, doublereal *y, doublereal *yh, 
	integer *nyh, doublereal *yh1, doublereal *ewt, doublereal *savf, 
	doublereal *acor, doublereal *wm, integer *iwm, S_fp df, U_fp djac, 
	doublereal *rpar, integer *ipar)
{
    /* System generated locals */
    integer yh_dim1, yh_offset, i__1, i__2;
    doublereal d__1, d__2, d__3;

    /* Local variables */
    static integer i__, j, m;
    static doublereal r__;
    static integer i1, jb;
    static doublereal rh, del, ddn;
    static integer ncf;
    static doublereal dsm, dup, dcon, delp, exdn, rhdn, told;
    static integer iret;
    static doublereal rhsm;
    static integer newq;
    static doublereal exsm, rhup, exup;
    extern /* Subroutine */ int dcfod_(integer *, doublereal *, doublereal *),
	     dpjac_(integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, integer *,
	     S_fp, U_fp, doublereal *, integer *);
    static integer iredo;
    extern /* Subroutine */ int dslvs_(doublereal *, integer *, doublereal *, 
	    doublereal *);
    extern doublereal dvnrms_(integer *, doublereal *, doublereal *);

/* ***BEGIN PROLOGUE  DSTOD */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to DDEBDF */
/* ***LIBRARY   SLATEC */
/* ***TYPE      DOUBLE PRECISION (STOD-S, DSTOD-D) */
/* ***AUTHOR  Watts, H. A., (SNLA) */
/* ***DESCRIPTION */

/*   DSTOD integrates a system of first order odes over one step in the */
/*   integrator package DDEBDF. */
/* ---------------------------------------------------------------------- */
/* DSTOD  performs one step of the integration of an initial value */
/* problem for a system of ordinary differential equations. */
/* Note.. DSTOD  is independent of the value of the iteration method */
/* indicator MITER, when this is .NE. 0, and hence is independent */
/* of the type of chord method used, or the Jacobian structure. */
/* Communication with DSTOD  is done with the following variables.. */

/* Y      = An array of length .GE. N used as the Y argument in */
/*          all calls to DF and DJAC. */
/* NEQ    = Integer array containing problem size in NEQ(1), and */
/*          passed as the NEQ argument in all calls to DF and DJAC. */
/* YH     = An NYH by LMAX array containing the dependent variables */
/*          and their approximate scaled derivatives, where */
/*          LMAX = MAXORD + 1.  YH(I,J+1) contains the approximate */
/*          J-th derivative of Y(I), scaled by H**J/FACTORIAL(J) */
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
/* WM,IWM = DOUBLE PRECISION and INTEGER work arrays associated with */
/*          matrix operations in chord iteration (MITER .NE. 0). */
/* DPJAC   = Name of routine to evaluate and preprocess Jacobian matrix */
/*          if a chord method is being used. */
/* DSLVS   = Name of routine to solve linear system in chord iteration. */
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

/* ***SEE ALSO  DDEBDF */
/* ***ROUTINES CALLED  DCFOD, DPJAC, DSLVS, DVNRMS */
/* ***COMMON BLOCKS    DDEBD1 */
/* ***REVISION HISTORY  (YYMMDD) */
/*   820301  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890911  Removed unnecessary intrinsics.  (WRB) */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900328  Added TYPE section.  (WRB) */
/*   910722  Updated AUTHOR section.  (ALS) */
/*   920422  Changed DIMENSION statement.  (WRB) */
/* ***END PROLOGUE  DSTOD */




/*     BEGIN BLOCK PERMITTING ...EXITS TO 690 */
/*        BEGIN BLOCK PERMITTING ...EXITS TO 60 */
/* ***FIRST EXECUTABLE STATEMENT  DSTOD */
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
    ddebd1_1.kflag = 0;
    told = ddebd1_1.tn;
    ncf = 0;
    if (ddebd1_1.jstart > 0) {
	goto L160;
    }
    if (ddebd1_1.jstart == -1) {
	goto L10;
    }
    if (ddebd1_1.jstart == -2) {
	goto L90;
    }
/*              --------------------------------------------------------- */
/*               ON THE FIRST CALL, THE ORDER IS SET TO 1, AND OTHER */
/*               VARIABLES ARE INITIALIZED.  RMAX IS THE MAXIMUM RATIO BY */
/*               WHICH H CAN BE INCREASED IN A SINGLE STEP.  IT IS */
/*               INITIALLY 1.E4 TO COMPENSATE FOR THE SMALL INITIAL H, */
/*               BUT THEN IS NORMALLY EQUAL TO 10.  IF A FAILURE OCCURS */
/*               (IN CORRECTOR CONVERGENCE OR ERROR TEST), RMAX IS SET AT */
/*               2 FOR THE NEXT INCREASE. */
/*              --------------------------------------------------------- */
    ddebd1_1.lmax = ddebd1_1.maxord + 1;
    ddebd1_1.nq = 1;
    ddebd1_1.l = 2;
    ddebd1_1.ialth = 2;
    ddebd1_1.rmax = 1e4;
    ddebd1_1.rc = 0.;
    ddebd1_1.el0 = 1.;
    ddebd1_1.crate = .7;
    delp = 0.;
    ddebd1_1.hold = ddebd1_1.h__;
    ddebd1_1.meo = ddebd1_1.meth;
    ddebd1_1.nstepj = 0;
    iret = 3;
    goto L50;
L10:
/*              BEGIN BLOCK PERMITTING ...EXITS TO 30 */
/*                 ------------------------------------------------------ */
/*                  THE FOLLOWING BLOCK HANDLES PRELIMINARIES NEEDED WHEN */
/*                  JSTART = -1.  IPUP IS SET TO MITER TO FORCE A MATRIX */
/*                  UPDATE.  IF AN ORDER INCREASE IS ABOUT TO BE */
/*                  CONSIDERED (IALTH = 1), IALTH IS RESET TO 2 TO */
/*                  POSTPONE CONSIDERATION ONE MORE STEP.  IF THE CALLER */
/*                  HAS CHANGED METH, DCFOD  IS CALLED TO RESET THE */
/*                  COEFFICIENTS OF THE METHOD.  IF THE CALLER HAS */
/*                  CHANGED MAXORD TO A VALUE LESS THAN THE CURRENT */
/*                  ORDER NQ, NQ IS REDUCED TO MAXORD, AND A NEW H CHOSEN */
/*                  ACCORDINGLY.  IF H IS TO BE CHANGED, YH MUST BE */
/*                  RESCALED.  IF H OR METH IS BEING CHANGED, IALTH IS */
/*                  RESET TO L = NQ + 1 TO PREVENT FURTHER CHANGES IN H */
/*                  FOR THAT MANY STEPS. */
/*                 ------------------------------------------------------ */
    ddebd1_1.ipup = ddebd1_1.miter;
    ddebd1_1.lmax = ddebd1_1.maxord + 1;
    if (ddebd1_1.ialth == 1) {
	ddebd1_1.ialth = 2;
    }
    if (ddebd1_1.meth == ddebd1_1.meo) {
	goto L20;
    }
    dcfod_(&ddebd1_1.meth, ddebd1_1.elco, ddebd1_1.tesco);
    ddebd1_1.meo = ddebd1_1.meth;
/*              ......EXIT */
    if (ddebd1_1.nq > ddebd1_1.maxord) {
	goto L30;
    }
    ddebd1_1.ialth = ddebd1_1.l;
    iret = 1;
/*        ............EXIT */
    goto L60;
L20:
    if (ddebd1_1.nq <= ddebd1_1.maxord) {
	goto L90;
    }
L30:
    ddebd1_1.nq = ddebd1_1.maxord;
    ddebd1_1.l = ddebd1_1.lmax;
    i__1 = ddebd1_1.l;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ddebd1_1.el[i__ - 1] = ddebd1_1.elco[i__ + ddebd1_1.nq * 13 - 14];
/* L40: */
    }
    ddebd1_1.nqnyh = ddebd1_1.nq * *nyh;
    ddebd1_1.rc = ddebd1_1.rc * ddebd1_1.el[0] / ddebd1_1.el0;
    ddebd1_1.el0 = ddebd1_1.el[0];
    ddebd1_1.conit = .5 / (ddebd1_1.nq + 2);
    ddn = dvnrms_(&ddebd1_1.n, &savf[1], &ewt[1]) / ddebd1_1.tesco[ddebd1_1.l 
	    * 3 - 3];
    exdn = 1. / ddebd1_1.l;
    rhdn = 1. / (pow_dd(&ddn, &exdn) * 1.3 + 1.3e-6);
    rh = min(rhdn,1.);
    iredo = 3;
    if (ddebd1_1.h__ == ddebd1_1.hold) {
	goto L660;
    }
/* Computing MIN */
    d__2 = rh, d__3 = (d__1 = ddebd1_1.h__ / ddebd1_1.hold, abs(d__1));
    rh = min(d__2,d__3);
    ddebd1_1.h__ = ddebd1_1.hold;
    goto L100;
L50:
/*           ------------------------------------------------------------ */
/*            DCFOD  IS CALLED TO GET ALL THE INTEGRATION COEFFICIENTS */
/*            FOR THE CURRENT METH.  THEN THE EL VECTOR AND RELATED */
/*            CONSTANTS ARE RESET WHENEVER THE ORDER NQ IS CHANGED, OR AT */
/*            THE START OF THE PROBLEM. */
/*           ------------------------------------------------------------ */
    dcfod_(&ddebd1_1.meth, ddebd1_1.elco, ddebd1_1.tesco);
L60:
L70:
/*           BEGIN BLOCK PERMITTING ...EXITS TO 680 */
    i__1 = ddebd1_1.l;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ddebd1_1.el[i__ - 1] = ddebd1_1.elco[i__ + ddebd1_1.nq * 13 - 14];
/* L80: */
    }
    ddebd1_1.nqnyh = ddebd1_1.nq * *nyh;
    ddebd1_1.rc = ddebd1_1.rc * ddebd1_1.el[0] / ddebd1_1.el0;
    ddebd1_1.el0 = ddebd1_1.el[0];
    ddebd1_1.conit = .5 / (ddebd1_1.nq + 2);
    switch (iret) {
	case 1:  goto L90;
	case 2:  goto L660;
	case 3:  goto L160;
    }
/*              --------------------------------------------------------- */
/*               IF H IS BEING CHANGED, THE H RATIO RH IS CHECKED AGAINST */
/*               RMAX, HMIN, AND HMXI, AND THE YH ARRAY RESCALED.  IALTH */
/*               IS SET TO L = NQ + 1 TO PREVENT A CHANGE OF H FOR THAT */
/*               MANY STEPS, UNLESS FORCED BY A CONVERGENCE OR ERROR TEST */
/*               FAILURE. */
/*              --------------------------------------------------------- */
L90:
    if (ddebd1_1.h__ == ddebd1_1.hold) {
	goto L160;
    }
    rh = ddebd1_1.h__ / ddebd1_1.hold;
    ddebd1_1.h__ = ddebd1_1.hold;
    iredo = 3;
L100:
L110:
    rh = min(rh,ddebd1_1.rmax);
/* Computing MAX */
    d__1 = 1., d__2 = abs(ddebd1_1.h__) * ddebd1_1.hmxi * rh;
    rh /= max(d__1,d__2);
    r__ = 1.;
    i__1 = ddebd1_1.l;
    for (j = 2; j <= i__1; ++j) {
	r__ *= rh;
	i__2 = ddebd1_1.n;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    yh[i__ + j * yh_dim1] *= r__;
/* L120: */
	}
/* L130: */
    }
    ddebd1_1.h__ *= rh;
    ddebd1_1.rc *= rh;
    ddebd1_1.ialth = ddebd1_1.l;
    if (iredo != 0) {
	goto L150;
    }
    ddebd1_1.rmax = 10.;
    r__ = 1. / ddebd1_1.tesco[ddebd1_1.nqu * 3 - 2];
    i__1 = ddebd1_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	acor[i__] *= r__;
/* L140: */
    }
/*     ...............EXIT */
    goto L690;
L150:
/*                 ------------------------------------------------------ */
/*                  THIS SECTION COMPUTES THE PREDICTED VALUES BY */
/*                  EFFECTIVELY MULTIPLYING THE YH ARRAY BY THE PASCAL */
/*                  TRIANGLE MATRIX.  RC IS THE RATIO OF NEW TO OLD */
/*                  VALUES OF THE COEFFICIENT  H*EL(1).  WHEN RC DIFFERS */
/*                  FROM 1 BY MORE THAN 30 PERCENT, IPUP IS SET TO MITER */
/*                  TO FORCE DPJAC TO BE CALLED, IF A JACOBIAN IS */
/*                  INVOLVED.  IN ANY CASE, DPJAC IS CALLED AT LEAST */
/*                  EVERY 20-TH STEP. */
/*                 ------------------------------------------------------ */
L160:
L170:
/*                    BEGIN BLOCK PERMITTING ...EXITS TO 610 */
/*                       BEGIN BLOCK PERMITTING ...EXITS TO 490 */
    if ((d__1 = ddebd1_1.rc - 1., abs(d__1)) > .3) {
	ddebd1_1.ipup = ddebd1_1.miter;
    }
    if (ddebd1_1.nst >= ddebd1_1.nstepj + 20) {
	ddebd1_1.ipup = ddebd1_1.miter;
    }
    ddebd1_1.tn += ddebd1_1.h__;
    i1 = ddebd1_1.nqnyh + 1;
    i__1 = ddebd1_1.nq;
    for (jb = 1; jb <= i__1; ++jb) {
	i1 -= *nyh;
	i__2 = ddebd1_1.nqnyh;
	for (i__ = i1; i__ <= i__2; ++i__) {
	    yh1[i__] += yh1[i__ + *nyh];
/* L180: */
	}
/* L190: */
    }
    ++ddebd1_1.ksteps;
/*                          --------------------------------------------- */
/*                           UP TO 3 CORRECTOR ITERATIONS ARE TAKEN.  A */
/*                           CONVERGENCE TEST IS MADE ON THE R.M.S. NORM */
/*                           OF EACH CORRECTION, WEIGHTED BY THE ERROR */
/*                           WEIGHT VECTOR EWT.  THE SUM OF THE */
/*                           CORRECTIONS IS ACCUMULATED IN THE VECTOR */
/*                           ACOR(I).  THE YH ARRAY IS NOT ALTERED IN THE */
/*                           CORRECTOR LOOP. */
/*                          --------------------------------------------- */
L200:
    m = 0;
    i__1 = ddebd1_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	y[i__] = yh[i__ + yh_dim1];
/* L210: */
    }
    (*df)(&ddebd1_1.tn, &y[1], &savf[1], &rpar[1], &ipar[1]);
    ++ddebd1_1.nfe;
    if (ddebd1_1.ipup <= 0) {
	goto L220;
    }
/*                                --------------------------------------- */
/*                                 IF INDICATED, THE MATRIX P = I - */
/*                                 H*EL(1)*J IS REEVALUATED AND */
/*                                 PREPROCESSED BEFORE STARTING THE */
/*                                 CORRECTOR ITERATION.  IPUP IS SET TO 0 */
/*                                 AS AN INDICATOR THAT THIS HAS BEEN */
/*                                 DONE. */
/*                                --------------------------------------- */
    ddebd1_1.ipup = 0;
    ddebd1_1.rc = 1.;
    ddebd1_1.nstepj = ddebd1_1.nst;
    ddebd1_1.crate = .7;
    dpjac_(neq, &y[1], &yh[yh_offset], nyh, &ewt[1], &acor[1], &savf[1], &wm[
	    1], &iwm[1], (S_fp)df, (U_fp)djac, &rpar[1], &ipar[1]);
/*                          ......EXIT */
    if (ddebd1_1.ier != 0) {
	goto L440;
    }
L220:
    i__1 = ddebd1_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	acor[i__] = 0.;
/* L230: */
    }
L240:
    if (ddebd1_1.miter != 0) {
	goto L270;
    }
/*                                   ------------------------------------ */
/*                                    IN THE CASE OF FUNCTIONAL */
/*                                    ITERATION, UPDATE Y DIRECTLY FROM */
/*                                    THE RESULT OF THE LAST FUNCTION */
/*                                    EVALUATION. */
/*                                   ------------------------------------ */
    i__1 = ddebd1_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	savf[i__] = ddebd1_1.h__ * savf[i__] - yh[i__ + (yh_dim1 << 1)];
	y[i__] = savf[i__] - acor[i__];
/* L250: */
    }
    del = dvnrms_(&ddebd1_1.n, &y[1], &ewt[1]);
    i__1 = ddebd1_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	y[i__] = yh[i__ + yh_dim1] + ddebd1_1.el[0] * savf[i__];
	acor[i__] = savf[i__];
/* L260: */
    }
    goto L300;
L270:
/*                                   ------------------------------------ */
/*                                    IN THE CASE OF THE CHORD METHOD, */
/*                                    COMPUTE THE CORRECTOR ERROR, AND */
/*                                    SOLVE THE LINEAR SYSTEM WITH THAT */
/*                                    AS RIGHT-HAND SIDE AND P AS */
/*                                    COEFFICIENT MATRIX. */
/*                                   ------------------------------------ */
    i__1 = ddebd1_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	y[i__] = ddebd1_1.h__ * savf[i__] - (yh[i__ + (yh_dim1 << 1)] + acor[
		i__]);
/* L280: */
    }
    dslvs_(&wm[1], &iwm[1], &y[1], &savf[1]);
/*                             ......EXIT */
    if (ddebd1_1.ier != 0) {
	goto L430;
    }
    del = dvnrms_(&ddebd1_1.n, &y[1], &ewt[1]);
    i__1 = ddebd1_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	acor[i__] += y[i__];
	y[i__] = yh[i__ + yh_dim1] + ddebd1_1.el[0] * acor[i__];
/* L290: */
    }
L300:
/*                                --------------------------------------- */
/*                                 TEST FOR CONVERGENCE.  IF M.GT.0, AN */
/*                                 ESTIMATE OF THE CONVERGENCE RATE */
/*                                 CONSTANT IS STORED IN CRATE, AND THIS */
/*                                 IS USED IN THE TEST. */
/*                                --------------------------------------- */
    if (m != 0) {
/* Computing MAX */
	d__1 = ddebd1_1.crate * .2, d__2 = del / delp;
	ddebd1_1.crate = max(d__1,d__2);
    }
/* Computing MIN */
    d__1 = 1., d__2 = ddebd1_1.crate * 1.5;
    dcon = del * min(d__1,d__2) / (ddebd1_1.tesco[ddebd1_1.nq * 3 - 2] * 
	    ddebd1_1.conit);
    if (dcon > 1.) {
	goto L420;
    }
/*                                   ------------------------------------ */
/*                                    THE CORRECTOR HAS CONVERGED.  IPUP */
/*                                    IS SET TO -1 IF MITER .NE. 0, TO */
/*                                    SIGNAL THAT THE JACOBIAN INVOLVED */
/*                                    MAY NEED UPDATING LATER.  THE LOCAL */
/*                                    ERROR TEST IS MADE AND CONTROL */
/*                                    PASSES TO STATEMENT 500 IF IT */
/*                                    FAILS. */
/*                                   ------------------------------------ */
    if (ddebd1_1.miter != 0) {
	ddebd1_1.ipup = -1;
    }
    if (m == 0) {
	dsm = del / ddebd1_1.tesco[ddebd1_1.nq * 3 - 2];
    }
    if (m > 0) {
	dsm = dvnrms_(&ddebd1_1.n, &acor[1], &ewt[1]) / ddebd1_1.tesco[
		ddebd1_1.nq * 3 - 2];
    }
    if (dsm > 1.) {
	goto L380;
    }
/*                                      BEGIN BLOCK */
/*                                      PERMITTING ...EXITS TO 360 */
/*                                         ------------------------------ */
/*                                          AFTER A SUCCESSFUL STEP, */
/*                                          UPDATE THE YH ARRAY. */
/*                                          CONSIDER CHANGING H IF IALTH */
/*                                          = 1.  OTHERWISE DECREASE */
/*                                          IALTH BY 1.  IF IALTH IS THEN */
/*                                          1 AND NQ .LT. MAXORD, THEN */
/*                                          ACOR IS SAVED FOR USE IN A */
/*                                          POSSIBLE ORDER INCREASE ON */
/*                                          THE NEXT STEP.  IF A CHANGE */
/*                                          IN H IS CONSIDERED, AN */
/*                                          INCREASE OR DECREASE IN ORDER */
/*                                          BY ONE IS CONSIDERED ALSO.  A */
/*                                          CHANGE IN H IS MADE ONLY IF */
/*                                          IT IS BY A FACTOR OF AT LEAST */
/*                                          1.1.  IF NOT, IALTH IS SET TO */
/*                                          3 TO PREVENT TESTING FOR THAT */
/*                                          MANY STEPS. */
/*                                         ------------------------------ */
    ddebd1_1.kflag = 0;
    iredo = 0;
    ++ddebd1_1.nst;
    ddebd1_1.hu = ddebd1_1.h__;
    ddebd1_1.nqu = ddebd1_1.nq;
    i__1 = ddebd1_1.l;
    for (j = 1; j <= i__1; ++j) {
	i__2 = ddebd1_1.n;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    yh[i__ + j * yh_dim1] += ddebd1_1.el[j - 1] * acor[i__];
/* L310: */
	}
/* L320: */
    }
    --ddebd1_1.ialth;
    if (ddebd1_1.ialth != 0) {
	goto L340;
    }
/*                                            --------------------------- */
/*                                             REGARDLESS OF THE SUCCESS */
/*                                             OR FAILURE OF THE STEP, */
/*                                             FACTORS RHDN, RHSM, AND */
/*                                             RHUP ARE COMPUTED, BY */
/*                                             WHICH H COULD BE */
/*                                             MULTIPLIED AT ORDER NQ - */
/*                                             1, ORDER NQ, OR ORDER NQ + */
/*                                             1, RESPECTIVELY.  IN THE */
/*                                             CASE OF FAILURE, RHUP = */
/*                                             0.0 TO AVOID AN ORDER */
/*                                             INCREASE.  THE LARGEST OF */
/*                                             THESE IS DETERMINED AND */
/*                                             THE NEW ORDER CHOSEN */
/*                                             ACCORDINGLY.  IF THE ORDER */
/*                                             IS TO BE INCREASED, WE */
/*                                             COMPUTE ONE ADDITIONAL */
/*                                             SCALED DERIVATIVE. */
/*                                            --------------------------- */
    rhup = 0.;
/*                       .....................EXIT */
    if (ddebd1_1.l == ddebd1_1.lmax) {
	goto L490;
    }
    i__1 = ddebd1_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	savf[i__] = acor[i__] - yh[i__ + ddebd1_1.lmax * yh_dim1];
/* L330: */
    }
    dup = dvnrms_(&ddebd1_1.n, &savf[1], &ewt[1]) / ddebd1_1.tesco[
	    ddebd1_1.nq * 3 - 1];
    exup = 1. / (ddebd1_1.l + 1);
    rhup = 1. / (pow_dd(&dup, &exup) * 1.4 + 1.4e-6);
/*                       .....................EXIT */
    goto L490;
L340:
/*                                      ...EXIT */
    if (ddebd1_1.ialth > 1) {
	goto L360;
    }
/*                                      ...EXIT */
    if (ddebd1_1.l == ddebd1_1.lmax) {
	goto L360;
    }
    i__1 = ddebd1_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	yh[i__ + ddebd1_1.lmax * yh_dim1] = acor[i__];
/* L350: */
    }
L360:
    r__ = 1. / ddebd1_1.tesco[ddebd1_1.nqu * 3 - 2];
    i__1 = ddebd1_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	acor[i__] *= r__;
/* L370: */
    }
/*     .................................EXIT */
    goto L690;
L380:
/*                                   ------------------------------------ */
/*                                    THE ERROR TEST FAILED.  KFLAG KEEPS */
/*                                    TRACK OF MULTIPLE FAILURES. */
/*                                    RESTORE TN AND THE YH ARRAY TO */
/*                                    THEIR PREVIOUS VALUES, AND PREPARE */
/*                                    TO TRY THE STEP AGAIN.  COMPUTE THE */
/*                                    OPTIMUM STEP SIZE FOR THIS OR ONE */
/*                                    LOWER ORDER.  AFTER 2 OR MORE */
/*                                    FAILURES, H IS FORCED TO DECREASE */
/*                                    BY A FACTOR OF 0.2 OR LESS. */
/*                                   ------------------------------------ */
    --ddebd1_1.kflag;
    ddebd1_1.tn = told;
    i1 = ddebd1_1.nqnyh + 1;
    i__1 = ddebd1_1.nq;
    for (jb = 1; jb <= i__1; ++jb) {
	i1 -= *nyh;
	i__2 = ddebd1_1.nqnyh;
	for (i__ = i1; i__ <= i__2; ++i__) {
	    yh1[i__] -= yh1[i__ + *nyh];
/* L390: */
	}
/* L400: */
    }
    ddebd1_1.rmax = 2.;
    if (abs(ddebd1_1.h__) > ddebd1_1.hmin * 1.00001) {
	goto L410;
    }
/*                                      --------------------------------- */
/*                                       ALL RETURNS ARE MADE THROUGH */
/*                                       THIS SECTION.  H IS SAVED IN */
/*                                       HOLD TO ALLOW THE CALLER TO */
/*                                       CHANGE H ON THE NEXT STEP. */
/*                                      --------------------------------- */
    ddebd1_1.kflag = -1;
/*     .................................EXIT */
    goto L690;
L410:
/*                    ...............EXIT */
    if (ddebd1_1.kflag <= -3) {
	goto L610;
    }
    iredo = 2;
    rhup = 0.;
/*                       ............EXIT */
    goto L490;
L420:
    ++m;
/*                             ...EXIT */
    if (m == 3) {
	goto L430;
    }
/*                             ...EXIT */
    if (m >= 2 && del > delp * 2.) {
	goto L430;
    }
    delp = del;
    (*df)(&ddebd1_1.tn, &y[1], &savf[1], &rpar[1], &ipar[1]);
    ++ddebd1_1.nfe;
    goto L240;
L430:
/*                             ------------------------------------------ */
/*                              THE CORRECTOR ITERATION FAILED TO */
/*                              CONVERGE IN 3 TRIES.  IF MITER .NE. 0 AND */
/*                              THE JACOBIAN IS OUT OF DATE, DPJAC IS */
/*                              CALLED FOR THE NEXT TRY.  OTHERWISE THE */
/*                              YH ARRAY IS RETRACTED TO ITS VALUES */
/*                              BEFORE PREDICTION, AND H IS REDUCED, IF */
/*                              POSSIBLE.  IF H CANNOT BE REDUCED OR 10 */
/*                              FAILURES HAVE OCCURRED, EXIT WITH KFLAG = */
/*                              -2. */
/*                             ------------------------------------------ */
/*                          ...EXIT */
    if (ddebd1_1.ipup == 0) {
	goto L440;
    }
    ddebd1_1.ipup = ddebd1_1.miter;
    goto L200;
L440:
    ddebd1_1.tn = told;
    ++ncf;
    ddebd1_1.rmax = 2.;
    i1 = ddebd1_1.nqnyh + 1;
    i__1 = ddebd1_1.nq;
    for (jb = 1; jb <= i__1; ++jb) {
	i1 -= *nyh;
	i__2 = ddebd1_1.nqnyh;
	for (i__ = i1; i__ <= i__2; ++i__) {
	    yh1[i__] -= yh1[i__ + *nyh];
/* L450: */
	}
/* L460: */
    }
    if (abs(ddebd1_1.h__) > ddebd1_1.hmin * 1.00001) {
	goto L470;
    }
    ddebd1_1.kflag = -2;
/*     ........................EXIT */
    goto L690;
L470:
    if (ncf != 10) {
	goto L480;
    }
    ddebd1_1.kflag = -2;
/*     ........................EXIT */
    goto L690;
L480:
    rh = .25;
    ddebd1_1.ipup = ddebd1_1.miter;
    iredo = 1;
/*                 .........EXIT */
    goto L650;
L490:
    exsm = 1. / ddebd1_1.l;
    rhsm = 1. / (pow_dd(&dsm, &exsm) * 1.2 + 1.2e-6);
    rhdn = 0.;
    if (ddebd1_1.nq == 1) {
	goto L500;
    }
    ddn = dvnrms_(&ddebd1_1.n, &yh[ddebd1_1.l * yh_dim1 + 1], &ewt[1]) / 
	    ddebd1_1.tesco[ddebd1_1.nq * 3 - 3];
    exdn = 1. / ddebd1_1.nq;
    rhdn = 1. / (pow_dd(&ddn, &exdn) * 1.3 + 1.3e-6);
L500:
    if (rhsm >= rhup) {
	goto L550;
    }
    if (rhup <= rhdn) {
	goto L540;
    }
    newq = ddebd1_1.l;
    rh = rhup;
    if (rh >= 1.1) {
	goto L520;
    }
    ddebd1_1.ialth = 3;
    r__ = 1. / ddebd1_1.tesco[ddebd1_1.nqu * 3 - 2];
    i__1 = ddebd1_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	acor[i__] *= r__;
/* L510: */
    }
/*     ...........................EXIT */
    goto L690;
L520:
    r__ = ddebd1_1.el[ddebd1_1.l - 1] / ddebd1_1.l;
    i__1 = ddebd1_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	yh[i__ + (newq + 1) * yh_dim1] = acor[i__] * r__;
/* L530: */
    }
    ddebd1_1.nq = newq;
    ddebd1_1.l = ddebd1_1.nq + 1;
    iret = 2;
/*           ..................EXIT */
    goto L680;
L540:
    goto L580;
L550:
    if (rhsm < rhdn) {
	goto L580;
    }
    newq = ddebd1_1.nq;
    rh = rhsm;
    if (ddebd1_1.kflag == 0 && rh < 1.1) {
	goto L560;
    }
    if (ddebd1_1.kflag <= -2) {
	rh = min(rh,.2);
    }
/*                             ------------------------------------------ */
/*                              IF THERE IS A CHANGE OF ORDER, RESET NQ, */
/*                              L, AND THE COEFFICIENTS.  IN ANY CASE H */
/*                              IS RESET ACCORDING TO RH AND THE YH ARRAY */
/*                              IS RESCALED.  THEN EXIT FROM 680 IF THE */
/*                              STEP WAS OK, OR REDO THE STEP OTHERWISE. */
/*                             ------------------------------------------ */
/*                 ............EXIT */
    if (newq == ddebd1_1.nq) {
	goto L650;
    }
    ddebd1_1.nq = newq;
    ddebd1_1.l = ddebd1_1.nq + 1;
    iret = 2;
/*           ..................EXIT */
    goto L680;
L560:
    ddebd1_1.ialth = 3;
    r__ = 1. / ddebd1_1.tesco[ddebd1_1.nqu * 3 - 2];
    i__1 = ddebd1_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	acor[i__] *= r__;
/* L570: */
    }
/*     .....................EXIT */
    goto L690;
L580:
    newq = ddebd1_1.nq - 1;
    rh = rhdn;
    if (ddebd1_1.kflag < 0 && rh > 1.) {
	rh = 1.;
    }
    if (ddebd1_1.kflag == 0 && rh < 1.1) {
	goto L590;
    }
    if (ddebd1_1.kflag <= -2) {
	rh = min(rh,.2);
    }
/*                          --------------------------------------------- */
/*                           IF THERE IS A CHANGE OF ORDER, RESET NQ, L, */
/*                           AND THE COEFFICIENTS.  IN ANY CASE H IS */
/*                           RESET ACCORDING TO RH AND THE YH ARRAY IS */
/*                           RESCALED.  THEN EXIT FROM 680 IF THE STEP */
/*                           WAS OK, OR REDO THE STEP OTHERWISE. */
/*                          --------------------------------------------- */
/*                 .........EXIT */
    if (newq == ddebd1_1.nq) {
	goto L650;
    }
    ddebd1_1.nq = newq;
    ddebd1_1.l = ddebd1_1.nq + 1;
    iret = 2;
/*           ...............EXIT */
    goto L680;
L590:
    ddebd1_1.ialth = 3;
    r__ = 1. / ddebd1_1.tesco[ddebd1_1.nqu * 3 - 2];
    i__1 = ddebd1_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	acor[i__] *= r__;
/* L600: */
    }
/*     ..................EXIT */
    goto L690;
L610:
/*                    --------------------------------------------------- */
/*                     CONTROL REACHES THIS SECTION IF 3 OR MORE FAILURES */
/*                     HAVE OCCURRED.  IF 10 FAILURES HAVE OCCURRED, EXIT */
/*                     WITH KFLAG = -1.  IT IS ASSUMED THAT THE */
/*                     DERIVATIVES THAT HAVE ACCUMULATED IN THE YH ARRAY */
/*                     HAVE ERRORS OF THE WRONG ORDER.  HENCE THE FIRST */
/*                     DERIVATIVE IS RECOMPUTED, AND THE ORDER IS SET TO */
/*                     1.  THEN H IS REDUCED BY A FACTOR OF 10, AND THE */
/*                     STEP IS RETRIED, UNTIL IT SUCCEEDS OR H REACHES */
/*                     HMIN. */
/*                    --------------------------------------------------- */
    if (ddebd1_1.kflag != -10) {
	goto L620;
    }
/*                       ------------------------------------------------ */
/*                        ALL RETURNS ARE MADE THROUGH THIS SECTION.  H */
/*                        IS SAVED IN HOLD TO ALLOW THE CALLER TO CHANGE */
/*                        H ON THE NEXT STEP. */
/*                       ------------------------------------------------ */
    ddebd1_1.kflag = -1;
/*     ..................EXIT */
    goto L690;
L620:
    rh = .1;
/* Computing MAX */
    d__1 = ddebd1_1.hmin / abs(ddebd1_1.h__);
    rh = max(d__1,rh);
    ddebd1_1.h__ *= rh;
    i__1 = ddebd1_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	y[i__] = yh[i__ + yh_dim1];
/* L630: */
    }
    (*df)(&ddebd1_1.tn, &y[1], &savf[1], &rpar[1], &ipar[1]);
    ++ddebd1_1.nfe;
    i__1 = ddebd1_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	yh[i__ + (yh_dim1 << 1)] = ddebd1_1.h__ * savf[i__];
/* L640: */
    }
    ddebd1_1.ipup = ddebd1_1.miter;
    ddebd1_1.ialth = 5;
/*              ......EXIT */
    if (ddebd1_1.nq != 1) {
	goto L670;
    }
    goto L170;
L650:
L660:
/* Computing MAX */
    d__1 = rh, d__2 = ddebd1_1.hmin / abs(ddebd1_1.h__);
    rh = max(d__1,d__2);
    goto L110;
L670:
    ddebd1_1.nq = 1;
    ddebd1_1.l = 2;
    iret = 3;
L680:
    goto L70;
L690:
    ddebd1_1.hold = ddebd1_1.h__;
    ddebd1_1.jstart = 1;
    return 0;
/*     ----------------------- END OF SUBROUTINE DSTOD */
/*     ----------------------- */
} /* dstod_ */

