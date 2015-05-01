/* steps.f -- translated by f2c (version 12.02.01).
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
static integer c__2 = 2;
static integer c__1 = 1;

/* DECK STEPS */
/* Subroutine */ int steps_(S_fp f, integer *neqn, real *y, real *x, real *
	h__, real *eps, real *wt, logical *start, real *hold, integer *k, 
	integer *kold, logical *crash, real *phi, real *p, real *yp, real *
	psi, real *alpha, real *beta, real *sig, real *v, real *w, real *g, 
	logical *phase1, integer *ns, logical *nornd, integer *ksteps, real *
	twou, real *fouru, real *xold, integer *kprev, integer *ivc, integer *
	iv, integer *kgi, real *gi, real *rpar, integer *ipar)
{
    /* Initialized data */

    static real two[13] = { 2.f,4.f,8.f,16.f,32.f,64.f,128.f,256.f,512.f,
	    1024.f,2048.f,4096.f,8192.f };
    static real gstr[13] = { .5f,.0833f,.0417f,.0264f,.0188f,.0143f,.0114f,
	    .00936f,.00789f,.00679f,.00592f,.00524f,.00468f };

    /* System generated locals */
    integer phi_dim1, phi_offset, i__1, i__2;
    real r__1, r__2, r__3;
    doublereal d__1, d__2;

    /* Local variables */
    static integer i__, j, l;
    static real r__, u;
    static integer iq, jv, im1, km1, km2, ip1, kp1, kp2;
    static real big, erk, err, tau, rho;
    static integer nsm2, nsp1, nsp2;
    static real absh, hnew;
    static integer knew;
    static real erkm1, erkm2, erkp1, temp1, temp2, temp3, temp4, temp5, temp6,
	     p5eps;
    static integer ifail;
    static real reali, round;
    extern doublereal r1mach_(integer *);
    static integer limit1, limit2;
    static real realns;
    extern /* Subroutine */ int hstart_(S_fp, integer *, real *, real *, real 
	    *, real *, real *, integer *, real *, real *, real *, real *, 
	    real *, real *, real *, integer *, real *);

/* ***BEGIN PROLOGUE  STEPS */
/* ***PURPOSE  Integrate a system of first order ordinary differential */
/*            equations one step. */
/* ***LIBRARY   SLATEC (DEPAC) */
/* ***CATEGORY  I1A1B */
/* ***TYPE      SINGLE PRECISION (STEPS-S, DSTEPS-D) */
/* ***KEYWORDS  ADAMS METHOD, DEPAC, INITIAL VALUE PROBLEMS, ODE, */
/*             ORDINARY DIFFERENTIAL EQUATIONS, PREDICTOR-CORRECTOR */
/* ***AUTHOR  Shampine, L. F., (SNLA) */
/*           Gordon, M. K., (SNLA) */
/*             MODIFIED BY H.A. WATTS */
/* ***DESCRIPTION */

/*   Written by L. F. Shampine and M. K. Gordon */

/*   Abstract */

/*   Subroutine  STEPS  is normally used indirectly through subroutine */
/*   DEABM .  Because  DEABM  suffices for most problems and is much */
/*   easier to use, using it should be considered before using  STEPS */
/*   alone. */

/*   Subroutine STEPS integrates a system of  NEQN  first order ordinary */
/*   differential equations one step, normally from X to X+H, using a */
/*   modified divided difference form of the Adams Pece formulas.  Local */
/*   extrapolation is used to improve absolute stability and accuracy. */
/*   The code adjusts its order and step size to control the local error */
/*   per unit step in a generalized sense.  Special devices are included */
/*   to control roundoff error and to detect when the user is requesting */
/*   too much accuracy. */

/*   This code is completely explained and documented in the text, */
/*   Computer Solution of Ordinary Differential Equations, The Initial */
/*   Value Problem  by L. F. Shampine and M. K. Gordon. */
/*   Further details on use of this code are available in "Solving */
/*   Ordinary Differential Equations with ODE, STEP, and INTRP", */
/*   by L. F. Shampine and M. K. Gordon, SLA-73-1060. */


/*   The parameters represent -- */
/*      F -- subroutine to evaluate derivatives */
/*      NEQN -- number of equations to be integrated */
/*      Y(*) -- solution vector at X */
/*      X -- independent variable */
/*      H -- appropriate step size for next step.  Normally determined by */
/*           code */
/*      EPS -- local error tolerance */
/*      WT(*) -- vector of weights for error criterion */
/*      START -- logical variable set .TRUE. for first step,  .FALSE. */
/*           otherwise */
/*      HOLD -- step size used for last successful step */
/*      K -- appropriate order for next step (determined by code) */
/*      KOLD -- order used for last successful step */
/*      CRASH -- logical variable set .TRUE. when no step can be taken, */
/*           .FALSE. otherwise. */
/*      YP(*) -- derivative of solution vector at  X  after successful */
/*           step */
/*      KSTEPS -- counter on attempted steps */
/*      TWOU -- 2.*U where U is machine unit roundoff quantity */
/*      FOURU -- 4.*U where U is machine unit roundoff quantity */
/*      RPAR,IPAR -- parameter arrays which you may choose to use */
/*            for communication between your program and subroutine F. */
/*            They are not altered or used by STEPS. */
/*   The variables X,XOLD,KOLD,KGI and IVC and the arrays Y,PHI,ALPHA,G, */
/*   W,P,IV and GI are required for the interpolation subroutine SINTRP. */
/*   The remaining variables and arrays are included in the call list */
/*   only to eliminate local retention of variables between calls. */

/*   Input to STEPS */

/*      First call -- */

/*   The user must provide storage in his calling program for all arrays */
/*   in the call list, namely */

/*     DIMENSION Y(NEQN),WT(NEQN),PHI(NEQN,16),P(NEQN),YP(NEQN),PSI(12), */
/*    1  ALPHA(12),BETA(12),SIG(13),V(12),W(12),G(13),GI(11),IV(10), */
/*    2  RPAR(*),IPAR(*) */

/*    **Note** */

/*   The user must also declare  START ,  CRASH ,  PHASE1  and  NORND */
/*   logical variables and  F  an EXTERNAL subroutine, supply the */
/*   subroutine  F(X,Y,YP)  to evaluate */
/*      DY(I)/DX = YP(I) = F(X,Y(1),Y(2),...,Y(NEQN)) */
/*   and initialize only the following parameters. */
/*      NEQN -- number of equations to be integrated */
/*      Y(*) -- vector of initial values of dependent variables */
/*      X -- initial value of the independent variable */
/*      H -- nominal step size indicating direction of integration */
/*           and maximum size of step.  Must be variable */
/*      EPS -- local error tolerance per step.  Must be variable */
/*      WT(*) -- vector of non-zero weights for error criterion */
/*      START -- .TRUE. */
/*      YP(*) -- vector of initial derivative values */
/*      KSTEPS -- set KSTEPS to zero */
/*      TWOU -- 2.*U where U is machine unit roundoff quantity */
/*      FOURU -- 4.*U where U is machine unit roundoff quantity */
/*   Define U to be the machine unit roundoff quantity by calling */
/*   the function routine  R1MACH,  U = R1MACH(4), or by */
/*   computing U so that U is the smallest positive number such */
/*   that 1.0+U .GT. 1.0. */

/*   STEPS  requires that the L2 norm of the vector with components */
/*   LOCAL ERROR(L)/WT(L)  be less than  EPS  for a successful step.  The */
/*   array  WT  allows the user to specify an error test appropriate */
/*   for his problem.  For example, */
/*      WT(L) = 1.0  specifies absolute error, */
/*            = ABS(Y(L))  error relative to the most recent value of the */
/*                 L-th component of the solution, */
/*            = ABS(YP(L))  error relative to the most recent value of */
/*                 the L-th component of the derivative, */
/*            = MAX(WT(L),ABS(Y(L)))  error relative to the largest */
/*                 magnitude of L-th component obtained so far, */
/*            = ABS(Y(L))*RELERR/EPS + ABSERR/EPS  specifies a mixed */
/*                 relative-absolute test where  RELERR  is relative */
/*                 error,  ABSERR  is absolute error and  EPS = */
/*                 MAX(RELERR,ABSERR) . */

/*      Subsequent calls -- */

/*   Subroutine  STEPS  is designed so that all information needed to */
/*   continue the integration, including the step size  H  and the order */
/*   K , is returned with each step.  With the exception of the step */
/*   size, the error tolerance, and the weights, none of the parameters */
/*   should be altered.  The array  WT  must be updated after each step */
/*   to maintain relative error tests like those above.  Normally the */
/*   integration is continued just beyond the desired endpoint and the */
/*   solution interpolated there with subroutine  SINTRP .  If it is */
/*   impossible to integrate beyond the endpoint, the step size may be */
/*   reduced to hit the endpoint since the code will not take a step */
/*   larger than the  H  input.  Changing the direction of integration, */
/*   i.e., the sign of  H , requires the user set  START = .TRUE. before */
/*   calling  STEPS  again.  This is the only situation in which  START */
/*   should be altered. */

/*   Output from STEPS */

/*      Successful Step -- */

/*   The subroutine returns after each successful step with  START  and */
/*   CRASH  set .FALSE. .  X  represents the independent variable */
/*   advanced one step of length  HOLD  from its value on input and  Y */
/*   the solution vector at the new value of  X .  All other parameters */
/*   represent information corresponding to the new  X  needed to */
/*   continue the integration. */

/*      Unsuccessful Step -- */

/*   When the error tolerance is too small for the machine precision, */
/*   the subroutine returns without taking a step and  CRASH = .TRUE. . */
/*   An appropriate step size and error tolerance for continuing are */
/*   estimated and all other information is restored as upon input */
/*   before returning.  To continue with the larger tolerance, the user */
/*   just calls the code again.  A restart is neither required nor */
/*   desirable. */

/* ***REFERENCES  L. F. Shampine and M. K. Gordon, Solving ordinary */
/*                 differential equations with ODE, STEP, and INTRP, */
/*                 Report SLA-73-1060, Sandia Laboratories, 1973. */
/* ***ROUTINES CALLED  HSTART, R1MACH */
/* ***REVISION HISTORY  (YYMMDD) */
/*   740101  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890831  Modified array declarations.  (WRB) */
/*   890831  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   920501  Reformatted the REFERENCES section.  (WRB) */
/* ***END PROLOGUE  STEPS */


    /* Parameter adjustments */
    phi_dim1 = *neqn;
    phi_offset = 1 + phi_dim1;
    phi -= phi_offset;
    --y;
    --wt;
    --p;
    --yp;
    --psi;
    --alpha;
    --beta;
    --sig;
    --v;
    --w;
    --g;
    --iv;
    --gi;
    --rpar;
    --ipar;

    /* Function Body */


/*       ***     BEGIN BLOCK 0     *** */
/*   CHECK IF STEP SIZE OR ERROR TOLERANCE IS TOO SMALL FOR MACHINE */
/*   PRECISION.  IF FIRST STEP, INITIALIZE PHI ARRAY AND ESTIMATE A */
/*   STARTING STEP SIZE. */
/*                   *** */

/*   IF STEP SIZE IS TOO SMALL, DETERMINE AN ACCEPTABLE ONE */

/* ***FIRST EXECUTABLE STATEMENT  STEPS */
    *crash = TRUE_;
    if (dabs(*h__) >= *fouru * dabs(*x)) {
	goto L5;
    }
    r__1 = *fouru * dabs(*x);
    *h__ = r_sign(&r__1, h__);
    return 0;
L5:
    p5eps = *eps * .5f;

/*   IF ERROR TOLERANCE IS TOO SMALL, INCREASE IT TO AN ACCEPTABLE VALUE */

    round = 0.f;
    i__1 = *neqn;
    for (l = 1; l <= i__1; ++l) {
/* L10: */
/* Computing 2nd power */
	r__1 = y[l] / wt[l];
	round += r__1 * r__1;
    }
    round = *twou * sqrt(round);
    if (p5eps >= round) {
	goto L15;
    }
    *eps = round * 2.f * (*fouru + 1.f);
    return 0;
L15:
    *crash = FALSE_;
    g[1] = 1.f;
    g[2] = .5f;
    sig[1] = 1.f;
    if (! (*start)) {
	goto L99;
    }

/*   INITIALIZE.  COMPUTE APPROPRIATE STEP SIZE FOR FIRST STEP */

/*     CALL F(X,Y,YP,RPAR,IPAR) */
/*     SUM = 0.0 */
    i__1 = *neqn;
    for (l = 1; l <= i__1; ++l) {
	phi[l + phi_dim1] = yp[l];
/* L20: */
	phi[l + (phi_dim1 << 1)] = 0.f;
    }
/* 20     SUM = SUM + (YP(L)/WT(L))**2 */
/*     SUM = SQRT(SUM) */
/*     ABSH = ABS(H) */
/*     IF(EPS .LT. 16.0*SUM*H*H) ABSH = 0.25*SQRT(EPS/SUM) */
/*     H = SIGN(MAX(ABSH,FOURU*ABS(X)),H) */

    u = r1mach_(&c__4);
    big = sqrt(r1mach_(&c__2));
    r__1 = *x + *h__;
    hstart_((S_fp)f, neqn, x, &r__1, &y[1], &yp[1], &wt[1], &c__1, &u, &big, &
	    phi[phi_dim1 * 3 + 1], &phi[(phi_dim1 << 2) + 1], &phi[phi_dim1 * 
	    5 + 1], &phi[phi_dim1 * 6 + 1], &rpar[1], &ipar[1], h__);

    *hold = 0.f;
    *k = 1;
    *kold = 0;
    *kprev = 0;
    *start = FALSE_;
    *phase1 = TRUE_;
    *nornd = TRUE_;
    if (p5eps > round * 100.f) {
	goto L99;
    }
    *nornd = FALSE_;
    i__1 = *neqn;
    for (l = 1; l <= i__1; ++l) {
/* L25: */
	phi[l + phi_dim1 * 15] = 0.f;
    }
L99:
    ifail = 0;
/*       ***     END BLOCK 0     *** */

/*       ***     BEGIN BLOCK 1     *** */
/*   COMPUTE COEFFICIENTS OF FORMULAS FOR THIS STEP.  AVOID COMPUTING */
/*   THOSE QUANTITIES NOT CHANGED WHEN STEP SIZE IS NOT CHANGED. */
/*                   *** */

L100:
    kp1 = *k + 1;
    kp2 = *k + 2;
    km1 = *k - 1;
    km2 = *k - 2;

/*   NS IS THE NUMBER OF STEPS TAKEN WITH SIZE H, INCLUDING THE CURRENT */
/*   ONE.  WHEN K.LT.NS, NO COEFFICIENTS CHANGE */

    if (*h__ != *hold) {
	*ns = 0;
    }
    if (*ns <= *kold) {
	++(*ns);
    }
    nsp1 = *ns + 1;
    if (*k < *ns) {
	goto L199;
    }

/*   COMPUTE THOSE COMPONENTS OF ALPHA(*),BETA(*),PSI(*),SIG(*) WHICH */
/*   ARE CHANGED */

    beta[*ns] = 1.f;
    realns = (real) (*ns);
    alpha[*ns] = 1.f / realns;
    temp1 = *h__ * realns;
    sig[nsp1] = 1.f;
    if (*k < nsp1) {
	goto L110;
    }
    i__1 = *k;
    for (i__ = nsp1; i__ <= i__1; ++i__) {
	im1 = i__ - 1;
	temp2 = psi[im1];
	psi[im1] = temp1;
	beta[i__] = beta[im1] * psi[im1] / temp2;
	temp1 = temp2 + *h__;
	alpha[i__] = *h__ / temp1;
	reali = (real) i__;
/* L105: */
	sig[i__ + 1] = reali * alpha[i__] * sig[i__];
    }
L110:
    psi[*k] = temp1;

/*   COMPUTE COEFFICIENTS G(*) */

/*   INITIALIZE V(*) AND SET W(*). */

    if (*ns > 1) {
	goto L120;
    }
    i__1 = *k;
    for (iq = 1; iq <= i__1; ++iq) {
	temp3 = (real) (iq * (iq + 1));
	v[iq] = 1.f / temp3;
/* L115: */
	w[iq] = v[iq];
    }
    *ivc = 0;
    *kgi = 0;
    if (*k == 1) {
	goto L140;
    }
    *kgi = 1;
    gi[1] = w[2];
    goto L140;

/*   IF ORDER WAS RAISED, UPDATE DIAGONAL PART OF V(*) */

L120:
    if (*k <= *kprev) {
	goto L130;
    }
    if (*ivc == 0) {
	goto L122;
    }
    jv = kp1 - iv[*ivc];
    --(*ivc);
    goto L123;
L122:
    jv = 1;
    temp4 = (real) (*k * kp1);
    v[*k] = 1.f / temp4;
    w[*k] = v[*k];
    if (*k != 2) {
	goto L123;
    }
    *kgi = 1;
    gi[1] = w[2];
L123:
    nsm2 = *ns - 2;
    if (nsm2 < jv) {
	goto L130;
    }
    i__1 = nsm2;
    for (j = jv; j <= i__1; ++j) {
	i__ = *k - j;
	v[i__] -= alpha[j + 1] * v[i__ + 1];
/* L125: */
	w[i__] = v[i__];
    }
    if (i__ != 2) {
	goto L130;
    }
    *kgi = *ns - 1;
    gi[*kgi] = w[2];

/*   UPDATE V(*) AND SET W(*) */

L130:
    limit1 = kp1 - *ns;
    temp5 = alpha[*ns];
    i__1 = limit1;
    for (iq = 1; iq <= i__1; ++iq) {
	v[iq] -= temp5 * v[iq + 1];
/* L135: */
	w[iq] = v[iq];
    }
    g[nsp1] = w[1];
    if (limit1 == 1) {
	goto L137;
    }
    *kgi = *ns;
    gi[*kgi] = w[2];
L137:
    w[limit1 + 1] = v[limit1 + 1];
    if (*k >= *kold) {
	goto L140;
    }
    ++(*ivc);
    iv[*ivc] = limit1 + 2;

/*   COMPUTE THE G(*) IN THE WORK VECTOR W(*) */

L140:
    nsp2 = *ns + 2;
    *kprev = *k;
    if (kp1 < nsp2) {
	goto L199;
    }
    i__1 = kp1;
    for (i__ = nsp2; i__ <= i__1; ++i__) {
	limit2 = kp2 - i__;
	temp6 = alpha[i__ - 1];
	i__2 = limit2;
	for (iq = 1; iq <= i__2; ++iq) {
/* L145: */
	    w[iq] -= temp6 * w[iq + 1];
	}
/* L150: */
	g[i__] = w[1];
    }
L199:
/*       ***     END BLOCK 1     *** */

/*       ***     BEGIN BLOCK 2     *** */
/*   PREDICT A SOLUTION P(*), EVALUATE DERIVATIVES USING PREDICTED */
/*   SOLUTION, ESTIMATE LOCAL ERROR AT ORDER K AND ERRORS AT ORDERS K, */
/*   K-1, K-2 AS IF CONSTANT STEP SIZE WERE USED. */
/*                   *** */

/*   INCREMENT COUNTER ON ATTEMPTED STEPS */

    ++(*ksteps);

/*   CHANGE PHI TO PHI STAR */

    if (*k < nsp1) {
	goto L215;
    }
    i__1 = *k;
    for (i__ = nsp1; i__ <= i__1; ++i__) {
	temp1 = beta[i__];
	i__2 = *neqn;
	for (l = 1; l <= i__2; ++l) {
/* L205: */
	    phi[l + i__ * phi_dim1] = temp1 * phi[l + i__ * phi_dim1];
	}
/* L210: */
    }

/*   PREDICT SOLUTION AND DIFFERENCES */

L215:
    i__1 = *neqn;
    for (l = 1; l <= i__1; ++l) {
	phi[l + kp2 * phi_dim1] = phi[l + kp1 * phi_dim1];
	phi[l + kp1 * phi_dim1] = 0.f;
/* L220: */
	p[l] = 0.f;
    }
    i__1 = *k;
    for (j = 1; j <= i__1; ++j) {
	i__ = kp1 - j;
	ip1 = i__ + 1;
	temp2 = g[i__];
	i__2 = *neqn;
	for (l = 1; l <= i__2; ++l) {
	    p[l] += temp2 * phi[l + i__ * phi_dim1];
/* L225: */
	    phi[l + i__ * phi_dim1] += phi[l + ip1 * phi_dim1];
	}
/* L230: */
    }
    if (*nornd) {
	goto L240;
    }
    i__1 = *neqn;
    for (l = 1; l <= i__1; ++l) {
	tau = *h__ * p[l] - phi[l + phi_dim1 * 15];
	p[l] = y[l] + tau;
/* L235: */
	phi[l + (phi_dim1 << 4)] = p[l] - y[l] - tau;
    }
    goto L250;
L240:
    i__1 = *neqn;
    for (l = 1; l <= i__1; ++l) {
/* L245: */
	p[l] = y[l] + *h__ * p[l];
    }
L250:
    *xold = *x;
    *x += *h__;
    absh = dabs(*h__);
    (*f)(x, &p[1], &yp[1], &rpar[1], &ipar[1]);

/*   ESTIMATE ERRORS AT ORDERS K,K-1,K-2 */

    erkm2 = 0.f;
    erkm1 = 0.f;
    erk = 0.f;
    i__1 = *neqn;
    for (l = 1; l <= i__1; ++l) {
	temp3 = 1.f / wt[l];
	temp4 = yp[l] - phi[l + phi_dim1];
	if (km2 < 0) {
	    goto L265;
	} else if (km2 == 0) {
	    goto L260;
	} else {
	    goto L255;
	}
L255:
/* Computing 2nd power */
	r__1 = (phi[l + km1 * phi_dim1] + temp4) * temp3;
	erkm2 += r__1 * r__1;
L260:
/* Computing 2nd power */
	r__1 = (phi[l + *k * phi_dim1] + temp4) * temp3;
	erkm1 += r__1 * r__1;
L265:
/* Computing 2nd power */
	r__1 = temp4 * temp3;
	erk += r__1 * r__1;
    }
    if (km2 < 0) {
	goto L280;
    } else if (km2 == 0) {
	goto L275;
    } else {
	goto L270;
    }
L270:
    erkm2 = absh * sig[km1] * gstr[km2 - 1] * sqrt(erkm2);
L275:
    erkm1 = absh * sig[*k] * gstr[km1 - 1] * sqrt(erkm1);
L280:
    temp5 = absh * sqrt(erk);
    err = temp5 * (g[*k] - g[kp1]);
    erk = temp5 * sig[kp1] * gstr[*k - 1];
    knew = *k;

/*   TEST IF ORDER SHOULD BE LOWERED */

    if (km2 < 0) {
	goto L299;
    } else if (km2 == 0) {
	goto L290;
    } else {
	goto L285;
    }
L285:
    if (dmax(erkm1,erkm2) <= erk) {
	knew = km1;
    }
    goto L299;
L290:
    if (erkm1 <= erk * .5f) {
	knew = km1;
    }

/*   TEST IF STEP SUCCESSFUL */

L299:
    if (err <= *eps) {
	goto L400;
    }
/*       ***     END BLOCK 2     *** */

/*       ***     BEGIN BLOCK 3     *** */
/*   THE STEP IS UNSUCCESSFUL.  RESTORE  X, PHI(*,*), PSI(*) . */
/*   IF THIRD CONSECUTIVE FAILURE, SET ORDER TO ONE.  IF STEP FAILS MORE */
/*   THAN THREE TIMES, CONSIDER AN OPTIMAL STEP SIZE.  DOUBLE ERROR */
/*   TOLERANCE AND RETURN IF ESTIMATED STEP SIZE IS TOO SMALL FOR MACHINE */
/*   PRECISION. */
/*                   *** */

/*   RESTORE X, PHI(*,*) AND PSI(*) */

    *phase1 = FALSE_;
    *x = *xold;
    i__1 = *k;
    for (i__ = 1; i__ <= i__1; ++i__) {
	temp1 = 1.f / beta[i__];
	ip1 = i__ + 1;
	i__2 = *neqn;
	for (l = 1; l <= i__2; ++l) {
/* L305: */
	    phi[l + i__ * phi_dim1] = temp1 * (phi[l + i__ * phi_dim1] - phi[
		    l + ip1 * phi_dim1]);
	}
/* L310: */
    }
    if (*k < 2) {
	goto L320;
    }
    i__1 = *k;
    for (i__ = 2; i__ <= i__1; ++i__) {
/* L315: */
	psi[i__ - 1] = psi[i__] - *h__;
    }

/*   ON THIRD FAILURE, SET ORDER TO ONE.  THEREAFTER, USE OPTIMAL STEP */
/*   SIZE */

L320:
    ++ifail;
    temp2 = .5f;
    if ((i__1 = ifail - 3) < 0) {
	goto L335;
    } else if (i__1 == 0) {
	goto L330;
    } else {
	goto L325;
    }
L325:
    if (p5eps < erk * .25f) {
	temp2 = sqrt(p5eps / erk);
    }
L330:
    knew = 1;
L335:
    *h__ = temp2 * *h__;
    *k = knew;
    *ns = 0;
    if (dabs(*h__) >= *fouru * dabs(*x)) {
	goto L340;
    }
    *crash = TRUE_;
    r__1 = *fouru * dabs(*x);
    *h__ = r_sign(&r__1, h__);
    *eps += *eps;
    return 0;
L340:
    goto L100;
/*       ***     END BLOCK 3     *** */

/*       ***     BEGIN BLOCK 4     *** */
/*   THE STEP IS SUCCESSFUL.  CORRECT THE PREDICTED SOLUTION, EVALUATE */
/*   THE DERIVATIVES USING THE CORRECTED SOLUTION AND UPDATE THE */
/*   DIFFERENCES.  DETERMINE BEST ORDER AND STEP SIZE FOR NEXT STEP. */
/*                   *** */
L400:
    *kold = *k;
    *hold = *h__;

/*   CORRECT AND EVALUATE */

    temp1 = *h__ * g[kp1];
    if (*nornd) {
	goto L410;
    }
    i__1 = *neqn;
    for (l = 1; l <= i__1; ++l) {
	temp3 = y[l];
	rho = temp1 * (yp[l] - phi[l + phi_dim1]) - phi[l + (phi_dim1 << 4)];
	y[l] = p[l] + rho;
	phi[l + phi_dim1 * 15] = y[l] - p[l] - rho;
/* L405: */
	p[l] = temp3;
    }
    goto L420;
L410:
    i__1 = *neqn;
    for (l = 1; l <= i__1; ++l) {
	temp3 = y[l];
	y[l] = p[l] + temp1 * (yp[l] - phi[l + phi_dim1]);
/* L415: */
	p[l] = temp3;
    }
L420:
    (*f)(x, &y[1], &yp[1], &rpar[1], &ipar[1]);

/*   UPDATE DIFFERENCES FOR NEXT STEP */

    i__1 = *neqn;
    for (l = 1; l <= i__1; ++l) {
	phi[l + kp1 * phi_dim1] = yp[l] - phi[l + phi_dim1];
/* L425: */
	phi[l + kp2 * phi_dim1] = phi[l + kp1 * phi_dim1] - phi[l + kp2 * 
		phi_dim1];
    }
    i__1 = *k;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = *neqn;
	for (l = 1; l <= i__2; ++l) {
/* L430: */
	    phi[l + i__ * phi_dim1] += phi[l + kp1 * phi_dim1];
	}
/* L435: */
    }

/*   ESTIMATE ERROR AT ORDER K+1 UNLESS: */
/*     IN FIRST PHASE WHEN ALWAYS RAISE ORDER, */
/*     ALREADY DECIDED TO LOWER ORDER, */
/*     STEP SIZE NOT CONSTANT SO ESTIMATE UNRELIABLE */

    erkp1 = 0.f;
    if (knew == km1 || *k == 12) {
	*phase1 = FALSE_;
    }
    if (*phase1) {
	goto L450;
    }
    if (knew == km1) {
	goto L455;
    }
    if (kp1 > *ns) {
	goto L460;
    }
    i__1 = *neqn;
    for (l = 1; l <= i__1; ++l) {
/* L440: */
/* Computing 2nd power */
	r__1 = phi[l + kp2 * phi_dim1] / wt[l];
	erkp1 += r__1 * r__1;
    }
    erkp1 = absh * gstr[kp1 - 1] * sqrt(erkp1);

/*   USING ESTIMATED ERROR AT ORDER K+1, DETERMINE APPROPRIATE ORDER */
/*   FOR NEXT STEP */

    if (*k > 1) {
	goto L445;
    }
    if (erkp1 >= erk * .5f) {
	goto L460;
    }
    goto L450;
L445:
    if (erkm1 <= dmin(erk,erkp1)) {
	goto L455;
    }
    if (erkp1 >= erk || *k == 12) {
	goto L460;
    }

/*   HERE ERKP1 .LT. ERK .LT. MAX(ERKM1,ERKM2) ELSE ORDER WOULD HAVE */
/*   BEEN LOWERED IN BLOCK 2.  THUS ORDER IS TO BE RAISED */

/*   RAISE ORDER */

L450:
    *k = kp1;
    erk = erkp1;
    goto L460;

/*   LOWER ORDER */

L455:
    *k = km1;
    erk = erkm1;

/*   WITH NEW ORDER DETERMINE APPROPRIATE STEP SIZE FOR NEXT STEP */

L460:
    hnew = *h__ + *h__;
    if (*phase1) {
	goto L465;
    }
    if (p5eps >= erk * two[*k]) {
	goto L465;
    }
    hnew = *h__;
    if (p5eps >= erk) {
	goto L465;
    }
    temp2 = (real) (*k + 1);
    d__1 = (doublereal) (p5eps / erk);
    d__2 = (doublereal) (1.f / temp2);
    r__ = pow_dd(&d__1, &d__2);
/* Computing MAX */
    r__1 = .5f, r__2 = dmin(.9f,r__);
    hnew = absh * dmax(r__1,r__2);
/* Computing MAX */
    r__2 = hnew, r__3 = *fouru * dabs(*x);
    r__1 = dmax(r__2,r__3);
    hnew = r_sign(&r__1, h__);
L465:
    *h__ = hnew;
    return 0;
/*       ***     END BLOCK 4     *** */
} /* steps_ */

