/* ddstp.f -- translated by f2c (version 12.02.01).
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

static integer c__1 = 1;
static integer c_n1 = -1;
static doublereal c_b22 = .2;

/* DECK DDSTP */
/* Subroutine */ int ddstp_(doublereal *eps, S_fp f, U_fp fa, doublereal *
	hmax, integer *impl, integer *ierror, U_fp jacobn, integer *matdim, 
	integer *maxord, integer *mint, integer *miter, integer *ml, integer *
	mu, integer *n, integer *nde, doublereal *ywt, doublereal *uround, 
	U_fp users, doublereal *avgh, doublereal *avgord, doublereal *h__, 
	doublereal *hused, integer *jtask, integer *mntold, integer *mtrold, 
	integer *nfe, integer *nje, integer *nqused, integer *nstep, 
	doublereal *t, doublereal *y, doublereal *yh, doublereal *a, logical *
	convrg, doublereal *dfdy, doublereal *el, doublereal *fac, doublereal 
	*hold, integer *ipvt, integer *jstate, integer *jstepl, integer *nq, 
	integer *nwait, doublereal *rc, doublereal *rmax, doublereal *save1, 
	doublereal *save2, doublereal *tq, doublereal *trend, integer *iswflg,
	 integer *mtrsv, integer *mxrdsv)
{
    /* Initialized data */

    static logical ier = FALSE_;

    /* System generated locals */
    integer a_dim1, a_offset, dfdy_dim1, dfdy_offset, yh_dim1, yh_offset, 
	    i__1, i__2;
    doublereal d__1, d__2, d__3, d__4;

    /* Local variables */
    static doublereal d__;
    static integer i__, j;
    static doublereal d1, hn, rh, hs, rh1, rh2, rh3, bnd;
    static integer nsv;
    static doublereal erdn, told;
    static integer iter;
    static doublereal erup;
    static integer ntry;
    extern doublereal dnrm2_(integer *, doublereal *, integer *);
    static doublereal y0nrm;
    static integer nfail;
    extern /* Subroutine */ int ddscl_(doublereal *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *), ddcor_(doublereal *, doublereal *, U_fp, 
	    doublereal *, integer *, integer *, integer *, integer *, integer 
	    *, integer *, integer *, integer *, integer *, integer *, 
	    doublereal *, U_fp, doublereal *, doublereal *, doublereal *, 
	    logical *, doublereal *, doublereal *, doublereal *, doublereal *,
	     integer *), ddpsc_(integer *, integer *, integer *, doublereal *)
	    , ddcst_(integer *, integer *, integer *, doublereal *, 
	    doublereal *);
    static doublereal denom;
    extern /* Subroutine */ int ddntl_(doublereal *, S_fp, U_fp, doublereal *,
	     doublereal *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, U_fp, doublereal *, 
	    doublereal *, doublereal *, integer *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, logical *, doublereal *,
	     doublereal *, logical *, integer *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, integer *, integer *), ddpst_(doublereal *, S_fp, 
	    U_fp, doublereal *, integer *, U_fp, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, doublereal 
	    *, doublereal *, U_fp, doublereal *, doublereal *, doublereal *, 
	    doublereal *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, logical *, integer *, doublereal *, integer *, 
	    doublereal *, integer *);
    static doublereal ctest, etest, numer;
    static logical evalfa, evaljc, switch__;

/* ***BEGIN PROLOGUE  DDSTP */
/* ***SUBSIDIARY */
/* ***PURPOSE  DDSTP performs one step of the integration of an initial */
/*            value problem for a system of ordinary differential */
/*            equations. */
/* ***LIBRARY   SLATEC (SDRIVE) */
/* ***TYPE      DOUBLE PRECISION (SDSTP-S, DDSTP-D, CDSTP-C) */
/* ***AUTHOR  Kahaner, D. K., (NIST) */
/*             National Institute of Standards and Technology */
/*             Gaithersburg, MD  20899 */
/*           Sutherland, C. D., (LANL) */
/*             Mail Stop D466 */
/*             Los Alamos National Laboratory */
/*             Los Alamos, NM  87545 */
/* ***DESCRIPTION */

/*  Communication with DDSTP is done with the following variables: */

/*    YH      An N by MAXORD+1 array containing the dependent variables */
/*              and their scaled derivatives.  MAXORD, the maximum order */
/*              used, is currently 12 for the Adams methods and 5 for the */
/*              Gear methods.  YH(I,J+1) contains the J-th derivative of */
/*              Y(I), scaled by H**J/factorial(J).  Only Y(I), */
/*              1 .LE. I .LE. N, need be set by the calling program on */
/*              the first entry.  The YH array should not be altered by */
/*              the calling program.  When referencing YH as a */
/*              2-dimensional array, use a column length of N, as this is */
/*              the value used in DDSTP. */
/*    DFDY    A block of locations used for partial derivatives if MITER */
/*              is not 0.  If MITER is 1 or 2 its length must be at least */
/*              N*N.  If MITER is 4 or 5 its length must be at least */
/*              (2*ML+MU+1)*N. */
/*    YWT     An array of N locations used in convergence and error tests */
/*    SAVE1 */
/*    SAVE2   Arrays of length N used for temporary storage. */
/*    IPVT    An integer array of length N used by the linear system */
/*              solvers for the storage of row interchange information. */
/*    A       A block of locations used to store the matrix A, when using */
/*              the implicit method.  If IMPL is 1, A is a MATDIM by N */
/*              array.  If MITER is 1 or 2 MATDIM is N, and if MITER is 4 */
/*              or 5 MATDIM is 2*ML+MU+1.  If IMPL is 2 its length is N. */
/*              If IMPL is 3, A is a MATDIM by NDE array. */
/*    JTASK   An integer used on input. */
/*              It has the following values and meanings: */
/*                 .EQ. 0  Perform the first step.  This value enables */
/*                         the subroutine to initialize itself. */
/*                .GT. 0  Take a new step continuing from the last. */
/*                         Assumes the last step was successful and */
/*                         user has not changed any parameters. */
/*                 .LT. 0  Take a new step with a new value of H and/or */
/*                         MINT and/or MITER. */
/*    JSTATE  A completion code with the following meanings: */
/*                1  The step was successful. */
/*                2  A solution could not be obtained with H .NE. 0. */
/*                3  A solution was not obtained in MXTRY attempts. */
/*                4  For IMPL .NE. 0, the matrix A is singular. */
/*              On a return with JSTATE .GT. 1, the values of T and */
/*              the YH array are as of the beginning of the last */
/*              step, and H is the last step size attempted. */

/* ***ROUTINES CALLED  DDCOR, DDCST, DDNTL, DDPSC, DDPST, DDSCL, DNRM2 */
/* ***REVISION HISTORY  (YYMMDD) */
/*   790601  DATE WRITTEN */
/*   900329  Initial submission to SLATEC. */
/* ***END PROLOGUE  DDSTP */
    /* Parameter adjustments */
    dfdy_dim1 = *matdim;
    dfdy_offset = 1 + dfdy_dim1;
    dfdy -= dfdy_offset;
    a_dim1 = *matdim;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    yh_dim1 = *n;
    yh_offset = 1 + yh_dim1;
    yh -= yh_offset;
    --ywt;
    --y;
    el -= 14;
    --fac;
    --ipvt;
    --save1;
    --save2;
    tq -= 4;

    /* Function Body */
/* ***FIRST EXECUTABLE STATEMENT  DDSTP */
    nsv = *n;
    bnd = 0.;
    switch__ = FALSE_;
    ntry = 0;
    told = *t;
    nfail = 0;
    if (*jtask <= 0) {
	ddntl_(eps, (S_fp)f, (U_fp)fa, hmax, hold, impl, jtask, matdim, 
		maxord, mint, miter, ml, mu, n, nde, &save1[1], t, uround, (
		U_fp)users, &y[1], &ywt[1], h__, mntold, mtrold, nfe, rc, &yh[
		yh_offset], &a[a_offset], convrg, &el[14], &fac[1], &ier, &
		ipvt[1], nq, nwait, &rh, rmax, &save2[1], &tq[4], trend, 
		iswflg, jstate);
	if (*n == 0) {
	    goto L440;
	}
	if (*h__ == 0.) {
	    goto L400;
	}
	if (ier) {
	    goto L420;
	}
    }
L100:
    ++ntry;
    if (ntry > 50) {
	goto L410;
    }
    *t += *h__;
    ddpsc_(&c__1, n, nq, &yh[yh_offset]);
    evaljc = ((d__1 = *rc - 1., abs(d__1)) > .3 || *nstep >= *jstepl + 10) && 
	    *miter != 0;
    evalfa = ! evaljc;

L110:
    iter = 0;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L115: */
	y[i__] = yh[i__ + yh_dim1];
    }
    (*f)(n, t, &y[1], &save2[1]);
    if (*n == 0) {
	*jstate = 6;
	goto L430;
    }
    ++(*nfe);
    if (evaljc || ier) {
	ddpst_(&el[14], (S_fp)f, (U_fp)fa, h__, impl, (U_fp)jacobn, matdim, 
		miter, ml, mu, n, nde, nq, &save2[1], t, (U_fp)users, &y[1], &
		yh[yh_offset], &ywt[1], uround, nfe, nje, &a[a_offset], &dfdy[
		dfdy_offset], &fac[1], &ier, &ipvt[1], &save1[1], iswflg, &
		bnd, jstate);
	if (*n == 0) {
	    goto L430;
	}
	if (ier) {
	    goto L160;
	}
	*convrg = FALSE_;
	*rc = 1.;
	*jstepl = *nstep;
    }
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L125: */
	save1[i__] = 0.;
    }
/*                      Up to MXITER corrector iterations are taken. */
/*                      Convergence is tested by requiring the r.m.s. */
/*                      norm of changes to be less than EPS.  The sum of */
/*                      the corrections is accumulated in the vector */
/*                      SAVE1(I).  It is approximately equal to the L-th */
/*                      derivative of Y multiplied by */
/*                      H**L/(factorial(L-1)*EL(L,NQ)), and is thus */
/*                      proportional to the actual errors to the lowest */
/*                      power of H present (H**L).  The YH array is not */
/*                      altered in the correction loop.  The norm of the */
/*                      iterate difference is stored in D.  If */
/*                      ITER .GT. 0, an estimate of the convergence rate */
/*                      constant is stored in TREND, and this is used in */
/*                      the convergence test. */

L130:
    ddcor_(&dfdy[dfdy_offset], &el[14], (U_fp)fa, h__, ierror, impl, &ipvt[1],
	     matdim, miter, ml, mu, n, nde, nq, t, (U_fp)users, &y[1], &yh[
	    yh_offset], &ywt[1], &evalfa, &save1[1], &save2[1], &a[a_offset], 
	    &d__, jstate);
    if (*n == 0) {
	goto L430;
    }
    if (*iswflg == 3 && *mint == 1) {
	if (iter == 0) {
	    numer = dnrm2_(n, &save1[1], &c__1);
	    i__1 = *n;
	    for (i__ = 1; i__ <= i__1; ++i__) {
/* L132: */
		dfdy[i__ * dfdy_dim1 + 1] = save1[i__];
	    }
	    y0nrm = dnrm2_(n, &yh[yh_offset], &c__1);
	} else {
	    denom = numer;
	    i__1 = *n;
	    for (i__ = 1; i__ <= i__1; ++i__) {
/* L134: */
		dfdy[i__ * dfdy_dim1 + 1] = save1[i__] - dfdy[i__ * dfdy_dim1 
			+ 1];
	    }
	    numer = dnrm2_(n, &dfdy[dfdy_offset], matdim);
	    if (el[*nq * 13 + 1] * numer <= *uround * 100. * y0nrm) {
		if (*rmax == 2.) {
		    switch__ = TRUE_;
		    goto L170;
		}
	    }
	    i__1 = *n;
	    for (i__ = 1; i__ <= i__1; ++i__) {
/* L136: */
		dfdy[i__ * dfdy_dim1 + 1] = save1[i__];
	    }
	    if (denom != 0.) {
/* Computing MAX */
		d__1 = bnd, d__2 = numer / (denom * abs(*h__) * el[*nq * 13 + 
			1]);
		bnd = max(d__1,d__2);
	    }
	}
    }
    if (iter > 0) {
/* Computing MAX */
	d__1 = *trend * .9, d__2 = d__ / d1;
	*trend = max(d__1,d__2);
    }
    d1 = d__;
/* Computing MIN */
    d__1 = *trend * 2.;
    ctest = min(d__1,1.) * d__;
    if (ctest <= *eps) {
	goto L170;
    }
    ++iter;
    if (iter < 3) {
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
/* L140: */
	    y[i__] = yh[i__ + yh_dim1] + el[*nq * 13 + 1] * save1[i__];
	}
	(*f)(n, t, &y[1], &save2[1]);
	if (*n == 0) {
	    *jstate = 6;
	    goto L430;
	}
	++(*nfe);
	goto L130;
    }
/*                     The corrector iteration failed to converge in */
/*                     MXITER tries.  If partials are involved but are */
/*                     not up to date, they are reevaluated for the next */
/*                     try.  Otherwise the YH array is retracted to its */
/*                     values before prediction, and H is reduced, if */
/*                     possible.  If not, a no-convergence exit is taken. */
    if (*convrg) {
	evaljc = TRUE_;
	evalfa = FALSE_;
	goto L110;
    }
L160:
    *t = told;
    ddpsc_(&c_n1, n, nq, &yh[yh_offset]);
    *nwait = *nq + 2;
    if (*jtask != 0 && *jtask != 2) {
	*rmax = 2.;
    }
    if (iter == 0) {
	rh = .3;
    } else {
	d__1 = *eps / ctest;
	rh = pow_dd(&d__1, &c_b22) * .9;
    }
    if (rh * *h__ == 0.) {
	goto L400;
    }
    ddscl_(hmax, n, nq, rmax, h__, rc, &rh, &yh[yh_offset]);
    goto L100;
/*                          The corrector has converged.  CONVRG is set */
/*                          to .TRUE. if partial derivatives were used, */
/*                          to indicate that they may need updating on */
/*                          subsequent steps.  The error test is made. */
L170:
    *convrg = *miter != 0;
    if (*ierror == 1 || *ierror == 5) {
	i__1 = *nde;
	for (i__ = 1; i__ <= i__1; ++i__) {
/* L180: */
	    save2[i__] = save1[i__] / ywt[i__];
	}
    } else {
	i__1 = *nde;
	for (i__ = 1; i__ <= i__1; ++i__) {
/* L185: */
/* Computing MAX */
	    d__2 = (d__1 = y[i__], abs(d__1)), d__3 = ywt[i__];
	    save2[i__] = save1[i__] / max(d__2,d__3);
	}
    }
    etest = dnrm2_(nde, &save2[1], &c__1) / (tq[*nq * 3 + 2] * sqrt((
	    doublereal) (*nde)));

/*                           The error test failed.  NFAIL keeps track of */
/*                           multiple failures.  Restore T and the YH */
/*                           array to their previous values, and prepare */
/*                           to try the step again.  Compute the optimum */
/*                           step size for this or one lower order. */
    if (etest > *eps) {
	*t = told;
	ddpsc_(&c_n1, n, nq, &yh[yh_offset]);
	++nfail;
	if (nfail < 3 || *nq == 1) {
	    if (*jtask != 0 && *jtask != 2) {
		*rmax = 2.;
	    }
	    d__1 = etest / *eps;
	    d__2 = 1. / (*nq + 1);
	    rh2 = 1. / (pow_dd(&d__1, &d__2) * 1.2);
	    if (*nq > 1) {
		if (*ierror == 1 || *ierror == 5) {
		    i__1 = *nde;
		    for (i__ = 1; i__ <= i__1; ++i__) {
/* L190: */
			save2[i__] = yh[i__ + (*nq + 1) * yh_dim1] / ywt[i__];
		    }
		} else {
		    i__1 = *nde;
		    for (i__ = 1; i__ <= i__1; ++i__) {
/* L195: */
/* Computing MAX */
			d__2 = (d__1 = y[i__], abs(d__1)), d__3 = ywt[i__];
			save2[i__] = yh[i__ + (*nq + 1) * yh_dim1] / max(d__2,
				d__3);
		    }
		}
		erdn = dnrm2_(nde, &save2[1], &c__1) / (tq[*nq * 3 + 1] * 
			sqrt((doublereal) (*nde)));
/* Computing MAX */
		d__3 = erdn / *eps;
		d__4 = 1. / *nq;
		d__1 = 1., d__2 = pow_dd(&d__3, &d__4) * 1.3;
		rh1 = 1. / max(d__1,d__2);
		if (rh2 < rh1) {
		    --(*nq);
		    *rc = *rc * el[*nq * 13 + 1] / el[(*nq + 1) * 13 + 1];
		    rh = rh1;
		} else {
		    rh = rh2;
		}
	    } else {
		rh = rh2;
	    }
	    *nwait = *nq + 2;
	    if (rh * *h__ == 0.) {
		goto L400;
	    }
	    ddscl_(hmax, n, nq, rmax, h__, rc, &rh, &yh[yh_offset]);
	    goto L100;
	}
/*                Control reaches this section if the error test has */
/*                failed MXFAIL or more times.  It is assumed that the */
/*                derivatives that have accumulated in the YH array have */
/*                errors of the wrong order.  Hence the first derivative */
/*                is recomputed, the order is set to 1, and the step is */
/*                retried. */
	nfail = 0;
	*jtask = 2;
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
/* L215: */
	    y[i__] = yh[i__ + yh_dim1];
	}
	ddntl_(eps, (S_fp)f, (U_fp)fa, hmax, hold, impl, jtask, matdim, 
		maxord, mint, miter, ml, mu, n, nde, &save1[1], t, uround, (
		U_fp)users, &y[1], &ywt[1], h__, mntold, mtrold, nfe, rc, &yh[
		yh_offset], &a[a_offset], convrg, &el[14], &fac[1], &ier, &
		ipvt[1], nq, nwait, &rh, rmax, &save2[1], &tq[4], trend, 
		iswflg, jstate);
	*rmax = 10.;
	if (*n == 0) {
	    goto L440;
	}
	if (*h__ == 0.) {
	    goto L400;
	}
	if (ier) {
	    goto L420;
	}
	goto L100;
    }
/*                          After a successful step, update the YH array. */
    ++(*nstep);
    *hused = *h__;
    *nqused = *nq;
    *avgh = ((*nstep - 1) * *avgh + *h__) / *nstep;
    *avgord = ((*nstep - 1) * *avgord + *nq) / *nstep;
    i__1 = *nq + 1;
    for (j = 1; j <= i__1; ++j) {
	i__2 = *n;
	for (i__ = 1; i__ <= i__2; ++i__) {
/* L230: */
	    yh[i__ + j * yh_dim1] += el[j + *nq * 13] * save1[i__];
	}
    }
    i__2 = *n;
    for (i__ = 1; i__ <= i__2; ++i__) {
/* L235: */
	y[i__] = yh[i__ + yh_dim1];
    }
/*                                          If ISWFLG is 3, consider */
/*                                          changing integration methods. */
    if (*iswflg == 3) {
	if (bnd != 0.) {
	    if (*mint == 1 && *nq <= 5) {
/* Computing MAX */
		d__3 = etest / *eps;
		d__4 = 1. / (*nq + 1);
		d__1 = *uround, d__2 = pow_dd(&d__3, &d__4);
		hn = abs(*h__) / max(d__1,d__2);
/* Computing MIN */
		d__1 = hn, d__2 = 1. / (el[*nq * 13 + 1] * 2. * bnd);
		hn = min(d__1,d__2);
/* Computing MAX */
		d__3 = etest / (*eps * el[*nq + 14]);
		d__4 = 1. / (*nq + 1);
		d__1 = *uround, d__2 = pow_dd(&d__3, &d__4);
		hs = abs(*h__) / max(d__1,d__2);
		if (hs > hn * 1.2) {
		    *mint = 2;
		    *mntold = *mint;
		    *miter = *mtrsv;
		    *mtrold = *miter;
		    *maxord = min(*mxrdsv,5);
		    *rc = 0.;
		    *rmax = 10.;
		    *trend = 1.;
		    ddcst_(maxord, mint, iswflg, &el[14], &tq[4]);
		    *nwait = *nq + 2;
		}
	    } else if (*mint == 2) {
/* Computing MAX */
		d__3 = etest / *eps;
		d__4 = 1. / (*nq + 1);
		d__1 = *uround, d__2 = pow_dd(&d__3, &d__4);
		hs = abs(*h__) / max(d__1,d__2);
/* Computing MAX */
		d__3 = etest * el[*nq + 14] / *eps;
		d__4 = 1. / (*nq + 1);
		d__1 = *uround, d__2 = pow_dd(&d__3, &d__4);
		hn = abs(*h__) / max(d__1,d__2);
/* Computing MIN */
		d__1 = hn, d__2 = 1. / (el[*nq * 13 + 1] * 2. * bnd);
		hn = min(d__1,d__2);
		if (hn >= hs) {
		    *mint = 1;
		    *mntold = *mint;
		    *miter = 0;
		    *mtrold = *miter;
		    *maxord = min(*mxrdsv,12);
		    *rmax = 10.;
		    *trend = 1.;
		    *convrg = FALSE_;
		    ddcst_(maxord, mint, iswflg, &el[14], &tq[4]);
		    *nwait = *nq + 2;
		}
	    }
	}
    }
    if (switch__) {
	*mint = 2;
	*mntold = *mint;
	*miter = *mtrsv;
	*mtrold = *miter;
	*maxord = min(*mxrdsv,5);
	*nq = min(*nq,*maxord);
	*rc = 0.;
	*rmax = 10.;
	*trend = 1.;
	ddcst_(maxord, mint, iswflg, &el[14], &tq[4]);
	*nwait = *nq + 2;
    }
/*                           Consider changing H if NWAIT = 1.  Otherwise */
/*                           decrease NWAIT by 1.  If NWAIT is then 1 and */
/*                           NQ.LT.MAXORD, then SAVE1 is saved for use in */
/*                           a possible order increase on the next step. */

    if (*jtask == 0 || *jtask == 2) {
/* Computing MAX */
	d__3 = etest / *eps;
	d__4 = 1. / (*nq + 1);
	d__1 = *uround, d__2 = pow_dd(&d__3, &d__4) * 1.2;
	rh = 1. / max(d__1,d__2);
	if (rh > 1.) {
	    ddscl_(hmax, n, nq, rmax, h__, rc, &rh, &yh[yh_offset]);
	}
    } else if (*nwait > 1) {
	--(*nwait);
	if (*nwait == 1 && *nq < *maxord) {
	    i__2 = *nde;
	    for (i__ = 1; i__ <= i__2; ++i__) {
/* L250: */
		yh[i__ + (*maxord + 1) * yh_dim1] = save1[i__];
	    }
	}
/*             If a change in H is considered, an increase or decrease in */
/*             order by one is considered also.  A change in H is made */
/*             only if it is by a factor of at least TRSHLD.  Factors */
/*             RH1, RH2, and RH3 are computed, by which H could be */
/*             multiplied at order NQ - 1, order NQ, or order NQ + 1, */
/*             respectively.  The largest of these is determined and the */
/*             new order chosen accordingly.  If the order is to be */
/*             increased, we compute one additional scaled derivative. */
/*             If there is a change of order, reset NQ and the */
/*             coefficients.  In any case H is reset according to RH and */
/*             the YH array is rescaled. */
    } else {
	if (*nq == 1) {
	    rh1 = 0.;
	} else {
	    if (*ierror == 1 || *ierror == 5) {
		i__2 = *nde;
		for (i__ = 1; i__ <= i__2; ++i__) {
/* L270: */
		    save2[i__] = yh[i__ + (*nq + 1) * yh_dim1] / ywt[i__];
		}
	    } else {
		i__2 = *nde;
		for (i__ = 1; i__ <= i__2; ++i__) {
/* L275: */
/* Computing MAX */
		    d__2 = (d__1 = y[i__], abs(d__1)), d__3 = ywt[i__];
		    save2[i__] = yh[i__ + (*nq + 1) * yh_dim1] / max(d__2,
			    d__3);
		}
	    }
	    erdn = dnrm2_(nde, &save2[1], &c__1) / (tq[*nq * 3 + 1] * sqrt((
		    doublereal) (*nde)));
/* Computing MAX */
	    d__3 = erdn / *eps;
	    d__4 = 1. / *nq;
	    d__1 = *uround, d__2 = pow_dd(&d__3, &d__4) * 1.3;
	    rh1 = 1. / max(d__1,d__2);
	}
/* Computing MAX */
	d__3 = etest / *eps;
	d__4 = 1. / (*nq + 1);
	d__1 = *uround, d__2 = pow_dd(&d__3, &d__4) * 1.2;
	rh2 = 1. / max(d__1,d__2);
	if (*nq == *maxord) {
	    rh3 = 0.;
	} else {
	    if (*ierror == 1 || *ierror == 5) {
		i__2 = *nde;
		for (i__ = 1; i__ <= i__2; ++i__) {
/* L290: */
		    save2[i__] = (save1[i__] - yh[i__ + (*maxord + 1) * 
			    yh_dim1]) / ywt[i__];
		}
	    } else {
		i__2 = *nde;
		for (i__ = 1; i__ <= i__2; ++i__) {
/* Computing MAX */
		    d__2 = (d__1 = y[i__], abs(d__1)), d__3 = ywt[i__];
		    save2[i__] = (save1[i__] - yh[i__ + (*maxord + 1) * 
			    yh_dim1]) / max(d__2,d__3);
/* L295: */
		}
	    }
	    erup = dnrm2_(nde, &save2[1], &c__1) / (tq[*nq * 3 + 3] * sqrt((
		    doublereal) (*nde)));
/* Computing MAX */
	    d__3 = erup / *eps;
	    d__4 = 1. / (*nq + 2);
	    d__1 = *uround, d__2 = pow_dd(&d__3, &d__4) * 1.4;
	    rh3 = 1. / max(d__1,d__2);
	}
	if (rh1 > rh2 && rh1 >= rh3) {
	    rh = rh1;
	    if (rh <= 1.) {
		goto L380;
	    }
	    --(*nq);
	    *rc = *rc * el[*nq * 13 + 1] / el[(*nq + 1) * 13 + 1];
	} else if (rh2 >= rh1 && rh2 >= rh3) {
	    rh = rh2;
	    if (rh <= 1.) {
		goto L380;
	    }
	} else {
	    rh = rh3;
	    if (rh <= 1.) {
		goto L380;
	    }
	    i__2 = *n;
	    for (i__ = 1; i__ <= i__2; ++i__) {
/* L360: */
		yh[i__ + (*nq + 2) * yh_dim1] = save1[i__] * el[*nq + 1 + *nq 
			* 13] / (*nq + 1);
	    }
	    ++(*nq);
	    *rc = *rc * el[*nq * 13 + 1] / el[(*nq - 1) * 13 + 1];
	}
	if (*iswflg == 3 && *mint == 1) {
	    if (bnd != 0.) {
/* Computing MIN */
		d__1 = rh, d__2 = 1. / (el[*nq * 13 + 1] * 2. * bnd * abs(*
			h__));
		rh = min(d__1,d__2);
	    }
	}
	ddscl_(hmax, n, nq, rmax, h__, rc, &rh, &yh[yh_offset]);
	*rmax = 10.;
L380:
	*nwait = *nq + 2;
    }
/*               All returns are made through this section.  H is saved */
/*               in HOLD to allow the caller to change H on the next step */
    *jstate = 1;
    *hold = *h__;
    return 0;

L400:
    *jstate = 2;
    *hold = *h__;
    i__2 = *n;
    for (i__ = 1; i__ <= i__2; ++i__) {
/* L405: */
	y[i__] = yh[i__ + yh_dim1];
    }
    return 0;

L410:
    *jstate = 3;
    *hold = *h__;
    return 0;

L420:
    *jstate = 4;
    *hold = *h__;
    return 0;

L430:
    *t = told;
    ddpsc_(&c_n1, &nsv, nq, &yh[yh_offset]);
    i__2 = nsv;
    for (i__ = 1; i__ <= i__2; ++i__) {
/* L435: */
	y[i__] = yh[i__ + yh_dim1];
    }
L440:
    *hold = *h__;
    return 0;
} /* ddstp_ */

