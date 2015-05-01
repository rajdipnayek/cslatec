/* dhstrt.f -- translated by f2c (version 12.02.01).
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

static doublereal c_b2 = .375;
static doublereal c_b20 = 10.;

/* DECK DHSTRT */
/* Subroutine */ int dhstrt_(S_fp df, integer *neq, doublereal *a, doublereal 
	*b, doublereal *y, doublereal *yprime, doublereal *etol, integer *
	morder, doublereal *small, doublereal *big, doublereal *spy, 
	doublereal *pv, doublereal *yp, doublereal *sf, doublereal *rpar, 
	integer *ipar, doublereal *h__)
{
    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1, d__2, d__3, d__4;

    /* Local variables */
    static integer j, k;
    static doublereal da;
    static integer lk;
    static doublereal dx, dy, fbnd, delf, dely, ydpb, tolp, dfdub, dfdxb, 
	    absdx, relper;
    extern doublereal dhvnrm_(doublereal *, integer *);
    static doublereal tolmin, srydpb, tolexp, tolsum;

/* ***BEGIN PROLOGUE  DHSTRT */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to DDEABM, DDEBDF and DDERKF */
/* ***LIBRARY   SLATEC */
/* ***TYPE      DOUBLE PRECISION (HSTART-S, DHSTRT-D) */
/* ***AUTHOR  Watts, H. A., (SNLA) */
/* ***DESCRIPTION */

/*   DHSTRT computes a starting step size to be used in solving initial */
/*   value problems in ordinary differential equations. */

/* ********************************************************************** */
/*  ABSTRACT */

/*     Subroutine DHSTRT computes a starting step size to be used by an */
/*     initial value method in solving ordinary differential equations. */
/*     It is based on an estimate of the local Lipschitz constant for the */
/*     differential equation   (lower bound on a norm of the Jacobian) , */
/*     a bound on the differential equation  (first derivative) , and */
/*     a bound on the partial derivative of the equation with respect to */
/*     the independent variable. */
/*     (all approximated near the initial point A) */

/*     Subroutine DHSTRT uses a function subprogram DHVNRM for computing */
/*     a vector norm. The maximum norm is presently utilized though it */
/*     can easily be replaced by any other vector norm. It is presumed */
/*     that any replacement norm routine would be carefully coded to */
/*     prevent unnecessary underflows or overflows from occurring, and */
/*     also, would not alter the vector or number of components. */

/* ********************************************************************** */
/*  On input you must provide the following */

/*      DF -- This is a subroutine of the form */
/*                               DF(X,U,UPRIME,RPAR,IPAR) */
/*             which defines the system of first order differential */
/*             equations to be solved. For the given values of X and the */
/*             vector  U(*)=(U(1),U(2),...,U(NEQ)) , the subroutine must */
/*             evaluate the NEQ components of the system of differential */
/*             equations  DU/DX=DF(X,U)  and store the derivatives in the */
/*             array UPRIME(*), that is,  UPRIME(I) = * DU(I)/DX *  for */
/*             equations I=1,...,NEQ. */

/*             Subroutine DF must not alter X or U(*). You must declare */
/*             the name DF in an external statement in your program that */
/*             calls DHSTRT. You must dimension U and UPRIME in DF. */

/*             RPAR and IPAR are DOUBLE PRECISION and INTEGER parameter */
/*             arrays which you can use for communication between your */
/*             program and subroutine DF. They are not used or altered by */
/*             DHSTRT. If you do not need RPAR or IPAR, ignore these */
/*             parameters by treating them as dummy arguments. If you do */
/*             choose to use them, dimension them in your program and in */
/*             DF as arrays of appropriate length. */

/*      NEQ -- This is the number of (first order) differential equations */
/*             to be integrated. */

/*      A -- This is the initial point of integration. */

/*      B -- This is a value of the independent variable used to define */
/*             the direction of integration. A reasonable choice is to */
/*             set  B  to the first point at which a solution is desired. */
/*             You can also use  B, if necessary, to restrict the length */
/*             of the first integration step because the algorithm will */
/*             not compute a starting step length which is bigger than */
/*             ABS(B-A), unless  B  has been chosen too close to  A. */
/*             (it is presumed that DHSTRT has been called with  B */
/*             different from  A  on the machine being used. Also see the */
/*             discussion about the parameter  SMALL.) */

/*      Y(*) -- This is the vector of initial values of the NEQ solution */
/*             components at the initial point  A. */

/*      YPRIME(*) -- This is the vector of derivatives of the NEQ */
/*             solution components at the initial point  A. */
/*             (defined by the differential equations in subroutine DF) */

/*      ETOL -- This is the vector of error tolerances corresponding to */
/*             the NEQ solution components. It is assumed that all */
/*             elements are positive. Following the first integration */
/*             step, the tolerances are expected to be used by the */
/*             integrator in an error test which roughly requires that */
/*                        ABS(LOCAL ERROR)  .LE.  ETOL */
/*             for each vector component. */

/*      MORDER -- This is the order of the formula which will be used by */
/*             the initial value method for taking the first integration */
/*             step. */

/*      SMALL -- This is a small positive machine dependent constant */
/*             which is used for protecting against computations with */
/*             numbers which are too small relative to the precision of */
/*             floating point arithmetic.  SMALL  should be set to */
/*             (approximately) the smallest positive DOUBLE PRECISION */
/*             number such that  (1.+SMALL) .GT. 1.  on the machine being */
/*             used. The quantity  SMALL**(3/8)  is used in computing */
/*             increments of variables for approximating derivatives by */
/*             differences.  Also the algorithm will not compute a */
/*             starting step length which is smaller than */
/*             100*SMALL*ABS(A). */

/*      BIG -- This is a large positive machine dependent constant which */
/*             is used for preventing machine overflows. A reasonable */
/*             choice is to set big to (approximately) the square root of */
/*             the largest DOUBLE PRECISION number which can be held in */
/*             the machine. */

/*      SPY(*),PV(*),YP(*),SF(*) -- These are DOUBLE PRECISION work */
/*             arrays of length NEQ which provide the routine with needed */
/*             storage space. */

/*      RPAR,IPAR -- These are parameter arrays, of DOUBLE PRECISION and */
/*             INTEGER type, respectively, which can be used for */
/*             communication between your program and the DF subroutine. */
/*             They are not used or altered by DHSTRT. */

/* ********************************************************************** */
/*  On Output  (after the return from DHSTRT), */

/*      H -- is an appropriate starting step size to be attempted by the */
/*             differential equation method. */

/*           All parameters in the call list remain unchanged except for */
/*           the working arrays SPY(*),PV(*),YP(*), and SF(*). */

/* ********************************************************************** */

/* ***SEE ALSO  DDEABM, DDEBDF, DDERKF */
/* ***ROUTINES CALLED  DHVNRM */
/* ***REVISION HISTORY  (YYMMDD) */
/*   820301  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890831  Modified array declarations.  (WRB) */
/*   890911  Removed unnecessary intrinsics.  (WRB) */
/*   891024  Changed references from DVNORM to DHVNRM.  (WRB) */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900328  Added TYPE section.  (WRB) */
/*   910722  Updated AUTHOR section.  (ALS) */
/* ***END PROLOGUE  DHSTRT */


/*     .................................................................. */

/*     BEGIN BLOCK PERMITTING ...EXITS TO 160 */
/* ***FIRST EXECUTABLE STATEMENT  DHSTRT */
    /* Parameter adjustments */
    --ipar;
    --rpar;
    --sf;
    --yp;
    --pv;
    --spy;
    --etol;
    --yprime;
    --y;

    /* Function Body */
    dx = *b - *a;
    absdx = abs(dx);
    relper = pow_dd(small, &c_b2);

/*        ............................................................... */

/*             COMPUTE AN APPROXIMATE BOUND (DFDXB) ON THE PARTIAL */
/*             DERIVATIVE OF THE EQUATION WITH RESPECT TO THE */
/*             INDEPENDENT VARIABLE. PROTECT AGAINST AN OVERFLOW. */
/*             ALSO COMPUTE A BOUND (FBND) ON THE FIRST DERIVATIVE */
/*             LOCALLY. */

/* Computing MAX */
/* Computing MIN */
    d__4 = relper * abs(*a);
    d__2 = min(d__4,absdx), d__3 = *small * 100. * abs(*a);
    d__1 = max(d__2,d__3);
    da = d_sign(&d__1, &dx);
    if (da == 0.) {
	da = relper * dx;
    }
    d__1 = *a + da;
    (*df)(&d__1, &y[1], &sf[1], &rpar[1], &ipar[1]);
    i__1 = *neq;
    for (j = 1; j <= i__1; ++j) {
	yp[j] = sf[j] - yprime[j];
/* L10: */
    }
    delf = dhvnrm_(&yp[1], neq);
    dfdxb = *big;
    if (delf < *big * abs(da)) {
	dfdxb = delf / abs(da);
    }
    fbnd = dhvnrm_(&sf[1], neq);

/*        ............................................................... */

/*             COMPUTE AN ESTIMATE (DFDUB) OF THE LOCAL LIPSCHITZ */
/*             CONSTANT FOR THE SYSTEM OF DIFFERENTIAL EQUATIONS. THIS */
/*             ALSO REPRESENTS AN ESTIMATE OF THE NORM OF THE JACOBIAN */
/*             LOCALLY.  THREE ITERATIONS (TWO WHEN NEQ=1) ARE USED TO */
/*             ESTIMATE THE LIPSCHITZ CONSTANT BY NUMERICAL DIFFERENCES. */
/*             THE FIRST PERTURBATION VECTOR IS BASED ON THE INITIAL */
/*             DERIVATIVES AND DIRECTION OF INTEGRATION. THE SECOND */
/*             PERTURBATION VECTOR IS FORMED USING ANOTHER EVALUATION OF */
/*             THE DIFFERENTIAL EQUATION.  THE THIRD PERTURBATION VECTOR */
/*             IS FORMED USING PERTURBATIONS BASED ONLY ON THE INITIAL */
/*             VALUES. COMPONENTS THAT ARE ZERO ARE ALWAYS CHANGED TO */
/*             NON-ZERO VALUES (EXCEPT ON THE FIRST ITERATION). WHEN */
/*             INFORMATION IS AVAILABLE, CARE IS TAKEN TO ENSURE THAT */
/*             COMPONENTS OF THE PERTURBATION VECTOR HAVE SIGNS WHICH ARE */
/*             CONSISTENT WITH THE SLOPES OF LOCAL SOLUTION CURVES. */
/*             ALSO CHOOSE THE LARGEST BOUND (FBND) FOR THE FIRST */
/*             DERIVATIVE. */

/*                               PERTURBATION VECTOR SIZE IS HELD */
/*                               CONSTANT FOR ALL ITERATIONS. COMPUTE */
/*                               THIS CHANGE FROM THE */
/*                                       SIZE OF THE VECTOR OF INITIAL */
/*                                       VALUES. */
    dely = relper * dhvnrm_(&y[1], neq);
    if (dely == 0.) {
	dely = relper;
    }
    dely = d_sign(&dely, &dx);
    delf = dhvnrm_(&yprime[1], neq);
    fbnd = max(fbnd,delf);
    if (delf == 0.) {
	goto L30;
    }
/*           USE INITIAL DERIVATIVES FOR FIRST PERTURBATION */
    i__1 = *neq;
    for (j = 1; j <= i__1; ++j) {
	spy[j] = yprime[j];
	yp[j] = yprime[j];
/* L20: */
    }
    goto L50;
L30:
/*           CANNOT HAVE A NULL PERTURBATION VECTOR */
    i__1 = *neq;
    for (j = 1; j <= i__1; ++j) {
	spy[j] = 0.;
	yp[j] = 1.;
/* L40: */
    }
    delf = dhvnrm_(&yp[1], neq);
L50:

    dfdub = 0.;
/* Computing MIN */
    i__1 = *neq + 1;
    lk = min(i__1,3);
    i__1 = lk;
    for (k = 1; k <= i__1; ++k) {
/*           DEFINE PERTURBED VECTOR OF INITIAL VALUES */
	i__2 = *neq;
	for (j = 1; j <= i__2; ++j) {
	    pv[j] = y[j] + dely * (yp[j] / delf);
/* L60: */
	}
	if (k == 2) {
	    goto L80;
	}
/*              EVALUATE DERIVATIVES ASSOCIATED WITH PERTURBED */
/*              VECTOR  AND  COMPUTE CORRESPONDING DIFFERENCES */
	(*df)(a, &pv[1], &yp[1], &rpar[1], &ipar[1]);
	i__2 = *neq;
	for (j = 1; j <= i__2; ++j) {
	    pv[j] = yp[j] - yprime[j];
/* L70: */
	}
	goto L100;
L80:
/*              USE A SHIFTED VALUE OF THE INDEPENDENT VARIABLE */
/*                                    IN COMPUTING ONE ESTIMATE */
	d__1 = *a + da;
	(*df)(&d__1, &pv[1], &yp[1], &rpar[1], &ipar[1]);
	i__2 = *neq;
	for (j = 1; j <= i__2; ++j) {
	    pv[j] = yp[j] - sf[j];
/* L90: */
	}
L100:
/*           CHOOSE LARGEST BOUNDS ON THE FIRST DERIVATIVE */
/*                          AND A LOCAL LIPSCHITZ CONSTANT */
/* Computing MAX */
	d__1 = fbnd, d__2 = dhvnrm_(&yp[1], neq);
	fbnd = max(d__1,d__2);
	delf = dhvnrm_(&pv[1], neq);
/*        ...EXIT */
	if (delf >= *big * abs(dely)) {
	    goto L150;
	}
/* Computing MAX */
	d__1 = dfdub, d__2 = delf / abs(dely);
	dfdub = max(d__1,d__2);
/*     ......EXIT */
	if (k == lk) {
	    goto L160;
	}
/*           CHOOSE NEXT PERTURBATION VECTOR */
	if (delf == 0.) {
	    delf = 1.;
	}
	i__2 = *neq;
	for (j = 1; j <= i__2; ++j) {
	    if (k == 2) {
		goto L110;
	    }
	    dy = (d__1 = pv[j], abs(d__1));
	    if (dy == 0.) {
		dy = delf;
	    }
	    goto L120;
L110:
	    dy = y[j];
	    if (dy == 0.) {
		dy = dely / relper;
	    }
L120:
	    if (spy[j] == 0.) {
		spy[j] = yp[j];
	    }
	    if (spy[j] != 0.) {
		dy = d_sign(&dy, &spy[j]);
	    }
	    yp[j] = dy;
/* L130: */
	}
	delf = dhvnrm_(&yp[1], neq);
/* L140: */
    }
L150:

/*        PROTECT AGAINST AN OVERFLOW */
    dfdub = *big;
L160:

/*     .................................................................. */

/*          COMPUTE A BOUND (YDPB) ON THE NORM OF THE SECOND DERIVATIVE */

    ydpb = dfdxb + dfdub * fbnd;

/*     .................................................................. */

/*          DEFINE THE TOLERANCE PARAMETER UPON WHICH THE STARTING STEP */
/*          SIZE IS TO BE BASED.  A VALUE IN THE MIDDLE OF THE ERROR */
/*          TOLERANCE RANGE IS SELECTED. */

    tolmin = *big;
    tolsum = 0.;
    i__1 = *neq;
    for (k = 1; k <= i__1; ++k) {
	tolexp = d_lg10(&etol[k]);
	tolmin = min(tolmin,tolexp);
	tolsum += tolexp;
/* L170: */
    }
    d__1 = (tolsum / *neq + tolmin) * .5 / (*morder + 1);
    tolp = pow_dd(&c_b20, &d__1);

/*     .................................................................. */

/*          COMPUTE A STARTING STEP SIZE BASED ON THE ABOVE FIRST AND */
/*          SECOND DERIVATIVE INFORMATION */

/*                            RESTRICT THE STEP LENGTH TO BE NOT BIGGER */
/*                            THAN ABS(B-A).   (UNLESS  B  IS TOO CLOSE */
/*                            TO  A) */
    *h__ = absdx;

    if (ydpb != 0. || fbnd != 0.) {
	goto L180;
    }

/*        BOTH FIRST DERIVATIVE TERM (FBND) AND SECOND */
/*                     DERIVATIVE TERM (YDPB) ARE ZERO */
    if (tolp < 1.) {
	*h__ = absdx * tolp;
    }
    goto L200;
L180:

    if (ydpb != 0.) {
	goto L190;
    }

/*        ONLY SECOND DERIVATIVE TERM (YDPB) IS ZERO */
    if (tolp < fbnd * absdx) {
	*h__ = tolp / fbnd;
    }
    goto L200;
L190:

/*        SECOND DERIVATIVE TERM (YDPB) IS NON-ZERO */
    srydpb = sqrt(ydpb * .5);
    if (tolp < srydpb * absdx) {
	*h__ = tolp / srydpb;
    }
L200:

/*     FURTHER RESTRICT THE STEP LENGTH TO BE NOT */
/*                               BIGGER THAN  1/DFDUB */
    if (*h__ * dfdub > 1.) {
	*h__ = 1. / dfdub;
    }

/*     FINALLY, RESTRICT THE STEP LENGTH TO BE NOT */
/*     SMALLER THAN  100*SMALL*ABS(A).  HOWEVER, IF */
/*     A=0. AND THE COMPUTED H UNDERFLOWED TO ZERO, */
/*     THE ALGORITHM RETURNS  SMALL*ABS(B)  FOR THE */
/*                                     STEP LENGTH. */
/* Computing MAX */
    d__1 = *h__, d__2 = *small * 100. * abs(*a);
    *h__ = max(d__1,d__2);
    if (*h__ == 0.) {
	*h__ = *small * abs(*b);
    }

/*     NOW SET DIRECTION OF INTEGRATION */
    *h__ = d_sign(h__, &dx);

    return 0;
} /* dhstrt_ */

