/* hstart.f -- translated by f2c (version 12.02.01).
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
static real c_b13 = 1.f;

/* DECK HSTART */
/* Subroutine */ int hstart_(S_fp f, integer *neq, real *a, real *b, real *y, 
	real *yprime, real *etol, integer *morder, real *small, real *big, 
	real *spy, real *pv, real *yp, real *sf, real *rpar, integer *ipar, 
	real *h__)
{
    /* System generated locals */
    integer i__1, i__2;
    real r__1, r__2, r__3, r__4;
    doublereal d__1, d__2;

    /* Local variables */
    static integer j, k;
    static real da;
    static integer lk;
    static real dx, dy, wtj, fbnd, delf, delx, dely, ydpb;
    static integer icase;
    static real dfdub, dfdxb, delxb, absdx;
    extern doublereal hvnrm_(real *, integer *);
    static real power, ynorm, relper, srydpb, ypnorm;

/* ***BEGIN PROLOGUE  HSTART */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to DEABM, DEBDF and DERKF */
/* ***LIBRARY   SLATEC */
/* ***TYPE      SINGLE PRECISION (HSTART-S, DHSTRT-D) */
/* ***AUTHOR  Watts, H. A., (SNLA) */
/* ***DESCRIPTION */

/*   HSTART computes a starting step size to be used in solving initial */
/*   value problems in ordinary differential equations. */
/* ********************************************************************** */
/*  Abstract */

/*     Subroutine HSTART computes a starting step size to be used by an */
/*     initial value method in solving ordinary differential equations. */
/*     It is based on an estimate of the local Lipschitz constant for the */
/*     differential equation (lower bound on a norm of the Jacobian), */
/*     a bound on the differential equation (first derivative), and */
/*     a bound on the partial derivative of the equation with respect to */
/*     the independent variable. */
/*     (All approximated near the initial point A.) */

/*     Subroutine HSTART uses a function subprogram HVNRM for computing */
/*     a vector norm.  The maximum norm is presently utilized though it */
/*     can easily be replaced by any other vector norm.  It is presumed */
/*     that any replacement norm routine would be carefully coded to */
/*     prevent unnecessary underflows or overflows from occurring, and */
/*     also, would not alter the vector or number of components. */

/* ********************************************************************** */
/*  On Input you must provide the following */

/*      F -- This is a subroutine of the form */
/*                               F(X,U,UPRIME,RPAR,IPAR) */
/*             which defines the system of first order differential */
/*             equations to be solved.  For the given values of X and the */
/*             vector  U(*)=(U(1),U(2),...,U(NEQ)) , the subroutine must */
/*             evaluate the NEQ components of the system of differential */
/*             equations  dU/DX=F(X,U)  and store the derivatives in the */
/*             array UPRIME(*), that is,  UPRIME(I) = * dU(I)/DX *  for */
/*             equations I=1,...,NEQ. */

/*             Subroutine F must not alter X or U(*).  You must declare */
/*             the name F in an EXTERNAL statement in your program that */
/*             calls HSTART.  You must dimension U and UPRIME in F. */

/*             RPAR and IPAR are real and integer parameter arrays which */
/*             you can use for communication between your program and */
/*             subroutine F.  They are not used or altered by HSTART.  If */
/*             you do not need RPAR or IPAR, ignore these parameters by */
/*             treating them as dummy arguments.  If you do choose to use */
/*             them, dimension them in your program and in F as arrays */
/*             of appropriate length. */

/*      NEQ -- This is the number of (first order) differential equations */
/*             to be integrated. */

/*      A -- This is the initial point of integration. */

/*      B -- This is a value of the independent variable used to define */
/*             the direction of integration.  A reasonable choice is to */
/*             set  B  to the first point at which a solution is desired. */
/*             You can also use  B, if necessary, to restrict the length */
/*             of the first integration step because the algorithm will */
/*             not compute a starting step length which is bigger than */
/*             ABS(B-A), unless  B  has been chosen too close to  A. */
/*             (It is presumed that HSTART has been called with  B */
/*             different from  A  on the machine being used.  Also see */
/*             the discussion about the parameter  SMALL.) */

/*      Y(*) -- This is the vector of initial values of the NEQ solution */
/*             components at the initial point  A. */

/*      YPRIME(*) -- This is the vector of derivatives of the NEQ */
/*             solution components at the initial point  A. */
/*             (defined by the differential equations in subroutine F) */

/*      ETOL -- This is the vector of error tolerances corresponding to */
/*             the NEQ solution components.  It is assumed that all */
/*             elements are positive.  Following the first integration */
/*             step, the tolerances are expected to be used by the */
/*             integrator in an error test which roughly requires that */
/*                        ABS(local error) .LE. ETOL */
/*             for each vector component. */

/*      MORDER -- This is the order of the formula which will be used by */
/*             the initial value method for taking the first integration */
/*             step. */

/*      SMALL -- This is a small positive machine dependent constant */
/*             which is used for protecting against computations with */
/*             numbers which are too small relative to the precision of */
/*             floating point arithmetic.  SMALL  should be set to */
/*             (approximately) the smallest positive real number such */
/*             that  (1.+SMALL) .GT. 1.  on the machine being used. the */
/*             quantity  SMALL**(3/8)  is used in computing increments of */
/*             variables for approximating derivatives by differences. */
/*             also the algorithm will not compute a starting step length */
/*             which is smaller than  100*SMALL*ABS(A). */

/*      BIG -- This is a large positive machine dependent constant which */
/*             is used for preventing machine overflows.  A reasonable */
/*             choice is to set big to (approximately) the square root of */
/*             the largest real number which can be held in the machine. */

/*      SPY(*),PV(*),YP(*),SF(*) -- These are real work arrays of length */
/*             NEQ which provide the routine with needed storage space. */

/*      RPAR,IPAR -- These are parameter arrays, of real and integer */
/*             type, respectively, which can be used for communication */
/*             between your program and the F subroutine.  They are not */
/*             used or altered by HSTART. */

/* ********************************************************************** */
/*  On Output  (after the return from HSTART), */

/*      H -- Is an appropriate starting step size to be attempted by the */
/*             differential equation method. */

/*           All parameters in the call list remain unchanged except for */
/*           the working arrays SPY(*),PV(*),YP(*) and SF(*). */

/* ********************************************************************** */

/* ***SEE ALSO  DEABM, DEBDF, DERKF */
/* ***ROUTINES CALLED  HVNRM */
/* ***REVISION HISTORY  (YYMMDD) */
/*   800501  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   891024  Changed references from VNORM to HVNRM.  (WRB) */
/*   891024  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900328  Added TYPE section.  (WRB) */
/*   910722  Updated AUTHOR section.  (ALS) */
/* ***END PROLOGUE  HSTART */


/* ....................................................................... */

/* ***FIRST EXECUTABLE STATEMENT  HSTART */
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
    absdx = dabs(dx);
    d__1 = (doublereal) (*small);
    relper = pow_dd(&d__1, &c_b2);
    ynorm = hvnrm_(&y[1], neq);

/* ....................................................................... */

/*     COMPUTE A WEIGHTED APPROXIMATE BOUND (DFDXB) ON THE PARTIAL */
/*     DERIVATIVE OF THE EQUATION WITH RESPECT TO THE */
/*     INDEPENDENT VARIABLE. PROTECT AGAINST AN OVERFLOW. ALSO */
/*     COMPUTE A WEIGHTED BOUND (FBND) ON THE FIRST DERIVATIVE LOCALLY. */

/* Computing MAX */
/* Computing MIN */
    r__4 = relper * dabs(*a);
    r__2 = dmin(r__4,absdx), r__3 = *small * 100.f * dabs(*a);
    r__1 = dmax(r__2,r__3);
    da = r_sign(&r__1, &dx);
    if (da == 0.f) {
	da = relper * dx;
    }
    r__1 = *a + da;
    (*f)(&r__1, &y[1], &sf[1], &rpar[1], &ipar[1]);

    if (*morder == 1) {
	goto L20;
    }
    power = 2.f / (*morder + 1);
    i__1 = *neq;
    for (j = 1; j <= i__1; ++j) {
	d__1 = (doublereal) etol[j];
	d__2 = (doublereal) power;
	wtj = pow_dd(&d__1, &d__2);
	spy[j] = sf[j] / wtj;
	yp[j] = yprime[j] / wtj;
/* L10: */
	pv[j] = spy[j] - yp[j];
    }
    goto L40;

L20:
    i__1 = *neq;
    for (j = 1; j <= i__1; ++j) {
	spy[j] = sf[j] / etol[j];
	yp[j] = yprime[j] / etol[j];
/* L30: */
	pv[j] = spy[j] - yp[j];
    }

L40:
    delf = hvnrm_(&pv[1], neq);
    dfdxb = *big;
    if (delf < *big * dabs(da)) {
	dfdxb = delf / dabs(da);
    }
    ypnorm = hvnrm_(&yp[1], neq);
/* Computing MAX */
    r__1 = hvnrm_(&spy[1], neq);
    fbnd = dmax(r__1,ypnorm);

/* ....................................................................... */

/*     COMPUTE AN ESTIMATE (DFDUB) OF THE LOCAL LIPSCHITZ CONSTANT FOR */
/*     THE SYSTEM OF DIFFERENTIAL EQUATIONS. THIS ALSO REPRESENTS AN */
/*     ESTIMATE OF THE NORM OF THE JACOBIAN LOCALLY. */
/*     THREE ITERATIONS (TWO WHEN NEQ=1) ARE USED TO ESTIMATE THE */
/*     LIPSCHITZ CONSTANT BY NUMERICAL DIFFERENCES. THE FIRST */
/*     PERTURBATION VECTOR IS BASED ON THE INITIAL DERIVATIVES AND */
/*     DIRECTION OF INTEGRATION. THE SECOND PERTURBATION VECTOR IS */
/*     FORMED USING ANOTHER EVALUATION OF THE DIFFERENTIAL EQUATION. */
/*     THE THIRD PERTURBATION VECTOR IS FORMED USING PERTURBATIONS BASED */
/*     ONLY ON THE INITIAL VALUES. COMPONENTS THAT ARE ZERO ARE ALWAYS */
/*     CHANGED TO NON-ZERO VALUES (EXCEPT ON THE FIRST ITERATION). WHEN */
/*     INFORMATION IS AVAILABLE, CARE IS TAKEN TO ENSURE THAT COMPONENTS */
/*     OF THE PERTURBATION VECTOR HAVE SIGNS WHICH ARE CONSISTENT WITH */
/*     THE SLOPES OF LOCAL SOLUTION CURVES. */
/*     ALSO CHOOSE THE LARGEST BOUND (FBND) FOR THE FIRST DERIVATIVE. */
/*     NO ATTEMPT IS MADE TO KEEP THE PERTURBATION VECTOR SIZE CONSTANT. */

    if (ypnorm == 0.f) {
	goto L60;
    }
/*                       USE INITIAL DERIVATIVES FOR FIRST PERTURBATION */
    icase = 1;
    i__1 = *neq;
    for (j = 1; j <= i__1; ++j) {
	spy[j] = yprime[j];
/* L50: */
	yp[j] = yprime[j];
    }
    goto L80;
/*                       CANNOT HAVE A NULL PERTURBATION VECTOR */
L60:
    icase = 2;
    i__1 = *neq;
    for (j = 1; j <= i__1; ++j) {
	spy[j] = yprime[j];
/* L70: */
	yp[j] = etol[j];
    }

L80:
    dfdub = 0.f;
/* Computing MIN */
    i__1 = *neq + 1;
    lk = min(i__1,3);
    i__1 = lk;
    for (k = 1; k <= i__1; ++k) {
/*                       SET YPNORM AND DELX */
	ypnorm = hvnrm_(&yp[1], neq);
	if (icase == 1 || icase == 3) {
	    goto L90;
	}
	delx = r_sign(&c_b13, &dx);
	goto L120;
/*                       TRY TO ENFORCE MEANINGFUL PERTURBATION VALUES */
L90:
	delx = dx;
	if (dabs(delx) * ypnorm >= relper * ynorm) {
	    goto L100;
	}
	delxb = *big;
	if (relper * ynorm < *big * ypnorm) {
	    delxb = relper * ynorm / ypnorm;
	}
	delx = r_sign(&delxb, &dx);
L100:
	i__2 = *neq;
	for (j = 1; j <= i__2; ++j) {
	    if ((r__1 = delx * yp[j], dabs(r__1)) > etol[j]) {
		r__2 = etol[j] / yp[j];
		delx = r_sign(&r__2, &dx);
	    }
/* L110: */
	}
/*                       DEFINE PERTURBED VECTOR OF INITIAL VALUES */
L120:
	i__2 = *neq;
	for (j = 1; j <= i__2; ++j) {
/* L130: */
	    pv[j] = y[j] + delx * yp[j];
	}
	if (k == 2) {
	    goto L150;
	}
/*                       EVALUATE DERIVATIVES ASSOCIATED WITH PERTURBED */
/*                       VECTOR  AND  COMPUTE CORRESPONDING DIFFERENCES */
	(*f)(a, &pv[1], &yp[1], &rpar[1], &ipar[1]);
	i__2 = *neq;
	for (j = 1; j <= i__2; ++j) {
/* L140: */
	    pv[j] = yp[j] - yprime[j];
	}
	goto L170;
/*                       USE A SHIFTED VALUE OF THE INDEPENDENT VARIABLE */
/*                                             IN COMPUTING ONE ESTIMATE */
L150:
	r__1 = *a + da;
	(*f)(&r__1, &pv[1], &yp[1], &rpar[1], &ipar[1]);
	i__2 = *neq;
	for (j = 1; j <= i__2; ++j) {
/* L160: */
	    pv[j] = yp[j] - sf[j];
	}
/*                       CHOOSE LARGEST BOUND ON THE WEIGHTED FIRST */
/*                                                   DERIVATIVE */
L170:
	if (*morder == 1) {
	    goto L190;
	}
	i__2 = *neq;
	for (j = 1; j <= i__2; ++j) {
/* L180: */
	    d__1 = (doublereal) etol[j];
	    d__2 = (doublereal) power;
	    yp[j] /= pow_dd(&d__1, &d__2);
	}
	goto L210;
L190:
	i__2 = *neq;
	for (j = 1; j <= i__2; ++j) {
/* L200: */
	    yp[j] /= etol[j];
	}
L210:
/* Computing MAX */
	r__1 = fbnd, r__2 = hvnrm_(&yp[1], neq);
	fbnd = dmax(r__1,r__2);
/*                       COMPUTE BOUND ON A LOCAL LIPSCHITZ CONSTANT */
	delf = hvnrm_(&pv[1], neq);
	if (delf == 0.f) {
	    goto L220;
	}
	dely = dabs(delx) * ypnorm;
	if (delf >= *big * dely) {
	    goto L270;
	}
/* Computing MAX */
	r__1 = dfdub, r__2 = delf / dely;
	dfdub = dmax(r__1,r__2);

L220:
	if (k == lk) {
	    goto L280;
	}
/*                       CHOOSE NEXT PERTURBATION VECTOR */
	i__2 = *neq;
	for (j = 1; j <= i__2; ++j) {
	    if (k == lk - 1) {
		goto L230;
	    }
	    icase = 3;
	    dy = (r__1 = pv[j], dabs(r__1));
	    if (dy == 0.f) {
/* Computing MAX */
		r__1 = delf, r__2 = etol[j];
		dy = dmax(r__1,r__2);
	    }
	    goto L240;
L230:
	    icase = 4;
/* Computing MAX */
	    r__2 = relper * (r__1 = y[j], dabs(r__1)), r__3 = etol[j];
	    dy = dmax(r__2,r__3);
L240:
	    if (spy[j] == 0.f) {
		spy[j] = yp[j];
	    }
	    if (spy[j] != 0.f) {
		dy = r_sign(&dy, &spy[j]);
	    }
/* L250: */
	    yp[j] = dy;
	}
/* L260: */
    }

/*                       PROTECT AGAINST AN OVERFLOW */
L270:
    dfdub = *big;

/* ....................................................................... */

/*     COMPUTE A BOUND (YDPB) ON THE NORM OF THE SECOND DERIVATIVE */

L280:
    ydpb = dfdxb + dfdub * fbnd;

/* ....................................................................... */

/*     COMPUTE A STARTING STEP SIZE BASED ON THE ABOVE FIRST AND SECOND */
/*     DERIVATIVE INFORMATION */

/*                       RESTRICT THE STEP LENGTH TO BE NOT BIGGER THAN */
/*                       ABS(B-A).   (UNLESS  B  IS TOO CLOSE TO  A) */
    *h__ = absdx;

    if (ydpb != 0.f || fbnd != 0.f) {
	goto L290;
    }

/*                       BOTH FIRST DERIVATIVE TERM (FBND) AND SECOND */
/*                                    DERIVATIVE TERM (YDPB) ARE ZERO */
    goto L310;

L290:
    if (ydpb != 0.f) {
	goto L300;
    }

/*                       ONLY SECOND DERIVATIVE TERM (YDPB) IS ZERO */
    if (1.f < fbnd * absdx) {
	*h__ = 1.f / fbnd;
    }
    goto L310;

/*                       SECOND DERIVATIVE TERM (YDPB) IS NON-ZERO */
L300:
    srydpb = sqrt(ydpb * .5f);
    if (1.f < srydpb * absdx) {
	*h__ = 1.f / srydpb;
    }

/*                       FURTHER RESTRICT THE STEP LENGTH TO BE NOT */
/*                                                 BIGGER THAN  1/DFDUB */
L310:
    if (*h__ * dfdub > 1.f) {
	*h__ = 1.f / dfdub;
    }

/*                       FINALLY, RESTRICT THE STEP LENGTH TO BE NOT */
/*                       SMALLER THAN  100*SMALL*ABS(A).  HOWEVER, IF */
/*                       A=0. AND THE COMPUTED H UNDERFLOWED TO ZERO, */
/*                       THE ALGORITHM RETURNS  SMALL*ABS(B)  FOR THE */
/*                                                       STEP LENGTH. */
/* Computing MAX */
    r__1 = *h__, r__2 = *small * 100.f * dabs(*a);
    *h__ = dmax(r__1,r__2);
    if (*h__ == 0.f) {
	*h__ = *small * dabs(*b);
    }

/*                       NOW SET DIRECTION OF INTEGRATION */
    *h__ = r_sign(h__, &dx);

    return 0;
} /* hstart_ */

