/* dqagie.f -- translated by f2c (version 12.02.01).
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
static doublereal c_b4 = 0.;
static doublereal c_b5 = 1.;
static integer c__1 = 1;
static integer c__2 = 2;

/* DECK DQAGIE */
/* Subroutine */ int dqagie_(D_fp f, doublereal *bound, integer *inf, 
	doublereal *epsabs, doublereal *epsrel, integer *limit, doublereal *
	result, doublereal *abserr, integer *neval, integer *ier, doublereal *
	alist__, doublereal *blist, doublereal *rlist, doublereal *elist, 
	integer *iord, integer *last)
{
    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1, d__2;

    /* Local variables */
    static integer k;
    static doublereal a1, a2, b1, b2;
    static integer id;
    static doublereal area, dres;
    static integer ksgn;
    static doublereal boun;
    static integer nres;
    static doublereal area1, area2, area12;
    extern /* Subroutine */ int dqelg_(integer *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, integer *), dqk15i_(D_fp, doublereal *
	    , integer *, doublereal *, doublereal *, doublereal *, doublereal 
	    *, doublereal *, doublereal *);
    static doublereal small, erro12;
    static integer ierro;
    static doublereal defab1, defab2;
    static integer ktmin, nrmax;
    static doublereal oflow, uflow;
    extern doublereal d1mach_(integer *);
    static logical noext;
    static integer iroff1, iroff2, iroff3;
    static doublereal res3la[3], error1, error2, rlist2[52];
    static integer numrl2;
    static doublereal defabs, epmach, erlarg, abseps, correc, errbnd, resabs;
    static integer jupbnd;
    static doublereal erlast, errmax;
    static integer maxerr;
    static doublereal reseps;
    static logical extrap;
    static doublereal ertest, errsum;
    extern /* Subroutine */ int dqpsrt_(integer *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, integer *);

/* ***BEGIN PROLOGUE  DQAGIE */
/* ***PURPOSE  The routine calculates an approximation result to a given */
/*            integral   I = Integral of F over (BOUND,+INFINITY) */
/*            or I = Integral of F over (-INFINITY,BOUND) */
/*            or I = Integral of F over (-INFINITY,+INFINITY), */
/*            hopefully satisfying following claim for accuracy */
/*            ABS(I-RESULT).LE.MAX(EPSABS,EPSREL*ABS(I)) */
/* ***LIBRARY   SLATEC (QUADPACK) */
/* ***CATEGORY  H2A3A1, H2A4A1 */
/* ***TYPE      DOUBLE PRECISION (QAGIE-S, DQAGIE-D) */
/* ***KEYWORDS  AUTOMATIC INTEGRATOR, EXTRAPOLATION, GENERAL-PURPOSE, */
/*             GLOBALLY ADAPTIVE, INFINITE INTERVALS, QUADPACK, */
/*             QUADRATURE, TRANSFORMATION */
/* ***AUTHOR  Piessens, Robert */
/*             Applied Mathematics and Programming Division */
/*             K. U. Leuven */
/*           de Doncker, Elise */
/*             Applied Mathematics and Programming Division */
/*             K. U. Leuven */
/* ***DESCRIPTION */

/* Integration over infinite intervals */
/* Standard fortran subroutine */

/*            F      - Double precision */
/*                     Function subprogram defining the integrand */
/*                     function F(X). The actual name for F needs to be */
/*                     declared E X T E R N A L in the driver program. */

/*            BOUND  - Double precision */
/*                     Finite bound of integration range */
/*                     (has no meaning if interval is doubly-infinite) */

/*            INF    - Double precision */
/*                     Indicating the kind of integration range involved */
/*                     INF = 1 corresponds to  (BOUND,+INFINITY), */
/*                     INF = -1            to  (-INFINITY,BOUND), */
/*                     INF = 2             to (-INFINITY,+INFINITY). */

/*            EPSABS - Double precision */
/*                     Absolute accuracy requested */
/*            EPSREL - Double precision */
/*                     Relative accuracy requested */
/*                     If  EPSABS.LE.0 */
/*                     and EPSREL.LT.MAX(50*REL.MACH.ACC.,0.5D-28), */
/*                     the routine will end with IER = 6. */

/*            LIMIT  - Integer */
/*                     Gives an upper bound on the number of subintervals */
/*                     in the partition of (A,B), LIMIT.GE.1 */

/*         ON RETURN */
/*            RESULT - Double precision */
/*                     Approximation to the integral */

/*            ABSERR - Double precision */
/*                     Estimate of the modulus of the absolute error, */
/*                     which should equal or exceed ABS(I-RESULT) */

/*            NEVAL  - Integer */
/*                     Number of integrand evaluations */

/*            IER    - Integer */
/*                     IER = 0 Normal and reliable termination of the */
/*                             routine. It is assumed that the requested */
/*                             accuracy has been achieved. */
/*                   - IER.GT.0 Abnormal termination of the routine. The */
/*                             estimates for result and error are less */
/*                             reliable. It is assumed that the requested */
/*                             accuracy has not been achieved. */
/*            ERROR MESSAGES */
/*                     IER = 1 Maximum number of subdivisions allowed */
/*                             has been achieved. One can allow more */
/*                             subdivisions by increasing the value of */
/*                             LIMIT (and taking the according dimension */
/*                             adjustments into account).  However, if */
/*                             this yields no improvement it is advised */
/*                             to analyze the integrand in order to */
/*                             determine the integration difficulties. */
/*                             If the position of a local difficulty can */
/*                             be determined (e.g. SINGULARITY, */
/*                             DISCONTINUITY within the interval) one */
/*                             will probably gain from splitting up the */
/*                             interval at this point and calling the */
/*                             integrator on the subranges. If possible, */
/*                             an appropriate special-purpose integrator */
/*                             should be used, which is designed for */
/*                             handling the type of difficulty involved. */
/*                         = 2 The occurrence of roundoff error is */
/*                             detected, which prevents the requested */
/*                             tolerance from being achieved. */
/*                             The error may be under-estimated. */
/*                         = 3 Extremely bad integrand behaviour occurs */
/*                             at some points of the integration */
/*                             interval. */
/*                         = 4 The algorithm does not converge. */
/*                             Roundoff error is detected in the */
/*                             extrapolation table. */
/*                             It is assumed that the requested tolerance */
/*                             cannot be achieved, and that the returned */
/*                             result is the best which can be obtained. */
/*                         = 5 The integral is probably divergent, or */
/*                             slowly convergent. It must be noted that */
/*                             divergence can occur with any other value */
/*                             of IER. */
/*                         = 6 The input is invalid, because */
/*                             (EPSABS.LE.0 and */
/*                              EPSREL.LT.MAX(50*REL.MACH.ACC.,0.5D-28), */
/*                             RESULT, ABSERR, NEVAL, LAST, RLIST(1), */
/*                             ELIST(1) and IORD(1) are set to zero. */
/*                             ALIST(1) and BLIST(1) are set to 0 */
/*                             and 1 respectively. */

/*            ALIST  - Double precision */
/*                     Vector of dimension at least LIMIT, the first */
/*                      LAST  elements of which are the left */
/*                     end points of the subintervals in the partition */
/*                     of the transformed integration range (0,1). */

/*            BLIST  - Double precision */
/*                     Vector of dimension at least LIMIT, the first */
/*                      LAST  elements of which are the right */
/*                     end points of the subintervals in the partition */
/*                     of the transformed integration range (0,1). */

/*            RLIST  - Double precision */
/*                     Vector of dimension at least LIMIT, the first */
/*                      LAST  elements of which are the integral */
/*                     approximations on the subintervals */

/*            ELIST  - Double precision */
/*                     Vector of dimension at least LIMIT,  the first */
/*                     LAST elements of which are the moduli of the */
/*                     absolute error estimates on the subintervals */

/*            IORD   - Integer */
/*                     Vector of dimension LIMIT, the first K */
/*                     elements of which are pointers to the */
/*                     error estimates over the subintervals, */
/*                     such that ELIST(IORD(1)), ..., ELIST(IORD(K)) */
/*                     form a decreasing sequence, with K = LAST */
/*                     If LAST.LE.(LIMIT/2+2), and K = LIMIT+1-LAST */
/*                     otherwise */

/*            LAST   - Integer */
/*                     Number of subintervals actually produced */
/*                     in the subdivision process */

/* ***REFERENCES  (NONE) */
/* ***ROUTINES CALLED  D1MACH, DQELG, DQK15I, DQPSRT */
/* ***REVISION HISTORY  (YYMMDD) */
/*   800101  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890831  Modified array declarations.  (WRB) */
/*   890831  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/* ***END PROLOGUE  DQAGIE */



/*            THE DIMENSION OF RLIST2 IS DETERMINED BY THE VALUE OF */
/*            LIMEXP IN SUBROUTINE DQELG. */


/*            LIST OF MAJOR VARIABLES */
/*            ----------------------- */

/*           ALIST     - LIST OF LEFT END POINTS OF ALL SUBINTERVALS */
/*                       CONSIDERED UP TO NOW */
/*           BLIST     - LIST OF RIGHT END POINTS OF ALL SUBINTERVALS */
/*                       CONSIDERED UP TO NOW */
/*           RLIST(I)  - APPROXIMATION TO THE INTEGRAL OVER */
/*                       (ALIST(I),BLIST(I)) */
/*           RLIST2    - ARRAY OF DIMENSION AT LEAST (LIMEXP+2), */
/*                       CONTAINING THE PART OF THE EPSILON TABLE */
/*                       WHICH IS STILL NEEDED FOR FURTHER COMPUTATIONS */
/*           ELIST(I)  - ERROR ESTIMATE APPLYING TO RLIST(I) */
/*           MAXERR    - POINTER TO THE INTERVAL WITH LARGEST ERROR */
/*                       ESTIMATE */
/*           ERRMAX    - ELIST(MAXERR) */
/*           ERLAST    - ERROR ON THE INTERVAL CURRENTLY SUBDIVIDED */
/*                       (BEFORE THAT SUBDIVISION HAS TAKEN PLACE) */
/*           AREA      - SUM OF THE INTEGRALS OVER THE SUBINTERVALS */
/*           ERRSUM    - SUM OF THE ERRORS OVER THE SUBINTERVALS */
/*           ERRBND    - REQUESTED ACCURACY MAX(EPSABS,EPSREL* */
/*                       ABS(RESULT)) */
/*           *****1    - VARIABLE FOR THE LEFT SUBINTERVAL */
/*           *****2    - VARIABLE FOR THE RIGHT SUBINTERVAL */
/*           LAST      - INDEX FOR SUBDIVISION */
/*           NRES      - NUMBER OF CALLS TO THE EXTRAPOLATION ROUTINE */
/*           NUMRL2    - NUMBER OF ELEMENTS CURRENTLY IN RLIST2. IF AN */
/*                       APPROPRIATE APPROXIMATION TO THE COMPOUNDED */
/*                       INTEGRAL HAS BEEN OBTAINED, IT IS PUT IN */
/*                       RLIST2(NUMRL2) AFTER NUMRL2 HAS BEEN INCREASED */
/*                       BY ONE. */
/*           SMALL     - LENGTH OF THE SMALLEST INTERVAL CONSIDERED UP */
/*                       TO NOW, MULTIPLIED BY 1.5 */
/*           ERLARG    - SUM OF THE ERRORS OVER THE INTERVALS LARGER */
/*                       THAN THE SMALLEST INTERVAL CONSIDERED UP TO NOW */
/*           EXTRAP    - LOGICAL VARIABLE DENOTING THAT THE ROUTINE */
/*                       IS ATTEMPTING TO PERFORM EXTRAPOLATION. I.E. */
/*                       BEFORE SUBDIVIDING THE SMALLEST INTERVAL WE */
/*                       TRY TO DECREASE THE VALUE OF ERLARG. */
/*           NOEXT     - LOGICAL VARIABLE DENOTING THAT EXTRAPOLATION */
/*                       IS NO LONGER ALLOWED (TRUE-VALUE) */

/*            MACHINE DEPENDENT CONSTANTS */
/*            --------------------------- */

/*           EPMACH IS THE LARGEST RELATIVE SPACING. */
/*           UFLOW IS THE SMALLEST POSITIVE MAGNITUDE. */
/*           OFLOW IS THE LARGEST POSITIVE MAGNITUDE. */

/* ***FIRST EXECUTABLE STATEMENT  DQAGIE */
    /* Parameter adjustments */
    --iord;
    --elist;
    --rlist;
    --blist;
    --alist__;

    /* Function Body */
    epmach = d1mach_(&c__4);

/*           TEST ON VALIDITY OF PARAMETERS */
/*           ----------------------------- */

    *ier = 0;
    *neval = 0;
    *last = 0;
    *result = 0.;
    *abserr = 0.;
    alist__[1] = 0.;
    blist[1] = 1.;
    rlist[1] = 0.;
    elist[1] = 0.;
    iord[1] = 0;
/* Computing MAX */
    d__1 = epmach * 50.;
    if (*epsabs <= 0. && *epsrel < max(d__1,5e-29)) {
	*ier = 6;
    }
    if (*ier == 6) {
	goto L999;
    }


/*           FIRST APPROXIMATION TO THE INTEGRAL */
/*           ----------------------------------- */

/*           DETERMINE THE INTERVAL TO BE MAPPED ONTO (0,1). */
/*           IF INF = 2 THE INTEGRAL IS COMPUTED AS I = I1+I2, WHERE */
/*           I1 = INTEGRAL OF F OVER (-INFINITY,0), */
/*           I2 = INTEGRAL OF F OVER (0,+INFINITY). */

    boun = *bound;
    if (*inf == 2) {
	boun = 0.;
    }
    dqk15i_((D_fp)f, &boun, inf, &c_b4, &c_b5, result, abserr, &defabs, &
	    resabs);

/*           TEST ON ACCURACY */

    *last = 1;
    rlist[1] = *result;
    elist[1] = *abserr;
    iord[1] = 1;
    dres = abs(*result);
/* Computing MAX */
    d__1 = *epsabs, d__2 = *epsrel * dres;
    errbnd = max(d__1,d__2);
    if (*abserr <= epmach * 100. * defabs && *abserr > errbnd) {
	*ier = 2;
    }
    if (*limit == 1) {
	*ier = 1;
    }
    if (*ier != 0 || *abserr <= errbnd && *abserr != resabs || *abserr == 0.) 
	    {
	goto L130;
    }

/*           INITIALIZATION */
/*           -------------- */

    uflow = d1mach_(&c__1);
    oflow = d1mach_(&c__2);
    rlist2[0] = *result;
    errmax = *abserr;
    maxerr = 1;
    area = *result;
    errsum = *abserr;
    *abserr = oflow;
    nrmax = 1;
    nres = 0;
    ktmin = 0;
    numrl2 = 2;
    extrap = FALSE_;
    noext = FALSE_;
    ierro = 0;
    iroff1 = 0;
    iroff2 = 0;
    iroff3 = 0;
    ksgn = -1;
    if (dres >= (1. - epmach * 50.) * defabs) {
	ksgn = 1;
    }

/*           MAIN DO-LOOP */
/*           ------------ */

    i__1 = *limit;
    for (*last = 2; *last <= i__1; ++(*last)) {

/*           BISECT THE SUBINTERVAL WITH NRMAX-TH LARGEST ERROR ESTIMATE. */

	a1 = alist__[maxerr];
	b1 = (alist__[maxerr] + blist[maxerr]) * .5;
	a2 = b1;
	b2 = blist[maxerr];
	erlast = errmax;
	dqk15i_((D_fp)f, &boun, inf, &a1, &b1, &area1, &error1, &resabs, &
		defab1);
	dqk15i_((D_fp)f, &boun, inf, &a2, &b2, &area2, &error2, &resabs, &
		defab2);

/*           IMPROVE PREVIOUS APPROXIMATIONS TO INTEGRAL */
/*           AND ERROR AND TEST FOR ACCURACY. */

	area12 = area1 + area2;
	erro12 = error1 + error2;
	errsum = errsum + erro12 - errmax;
	area = area + area12 - rlist[maxerr];
	if (defab1 == error1 || defab2 == error2) {
	    goto L15;
	}
	if ((d__1 = rlist[maxerr] - area12, abs(d__1)) > abs(area12) * 1e-5 ||
		 erro12 < errmax * .99) {
	    goto L10;
	}
	if (extrap) {
	    ++iroff2;
	}
	if (! extrap) {
	    ++iroff1;
	}
L10:
	if (*last > 10 && erro12 > errmax) {
	    ++iroff3;
	}
L15:
	rlist[maxerr] = area1;
	rlist[*last] = area2;
/* Computing MAX */
	d__1 = *epsabs, d__2 = *epsrel * abs(area);
	errbnd = max(d__1,d__2);

/*           TEST FOR ROUNDOFF ERROR AND EVENTUALLY SET ERROR FLAG. */

	if (iroff1 + iroff2 >= 10 || iroff3 >= 20) {
	    *ier = 2;
	}
	if (iroff2 >= 5) {
	    ierro = 3;
	}

/*           SET ERROR FLAG IN THE CASE THAT THE NUMBER OF */
/*           SUBINTERVALS EQUALS LIMIT. */

	if (*last == *limit) {
	    *ier = 1;
	}

/*           SET ERROR FLAG IN THE CASE OF BAD INTEGRAND BEHAVIOUR */
/*           AT SOME POINTS OF THE INTEGRATION RANGE. */

/* Computing MAX */
	d__1 = abs(a1), d__2 = abs(b2);
	if (max(d__1,d__2) <= (epmach * 100. + 1.) * (abs(a2) + uflow * 1e3)) 
		{
	    *ier = 4;
	}

/*           APPEND THE NEWLY-CREATED INTERVALS TO THE LIST. */

	if (error2 > error1) {
	    goto L20;
	}
	alist__[*last] = a2;
	blist[maxerr] = b1;
	blist[*last] = b2;
	elist[maxerr] = error1;
	elist[*last] = error2;
	goto L30;
L20:
	alist__[maxerr] = a2;
	alist__[*last] = a1;
	blist[*last] = b1;
	rlist[maxerr] = area2;
	rlist[*last] = area1;
	elist[maxerr] = error2;
	elist[*last] = error1;

/*           CALL SUBROUTINE DQPSRT TO MAINTAIN THE DESCENDING ORDERING */
/*           IN THE LIST OF ERROR ESTIMATES AND SELECT THE SUBINTERVAL */
/*           WITH NRMAX-TH LARGEST ERROR ESTIMATE (TO BE BISECTED NEXT). */

L30:
	dqpsrt_(limit, last, &maxerr, &errmax, &elist[1], &iord[1], &nrmax);
	if (errsum <= errbnd) {
	    goto L115;
	}
	if (*ier != 0) {
	    goto L100;
	}
	if (*last == 2) {
	    goto L80;
	}
	if (noext) {
	    goto L90;
	}
	erlarg -= erlast;
	if ((d__1 = b1 - a1, abs(d__1)) > small) {
	    erlarg += erro12;
	}
	if (extrap) {
	    goto L40;
	}

/*           TEST WHETHER THE INTERVAL TO BE BISECTED NEXT IS THE */
/*           SMALLEST INTERVAL. */

	if ((d__1 = blist[maxerr] - alist__[maxerr], abs(d__1)) > small) {
	    goto L90;
	}
	extrap = TRUE_;
	nrmax = 2;
L40:
	if (ierro == 3 || erlarg <= ertest) {
	    goto L60;
	}

/*           THE SMALLEST INTERVAL HAS THE LARGEST ERROR. */
/*           BEFORE BISECTING DECREASE THE SUM OF THE ERRORS OVER THE */
/*           LARGER INTERVALS (ERLARG) AND PERFORM EXTRAPOLATION. */

	id = nrmax;
	jupbnd = *last;
	if (*last > *limit / 2 + 2) {
	    jupbnd = *limit + 3 - *last;
	}
	i__2 = jupbnd;
	for (k = id; k <= i__2; ++k) {
	    maxerr = iord[nrmax];
	    errmax = elist[maxerr];
	    if ((d__1 = blist[maxerr] - alist__[maxerr], abs(d__1)) > small) {
		goto L90;
	    }
	    ++nrmax;
/* L50: */
	}

/*           PERFORM EXTRAPOLATION. */

L60:
	++numrl2;
	rlist2[numrl2 - 1] = area;
	dqelg_(&numrl2, rlist2, &reseps, &abseps, res3la, &nres);
	++ktmin;
	if (ktmin > 5 && *abserr < errsum * .001) {
	    *ier = 5;
	}
	if (abseps >= *abserr) {
	    goto L70;
	}
	ktmin = 0;
	*abserr = abseps;
	*result = reseps;
	correc = erlarg;
/* Computing MAX */
	d__1 = *epsabs, d__2 = *epsrel * abs(reseps);
	ertest = max(d__1,d__2);
	if (*abserr <= ertest) {
	    goto L100;
	}

/*            PREPARE BISECTION OF THE SMALLEST INTERVAL. */

L70:
	if (numrl2 == 1) {
	    noext = TRUE_;
	}
	if (*ier == 5) {
	    goto L100;
	}
	maxerr = iord[1];
	errmax = elist[maxerr];
	nrmax = 1;
	extrap = FALSE_;
	small *= .5;
	erlarg = errsum;
	goto L90;
L80:
	small = .375;
	erlarg = errsum;
	ertest = errbnd;
	rlist2[1] = area;
L90:
	;
    }

/*           SET FINAL RESULT AND ERROR ESTIMATE. */
/*           ------------------------------------ */

L100:
    if (*abserr == oflow) {
	goto L115;
    }
    if (*ier + ierro == 0) {
	goto L110;
    }
    if (ierro == 3) {
	*abserr += correc;
    }
    if (*ier == 0) {
	*ier = 3;
    }
    if (*result != 0. && area != 0.) {
	goto L105;
    }
    if (*abserr > errsum) {
	goto L115;
    }
    if (area == 0.) {
	goto L130;
    }
    goto L110;
L105:
    if (*abserr / abs(*result) > errsum / abs(area)) {
	goto L115;
    }

/*           TEST ON DIVERGENCE */

L110:
/* Computing MAX */
    d__1 = abs(*result), d__2 = abs(area);
    if (ksgn == -1 && max(d__1,d__2) <= defabs * .01) {
	goto L130;
    }
    if (.01 > *result / area || *result / area > 100. || errsum > abs(area)) {
	*ier = 6;
    }
    goto L130;

/*           COMPUTE GLOBAL INTEGRAL SUM. */

L115:
    *result = 0.;
    i__1 = *last;
    for (k = 1; k <= i__1; ++k) {
	*result += rlist[k];
/* L120: */
    }
    *abserr = errsum;
L130:
    *neval = *last * 30 - 15;
    if (*inf == 2) {
	*neval <<= 1;
    }
    if (*ier > 2) {
	--(*ier);
    }
L999:
    return 0;
} /* dqagie_ */

