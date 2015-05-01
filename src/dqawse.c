/* dqawse.f -- translated by f2c (version 12.02.01).
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
static integer c__1 = 1;

/* DECK DQAWSE */
/* Subroutine */ int dqawse_(D_fp f, doublereal *a, doublereal *b, doublereal 
	*alfa, doublereal *beta, integer *integr, doublereal *epsabs, 
	doublereal *epsrel, integer *limit, doublereal *result, doublereal *
	abserr, integer *neval, integer *ier, doublereal *alist__, doublereal 
	*blist, doublereal *rlist, doublereal *elist, integer *iord, integer *
	last)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2;

    /* Local variables */
    static integer k;
    static doublereal a1, a2, b1, b2, rg[25], rh[25], ri[25], rj[25];
    static integer nev;
    static doublereal area, area1, area2, area12;
    extern /* Subroutine */ int dqc25s_(D_fp, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, integer *, integer *);
    static doublereal erro12;
    static integer nrmax;
    static doublereal uflow;
    extern doublereal d1mach_(integer *);
    static integer iroff1, iroff2;
    static doublereal resas1, resas2, error1, error2, epmach, errbnd, centre;
    extern /* Subroutine */ int dqmomo_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, integer *)
	    ;
    static doublereal errmax;
    static integer maxerr;
    static doublereal errsum;
    extern /* Subroutine */ int dqpsrt_(integer *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, integer *);

/* ***BEGIN PROLOGUE  DQAWSE */
/* ***PURPOSE  The routine calculates an approximation result to a given */
/*            definite integral I = Integral of F*W over (A,B), */
/*            (where W shows a singular behaviour at the end points, */
/*            see parameter INTEGR). */
/*            Hopefully satisfying following claim for accuracy */
/*            ABS(I-RESULT).LE.MAX(EPSABS,EPSREL*ABS(I)). */
/* ***LIBRARY   SLATEC (QUADPACK) */
/* ***CATEGORY  H2A2A1 */
/* ***TYPE      DOUBLE PRECISION (QAWSE-S, DQAWSE-D) */
/* ***KEYWORDS  ALGEBRAIC-LOGARITHMIC END POINT SINGULARITIES, */
/*             AUTOMATIC INTEGRATOR, CLENSHAW-CURTIS METHOD, QUADPACK, */
/*             QUADRATURE, SPECIAL-PURPOSE */
/* ***AUTHOR  Piessens, Robert */
/*             Applied Mathematics and Programming Division */
/*             K. U. Leuven */
/*           de Doncker, Elise */
/*             Applied Mathematics and Programming Division */
/*             K. U. Leuven */
/* ***DESCRIPTION */

/*        Integration of functions having algebraico-logarithmic */
/*        end point singularities */
/*        Standard fortran subroutine */
/*        Double precision version */

/*        PARAMETERS */
/*         ON ENTRY */
/*            F      - Double precision */
/*                     Function subprogram defining the integrand */
/*                     function F(X). The actual name for F needs to be */
/*                     declared E X T E R N A L in the driver program. */

/*            A      - Double precision */
/*                     Lower limit of integration */

/*            B      - Double precision */
/*                     Upper limit of integration, B.GT.A */
/*                     If B.LE.A, the routine will end with IER = 6. */

/*            ALFA   - Double precision */
/*                     Parameter in the WEIGHT function, ALFA.GT.(-1) */
/*                     If ALFA.LE.(-1), the routine will end with */
/*                     IER = 6. */

/*            BETA   - Double precision */
/*                     Parameter in the WEIGHT function, BETA.GT.(-1) */
/*                     If BETA.LE.(-1), the routine will end with */
/*                     IER = 6. */

/*            INTEGR - Integer */
/*                     Indicates which WEIGHT function is to be used */
/*                     = 1  (X-A)**ALFA*(B-X)**BETA */
/*                     = 2  (X-A)**ALFA*(B-X)**BETA*LOG(X-A) */
/*                     = 3  (X-A)**ALFA*(B-X)**BETA*LOG(B-X) */
/*                     = 4  (X-A)**ALFA*(B-X)**BETA*LOG(X-A)*LOG(B-X) */
/*                     If INTEGR.LT.1 or INTEGR.GT.4, the routine */
/*                     will end with IER = 6. */

/*            EPSABS - Double precision */
/*                     Absolute accuracy requested */
/*            EPSREL - Double precision */
/*                     Relative accuracy requested */
/*                     If  EPSABS.LE.0 */
/*                     and EPSREL.LT.MAX(50*REL.MACH.ACC.,0.5D-28), */
/*                     the routine will end with IER = 6. */

/*            LIMIT  - Integer */
/*                     Gives an upper bound on the number of subintervals */
/*                     in the partition of (A,B), LIMIT.GE.2 */
/*                     If LIMIT.LT.2, the routine will end with IER = 6. */

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
/*                     IER.GT.0 Abnormal termination of the routine */
/*                             the estimates for the integral and error */
/*                             are less reliable. It is assumed that the */
/*                             requested accuracy has not been achieved. */
/*            ERROR MESSAGES */
/*                         = 1 Maximum number of subdivisions allowed */
/*                             has been achieved. One can allow more */
/*                             subdivisions by increasing the value of */
/*                             LIMIT. However, if this yields no */
/*                             improvement, it is advised to analyze the */
/*                             integrand in order to determine the */
/*                             integration difficulties which prevent the */
/*                             requested tolerance from being achieved. */
/*                             In case of a jump DISCONTINUITY or a local */
/*                             SINGULARITY of algebraico-logarithmic type */
/*                             at one or more interior points of the */
/*                             integration range, one should proceed by */
/*                             splitting up the interval at these */
/*                             points and calling the integrator on the */
/*                             subranges. */
/*                         = 2 The occurrence of roundoff error is */
/*                             detected, which prevents the requested */
/*                             tolerance from being achieved. */
/*                         = 3 Extremely bad integrand behaviour occurs */
/*                             at some points of the integration */
/*                             interval. */
/*                         = 6 The input is invalid, because */
/*                             B.LE.A or ALFA.LE.(-1) or BETA.LE.(-1), or */
/*                             INTEGR.LT.1 or INTEGR.GT.4, or */
/*                             (EPSABS.LE.0 and */
/*                              EPSREL.LT.MAX(50*REL.MACH.ACC.,0.5D-28), */
/*                             or LIMIT.LT.2. */
/*                             RESULT, ABSERR, NEVAL, RLIST(1), ELIST(1), */
/*                             IORD(1) and LAST are set to zero. ALIST(1) */
/*                             and BLIST(1) are set to A and B */
/*                             respectively. */

/*            ALIST  - Double precision */
/*                     Vector of dimension at least LIMIT, the first */
/*                      LAST  elements of which are the left */
/*                     end points of the subintervals in the partition */
/*                     of the given integration range (A,B) */

/*            BLIST  - Double precision */
/*                     Vector of dimension at least LIMIT, the first */
/*                      LAST  elements of which are the right */
/*                     end points of the subintervals in the partition */
/*                     of the given integration range (A,B) */

/*            RLIST  - Double precision */
/*                     Vector of dimension at least LIMIT, the first */
/*                      LAST  elements of which are the integral */
/*                     approximations on the subintervals */

/*            ELIST  - Double precision */
/*                     Vector of dimension at least LIMIT, the first */
/*                      LAST  elements of which are the moduli of the */
/*                     absolute error estimates on the subintervals */

/*            IORD   - Integer */
/*                     Vector of dimension at least LIMIT, the first K */
/*                     of which are pointers to the error */
/*                     estimates over the subintervals, so that */
/*                     ELIST(IORD(1)), ..., ELIST(IORD(K)) with K = LAST */
/*                     If LAST.LE.(LIMIT/2+2), and K = LIMIT+1-LAST */
/*                     otherwise form a decreasing sequence */

/*            LAST   - Integer */
/*                     Number of subintervals actually produced in */
/*                     the subdivision process */

/* ***REFERENCES  (NONE) */
/* ***ROUTINES CALLED  D1MACH, DQC25S, DQMOMO, DQPSRT */
/* ***REVISION HISTORY  (YYMMDD) */
/*   800101  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890831  Modified array declarations.  (WRB) */
/*   890831  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/* ***END PROLOGUE  DQAWSE */




/*            LIST OF MAJOR VARIABLES */
/*            ----------------------- */

/*           ALIST     - LIST OF LEFT END POINTS OF ALL SUBINTERVALS */
/*                       CONSIDERED UP TO NOW */
/*           BLIST     - LIST OF RIGHT END POINTS OF ALL SUBINTERVALS */
/*                       CONSIDERED UP TO NOW */
/*           RLIST(I)  - APPROXIMATION TO THE INTEGRAL OVER */
/*                       (ALIST(I),BLIST(I)) */
/*           ELIST(I)  - ERROR ESTIMATE APPLYING TO RLIST(I) */
/*           MAXERR    - POINTER TO THE INTERVAL WITH LARGEST */
/*                       ERROR ESTIMATE */
/*           ERRMAX    - ELIST(MAXERR) */
/*           AREA      - SUM OF THE INTEGRALS OVER THE SUBINTERVALS */
/*           ERRSUM    - SUM OF THE ERRORS OVER THE SUBINTERVALS */
/*           ERRBND    - REQUESTED ACCURACY MAX(EPSABS,EPSREL* */
/*                       ABS(RESULT)) */
/*           *****1    - VARIABLE FOR THE LEFT SUBINTERVAL */
/*           *****2    - VARIABLE FOR THE RIGHT SUBINTERVAL */
/*           LAST      - INDEX FOR SUBDIVISION */


/*            MACHINE DEPENDENT CONSTANTS */
/*            --------------------------- */

/*           EPMACH IS THE LARGEST RELATIVE SPACING. */
/*           UFLOW IS THE SMALLEST POSITIVE MAGNITUDE. */

/* ***FIRST EXECUTABLE STATEMENT  DQAWSE */
    /* Parameter adjustments */
    --iord;
    --elist;
    --rlist;
    --blist;
    --alist__;

    /* Function Body */
    epmach = d1mach_(&c__4);
    uflow = d1mach_(&c__1);

/*           TEST ON VALIDITY OF PARAMETERS */
/*           ------------------------------ */

    *ier = 6;
    *neval = 0;
    *last = 0;
    rlist[1] = 0.;
    elist[1] = 0.;
    iord[1] = 0;
    *result = 0.;
    *abserr = 0.;
/* Computing MAX */
    d__1 = epmach * 50.;
    if (*b <= *a || *epsabs == 0. && *epsrel < max(d__1,5e-29) || *alfa <= 
	    -1. || *beta <= -1. || *integr < 1 || *integr > 4 || *limit < 2) {
	goto L999;
    }
    *ier = 0;

/*           COMPUTE THE MODIFIED CHEBYSHEV MOMENTS. */

    dqmomo_(alfa, beta, ri, rj, rg, rh, integr);

/*           INTEGRATE OVER THE INTERVALS (A,(A+B)/2) AND ((A+B)/2,B). */

    centre = (*b + *a) * .5;
    dqc25s_((D_fp)f, a, b, a, &centre, alfa, beta, ri, rj, rg, rh, &area1, &
	    error1, &resas1, integr, &nev);
    *neval = nev;
    dqc25s_((D_fp)f, a, b, &centre, b, alfa, beta, ri, rj, rg, rh, &area2, &
	    error2, &resas2, integr, &nev);
    *last = 2;
    *neval += nev;
    *result = area1 + area2;
    *abserr = error1 + error2;

/*           TEST ON ACCURACY. */

/* Computing MAX */
    d__1 = *epsabs, d__2 = *epsrel * abs(*result);
    errbnd = max(d__1,d__2);

/*           INITIALIZATION */
/*           -------------- */

    if (error2 > error1) {
	goto L10;
    }
    alist__[1] = *a;
    alist__[2] = centre;
    blist[1] = centre;
    blist[2] = *b;
    rlist[1] = area1;
    rlist[2] = area2;
    elist[1] = error1;
    elist[2] = error2;
    goto L20;
L10:
    alist__[1] = centre;
    alist__[2] = *a;
    blist[1] = *b;
    blist[2] = centre;
    rlist[1] = area2;
    rlist[2] = area1;
    elist[1] = error2;
    elist[2] = error1;
L20:
    iord[1] = 1;
    iord[2] = 2;
    if (*limit == 2) {
	*ier = 1;
    }
    if (*abserr <= errbnd || *ier == 1) {
	goto L999;
    }
    errmax = elist[1];
    maxerr = 1;
    nrmax = 1;
    area = *result;
    errsum = *abserr;
    iroff1 = 0;
    iroff2 = 0;

/*            MAIN DO-LOOP */
/*            ------------ */

    i__1 = *limit;
    for (*last = 3; *last <= i__1; ++(*last)) {

/*           BISECT THE SUBINTERVAL WITH LARGEST ERROR ESTIMATE. */

	a1 = alist__[maxerr];
	b1 = (alist__[maxerr] + blist[maxerr]) * .5;
	a2 = b1;
	b2 = blist[maxerr];

	dqc25s_((D_fp)f, a, b, &a1, &b1, alfa, beta, ri, rj, rg, rh, &area1, &
		error1, &resas1, integr, &nev);
	*neval += nev;
	dqc25s_((D_fp)f, a, b, &a2, &b2, alfa, beta, ri, rj, rg, rh, &area2, &
		error2, &resas2, integr, &nev);
	*neval += nev;

/*           IMPROVE PREVIOUS APPROXIMATIONS INTEGRAL AND ERROR */
/*           AND TEST FOR ACCURACY. */

	area12 = area1 + area2;
	erro12 = error1 + error2;
	errsum = errsum + erro12 - errmax;
	area = area + area12 - rlist[maxerr];
	if (*a == a1 || *b == b2) {
	    goto L30;
	}
	if (resas1 == error1 || resas2 == error2) {
	    goto L30;
	}

/*           TEST FOR ROUNDOFF ERROR. */

	if ((d__1 = rlist[maxerr] - area12, abs(d__1)) < abs(area12) * 1e-5 &&
		 erro12 >= errmax * .99) {
	    ++iroff1;
	}
	if (*last > 10 && erro12 > errmax) {
	    ++iroff2;
	}
L30:
	rlist[maxerr] = area1;
	rlist[*last] = area2;

/*           TEST ON ACCURACY. */

/* Computing MAX */
	d__1 = *epsabs, d__2 = *epsrel * abs(area);
	errbnd = max(d__1,d__2);
	if (errsum <= errbnd) {
	    goto L35;
	}

/*           SET ERROR FLAG IN THE CASE THAT THE NUMBER OF INTERVAL */
/*           BISECTIONS EXCEEDS LIMIT. */

	if (*last == *limit) {
	    *ier = 1;
	}


/*           SET ERROR FLAG IN THE CASE OF ROUNDOFF ERROR. */

	if (iroff1 >= 6 || iroff2 >= 20) {
	    *ier = 2;
	}

/*           SET ERROR FLAG IN THE CASE OF BAD INTEGRAND BEHAVIOUR */
/*           AT INTERIOR POINTS OF INTEGRATION RANGE. */

/* Computing MAX */
	d__1 = abs(a1), d__2 = abs(b2);
	if (max(d__1,d__2) <= (epmach * 100. + 1.) * (abs(a2) + uflow * 1e3)) 
		{
	    *ier = 3;
	}

/*           APPEND THE NEWLY-CREATED INTERVALS TO THE LIST. */

L35:
	if (error2 > error1) {
	    goto L40;
	}
	alist__[*last] = a2;
	blist[maxerr] = b1;
	blist[*last] = b2;
	elist[maxerr] = error1;
	elist[*last] = error2;
	goto L50;
L40:
	alist__[maxerr] = a2;
	alist__[*last] = a1;
	blist[*last] = b1;
	rlist[maxerr] = area2;
	rlist[*last] = area1;
	elist[maxerr] = error2;
	elist[*last] = error1;

/*           CALL SUBROUTINE DQPSRT TO MAINTAIN THE DESCENDING ORDERING */
/*           IN THE LIST OF ERROR ESTIMATES AND SELECT THE SUBINTERVAL */
/*           WITH LARGEST ERROR ESTIMATE (TO BE BISECTED NEXT). */

L50:
	dqpsrt_(limit, last, &maxerr, &errmax, &elist[1], &iord[1], &nrmax);
/* ***JUMP OUT OF DO-LOOP */
	if (*ier != 0 || errsum <= errbnd) {
	    goto L70;
	}
/* L60: */
    }

/*           COMPUTE FINAL RESULT. */
/*           --------------------- */

L70:
    *result = 0.;
    i__1 = *last;
    for (k = 1; k <= i__1; ++k) {
	*result += rlist[k];
/* L80: */
    }
    *abserr = errsum;
L999:
    return 0;
} /* dqawse_ */

