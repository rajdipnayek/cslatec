/* dqawce.f -- translated by f2c (version 12.02.01).
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

/* DECK DQAWCE */
/* Subroutine */ int dqawce_(D_fp f, doublereal *a, doublereal *b, doublereal 
	*c__, doublereal *epsabs, doublereal *epsrel, integer *limit, 
	doublereal *result, doublereal *abserr, integer *neval, integer *ier, 
	doublereal *alist__, doublereal *blist, doublereal *rlist, doublereal 
	*elist, integer *iord, integer *last)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2;

    /* Local variables */
    static integer k;
    static doublereal a1, a2, b1, b2, aa, bb;
    static integer nev;
    static doublereal area, area1, area2, area12;
    extern /* Subroutine */ int dqc25c_(D_fp, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, integer *, integer *);
    static doublereal erro12;
    static integer krule, nrmax;
    static doublereal uflow;
    extern doublereal d1mach_(integer *);
    static integer iroff1, iroff2;
    static doublereal error1, error2, epmach, errbnd, errmax;
    static integer maxerr;
    static doublereal errsum;
    extern /* Subroutine */ int dqpsrt_(integer *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, integer *);

/* ***BEGIN PROLOGUE  DQAWCE */
/* ***PURPOSE  The routine calculates an approximation result to a */
/*            CAUCHY PRINCIPAL VALUE I = Integral of F*W over (A,B) */
/*            (W(X) = 1/(X-C), (C.NE.A, C.NE.B), hopefully satisfying */
/*            following claim for accuracy */
/*            ABS(I-RESULT).LE.MAX(EPSABS,EPSREL*ABS(I)) */
/* ***LIBRARY   SLATEC (QUADPACK) */
/* ***CATEGORY  H2A2A1, J4 */
/* ***TYPE      DOUBLE PRECISION (QAWCE-S, DQAWCE-D) */
/* ***KEYWORDS  AUTOMATIC INTEGRATOR, CAUCHY PRINCIPAL VALUE, */
/*             CLENSHAW-CURTIS METHOD, QUADPACK, QUADRATURE, */
/*             SPECIAL-PURPOSE */
/* ***AUTHOR  Piessens, Robert */
/*             Applied Mathematics and Programming Division */
/*             K. U. Leuven */
/*           de Doncker, Elise */
/*             Applied Mathematics and Programming Division */
/*             K. U. Leuven */
/* ***DESCRIPTION */

/*        Computation of a CAUCHY PRINCIPAL VALUE */
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
/*                     Upper limit of integration */

/*            C      - Double precision */
/*                     Parameter in the WEIGHT function, C.NE.A, C.NE.B */
/*                     If C = A OR C = B, the routine will end with */
/*                     IER = 6. */

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
/*                     IER.GT.0 Abnormal termination of the routine */
/*                             the estimates for integral and error are */
/*                             less reliable. It is assumed that the */
/*                             requested accuracy has not been achieved. */
/*            ERROR MESSAGES */
/*                     IER = 1 Maximum number of subdivisions allowed */
/*                             has been achieved. One can allow more sub- */
/*                             divisions by increasing the value of */
/*                             LIMIT. However, if this yields no */
/*                             improvement it is advised to analyze the */
/*                             the integrand, in order to determine the */
/*                             the integration difficulties. If the */
/*                             position of a local difficulty can be */
/*                             determined (e.g. SINGULARITY, */
/*                             DISCONTINUITY within the interval) one */
/*                             will probably gain from splitting up the */
/*                             interval at this point and calling */
/*                             appropriate integrators on the subranges. */
/*                         = 2 The occurrence of roundoff error is detec- */
/*                             ted, which prevents the requested */
/*                             tolerance from being achieved. */
/*                         = 3 Extremely bad integrand behaviour */
/*                             occurs at some interior points of */
/*                             the integration interval. */
/*                         = 6 The input is invalid, because */
/*                             C = A or C = B or */
/*                             (EPSABS.LE.0 and */
/*                              EPSREL.LT.MAX(50*REL.MACH.ACC.,0.5D-28)) */
/*                             or LIMIT.LT.1. */
/*                             RESULT, ABSERR, NEVAL, RLIST(1), ELIST(1), */
/*                             IORD(1) and LAST are set to zero. ALIST(1) */
/*                             and BLIST(1) are set to A and B */
/*                             respectively. */

/*            ALIST   - Double precision */
/*                      Vector of dimension at least LIMIT, the first */
/*                       LAST  elements of which are the left */
/*                      end points of the subintervals in the partition */
/*                      of the given integration range (A,B) */

/*            BLIST   - Double precision */
/*                      Vector of dimension at least LIMIT, the first */
/*                       LAST  elements of which are the right */
/*                      end points of the subintervals in the partition */
/*                      of the given integration range (A,B) */

/*            RLIST   - Double precision */
/*                      Vector of dimension at least LIMIT, the first */
/*                       LAST  elements of which are the integral */
/*                      approximations on the subintervals */

/*            ELIST   - Double precision */
/*                      Vector of dimension LIMIT, the first  LAST */
/*                      elements of which are the moduli of the absolute */
/*                      error estimates on the subintervals */

/*            IORD    - Integer */
/*                      Vector of dimension at least LIMIT, the first K */
/*                      elements of which are pointers to the error */
/*                      estimates over the subintervals, so that */
/*                      ELIST(IORD(1)), ..., ELIST(IORD(K)) with K = LAST */
/*                      If LAST.LE.(LIMIT/2+2), and K = LIMIT+1-LAST */
/*                      otherwise, form a decreasing sequence */

/*            LAST    - Integer */
/*                      Number of subintervals actually produced in */
/*                      the subdivision process */

/* ***REFERENCES  (NONE) */
/* ***ROUTINES CALLED  D1MACH, DQC25C, DQPSRT */
/* ***REVISION HISTORY  (YYMMDD) */
/*   800101  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890831  Modified array declarations.  (WRB) */
/*   890831  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/* ***END PROLOGUE  DQAWCE */




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

/* ***FIRST EXECUTABLE STATEMENT  DQAWCE */
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
    alist__[1] = *a;
    blist[1] = *b;
    rlist[1] = 0.;
    elist[1] = 0.;
    iord[1] = 0;
    *result = 0.;
    *abserr = 0.;
/* Computing MAX */
    d__1 = epmach * 50.;
    if (*c__ == *a || *c__ == *b || *epsabs <= 0. && *epsrel < max(d__1,5e-29)
	    ) {
	goto L999;
    }

/*           FIRST APPROXIMATION TO THE INTEGRAL */
/*           ----------------------------------- */

    aa = *a;
    bb = *b;
    if (*a <= *b) {
	goto L10;
    }
    aa = *b;
    bb = *a;
L10:
    *ier = 0;
    krule = 1;
    dqc25c_((D_fp)f, &aa, &bb, c__, result, abserr, &krule, neval);
    *last = 1;
    rlist[1] = *result;
    elist[1] = *abserr;
    iord[1] = 1;
    alist__[1] = *a;
    blist[1] = *b;

/*           TEST ON ACCURACY */

/* Computing MAX */
    d__1 = *epsabs, d__2 = *epsrel * abs(*result);
    errbnd = max(d__1,d__2);
    if (*limit == 1) {
	*ier = 1;
    }
/* Computing MIN */
    d__1 = abs(*result) * .01;
    if (*abserr < min(d__1,errbnd) || *ier == 1) {
	goto L70;
    }

/*           INITIALIZATION */
/*           -------------- */

    alist__[1] = aa;
    blist[1] = bb;
    rlist[1] = *result;
    errmax = *abserr;
    maxerr = 1;
    area = *result;
    errsum = *abserr;
    nrmax = 1;
    iroff1 = 0;
    iroff2 = 0;

/*           MAIN DO-LOOP */
/*           ------------ */

    i__1 = *limit;
    for (*last = 2; *last <= i__1; ++(*last)) {

/*           BISECT THE SUBINTERVAL WITH NRMAX-TH LARGEST */
/*           ERROR ESTIMATE. */

	a1 = alist__[maxerr];
	b1 = (alist__[maxerr] + blist[maxerr]) * .5;
	b2 = blist[maxerr];
	if (*c__ <= b1 && *c__ > a1) {
	    b1 = (*c__ + b2) * .5;
	}
	if (*c__ > b1 && *c__ < b2) {
	    b1 = (a1 + *c__) * .5;
	}
	a2 = b1;
	krule = 2;
	dqc25c_((D_fp)f, &a1, &b1, c__, &area1, &error1, &krule, &nev);
	*neval += nev;
	dqc25c_((D_fp)f, &a2, &b2, c__, &area2, &error2, &krule, &nev);
	*neval += nev;

/*           IMPROVE PREVIOUS APPROXIMATIONS TO INTEGRAL */
/*           AND ERROR AND TEST FOR ACCURACY. */

	area12 = area1 + area2;
	erro12 = error1 + error2;
	errsum = errsum + erro12 - errmax;
	area = area + area12 - rlist[maxerr];
	if ((d__1 = rlist[maxerr] - area12, abs(d__1)) < abs(area12) * 1e-5 &&
		 erro12 >= errmax * .99 && krule == 0) {
	    ++iroff1;
	}
	if (*last > 10 && erro12 > errmax && krule == 0) {
	    ++iroff2;
	}
	rlist[maxerr] = area1;
	rlist[*last] = area2;
/* Computing MAX */
	d__1 = *epsabs, d__2 = *epsrel * abs(area);
	errbnd = max(d__1,d__2);
	if (errsum <= errbnd) {
	    goto L15;
	}

/*           TEST FOR ROUNDOFF ERROR AND EVENTUALLY SET ERROR FLAG. */

	if (iroff1 >= 6 && iroff2 > 20) {
	    *ier = 2;
	}

/*           SET ERROR FLAG IN THE CASE THAT NUMBER OF INTERVAL */
/*           BISECTIONS EXCEEDS LIMIT. */

	if (*last == *limit) {
	    *ier = 1;
	}

/*           SET ERROR FLAG IN THE CASE OF BAD INTEGRAND BEHAVIOUR */
/*           AT A POINT OF THE INTEGRATION RANGE. */

/* Computing MAX */
	d__1 = abs(a1), d__2 = abs(b2);
	if (max(d__1,d__2) <= (epmach * 100. + 1.) * (abs(a2) + uflow * 1e3)) 
		{
	    *ier = 3;
	}

/*           APPEND THE NEWLY-CREATED INTERVALS TO THE LIST. */

L15:
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
/* ***JUMP OUT OF DO-LOOP */
	if (*ier != 0 || errsum <= errbnd) {
	    goto L50;
	}
/* L40: */
    }

/*           COMPUTE FINAL RESULT. */
/*           --------------------- */

L50:
    *result = 0.;
    i__1 = *last;
    for (k = 1; k <= i__1; ++k) {
	*result += rlist[k];
/* L60: */
    }
    *abserr = errsum;
L70:
    if (aa == *b) {
	*result = -(*result);
    }
L999:
    return 0;
} /* dqawce_ */

