/* qage.f -- translated by f2c (version 12.02.01).
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

/* DECK QAGE */
/* Subroutine */ int qage_(E_fp f, real *a, real *b, real *epsabs, real *
	epsrel, integer *key, integer *limit, real *result, real *abserr, 
	integer *neval, integer *ier, real *alist__, real *blist, real *rlist,
	 real *elist, integer *iord, integer *last)
{
    /* System generated locals */
    integer i__1;
    real r__1, r__2;

    /* Local variables */
    static integer k;
    static real a1, a2, b1, b2;
    extern /* Subroutine */ int qk21_(E_fp, real *, real *, real *, real *, 
	    real *, real *), qk31_(E_fp, real *, real *, real *, real *, real 
	    *, real *), qk15_(E_fp, real *, real *, real *, real *, real *, 
	    real *), qk41_(E_fp, real *, real *, real *, real *, real *, real 
	    *), qk51_(E_fp, real *, real *, real *, real *, real *, real *), 
	    qk61_(E_fp, real *, real *, real *, real *, real *, real *);
    static real area;
    static integer keyf;
    static real area1, area2, area12, erro12, defab1, defab2;
    static integer nrmax;
    static real uflow;
    extern /* Subroutine */ int qpsrt_(integer *, integer *, integer *, real *
	    , real *, integer *, integer *);
    extern doublereal r1mach_(integer *);
    static integer iroff1, iroff2;
    static real error1, error2, defabs, epmach, errbnd, resabs, errmax;
    static integer maxerr;
    static real errsum;

/* ***BEGIN PROLOGUE  QAGE */
/* ***PURPOSE  The routine calculates an approximation result to a given */
/*            definite integral   I = Integral of F over (A,B), */
/*            hopefully satisfying following claim for accuracy */
/*            ABS(I-RESLT).LE.MAX(EPSABS,EPSREL*ABS(I)). */
/* ***LIBRARY   SLATEC (QUADPACK) */
/* ***CATEGORY  H2A1A1 */
/* ***TYPE      SINGLE PRECISION (QAGE-S, DQAGE-D) */
/* ***KEYWORDS  AUTOMATIC INTEGRATOR, GAUSS-KRONROD RULES, */
/*             GENERAL-PURPOSE, GLOBALLY ADAPTIVE, INTEGRAND EXAMINATOR, */
/*             QUADPACK, QUADRATURE */
/* ***AUTHOR  Piessens, Robert */
/*             Applied Mathematics and Programming Division */
/*             K. U. Leuven */
/*           de Doncker, Elise */
/*             Applied Mathematics and Programming Division */
/*             K. U. Leuven */
/* ***DESCRIPTION */

/*        Computation of a definite integral */
/*        Standard fortran subroutine */
/*        Real version */

/*        PARAMETERS */
/*         ON ENTRY */
/*            F      - Real */
/*                     Function subprogram defining the integrand */
/*                     function F(X). The actual name for F needs to be */
/*                     declared E X T E R N A L in the driver program. */

/*            A      - Real */
/*                     Lower limit of integration */

/*            B      - Real */
/*                     Upper limit of integration */

/*            EPSABS - Real */
/*                     Absolute accuracy requested */
/*            EPSREL - Real */
/*                     Relative accuracy requested */
/*                     If  EPSABS.LE.0 */
/*                     and EPSREL.LT.MAX(50*REL.MACH.ACC.,0.5D-28), */
/*                     the routine will end with IER = 6. */

/*            KEY    - Integer */
/*                     Key for choice of local integration rule */
/*                     A Gauss-Kronrod pair is used with */
/*                          7 - 15 points if KEY.LT.2, */
/*                         10 - 21 points if KEY = 2, */
/*                         15 - 31 points if KEY = 3, */
/*                         20 - 41 points if KEY = 4, */
/*                         25 - 51 points if KEY = 5, */
/*                         30 - 61 points if KEY.GT.5. */

/*            LIMIT  - Integer */
/*                     Gives an upper bound on the number of subintervals */
/*                     in the partition of (A,B), LIMIT.GE.1. */

/*         ON RETURN */
/*            RESULT - Real */
/*                     Approximation to the integral */

/*            ABSERR - Real */
/*                     Estimate of the modulus of the absolute error, */
/*                     which should equal or exceed ABS(I-RESULT) */

/*            NEVAL  - Integer */
/*                     Number of integrand evaluations */

/*            IER    - Integer */
/*                     IER = 0 Normal and reliable termination of the */
/*                             routine. It is assumed that the requested */
/*                             accuracy has been achieved. */
/*                     IER.GT.0 Abnormal termination of the routine */
/*                             The estimates for result and error are */
/*                             less reliable. It is assumed that the */
/*                             requested accuracy has not been achieved. */
/*            ERROR MESSAGES */
/*                     IER = 1 Maximum number of subdivisions allowed */
/*                             has been achieved. One can allow more */
/*                             subdivisions by increasing the value */
/*                             of LIMIT. */
/*                             However, if this yields no improvement it */
/*                             is rather advised to analyze the integrand */
/*                             in order to determine the integration */
/*                             difficulties. If the position of a local */
/*                             difficulty can be determined(e.g. */
/*                             SINGULARITY, DISCONTINUITY within the */
/*                             interval) one will probably gain from */
/*                             splitting up the interval at this point */
/*                             and calling the integrator on the */
/*                             subranges. If possible, an appropriate */
/*                             special-purpose integrator should be used */
/*                             which is designed for handling the type of */
/*                             difficulty involved. */
/*                         = 2 The occurrence of roundoff error is */
/*                             detected, which prevents the requested */
/*                             tolerance from being achieved. */
/*                         = 3 Extremely bad integrand behaviour occurs */
/*                             at some points of the integration */
/*                             interval. */
/*                         = 6 The input is invalid, because */
/*                             (EPSABS.LE.0 and */
/*                              EPSREL.LT.MAX(50*REL.MACH.ACC.,0.5D-28), */
/*                             RESULT, ABSERR, NEVAL, LAST, RLIST(1) , */
/*                             ELIST(1) and IORD(1) are set to zero. */
/*                             ALIST(1) and BLIST(1) are set to A and B */
/*                             respectively. */

/*            ALIST   - Real */
/*                      Vector of dimension at least LIMIT, the first */
/*                       LAST  elements of which are the left */
/*                      end points of the subintervals in the partition */
/*                      of the given integration range (A,B) */

/*            BLIST   - Real */
/*                      Vector of dimension at least LIMIT, the first */
/*                       LAST  elements of which are the right */
/*                      end points of the subintervals in the partition */
/*                      of the given integration range (A,B) */

/*            RLIST   - Real */
/*                      Vector of dimension at least LIMIT, the first */
/*                       LAST  elements of which are the */
/*                      integral approximations on the subintervals */

/*            ELIST   - Real */
/*                      Vector of dimension at least LIMIT, the first */
/*                       LAST  elements of which are the moduli of the */
/*                      absolute error estimates on the subintervals */

/*            IORD    - Integer */
/*                      Vector of dimension at least LIMIT, the first K */
/*                      elements of which are pointers to the */
/*                      error estimates over the subintervals, */
/*                      such that ELIST(IORD(1)), ..., */
/*                      ELIST(IORD(K)) form a decreasing sequence, */
/*                      with K = LAST if LAST.LE.(LIMIT/2+2), and */
/*                      K = LIMIT+1-LAST otherwise */

/*            LAST    - Integer */
/*                      Number of subintervals actually produced in the */
/*                      subdivision process */

/* ***REFERENCES  (NONE) */
/* ***ROUTINES CALLED  QK15, QK21, QK31, QK41, QK51, QK61, QPSRT, R1MACH */
/* ***REVISION HISTORY  (YYMMDD) */
/*   800101  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890831  Modified array declarations.  (WRB) */
/*   890831  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/* ***END PROLOGUE  QAGE */




/*            LIST OF MAJOR VARIABLES */
/*            ----------------------- */

/*           ALIST     - LIST OF LEFT END POINTS OF ALL SUBINTERVALS */
/*                       CONSIDERED UP TO NOW */
/*           BLIST     - LIST OF RIGHT END POINTS OF ALL SUBINTERVALS */
/*                       CONSIDERED UP TO NOW */
/*           RLIST(I)  - APPROXIMATION TO THE INTEGRAL OVER */
/*                      (ALIST(I),BLIST(I)) */
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


/*           MACHINE DEPENDENT CONSTANTS */
/*           --------------------------- */

/*           EPMACH  IS THE LARGEST RELATIVE SPACING. */
/*           UFLOW  IS THE SMALLEST POSITIVE MAGNITUDE. */

/* ***FIRST EXECUTABLE STATEMENT  QAGE */
    /* Parameter adjustments */
    --iord;
    --elist;
    --rlist;
    --blist;
    --alist__;

    /* Function Body */
    epmach = r1mach_(&c__4);
    uflow = r1mach_(&c__1);

/*           TEST ON VALIDITY OF PARAMETERS */
/*           ------------------------------ */

    *ier = 0;
    *neval = 0;
    *last = 0;
    *result = 0.f;
    *abserr = 0.f;
    alist__[1] = *a;
    blist[1] = *b;
    rlist[1] = 0.f;
    elist[1] = 0.f;
    iord[1] = 0;
/* Computing MAX */
    r__1 = epmach * 50.f;
    if (*epsabs <= 0.f && *epsrel < dmax(r__1,5e-15f)) {
	*ier = 6;
    }
    if (*ier == 6) {
	goto L999;
    }

/*           FIRST APPROXIMATION TO THE INTEGRAL */
/*           ----------------------------------- */

    keyf = *key;
    if (*key <= 0) {
	keyf = 1;
    }
    if (*key >= 7) {
	keyf = 6;
    }
    *neval = 0;
    if (keyf == 1) {
	qk15_((E_fp)f, a, b, result, abserr, &defabs, &resabs);
    }
    if (keyf == 2) {
	qk21_((E_fp)f, a, b, result, abserr, &defabs, &resabs);
    }
    if (keyf == 3) {
	qk31_((E_fp)f, a, b, result, abserr, &defabs, &resabs);
    }
    if (keyf == 4) {
	qk41_((E_fp)f, a, b, result, abserr, &defabs, &resabs);
    }
    if (keyf == 5) {
	qk51_((E_fp)f, a, b, result, abserr, &defabs, &resabs);
    }
    if (keyf == 6) {
	qk61_((E_fp)f, a, b, result, abserr, &defabs, &resabs);
    }
    *last = 1;
    rlist[1] = *result;
    elist[1] = *abserr;
    iord[1] = 1;

/*           TEST ON ACCURACY. */

/* Computing MAX */
    r__1 = *epsabs, r__2 = *epsrel * dabs(*result);
    errbnd = dmax(r__1,r__2);
    if (*abserr <= epmach * 50.f * defabs && *abserr > errbnd) {
	*ier = 2;
    }
    if (*limit == 1) {
	*ier = 1;
    }
    if (*ier != 0 || *abserr <= errbnd && *abserr != resabs || *abserr == 0.f)
	     {
	goto L60;
    }

/*           INITIALIZATION */
/*           -------------- */


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

/*           BISECT THE SUBINTERVAL WITH THE LARGEST ERROR ESTIMATE. */

	a1 = alist__[maxerr];
	b1 = (alist__[maxerr] + blist[maxerr]) * .5f;
	a2 = b1;
	b2 = blist[maxerr];
	if (keyf == 1) {
	    qk15_((E_fp)f, &a1, &b1, &area1, &error1, &resabs, &defab1);
	}
	if (keyf == 2) {
	    qk21_((E_fp)f, &a1, &b1, &area1, &error1, &resabs, &defab1);
	}
	if (keyf == 3) {
	    qk31_((E_fp)f, &a1, &b1, &area1, &error1, &resabs, &defab1);
	}
	if (keyf == 4) {
	    qk41_((E_fp)f, &a1, &b1, &area1, &error1, &resabs, &defab1);
	}
	if (keyf == 5) {
	    qk51_((E_fp)f, &a1, &b1, &area1, &error1, &resabs, &defab1);
	}
	if (keyf == 6) {
	    qk61_((E_fp)f, &a1, &b1, &area1, &error1, &resabs, &defab1);
	}
	if (keyf == 1) {
	    qk15_((E_fp)f, &a2, &b2, &area2, &error2, &resabs, &defab2);
	}
	if (keyf == 2) {
	    qk21_((E_fp)f, &a2, &b2, &area2, &error2, &resabs, &defab2);
	}
	if (keyf == 3) {
	    qk31_((E_fp)f, &a2, &b2, &area2, &error2, &resabs, &defab2);
	}
	if (keyf == 4) {
	    qk41_((E_fp)f, &a2, &b2, &area2, &error2, &resabs, &defab2);
	}
	if (keyf == 5) {
	    qk51_((E_fp)f, &a2, &b2, &area2, &error2, &resabs, &defab2);
	}
	if (keyf == 6) {
	    qk61_((E_fp)f, &a2, &b2, &area2, &error2, &resabs, &defab2);
	}

/*           IMPROVE PREVIOUS APPROXIMATIONS TO INTEGRAL */
/*           AND ERROR AND TEST FOR ACCURACY. */

	++(*neval);
	area12 = area1 + area2;
	erro12 = error1 + error2;
	errsum = errsum + erro12 - errmax;
	area = area + area12 - rlist[maxerr];
	if (defab1 == error1 || defab2 == error2) {
	    goto L5;
	}
	if ((r__1 = rlist[maxerr] - area12, dabs(r__1)) <= dabs(area12) * 
		1e-5f && erro12 >= errmax * .99f) {
	    ++iroff1;
	}
	if (*last > 10 && erro12 > errmax) {
	    ++iroff2;
	}
L5:
	rlist[maxerr] = area1;
	rlist[*last] = area2;
/* Computing MAX */
	r__1 = *epsabs, r__2 = *epsrel * dabs(area);
	errbnd = dmax(r__1,r__2);
	if (errsum <= errbnd) {
	    goto L8;
	}

/*           TEST FOR ROUNDOFF ERROR AND EVENTUALLY */
/*           SET ERROR FLAG. */

	if (iroff1 >= 6 || iroff2 >= 20) {
	    *ier = 2;
	}

/*           SET ERROR FLAG IN THE CASE THAT THE NUMBER OF */
/*           SUBINTERVALS EQUALS LIMIT. */

	if (*last == *limit) {
	    *ier = 1;
	}

/*           SET ERROR FLAG IN THE CASE OF BAD INTEGRAND BEHAVIOUR */
/*           AT A POINT OF THE INTEGRATION RANGE. */

/* Computing MAX */
	r__1 = dabs(a1), r__2 = dabs(b2);
	if (dmax(r__1,r__2) <= (epmach * 100.f + 1.f) * (dabs(a2) + uflow * 
		1e3f)) {
	    *ier = 3;
	}

/*           APPEND THE NEWLY-CREATED INTERVALS TO THE LIST. */

L8:
	if (error2 > error1) {
	    goto L10;
	}
	alist__[*last] = a2;
	blist[maxerr] = b1;
	blist[*last] = b2;
	elist[maxerr] = error1;
	elist[*last] = error2;
	goto L20;
L10:
	alist__[maxerr] = a2;
	alist__[*last] = a1;
	blist[*last] = b1;
	rlist[maxerr] = area2;
	rlist[*last] = area1;
	elist[maxerr] = error2;
	elist[*last] = error1;

/*           CALL SUBROUTINE QPSRT TO MAINTAIN THE DESCENDING ORDERING */
/*           IN THE LIST OF ERROR ESTIMATES AND SELECT THE */
/*           SUBINTERVAL WITH THE LARGEST ERROR ESTIMATE (TO BE */
/*           BISECTED NEXT). */

L20:
	qpsrt_(limit, last, &maxerr, &errmax, &elist[1], &iord[1], &nrmax);
/* ***JUMP OUT OF DO-LOOP */
	if (*ier != 0 || errsum <= errbnd) {
	    goto L40;
	}
/* L30: */
    }

/*           COMPUTE FINAL RESULT. */
/*           --------------------- */

L40:
    *result = 0.f;
    i__1 = *last;
    for (k = 1; k <= i__1; ++k) {
	*result += rlist[k];
/* L50: */
    }
    *abserr = errsum;
L60:
    if (keyf != 1) {
	*neval = (keyf * 10 + 1) * ((*neval << 1) + 1);
    }
    if (keyf == 1) {
	*neval = *neval * 30 + 15;
    }
L999:
    return 0;
} /* qage_ */

