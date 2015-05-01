/* qagpe.f -- translated by f2c (version 12.02.01).
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
static integer c__2 = 2;

/* DECK QAGPE */
/* Subroutine */ int qagpe_(E_fp f, real *a, real *b, integer *npts2, real *
	points, real *epsabs, real *epsrel, integer *limit, real *result, 
	real *abserr, integer *neval, integer *ier, real *alist__, real *
	blist, real *rlist, real *elist, real *pts, integer *iord, integer *
	level, integer *ndin, integer *last)
{
    /* System generated locals */
    integer i__1, i__2;
    real r__1, r__2;

    /* Local variables */
    static integer i__, j, k;
    static real a1, a2, b1, b2;
    static integer id, ip1;
    extern /* Subroutine */ int qk21_(E_fp, real *, real *, real *, real *, 
	    real *, real *);
    static integer ind1, ind2;
    static real area;
    extern /* Subroutine */ int qelg_(integer *, real *, real *, real *, real 
	    *, integer *);
    static real resa, dres, sign;
    static integer ksgn;
    static real temp;
    static integer nres, nint, jlow, npts;
    static real area1, area2, area12, erro12;
    static integer ierro;
    static real defab1, defab2;
    static integer ktmin, nrmax;
    static real oflow, uflow;
    static logical noext;
    extern /* Subroutine */ int qpsrt_(integer *, integer *, integer *, real *
	    , real *, integer *, integer *);
    extern doublereal r1mach_(integer *);
    static integer iroff1, iroff2, iroff3;
    static real res3la[3];
    static integer nintp1;
    static real error1, error2, rlist2[52];
    static integer numrl2;
    static real defabs, epmach, erlarg, abseps, correc, errbnd, resabs;
    static integer jupbnd;
    static real erlast;
    static integer levmax;
    static real errmax;
    static integer maxerr, levcur;
    static real reseps;
    static logical extrap;
    static real ertest, errsum;

/* ***BEGIN PROLOGUE  QAGPE */
/* ***PURPOSE  Approximate a given definite integral I = Integral of F */
/*            over (A,B), hopefully satisfying the accuracy claim: */
/*                  ABS(I-RESULT).LE.MAX(EPSABS,EPSREL*ABS(I)). */
/*            Break points of the integration interval, where local */
/*            difficulties of the integrand may occur (e.g. singularities */
/*            or discontinuities) are provided by the user. */
/* ***LIBRARY   SLATEC (QUADPACK) */
/* ***CATEGORY  H2A2A1 */
/* ***TYPE      SINGLE PRECISION (QAGPE-S, DQAGPE-D) */
/* ***KEYWORDS  AUTOMATIC INTEGRATOR, EXTRAPOLATION, GENERAL-PURPOSE, */
/*             GLOBALLY ADAPTIVE, QUADPACK, QUADRATURE, */
/*             SINGULARITIES AT USER SPECIFIED POINTS */
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

/*            NPTS2  - Integer */
/*                     Number equal to two more than the number of */
/*                     user-supplied break points within the integration */
/*                     range, NPTS2.GE.2. */
/*                     If NPTS2.LT.2, the routine will end with IER = 6. */

/*            POINTS - Real */
/*                     Vector of dimension NPTS2, the first (NPTS2-2) */
/*                     elements of which are the user provided break */
/*                     POINTS. If these POINTS do not constitute an */
/*                     ascending sequence there will be an automatic */
/*                     sorting. */

/*            EPSABS - Real */
/*                     Absolute accuracy requested */
/*            EPSREL - Real */
/*                     Relative accuracy requested */
/*                     If  EPSABS.LE.0 */
/*                     and EPSREL.LT.MAX(50*REL.MACH.ACC.,0.5D-28), */
/*                     the routine will end with IER = 6. */

/*            LIMIT  - Integer */
/*                     Gives an upper bound on the number of subintervals */
/*                     in the partition of (A,B), LIMIT.GE.NPTS2 */
/*                     If LIMIT.LT.NPTS2, the routine will end with */
/*                     IER = 6. */

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
/*                     IER.GT.0 Abnormal termination of the routine. */
/*                             The estimates for integral and error are */
/*                             less reliable. It is assumed that the */
/*                             requested accuracy has not been achieved. */
/*            ERROR MESSAGES */
/*                     IER = 1 Maximum number of subdivisions allowed */
/*                             has been achieved. One can allow more */
/*                             subdivisions by increasing the value of */
/*                             LIMIT (and taking the according dimension */
/*                             adjustments into account). However, if */
/*                             this yields no improvement it is advised */
/*                             to analyze the integrand in order to */
/*                             determine the integration difficulties. If */
/*                             the position of a local difficulty can be */
/*                             determined (i.e. SINGULARITY, */
/*                             DISCONTINUITY within the interval), it */
/*                             should be supplied to the routine as an */
/*                             element of the vector points. If necessary */
/*                             an appropriate special-purpose integrator */
/*                             must be used, which is designed for */
/*                             handling the type of difficulty involved. */
/*                         = 2 The occurrence of roundoff error is */
/*                             detected, which prevents the requested */
/*                             tolerance from being achieved. */
/*                             The error may be under-estimated. */
/*                         = 3 Extremely bad integrand behaviour occurs */
/*                             At some points of the integration */
/*                             interval. */
/*                         = 4 The algorithm does not converge. */
/*                             Roundoff error is detected in the */
/*                             extrapolation table. It is presumed that */
/*                             the requested tolerance cannot be */
/*                             achieved, and that the returned result is */
/*                             the best which can be obtained. */
/*                         = 5 The integral is probably divergent, or */
/*                             slowly convergent. It must be noted that */
/*                             divergence can occur with any other value */
/*                             of IER.GT.0. */
/*                         = 6 The input is invalid because */
/*                             NPTS2.LT.2 or */
/*                             Break points are specified outside */
/*                             the integration range or */
/*                             (EPSABS.LE.0 and */
/*                              EPSREL.LT.MAX(50*REL.MACH.ACC.,0.5D-28)) */
/*                             or LIMIT.LT.NPTS2. */
/*                             RESULT, ABSERR, NEVAL, LAST, RLIST(1), */
/*                             and ELIST(1) are set to zero. ALIST(1) and */
/*                             BLIST(1) are set to A and B respectively. */

/*            ALIST  - Real */
/*                     Vector of dimension at least LIMIT, the first */
/*                      LAST  elements of which are the left end points */
/*                     of the subintervals in the partition of the given */
/*                     integration range (A,B) */

/*            BLIST  - Real */
/*                     Vector of dimension at least LIMIT, the first */
/*                      LAST  elements of which are the right end points */
/*                     of the subintervals in the partition of the given */
/*                     integration range (A,B) */

/*            RLIST  - Real */
/*                     Vector of dimension at least LIMIT, the first */
/*                      LAST  elements of which are the integral */
/*                     approximations on the subintervals */

/*            ELIST  - Real */
/*                     Vector of dimension at least LIMIT, the first */
/*                      LAST  elements of which are the moduli of the */
/*                     absolute error estimates on the subintervals */

/*            PTS    - Real */
/*                     Vector of dimension at least NPTS2, containing the */
/*                     integration limits and the break points of the */
/*                     interval in ascending sequence. */

/*            LEVEL  - Integer */
/*                     Vector of dimension at least LIMIT, containing the */
/*                     subdivision levels of the subinterval, i.e. if */
/*                     (AA,BB) is a subinterval of (P1,P2) where P1 as */
/*                     well as P2 is a user-provided break point or */
/*                     integration limit, then (AA,BB) has level L if */
/*                     ABS(BB-AA) = ABS(P2-P1)*2**(-L). */

/*            NDIN   - Integer */
/*                     Vector of dimension at least NPTS2, after first */
/*                     integration over the intervals (PTS(I)),PTS(I+1), */
/*                     I = 0,1, ..., NPTS2-2, the error estimates over */
/*                     some of the intervals may have been increased */
/*                     artificially, in order to put their subdivision */
/*                     forward. If this happens for the subinterval */
/*                     numbered K, NDIN(K) is put to 1, otherwise */
/*                     NDIN(K) = 0. */

/*            IORD   - Integer */
/*                     Vector of dimension at least LIMIT, the first K */
/*                     elements of which are pointers to the */
/*                     error estimates over the subintervals, */
/*                     such that ELIST(IORD(1)), ..., ELIST(IORD(K)) */
/*                     form a decreasing sequence, with K = LAST */
/*                     If LAST.LE.(LIMIT/2+2), and K = LIMIT+1-LAST */
/*                     otherwise */

/*            LAST   - Integer */
/*                     Number of subintervals actually produced in the */
/*                     subdivisions process */

/* ***REFERENCES  (NONE) */
/* ***ROUTINES CALLED  QELG, QK21, QPSRT, R1MACH */
/* ***REVISION HISTORY  (YYMMDD) */
/*   800101  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890831  Modified array declarations.  (WRB) */
/*   890831  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/* ***END PROLOGUE  QAGPE */




/*            THE DIMENSION OF RLIST2 IS DETERMINED BY THE VALUE OF */
/*            LIMEXP IN SUBROUTINE EPSALG (RLIST2 SHOULD BE OF DIMENSION */
/*            (LIMEXP+2) AT LEAST). */


/*            LIST OF MAJOR VARIABLES */
/*            ----------------------- */

/*           ALIST     - LIST OF LEFT END POINTS OF ALL SUBINTERVALS */
/*                       CONSIDERED UP TO NOW */
/*           BLIST     - LIST OF RIGHT END POINTS OF ALL SUBINTERVALS */
/*                       CONSIDERED UP TO NOW */
/*           RLIST(I)  - APPROXIMATION TO THE INTEGRAL OVER */
/*                       (ALIST(I),BLIST(I)) */
/*           RLIST2    - ARRAY OF DIMENSION AT LEAST LIMEXP+2 */
/*                       CONTAINING THE PART OF THE EPSILON TABLE WHICH */
/*                       IS STILL NEEDED FOR FURTHER COMPUTATIONS */
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
/*           NUMRL2    - NUMBER OF ELEMENTS IN RLIST2. IF AN */
/*                       APPROPRIATE APPROXIMATION TO THE COMPOUNDED */
/*                       INTEGRAL HAS BEEN OBTAINED, IT IS PUT IN */
/*                       RLIST2(NUMRL2) AFTER NUMRL2 HAS BEEN INCREASED */
/*                       BY ONE. */
/*           ERLARG    - SUM OF THE ERRORS OVER THE INTERVALS LARGER */
/*                       THAN THE SMALLEST INTERVAL CONSIDERED UP TO NOW */
/*           EXTRAP    - LOGICAL VARIABLE DENOTING THAT THE ROUTINE */
/*                       IS ATTEMPTING TO PERFORM EXTRAPOLATION. I.E. */
/*                       BEFORE SUBDIVIDING THE SMALLEST INTERVAL WE */
/*                       TRY TO DECREASE THE VALUE OF ERLARG. */
/*           NOEXT     - LOGICAL VARIABLE DENOTING THAT EXTRAPOLATION IS */
/*                       NO LONGER ALLOWED (TRUE-VALUE) */

/*            MACHINE DEPENDENT CONSTANTS */
/*            --------------------------- */

/*           EPMACH IS THE LARGEST RELATIVE SPACING. */
/*           UFLOW IS THE SMALLEST POSITIVE MAGNITUDE. */
/*           OFLOW IS THE LARGEST POSITIVE MAGNITUDE. */

/* ***FIRST EXECUTABLE STATEMENT  QAGPE */
    /* Parameter adjustments */
    --ndin;
    --level;
    --iord;
    --pts;
    --elist;
    --rlist;
    --blist;
    --alist__;
    --points;

    /* Function Body */
    epmach = r1mach_(&c__4);

/*            TEST ON VALIDITY OF PARAMETERS */
/*            ----------------------------- */

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
    level[1] = 0;
    npts = *npts2 - 2;
/* Computing MAX */
    r__1 = epmach * 50.f;
    if (*npts2 < 2 || *limit <= npts || *epsabs <= 0.f && *epsrel < dmax(r__1,
	    5e-15f)) {
	*ier = 6;
    }
    if (*ier == 6) {
	goto L999;
    }

/*            IF ANY BREAK POINTS ARE PROVIDED, SORT THEM INTO AN */
/*            ASCENDING SEQUENCE. */

    sign = 1.f;
    if (*a > *b) {
	sign = -1.f;
    }
    pts[1] = dmin(*a,*b);
    if (npts == 0) {
	goto L15;
    }
    i__1 = npts;
    for (i__ = 1; i__ <= i__1; ++i__) {
	pts[i__ + 1] = points[i__];
/* L10: */
    }
L15:
    pts[npts + 2] = dmax(*a,*b);
    nint = npts + 1;
    a1 = pts[1];
    if (npts == 0) {
	goto L40;
    }
    nintp1 = nint + 1;
    i__1 = nint;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ip1 = i__ + 1;
	i__2 = nintp1;
	for (j = ip1; j <= i__2; ++j) {
	    if (pts[i__] <= pts[j]) {
		goto L20;
	    }
	    temp = pts[i__];
	    pts[i__] = pts[j];
	    pts[j] = temp;
L20:
	    ;
	}
    }
    if (pts[1] != dmin(*a,*b) || pts[nintp1] != dmax(*a,*b)) {
	*ier = 6;
    }
    if (*ier == 6) {
	goto L999;
    }

/*            COMPUTE FIRST INTEGRAL AND ERROR APPROXIMATIONS. */
/*            ------------------------------------------------ */

L40:
    resabs = 0.f;
    i__2 = nint;
    for (i__ = 1; i__ <= i__2; ++i__) {
	b1 = pts[i__ + 1];
	qk21_((E_fp)f, &a1, &b1, &area1, &error1, &defabs, &resa);
	*abserr += error1;
	*result += area1;
	ndin[i__] = 0;
	if (error1 == resa && error1 != 0.f) {
	    ndin[i__] = 1;
	}
	resabs += defabs;
	level[i__] = 0;
	elist[i__] = error1;
	alist__[i__] = a1;
	blist[i__] = b1;
	rlist[i__] = area1;
	iord[i__] = i__;
	a1 = b1;
/* L50: */
    }
    errsum = 0.f;
    i__2 = nint;
    for (i__ = 1; i__ <= i__2; ++i__) {
	if (ndin[i__] == 1) {
	    elist[i__] = *abserr;
	}
	errsum += elist[i__];
/* L55: */
    }

/*           TEST ON ACCURACY. */

    *last = nint;
    *neval = nint * 21;
    dres = dabs(*result);
/* Computing MAX */
    r__1 = *epsabs, r__2 = *epsrel * dres;
    errbnd = dmax(r__1,r__2);
    if (*abserr <= epmach * 100.f * resabs && *abserr > errbnd) {
	*ier = 2;
    }
    if (nint == 1) {
	goto L80;
    }
    i__2 = npts;
    for (i__ = 1; i__ <= i__2; ++i__) {
	jlow = i__ + 1;
	ind1 = iord[i__];
	i__1 = nint;
	for (j = jlow; j <= i__1; ++j) {
	    ind2 = iord[j];
	    if (elist[ind1] > elist[ind2]) {
		goto L60;
	    }
	    ind1 = ind2;
	    k = j;
L60:
	    ;
	}
	if (ind1 == iord[i__]) {
	    goto L70;
	}
	iord[k] = iord[i__];
	iord[i__] = ind1;
L70:
	;
    }
    if (*limit < *npts2) {
	*ier = 1;
    }
L80:
    if (*ier != 0 || *abserr <= errbnd) {
	goto L999;
    }

/*           INITIALIZATION */
/*           -------------- */

    rlist2[0] = *result;
    maxerr = iord[1];
    errmax = elist[maxerr];
    area = *result;
    nrmax = 1;
    nres = 0;
    numrl2 = 1;
    ktmin = 0;
    extrap = FALSE_;
    noext = FALSE_;
    erlarg = errsum;
    ertest = errbnd;
    levmax = 1;
    iroff1 = 0;
    iroff2 = 0;
    iroff3 = 0;
    ierro = 0;
    uflow = r1mach_(&c__1);
    oflow = r1mach_(&c__2);
    *abserr = oflow;
    ksgn = -1;
    if (dres >= (1.f - epmach * 50.f) * resabs) {
	ksgn = 1;
    }

/*           MAIN DO-LOOP */
/*           ------------ */

    i__2 = *limit;
    for (*last = *npts2; *last <= i__2; ++(*last)) {

/*           BISECT THE SUBINTERVAL WITH THE NRMAX-TH LARGEST */
/*           ERROR ESTIMATE. */

	levcur = level[maxerr] + 1;
	a1 = alist__[maxerr];
	b1 = (alist__[maxerr] + blist[maxerr]) * .5f;
	a2 = b1;
	b2 = blist[maxerr];
	erlast = errmax;
	qk21_((E_fp)f, &a1, &b1, &area1, &error1, &resa, &defab1);
	qk21_((E_fp)f, &a2, &b2, &area2, &error2, &resa, &defab2);

/*           IMPROVE PREVIOUS APPROXIMATIONS TO INTEGRAL */
/*           AND ERROR AND TEST FOR ACCURACY. */

	*neval += 42;
	area12 = area1 + area2;
	erro12 = error1 + error2;
	errsum = errsum + erro12 - errmax;
	area = area + area12 - rlist[maxerr];
	if (defab1 == error1 || defab2 == error2) {
	    goto L95;
	}
	if ((r__1 = rlist[maxerr] - area12, dabs(r__1)) > dabs(area12) * 
		1e-5f || erro12 < errmax * .99f) {
	    goto L90;
	}
	if (extrap) {
	    ++iroff2;
	}
	if (! extrap) {
	    ++iroff1;
	}
L90:
	if (*last > 10 && erro12 > errmax) {
	    ++iroff3;
	}
L95:
	level[maxerr] = levcur;
	level[*last] = levcur;
	rlist[maxerr] = area1;
	rlist[*last] = area2;
/* Computing MAX */
	r__1 = *epsabs, r__2 = *epsrel * dabs(area);
	errbnd = dmax(r__1,r__2);

/*           TEST FOR ROUNDOFF ERROR AND EVENTUALLY */
/*           SET ERROR FLAG. */

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
/*           AT A POINT OF THE INTEGRATION RANGE */

/* Computing MAX */
	r__1 = dabs(a1), r__2 = dabs(b2);
	if (dmax(r__1,r__2) <= (epmach * 100.f + 1.f) * (dabs(a2) + uflow * 
		1e3f)) {
	    *ier = 4;
	}

/*           APPEND THE NEWLY-CREATED INTERVALS TO THE LIST. */

	if (error2 > error1) {
	    goto L100;
	}
	alist__[*last] = a2;
	blist[maxerr] = b1;
	blist[*last] = b2;
	elist[maxerr] = error1;
	elist[*last] = error2;
	goto L110;
L100:
	alist__[maxerr] = a2;
	alist__[*last] = a1;
	blist[*last] = b1;
	rlist[maxerr] = area2;
	rlist[*last] = area1;
	elist[maxerr] = error2;
	elist[*last] = error1;

/*           CALL SUBROUTINE QPSRT TO MAINTAIN THE DESCENDING ORDERING */
/*           IN THE LIST OF ERROR ESTIMATES AND SELECT THE */
/*           SUBINTERVAL WITH NRMAX-TH LARGEST ERROR ESTIMATE (TO BE */
/*           BISECTED NEXT). */

L110:
	qpsrt_(limit, last, &maxerr, &errmax, &elist[1], &iord[1], &nrmax);
/* ***JUMP OUT OF DO-LOOP */
	if (errsum <= errbnd) {
	    goto L190;
	}
/* ***JUMP OUT OF DO-LOOP */
	if (*ier != 0) {
	    goto L170;
	}
	if (noext) {
	    goto L160;
	}
	erlarg -= erlast;
	if (levcur + 1 <= levmax) {
	    erlarg += erro12;
	}
	if (extrap) {
	    goto L120;
	}

/*           TEST WHETHER THE INTERVAL TO BE BISECTED NEXT IS THE */
/*           SMALLEST INTERVAL. */

	if (level[maxerr] + 1 <= levmax) {
	    goto L160;
	}
	extrap = TRUE_;
	nrmax = 2;
L120:
	if (ierro == 3 || erlarg <= ertest) {
	    goto L140;
	}

/*           THE SMALLEST INTERVAL HAS THE LARGEST ERROR. */
/*           BEFORE BISECTING DECREASE THE SUM OF THE ERRORS */
/*           OVER THE LARGER INTERVALS (ERLARG) AND PERFORM */
/*           EXTRAPOLATION. */

	id = nrmax;
	jupbnd = *last;
	if (*last > *limit / 2 + 2) {
	    jupbnd = *limit + 3 - *last;
	}
	i__1 = jupbnd;
	for (k = id; k <= i__1; ++k) {
	    maxerr = iord[nrmax];
	    errmax = elist[maxerr];
/* ***JUMP OUT OF DO-LOOP */
	    if (level[maxerr] + 1 <= levmax) {
		goto L160;
	    }
	    ++nrmax;
/* L130: */
	}

/*           PERFORM EXTRAPOLATION. */

L140:
	++numrl2;
	rlist2[numrl2 - 1] = area;
	if (numrl2 <= 2) {
	    goto L155;
	}
	qelg_(&numrl2, rlist2, &reseps, &abseps, res3la, &nres);
	++ktmin;
	if (ktmin > 5 && *abserr < errsum * .001f) {
	    *ier = 5;
	}
	if (abseps >= *abserr) {
	    goto L150;
	}
	ktmin = 0;
	*abserr = abseps;
	*result = reseps;
	correc = erlarg;
/* Computing MAX */
	r__1 = *epsabs, r__2 = *epsrel * dabs(reseps);
	ertest = dmax(r__1,r__2);
/* ***JUMP OUT OF DO-LOOP */
	if (*abserr < ertest) {
	    goto L170;
	}

/*           PREPARE BISECTION OF THE SMALLEST INTERVAL. */

L150:
	if (numrl2 == 1) {
	    noext = TRUE_;
	}
	if (*ier >= 5) {
	    goto L170;
	}
L155:
	maxerr = iord[1];
	errmax = elist[maxerr];
	nrmax = 1;
	extrap = FALSE_;
	++levmax;
	erlarg = errsum;
L160:
	;
    }

/*           SET THE FINAL RESULT. */
/*           --------------------- */


L170:
    if (*abserr == oflow) {
	goto L190;
    }
    if (*ier + ierro == 0) {
	goto L180;
    }
    if (ierro == 3) {
	*abserr += correc;
    }
    if (*ier == 0) {
	*ier = 3;
    }
    if (*result != 0.f && area != 0.f) {
	goto L175;
    }
    if (*abserr > errsum) {
	goto L190;
    }
    if (area == 0.f) {
	goto L210;
    }
    goto L180;
L175:
    if (*abserr / dabs(*result) > errsum / dabs(area)) {
	goto L190;
    }

/*           TEST ON DIVERGENCE. */

L180:
/* Computing MAX */
    r__1 = dabs(*result), r__2 = dabs(area);
    if (ksgn == -1 && dmax(r__1,r__2) <= defabs * .01f) {
	goto L210;
    }
    if (.01f > *result / area || *result / area > 100.f || errsum > dabs(area)
	    ) {
	*ier = 6;
    }
    goto L210;

/*           COMPUTE GLOBAL INTEGRAL SUM. */

L190:
    *result = 0.f;
    i__2 = *last;
    for (k = 1; k <= i__2; ++k) {
	*result += rlist[k];
/* L200: */
    }
    *abserr = errsum;
L210:
    if (*ier > 2) {
	--(*ier);
    }
    *result *= sign;
L999:
    return 0;
} /* qagpe_ */

