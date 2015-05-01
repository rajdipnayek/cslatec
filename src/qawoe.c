/* qawoe.f -- translated by f2c (version 12.02.01).
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
static integer c__0 = 0;
static integer c__1 = 1;
static integer c__2 = 2;

/* DECK QAWOE */
/* Subroutine */ int qawoe_(E_fp f, real *a, real *b, real *omega, integer *
	integr, real *epsabs, real *epsrel, integer *limit, integer *icall, 
	integer *maxp1, real *result, real *abserr, integer *neval, integer *
	ier, integer *last, real *alist__, real *blist, real *rlist, real *
	elist, integer *iord, integer *nnlog, integer *momcom, real *chebmo)
{
    /* System generated locals */
    integer chebmo_dim1, chebmo_offset, i__1, i__2;
    real r__1, r__2;

    /* Local variables */
    static integer k;
    static real a1, a2, b1, b2;
    static integer id, nev;
    static real area;
    extern /* Subroutine */ int qc25f_(E_fp, real *, real *, real *, integer *
	    , integer *, integer *, integer *, real *, real *, integer *, 
	    real *, real *, integer *, real *), qelg_(integer *, real *, real 
	    *, real *, real *, integer *);
    static real dres;
    static integer ksgn, nres;
    static real area1, area2, area12, small, erro12, defab1, defab2, width;
    static integer ierro;
    static real oflow;
    static integer ktmin, nrmax, nrmom;
    static real uflow;
    static logical noext;
    extern /* Subroutine */ int qpsrt_(integer *, integer *, integer *, real *
	    , real *, integer *, integer *);
    extern doublereal r1mach_(integer *);
    static integer iroff1, iroff2, iroff3;
    static real res3la[3], error1, error2, rlist2[52];
    static integer numrl2;
    static real defabs, domega, epmach, erlarg, abseps, correc, errbnd, 
	    resabs;
    static integer jupbnd;
    static logical extall;
    static real erlast, errmax;
    static integer maxerr;
    static real reseps;
    static logical extrap;
    static real ertest, errsum;

/* ***BEGIN PROLOGUE  QAWOE */
/* ***PURPOSE  Calculate an approximation to a given definite integral */
/*               I = Integral of F(X)*W(X) over (A,B), where */
/*                  W(X) = COS(OMEGA*X) */
/*               or W(X) = SIN(OMEGA*X), */
/*            hopefully satisfying the following claim for accuracy */
/*               ABS(I-RESULT).LE.MAX(EPSABS,EPSREL*ABS(I)). */
/* ***LIBRARY   SLATEC (QUADPACK) */
/* ***CATEGORY  H2A2A1 */
/* ***TYPE      SINGLE PRECISION (QAWOE-S, DQAWOE-D) */
/* ***KEYWORDS  AUTOMATIC INTEGRATOR, CLENSHAW-CURTIS METHOD, */
/*             EXTRAPOLATION, GLOBALLY ADAPTIVE, */
/*             INTEGRAND WITH OSCILLATORY COS OR SIN FACTOR, QUADPACK, */
/*             QUADRATURE, SPECIAL-PURPOSE */
/* ***AUTHOR  Piessens, Robert */
/*             Applied Mathematics and Programming Division */
/*             K. U. Leuven */
/*           de Doncker, Elise */
/*             Applied Mathematics and Programming Division */
/*             K. U. Leuven */
/* ***DESCRIPTION */

/*        Computation of Oscillatory integrals */
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

/*            OMEGA  - Real */
/*                     Parameter in the integrand weight function */

/*            INTEGR - Integer */
/*                     Indicates which of the WEIGHT functions is to be */
/*                     used */
/*                     INTEGR = 1      W(X) = COS(OMEGA*X) */
/*                     INTEGR = 2      W(X) = SIN(OMEGA*X) */
/*                     If INTEGR.NE.1 and INTEGR.NE.2, the routine */
/*                     will end with IER = 6. */

/*            EPSABS - Real */
/*                     Absolute accuracy requested */
/*            EPSREL - Real */
/*                     Relative accuracy requested */
/*                     If  EPSABS.LE.0 */
/*                     and EPSREL.LT.MAX(50*REL.MACH.ACC.,0.5D-28), */
/*                     the routine will end with IER = 6. */

/*            LIMIT  - Integer */
/*                     Gives an upper bound on the number of subdivisions */
/*                     in the partition of (A,B), LIMIT.GE.1. */

/*            ICALL  - Integer */
/*                     If QAWOE is to be used only once, ICALL must */
/*                     be set to 1.  Assume that during this call, the */
/*                     Chebyshev moments (for CLENSHAW-CURTIS integration */
/*                     of degree 24) have been computed for intervals of */
/*                     lengths (ABS(B-A))*2**(-L), L=0,1,2,...MOMCOM-1. */
/*                     If ICALL.GT.1 this means that QAWOE has been */
/*                     called twice or more on intervals of the same */
/*                     length ABS(B-A). The Chebyshev moments already */
/*                     computed are then re-used in subsequent calls. */
/*                     If ICALL.LT.1, the routine will end with IER = 6. */

/*            MAXP1  - Integer */
/*                     Gives an upper bound on the number of Chebyshev */
/*                     moments which can be stored, i.e. for the */
/*                     intervals of lengths ABS(B-A)*2**(-L), */
/*                     L=0,1, ..., MAXP1-2, MAXP1.GE.1. */
/*                     If MAXP1.LT.1, the routine will end with IER = 6. */

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
/*                             routine. It is assumed that the */
/*                             requested accuracy has been achieved. */
/*                   - IER.GT.0 Abnormal termination of the routine. */
/*                             The estimates for integral and error are */
/*                             less reliable. It is assumed that the */
/*                             requested accuracy has not been achieved. */
/*            ERROR MESSAGES */
/*                     IER = 1 Maximum number of subdivisions allowed */
/*                             has been achieved. One can allow more */
/*                             subdivisions by increasing the value of */
/*                             LIMIT (and taking according dimension */
/*                             adjustments into account). However, if */
/*                             this yields no improvement it is advised */
/*                             to analyze the integrand, in order to */
/*                             determine the integration difficulties. */
/*                             If the position of a local difficulty can */
/*                             be determined (e.g. SINGULARITY, */
/*                             DISCONTINUITY within the interval) one */
/*                             will probably gain from splitting up the */
/*                             interval at this point and calling the */
/*                             integrator on the subranges. If possible, */
/*                             an appropriate special-purpose integrator */
/*                             should be used which is designed for */
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
/*                             It is presumed that the requested */
/*                             tolerance cannot be achieved due to */
/*                             roundoff in the extrapolation table, */
/*                             and that the returned result is the */
/*                             best which can be obtained. */
/*                         = 5 The integral is probably divergent, or */
/*                             slowly convergent. It must be noted that */
/*                             divergence can occur with any other value */
/*                             of IER.GT.0. */
/*                         = 6 The input is invalid, because */
/*                             (EPSABS.LE.0 and */
/*                              EPSREL.LT.MAX(50*REL.MACH.ACC.,0.5D-28)) */
/*                             or (INTEGR.NE.1 and INTEGR.NE.2) or */
/*                             ICALL.LT.1 or MAXP1.LT.1. */
/*                             RESULT, ABSERR, NEVAL, LAST, RLIST(1), */
/*                             ELIST(1), IORD(1) and NNLOG(1) are set */
/*                             to ZERO. ALIST(1) and BLIST(1) are set */
/*                             to A and B respectively. */

/*            LAST  -  Integer */
/*                     On return, LAST equals the number of */
/*                     subintervals produces in the subdivision */
/*                     process, which determines the number of */
/*                     significant elements actually in the */
/*                     WORK ARRAYS. */
/*            ALIST  - Real */
/*                     Vector of dimension at least LIMIT, the first */
/*                      LAST  elements of which are the left */
/*                     end points of the subintervals in the partition */
/*                     of the given integration range (A,B) */

/*            BLIST  - Real */
/*                     Vector of dimension at least LIMIT, the first */
/*                      LAST  elements of which are the right */
/*                     end points of the subintervals in the partition */
/*                     of the given integration range (A,B) */

/*            RLIST  - Real */
/*                     Vector of dimension at least LIMIT, the first */
/*                      LAST  elements of which are the integral */
/*                     approximations on the subintervals */

/*            ELIST  - Real */
/*                     Vector of dimension at least LIMIT, the first */
/*                      LAST  elements of which are the moduli of the */
/*                     absolute error estimates on the subintervals */

/*            IORD   - Integer */
/*                     Vector of dimension at least LIMIT, the first K */
/*                     elements of which are pointers to the error */
/*                     estimates over the subintervals, */
/*                     such that ELIST(IORD(1)), ..., */
/*                     ELIST(IORD(K)) form a decreasing sequence, with */
/*                     K = LAST if LAST.LE.(LIMIT/2+2), and */
/*                     K = LIMIT+1-LAST otherwise. */

/*            NNLOG  - Integer */
/*                     Vector of dimension at least LIMIT, containing the */
/*                     subdivision levels of the subintervals, i.e. */
/*                     IWORK(I) = L means that the subinterval */
/*                     numbered I is of length ABS(B-A)*2**(1-L) */

/*         ON ENTRY AND RETURN */
/*            MOMCOM - Integer */
/*                     Indicating that the Chebyshev moments */
/*                     have been computed for intervals of lengths */
/*                     (ABS(B-A))*2**(-L), L=0,1,2, ..., MOMCOM-1, */
/*                     MOMCOM.LT.MAXP1 */

/*            CHEBMO - Real */
/*                     Array of dimension (MAXP1,25) containing the */
/*                     Chebyshev moments */

/* ***REFERENCES  (NONE) */
/* ***ROUTINES CALLED  QC25F, QELG, QPSRT, R1MACH */
/* ***REVISION HISTORY  (YYMMDD) */
/*   800101  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890831  Modified array declarations.  (WRB) */
/*   890831  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/* ***END PROLOGUE  QAWOE */




/*            THE DIMENSION OF RLIST2 IS DETERMINED BY  THE VALUE OF */
/*            LIMEXP IN SUBROUTINE QELG (RLIST2 SHOULD BE OF */
/*            DIMENSION (LIMEXP+2) AT LEAST). */

/*            LIST OF MAJOR VARIABLES */
/*            ----------------------- */

/*           ALIST     - LIST OF LEFT END POINTS OF ALL SUBINTERVALS */
/*                       CONSIDERED UP TO NOW */
/*           BLIST     - LIST OF RIGHT END POINTS OF ALL SUBINTERVALS */
/*                       CONSIDERED UP TO NOW */
/*           RLIST(I)  - APPROXIMATION TO THE INTEGRAL OVER */
/*                       (ALIST(I),BLIST(I)) */
/*           RLIST2    - ARRAY OF DIMENSION AT LEAST LIMEXP+2 */
/*                       CONTAINING THE PART OF THE EPSILON TABLE */
/*                       WHICH IS STILL NEEDED FOR FURTHER COMPUTATIONS */
/*           ELIST(I)  - ERROR ESTIMATE APPLYING TO RLIST(I) */
/*           MAXERR    - POINTER TO THE INTERVAL WITH LARGEST */
/*                       ERROR ESTIMATE */
/*           ERRMAX    - ELIST(MAXERR) */
/*           ERLAST    - ERROR ON THE INTERVAL CURRENTLY SUBDIVIDED */
/*           AREA      - SUM OF THE INTEGRALS OVER THE SUBINTERVALS */
/*           ERRSUM    - SUM OF THE ERRORS OVER THE SUBINTERVALS */
/*           ERRBND    - REQUESTED ACCURACY MAX(EPSABS,EPSREL* */
/*                       ABS(RESULT)) */
/*           *****1    - VARIABLE FOR THE LEFT SUBINTERVAL */
/*           *****2    - VARIABLE FOR THE RIGHT SUBINTERVAL */
/*           LAST      - INDEX FOR SUBDIVISION */
/*           NRES      - NUMBER OF CALLS TO THE EXTRAPOLATION ROUTINE */
/*           NUMRL2    - NUMBER OF ELEMENTS IN RLIST2. IF AN APPROPRIATE */
/*                       APPROXIMATION TO THE COMPOUNDED INTEGRAL HAS */
/*                       BEEN OBTAINED IT IS PUT IN RLIST2(NUMRL2) AFTER */
/*                       NUMRL2 HAS BEEN INCREASED BY ONE */
/*           SMALL     - LENGTH OF THE SMALLEST INTERVAL CONSIDERED */
/*                       UP TO NOW, MULTIPLIED BY 1.5 */
/*           ERLARG    - SUM OF THE ERRORS OVER THE INTERVALS LARGER */
/*                       THAN THE SMALLEST INTERVAL CONSIDERED UP TO NOW */
/*           EXTRAP    - LOGICAL VARIABLE DENOTING THAT THE ROUTINE IS */
/*                       ATTEMPTING TO PERFORM EXTRAPOLATION, I.E. BEFORE */
/*                       SUBDIVIDING THE SMALLEST INTERVAL WE TRY TO */
/*                       DECREASE THE VALUE OF ERLARG */
/*           NOEXT     - LOGICAL VARIABLE DENOTING THAT EXTRAPOLATION */
/*                       IS NO LONGER ALLOWED (TRUE VALUE) */

/*            MACHINE DEPENDENT CONSTANTS */
/*            --------------------------- */

/*           EPMACH IS THE LARGEST RELATIVE SPACING. */
/*           UFLOW IS THE SMALLEST POSITIVE MAGNITUDE. */
/*           OFLOW IS THE LARGEST POSITIVE MAGNITUDE. */

/* ***FIRST EXECUTABLE STATEMENT  QAWOE */
    /* Parameter adjustments */
    chebmo_dim1 = *maxp1;
    chebmo_offset = 1 + chebmo_dim1;
    chebmo -= chebmo_offset;
    --alist__;
    --blist;
    --rlist;
    --elist;
    --iord;
    --nnlog;

    /* Function Body */
    epmach = r1mach_(&c__4);

/*         TEST ON VALIDITY OF PARAMETERS */
/*         ------------------------------ */

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
    nnlog[1] = 0;
/* Computing MAX */
    r__1 = epmach * 50.f;
    if (*integr != 1 && *integr != 2 || *epsabs <= 0.f && *epsrel < dmax(r__1,
	    5e-15f) || *icall < 1 || *maxp1 < 1) {
	*ier = 6;
    }
    if (*ier == 6) {
	goto L999;
    }

/*           FIRST APPROXIMATION TO THE INTEGRAL */
/*           ----------------------------------- */

    domega = dabs(*omega);
    nrmom = 0;
    if (*icall > 1) {
	goto L5;
    }
    *momcom = 0;
L5:
    qc25f_((E_fp)f, a, b, &domega, integr, &nrmom, maxp1, &c__0, result, 
	    abserr, neval, &defabs, &resabs, momcom, &chebmo[chebmo_offset]);

/*           TEST ON ACCURACY. */

    dres = dabs(*result);
/* Computing MAX */
    r__1 = *epsabs, r__2 = *epsrel * dres;
    errbnd = dmax(r__1,r__2);
    rlist[1] = *result;
    elist[1] = *abserr;
    iord[1] = 1;
    if (*abserr <= epmach * 100.f * defabs && *abserr > errbnd) {
	*ier = 2;
    }
    if (*limit == 1) {
	*ier = 1;
    }
    if (*ier != 0 || *abserr <= errbnd) {
	goto L200;
    }

/*           INITIALIZATIONS */
/*           --------------- */

    uflow = r1mach_(&c__1);
    oflow = r1mach_(&c__2);
    errmax = *abserr;
    maxerr = 1;
    area = *result;
    errsum = *abserr;
    *abserr = oflow;
    nrmax = 1;
    extrap = FALSE_;
    noext = FALSE_;
    ierro = 0;
    iroff1 = 0;
    iroff2 = 0;
    iroff3 = 0;
    ktmin = 0;
    small = (r__1 = *b - *a, dabs(r__1)) * .75f;
    nres = 0;
    numrl2 = 0;
    extall = FALSE_;
    if ((r__1 = *b - *a, dabs(r__1)) * .5f * domega > 2.f) {
	goto L10;
    }
    numrl2 = 1;
    extall = TRUE_;
    rlist2[0] = *result;
L10:
    if ((r__1 = *b - *a, dabs(r__1)) * .25f * domega <= 2.f) {
	extall = TRUE_;
    }
    ksgn = -1;
    if (dres >= (1.f - epmach * 50.f) * defabs) {
	ksgn = 1;
    }

/*           MAIN DO-LOOP */
/*           ------------ */

    i__1 = *limit;
    for (*last = 2; *last <= i__1; ++(*last)) {

/*           BISECT THE SUBINTERVAL WITH THE NRMAX-TH LARGEST */
/*           ERROR ESTIMATE. */

	nrmom = nnlog[maxerr] + 1;
	a1 = alist__[maxerr];
	b1 = (alist__[maxerr] + blist[maxerr]) * .5f;
	a2 = b1;
	b2 = blist[maxerr];
	erlast = errmax;
	qc25f_((E_fp)f, &a1, &b1, &domega, integr, &nrmom, maxp1, &c__0, &
		area1, &error1, &nev, &resabs, &defab1, momcom, &chebmo[
		chebmo_offset]);
	*neval += nev;
	qc25f_((E_fp)f, &a2, &b2, &domega, integr, &nrmom, maxp1, &c__1, &
		area2, &error2, &nev, &resabs, &defab2, momcom, &chebmo[
		chebmo_offset]);
	*neval += nev;

/*           IMPROVE PREVIOUS APPROXIMATIONS TO INTEGRAL */
/*           AND ERROR AND TEST FOR ACCURACY. */

	area12 = area1 + area2;
	erro12 = error1 + error2;
	errsum = errsum + erro12 - errmax;
	area = area + area12 - rlist[maxerr];
	if (defab1 == error1 || defab2 == error2) {
	    goto L25;
	}
	if ((r__1 = rlist[maxerr] - area12, dabs(r__1)) > dabs(area12) * 
		1e-5f || erro12 < errmax * .99f) {
	    goto L20;
	}
	if (extrap) {
	    ++iroff2;
	}
	if (! extrap) {
	    ++iroff1;
	}
L20:
	if (*last > 10 && erro12 > errmax) {
	    ++iroff3;
	}
L25:
	rlist[maxerr] = area1;
	rlist[*last] = area2;
	nnlog[maxerr] = nrmom;
	nnlog[*last] = nrmom;
/* Computing MAX */
	r__1 = *epsabs, r__2 = *epsrel * dabs(area);
	errbnd = dmax(r__1,r__2);

/*           TEST FOR ROUNDOFF ERROR AND EVENTUALLY */
/*           SET ERROR FLAG */

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
/*           AT A POINT OF THE INTEGRATION RANGE. */

/* Computing MAX */
	r__1 = dabs(a1), r__2 = dabs(b2);
	if (dmax(r__1,r__2) <= (epmach * 100.f + 1.f) * (dabs(a2) + uflow * 
		1e3f)) {
	    *ier = 4;
	}

/*           APPEND THE NEWLY-CREATED INTERVALS TO THE LIST. */

	if (error2 > error1) {
	    goto L30;
	}
	alist__[*last] = a2;
	blist[maxerr] = b1;
	blist[*last] = b2;
	elist[maxerr] = error1;
	elist[*last] = error2;
	goto L40;
L30:
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

L40:
	qpsrt_(limit, last, &maxerr, &errmax, &elist[1], &iord[1], &nrmax);
/* ***JUMP OUT OF DO-LOOP */
	if (errsum <= errbnd) {
	    goto L170;
	}
	if (*ier != 0) {
	    goto L150;
	}
	if (*last == 2 && extall) {
	    goto L120;
	}
	if (noext) {
	    goto L140;
	}
	if (! extall) {
	    goto L50;
	}
	erlarg -= erlast;
	if ((r__1 = b1 - a1, dabs(r__1)) > small) {
	    erlarg += erro12;
	}
	if (extrap) {
	    goto L70;
	}

/*           TEST WHETHER THE INTERVAL TO BE BISECTED NEXT IS THE */
/*           SMALLEST INTERVAL. */

L50:
	width = (r__1 = blist[maxerr] - alist__[maxerr], dabs(r__1));
	if (width > small) {
	    goto L140;
	}
	if (extall) {
	    goto L60;
	}

/*           TEST WHETHER WE CAN START WITH THE EXTRAPOLATION */
/*           PROCEDURE (WE DO THIS IF WE INTEGRATE OVER THE */
/*           NEXT INTERVAL WITH USE OF A GAUSS-KRONROD RULE - SEE */
/*           SUBROUTINE QC25F). */

	small *= .5f;
	if (width * .25f * domega > 2.f) {
	    goto L140;
	}
	extall = TRUE_;
	goto L130;
L60:
	extrap = TRUE_;
	nrmax = 2;
L70:
	if (ierro == 3 || erlarg <= ertest) {
	    goto L90;
	}

/*           THE SMALLEST INTERVAL HAS THE LARGEST ERROR. */
/*           BEFORE BISECTING DECREASE THE SUM OF THE ERRORS */
/*           OVER THE LARGER INTERVALS (ERLARG) AND PERFORM */
/*           EXTRAPOLATION. */

	jupbnd = *last;
	if (*last > *limit / 2 + 2) {
	    jupbnd = *limit + 3 - *last;
	}
	id = nrmax;
	i__2 = jupbnd;
	for (k = id; k <= i__2; ++k) {
	    maxerr = iord[nrmax];
	    errmax = elist[maxerr];
	    if ((r__1 = blist[maxerr] - alist__[maxerr], dabs(r__1)) > small) 
		    {
		goto L140;
	    }
	    ++nrmax;
/* L80: */
	}

/*           PERFORM EXTRAPOLATION. */

L90:
	++numrl2;
	rlist2[numrl2 - 1] = area;
	if (numrl2 < 3) {
	    goto L110;
	}
	qelg_(&numrl2, rlist2, &reseps, &abseps, res3la, &nres);
	++ktmin;
	if (ktmin > 5 && *abserr < errsum * .001f) {
	    *ier = 5;
	}
	if (abseps >= *abserr) {
	    goto L100;
	}
	ktmin = 0;
	*abserr = abseps;
	*result = reseps;
	correc = erlarg;
/* Computing MAX */
	r__1 = *epsabs, r__2 = *epsrel * dabs(reseps);
	ertest = dmax(r__1,r__2);
/* ***JUMP OUT OF DO-LOOP */
	if (*abserr <= ertest) {
	    goto L150;
	}

/*           PREPARE BISECTION OF THE SMALLEST INTERVAL. */

L100:
	if (numrl2 == 1) {
	    noext = TRUE_;
	}
	if (*ier == 5) {
	    goto L150;
	}
L110:
	maxerr = iord[1];
	errmax = elist[maxerr];
	nrmax = 1;
	extrap = FALSE_;
	small *= .5f;
	erlarg = errsum;
	goto L140;
L120:
	small *= .5f;
	++numrl2;
	rlist2[numrl2 - 1] = area;
L130:
	ertest = errbnd;
	erlarg = errsum;
L140:
	;
    }

/*           SET THE FINAL RESULT. */
/*           --------------------- */

L150:
    if (*abserr == oflow || nres == 0) {
	goto L170;
    }
    if (*ier + ierro == 0) {
	goto L165;
    }
    if (ierro == 3) {
	*abserr += correc;
    }
    if (*ier == 0) {
	*ier = 3;
    }
    if (*result != 0.f && area != 0.f) {
	goto L160;
    }
    if (*abserr > errsum) {
	goto L170;
    }
    if (area == 0.f) {
	goto L190;
    }
    goto L165;
L160:
    if (*abserr / dabs(*result) > errsum / dabs(area)) {
	goto L170;
    }

/*           TEST ON DIVERGENCE. */

L165:
/* Computing MAX */
    r__1 = dabs(*result), r__2 = dabs(area);
    if (ksgn == -1 && dmax(r__1,r__2) <= defabs * .01f) {
	goto L190;
    }
    if (.01f > *result / area || *result / area > 100.f || errsum >= dabs(
	    area)) {
	*ier = 6;
    }
    goto L190;

/*           COMPUTE GLOBAL INTEGRAL SUM. */

L170:
    *result = 0.f;
    i__1 = *last;
    for (k = 1; k <= i__1; ++k) {
	*result += rlist[k];
/* L180: */
    }
    *abserr = errsum;
L190:
    if (*ier > 2) {
	--(*ier);
    }
L200:
    if (*integr == 2 && *omega < 0.f) {
	*result = -(*result);
    }
L999:
    return 0;
} /* qawoe_ */

