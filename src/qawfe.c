/* qawfe.f -- translated by f2c (version 12.02.01).
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
static real c_b5 = 0.f;

/* DECK QAWFE */
/* Subroutine */ int qawfe_(U_fp f, real *a, real *omega, integer *integr, 
	real *epsabs, integer *limlst, integer *limit, integer *maxp1, real *
	result, real *abserr, integer *neval, integer *ier, real *rslst, real 
	*erlst, integer *ierlst, integer *lst, real *alist__, real *blist, 
	real *rlist, real *elist, integer *iord, integer *nnlog, real *chebmo)
{
    /* Initialized data */

    static real p = .9f;
    static real pi = 3.1415926535897932f;

    /* System generated locals */
    integer chebmo_dim1, chebmo_offset, i__1;
    real r__1, r__2;

    /* Local variables */
    static integer l;
    static real c1, c2, p1, dl, ep;
    static integer ll;
    static real drl, eps;
    static integer nev;
    static real fact, epsa;
    extern /* Subroutine */ int qelg_(integer *, real *, real *, real *, real 
	    *, integer *);
    static integer last, nres;
    static real psum[52];
    extern /* Subroutine */ int qagie_(U_fp, real *, integer *, real *, real *
	    , integer *, real *, real *, integer *, integer *, real *, real *,
	     real *, real *, integer *, integer *);
    static real cycle;
    extern /* Subroutine */ int qawoe_(U_fp, real *, real *, real *, integer *
	    , real *, real *, integer *, integer *, integer *, real *, real *,
	     integer *, integer *, integer *, real *, real *, real *, real *, 
	    integer *, integer *, integer *, real *);
    static integer ktmin;
    static real uflow;
    extern doublereal r1mach_(integer *);
    static real res3la[3];
    static integer numrl2;
    static real abseps, correc;
    static integer momcom;
    static real reseps, errsum;

/* ***BEGIN PROLOGUE  QAWFE */
/* ***PURPOSE  The routine calculates an approximation result to a */
/*            given Fourier integral */
/*            I = Integral of F(X)*W(X) over (A,INFINITY) */
/*             where W(X) = COS(OMEGA*X) or W(X) = SIN(OMEGA*X), */
/*            hopefully satisfying following claim for accuracy */
/*            ABS(I-RESULT).LE.EPSABS. */
/* ***LIBRARY   SLATEC (QUADPACK) */
/* ***CATEGORY  H2A3A1 */
/* ***TYPE      SINGLE PRECISION (QAWFE-S, DQAWFE-D) */
/* ***KEYWORDS  AUTOMATIC INTEGRATOR, CONVERGENCE ACCELERATION, */
/*             FOURIER INTEGRALS, INTEGRATION BETWEEN ZEROS, QUADPACK, */
/*             QUADRATURE, SPECIAL-PURPOSE INTEGRAL */
/* ***AUTHOR  Piessens, Robert */
/*             Applied Mathematics and Programming Division */
/*             K. U. Leuven */
/*           de Doncker, Elise */
/*             Applied Mathematics and Programming Division */
/*             K. U. Leuven */
/* ***DESCRIPTION */

/*        Computation of Fourier integrals */
/*        Standard fortran subroutine */
/*        Real version */

/*        PARAMETERS */
/*         ON ENTRY */
/*            F      - Real */
/*                     Function subprogram defining the integrand */
/*                     Function F(X). The actual name for F needs to */
/*                     be declared E X T E R N A L in the driver program. */

/*            A      - Real */
/*                     Lower limit of integration */

/*            OMEGA  - Real */
/*                     Parameter in the WEIGHT function */

/*            INTEGR - Integer */
/*                     Indicates which WEIGHT function is used */
/*                     INTEGR = 1      W(X) = COS(OMEGA*X) */
/*                     INTEGR = 2      W(X) = SIN(OMEGA*X) */
/*                     If INTEGR.NE.1.AND.INTEGR.NE.2, the routine will */
/*                     end with IER = 6. */

/*            EPSABS - Real */
/*                     absolute accuracy requested, EPSABS.GT.0 */
/*                     If EPSABS.LE.0, the routine will end with IER = 6. */

/*            LIMLST - Integer */
/*                     LIMLST gives an upper bound on the number of */
/*                     cycles, LIMLST.GE.1. */
/*                     If LIMLST.LT.3, the routine will end with IER = 6. */

/*            LIMIT  - Integer */
/*                     Gives an upper bound on the number of subintervals */
/*                     allowed in the partition of each cycle, LIMIT.GE.1 */
/*                     each cycle, LIMIT.GE.1. */

/*            MAXP1  - Integer */
/*                     Gives an upper bound on the number of */
/*                     Chebyshev moments which can be stored, I.E. */
/*                     for the intervals of lengths ABS(B-A)*2**(-L), */
/*                     L=0,1, ..., MAXP1-2, MAXP1.GE.1 */

/*         ON RETURN */
/*            RESULT - Real */
/*                     Approximation to the integral X */

/*            ABSERR - Real */
/*                     Estimate of the modulus of the absolute error, */
/*                     which should equal or exceed ABS(I-RESULT) */

/*            NEVAL  - Integer */
/*                     Number of integrand evaluations */

/*            IER    - IER = 0 Normal and reliable termination of */
/*                             the routine. It is assumed that the */
/*                             requested accuracy has been achieved. */
/*                     IER.GT.0 Abnormal termination of the routine. The */
/*                             estimates for integral and error are less */
/*                             reliable. It is assumed that the requested */
/*                             accuracy has not been achieved. */
/*            ERROR MESSAGES */
/*                    If OMEGA.NE.0 */
/*                     IER = 1 Maximum number of  cycles  allowed */
/*                             Has been achieved., i.e. of subintervals */
/*                             (A+(K-1)C,A+KC) where */
/*                             C = (2*INT(ABS(OMEGA))+1)*PI/ABS(OMEGA), */
/*                             for K = 1, 2, ..., LST. */
/*                             One can allow more cycles by increasing */
/*                             the value of LIMLST (and taking the */
/*                             according dimension adjustments into */
/*                             account). */
/*                             Examine the array IWORK which contains */
/*                             the error flags on the cycles, in order to */
/*                             look for eventual local integration */
/*                             difficulties. If the position of a local */
/*                             difficulty can be determined (e.g. */
/*                             SINGULARITY, DISCONTINUITY within the */
/*                             interval) one will probably gain from */
/*                             splitting up the interval at this point */
/*                             and calling appropriate integrators on */
/*                             the subranges. */
/*                         = 4 The extrapolation table constructed for */
/*                             convergence acceleration of the series */
/*                             formed by the integral contributions over */
/*                             the cycles, does not converge to within */
/*                             the requested accuracy. As in the case of */
/*                             IER = 1, it is advised to examine the */
/*                             array IWORK which contains the error */
/*                             flags on the cycles. */
/*                         = 6 The input is invalid because */
/*                             (INTEGR.NE.1 AND INTEGR.NE.2) or */
/*                              EPSABS.LE.0 or LIMLST.LT.3. */
/*                              RESULT, ABSERR, NEVAL, LST are set */
/*                              to zero. */
/*                         = 7 Bad integrand behaviour occurs within one */
/*                             or more of the cycles. Location and type */
/*                             of the difficulty involved can be */
/*                             determined from the vector IERLST. Here */
/*                             LST is the number of cycles actually */
/*                             needed (see below). */
/*                             IERLST(K) = 1 The maximum number of */
/*                                           subdivisions (= LIMIT) has */
/*                                           been achieved on the K th */
/*                                           cycle. */
/*                                       = 2 Occurrence of roundoff error */
/*                                           is detected and prevents the */
/*                                           tolerance imposed on the */
/*                                           K th cycle, from being */
/*                                           achieved. */
/*                                       = 3 Extremely bad integrand */
/*                                           behaviour occurs at some */
/*                                           points of the K th cycle. */
/*                                       = 4 The integration procedure */
/*                                           over the K th cycle does */
/*                                           not converge (to within the */
/*                                           required accuracy) due to */
/*                                           roundoff in the */
/*                                           extrapolation procedure */
/*                                           invoked on this cycle. It */
/*                                           is assumed that the result */
/*                                           on this interval is the */
/*                                           best which can be obtained. */
/*                                       = 5 The integral over the K th */
/*                                           cycle is probably divergent */
/*                                           or slowly convergent. It */
/*                                           must be noted that */
/*                                           divergence can occur with */
/*                                           any other value of */
/*                                           IERLST(K). */
/*                    If OMEGA = 0 and INTEGR = 1, */
/*                    The integral is calculated by means of DQAGIE */
/*                    and IER = IERLST(1) (with meaning as described */
/*                    for IERLST(K), K = 1). */

/*            RSLST  - Real */
/*                     Vector of dimension at least LIMLST */
/*                     RSLST(K) contains the integral contribution */
/*                     over the interval (A+(K-1)C,A+KC) where */
/*                     C = (2*INT(ABS(OMEGA))+1)*PI/ABS(OMEGA), */
/*                     K = 1, 2, ..., LST. */
/*                     Note that, if OMEGA = 0, RSLST(1) contains */
/*                     the value of the integral over (A,INFINITY). */

/*            ERLST  - Real */
/*                     Vector of dimension at least LIMLST */
/*                     ERLST(K) contains the error estimate corresponding */
/*                     with RSLST(K). */

/*            IERLST - Integer */
/*                     Vector of dimension at least LIMLST */
/*                     IERLST(K) contains the error flag corresponding */
/*                     with RSLST(K). For the meaning of the local error */
/*                     flags see description of output parameter IER. */

/*            LST    - Integer */
/*                     Number of subintervals needed for the integration */
/*                     If OMEGA = 0 then LST is set to 1. */

/*            ALIST, BLIST, RLIST, ELIST - Real */
/*                     vector of dimension at least LIMIT, */

/*            IORD, NNLOG - Integer */
/*                     Vector of dimension at least LIMIT, providing */
/*                     space for the quantities needed in the subdivision */
/*                     process of each cycle */

/*            CHEBMO - Real */
/*                     Array of dimension at least (MAXP1,25), providing */
/*                     space for the Chebyshev moments needed within the */
/*                     cycles */

/* ***REFERENCES  (NONE) */
/* ***ROUTINES CALLED  QAGIE, QAWOE, QELG, R1MACH */
/* ***REVISION HISTORY  (YYMMDD) */
/*   800101  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890831  Modified array declarations.  (WRB) */
/*   891009  Removed unreferenced variable.  (WRB) */
/*   891009  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/* ***END PROLOGUE  QAWFE */





/*            THE DIMENSION OF  PSUM  IS DETERMINED BY THE VALUE OF */
/*            LIMEXP IN SUBROUTINE QELG (PSUM MUST BE */
/*            OF DIMENSION (LIMEXP+2) AT LEAST). */

/*           LIST OF MAJOR VARIABLES */
/*           ----------------------- */

/*           C1, C2    - END POINTS OF SUBINTERVAL (OF LENGTH */
/*                       CYCLE) */
/*           CYCLE     - (2*INT(ABS(OMEGA))+1)*PI/ABS(OMEGA) */
/*           PSUM      - VECTOR OF DIMENSION AT LEAST (LIMEXP+2) */
/*                       (SEE ROUTINE QELG) */
/*                       PSUM CONTAINS THE PART OF THE EPSILON */
/*                       TABLE WHICH IS STILL NEEDED FOR FURTHER */
/*                       COMPUTATIONS. */
/*                       EACH ELEMENT OF PSUM IS A PARTIAL SUM OF */
/*                       THE SERIES WHICH SHOULD SUM TO THE VALUE OF */
/*                       THE INTEGRAL. */
/*           ERRSUM    - SUM OF ERROR ESTIMATES OVER THE */
/*                       SUBINTERVALS, CALCULATED CUMULATIVELY */
/*           EPSA      - ABSOLUTE TOLERANCE REQUESTED OVER CURRENT */
/*                       SUBINTERVAL */
/*           CHEBMO    - ARRAY CONTAINING THE MODIFIED CHEBYSHEV */
/*                       MOMENTS (SEE ALSO ROUTINE QC25F) */

    /* Parameter adjustments */
    chebmo_dim1 = *maxp1;
    chebmo_offset = 1 + chebmo_dim1;
    chebmo -= chebmo_offset;
    --rslst;
    --erlst;
    --ierlst;
    --alist__;
    --blist;
    --rlist;
    --elist;
    --iord;
    --nnlog;

    /* Function Body */

/*           TEST ON VALIDITY OF PARAMETERS */
/*           ------------------------------ */

/* ***FIRST EXECUTABLE STATEMENT  QAWFE */
    *result = 0.f;
    *abserr = 0.f;
    *neval = 0;
    *lst = 0;
    *ier = 0;
    if (*integr != 1 && *integr != 2 || *epsabs <= 0.f || *limlst < 3) {
	*ier = 6;
    }
    if (*ier == 6) {
	goto L999;
    }
    if (*omega != 0.f) {
	goto L10;
    }

/*           INTEGRATION BY QAGIE IF OMEGA IS ZERO */
/*           -------------------------------------- */

    if (*integr == 1) {
	qagie_((U_fp)f, a, &c__1, epsabs, &c_b5, limit, result, abserr, neval,
		 ier, &alist__[1], &blist[1], &rlist[1], &elist[1], &iord[1], 
		&last);
    }
    rslst[1] = *result;
    erlst[1] = *abserr;
    ierlst[1] = *ier;
    *lst = 1;
    goto L999;

/*           INITIALIZATIONS */
/*           --------------- */

L10:
    l = dabs(*omega);
    dl = (real) ((l << 1) + 1);
    cycle = dl * pi / dabs(*omega);
    *ier = 0;
    ktmin = 0;
    *neval = 0;
    numrl2 = 0;
    nres = 0;
    c1 = *a;
    c2 = cycle + *a;
    p1 = 1.f - p;
    eps = *epsabs;
    uflow = r1mach_(&c__1);
    if (*epsabs > uflow / p1) {
	eps = *epsabs * p1;
    }
    ep = eps;
    fact = 1.f;
    correc = 0.f;
    *abserr = 0.f;
    errsum = 0.f;

/*           MAIN DO-LOOP */
/*           ------------ */

    i__1 = *limlst;
    for (*lst = 1; *lst <= i__1; ++(*lst)) {

/*           INTEGRATE OVER CURRENT SUBINTERVAL. */

	epsa = eps * fact;
	qawoe_((U_fp)f, &c1, &c2, omega, integr, &epsa, &c_b5, limit, lst, 
		maxp1, &rslst[*lst], &erlst[*lst], &nev, &ierlst[*lst], &last,
		 &alist__[1], &blist[1], &rlist[1], &elist[1], &iord[1], &
		nnlog[1], &momcom, &chebmo[chebmo_offset]);
	*neval += nev;
	fact *= p;
	errsum += erlst[*lst];
	drl = (r__1 = rslst[*lst], dabs(r__1)) * 50.f;

/*           TEST ON ACCURACY WITH PARTIAL SUM */

	if (errsum + drl <= *epsabs && *lst >= 6) {
	    goto L80;
	}
/* Computing MAX */
	r__1 = correc, r__2 = erlst[*lst];
	correc = dmax(r__1,r__2);
	if (ierlst[*lst] != 0) {
/* Computing MAX */
	    r__1 = ep, r__2 = correc * p1;
	    eps = dmax(r__1,r__2);
	}
	if (ierlst[*lst] != 0) {
	    *ier = 7;
	}
	if (*ier == 7 && errsum + drl <= correc * 10.f && *lst > 5) {
	    goto L80;
	}
	++numrl2;
	if (*lst > 1) {
	    goto L20;
	}
	psum[0] = rslst[1];
	goto L40;
L20:
	psum[numrl2 - 1] = psum[ll - 1] + rslst[*lst];
	if (*lst == 2) {
	    goto L40;
	}

/*           TEST ON MAXIMUM NUMBER OF SUBINTERVALS */

	if (*lst == *limlst) {
	    *ier = 1;
	}

/*           PERFORM NEW EXTRAPOLATION */

	qelg_(&numrl2, psum, &reseps, &abseps, res3la, &nres);

/*           TEST WHETHER EXTRAPOLATED RESULT IS INFLUENCED BY */
/*           ROUNDOFF */

	++ktmin;
	if (ktmin >= 15 && *abserr <= (errsum + drl) * .001f) {
	    *ier = 4;
	}
	if (abseps > *abserr && *lst != 3) {
	    goto L30;
	}
	*abserr = abseps;
	*result = reseps;
	ktmin = 0;

/*           IF IER IS NOT 0, CHECK WHETHER DIRECT RESULT (PARTIAL */
/*           SUM) OR EXTRAPOLATED RESULT YIELDS THE BEST INTEGRAL */
/*           APPROXIMATION */

	if (*abserr + correc * 10.f <= *epsabs || *abserr <= *epsabs && 
		correc * 10.f >= *epsabs) {
	    goto L60;
	}
L30:
	if (*ier != 0 && *ier != 7) {
	    goto L60;
	}
L40:
	ll = numrl2;
	c1 = c2;
	c2 += cycle;
/* L50: */
    }

/*         SET FINAL RESULT AND ERROR ESTIMATE */
/*         ----------------------------------- */

L60:
    *abserr += correc * 10.f;
    if (*ier == 0) {
	goto L999;
    }
    if (*result != 0.f && psum[numrl2 - 1] != 0.f) {
	goto L70;
    }
    if (*abserr > errsum) {
	goto L80;
    }
    if (psum[numrl2 - 1] == 0.f) {
	goto L999;
    }
L70:
    if (*abserr / dabs(*result) > (errsum + drl) / (r__1 = psum[numrl2 - 1], 
	    dabs(r__1))) {
	goto L80;
    }
    if (*ier >= 1 && *ier != 7) {
	*abserr += drl;
    }
    goto L999;
L80:
    *result = psum[numrl2 - 1];
    *abserr = errsum + drl;
L999:
    return 0;
} /* qawfe_ */

