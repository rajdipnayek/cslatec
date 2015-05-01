/* dqagp.f -- translated by f2c (version 12.02.01).
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

/* DECK DQAGP */
/* Subroutine */ int dqagp_(D_fp f, doublereal *a, doublereal *b, integer *
	npts2, doublereal *points, doublereal *epsabs, doublereal *epsrel, 
	doublereal *result, doublereal *abserr, integer *neval, integer *ier, 
	integer *leniw, integer *lenw, integer *last, integer *iwork, 
	doublereal *work)
{
    static integer l1, l2, l3, l4, lvl, limit;
    extern /* Subroutine */ int dqagpe_(D_fp, doublereal *, doublereal *, 
	    integer *, doublereal *, doublereal *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, integer *,
	     integer *, integer *, integer *), xermsg_(char *, char *, char *,
	     integer *, integer *, ftnlen, ftnlen, ftnlen);

/* ***BEGIN PROLOGUE  DQAGP */
/* ***PURPOSE  The routine calculates an approximation result to a given */
/*            definite integral I = Integral of F over (A,B), */
/*            hopefully satisfying following claim for accuracy */
/*            break points of the integration interval, where local */
/*            difficulties of the integrand may occur (e.g. */
/*            SINGULARITIES, DISCONTINUITIES), are provided by the user. */
/* ***LIBRARY   SLATEC (QUADPACK) */
/* ***CATEGORY  H2A2A1 */
/* ***TYPE      DOUBLE PRECISION (QAGP-S, DQAGP-D) */
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
/*        Double precision version */

/*        PARAMETERS */
/*         ON ENTRY */
/*            F      - Double precision */
/*                     Function subprogram defining the integrand */
/*                     Function F(X). The actual name for F needs to be */
/*                     declared E X T E R N A L in the driver program. */

/*            A      - Double precision */
/*                     Lower limit of integration */

/*            B      - Double precision */
/*                     Upper limit of integration */

/*            NPTS2  - Integer */
/*                     Number equal to two more than the number of */
/*                     user-supplied break points within the integration */
/*                     range, NPTS.GE.2. */
/*                     If NPTS2.LT.2, The routine will end with IER = 6. */

/*            POINTS - Double precision */
/*                     Vector of dimension NPTS2, the first (NPTS2-2) */
/*                     elements of which are the user provided break */
/*                     points. If these points do not constitute an */
/*                     ascending sequence there will be an automatic */
/*                     sorting. */

/*            EPSABS - Double precision */
/*                     Absolute accuracy requested */
/*            EPSREL - Double precision */
/*                     Relative accuracy requested */
/*                     If  EPSABS.LE.0 */
/*                     And EPSREL.LT.MAX(50*REL.MACH.ACC.,0.5D-28), */
/*                     The routine will end with IER = 6. */

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
/*                     IER.GT.0 Abnormal termination of the routine. */
/*                             The estimates for integral and error are */
/*                             less reliable. it is assumed that the */
/*                             requested accuracy has not been achieved. */
/*            ERROR MESSAGES */
/*                     IER = 1 Maximum number of subdivisions allowed */
/*                             has been achieved. one can allow more */
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
/*                             at some points of the integration */
/*                             interval. */
/*                         = 4 The algorithm does not converge. */
/*                             roundoff error is detected in the */
/*                             extrapolation table. */
/*                             It is presumed that the requested */
/*                             tolerance cannot be achieved, and that */
/*                             the returned RESULT is the best which */
/*                             can be obtained. */
/*                         = 5 The integral is probably divergent, or */
/*                             slowly convergent. it must be noted that */
/*                             divergence can occur with any other value */
/*                             of IER.GT.0. */
/*                         = 6 The input is invalid because */
/*                             NPTS2.LT.2 or */
/*                             break points are specified outside */
/*                             the integration range or */
/*                             (EPSABS.LE.0 and */
/*                              EPSREL.LT.MAX(50*REL.MACH.ACC.,0.5D-28)) */
/*                             RESULT, ABSERR, NEVAL, LAST are set to */
/*                             zero.  Except when LENIW or LENW or NPTS2 */
/*                             is invalid, IWORK(1), IWORK(LIMIT+1), */
/*                             WORK(LIMIT*2+1) and WORK(LIMIT*3+1) */
/*                             are set to zero. */
/*                             WORK(1) is set to A and WORK(LIMIT+1) */
/*                             to B (where LIMIT = (LENIW-NPTS2)/2). */

/*         DIMENSIONING PARAMETERS */
/*            LENIW - Integer */
/*                    Dimensioning parameter for IWORK */
/*                    LENIW determines LIMIT = (LENIW-NPTS2)/2, */
/*                    which is the maximum number of subintervals in the */
/*                    partition of the given integration interval (A,B), */
/*                    LENIW.GE.(3*NPTS2-2). */
/*                    If LENIW.LT.(3*NPTS2-2), the routine will end with */
/*                    IER = 6. */

/*            LENW  - Integer */
/*                    Dimensioning parameter for WORK */
/*                    LENW must be at least LENIW*2-NPTS2. */
/*                    If LENW.LT.LENIW*2-NPTS2, the routine will end */
/*                    with IER = 6. */

/*            LAST  - Integer */
/*                    On return, LAST equals the number of subintervals */
/*                    produced in the subdivision process, which */
/*                    determines the number of significant elements */
/*                    actually in the WORK ARRAYS. */

/*         WORK ARRAYS */
/*            IWORK - Integer */
/*                    Vector of dimension at least LENIW. on return, */
/*                    the first K elements of which contain */
/*                    pointers to the error estimates over the */
/*                    subintervals, such that WORK(LIMIT*3+IWORK(1)),..., */
/*                    WORK(LIMIT*3+IWORK(K)) form a decreasing */
/*                    sequence, with K = LAST if LAST.LE.(LIMIT/2+2), and */
/*                    K = LIMIT+1-LAST otherwise */
/*                    IWORK(LIMIT+1), ...,IWORK(LIMIT+LAST) Contain the */
/*                     subdivision levels of the subintervals, i.e. */
/*                     if (AA,BB) is a subinterval of (P1,P2) */
/*                     where P1 as well as P2 is a user-provided */
/*                     break point or integration LIMIT, then (AA,BB) has */
/*                     level L if ABS(BB-AA) = ABS(P2-P1)*2**(-L), */
/*                    IWORK(LIMIT*2+1), ..., IWORK(LIMIT*2+NPTS2) have */
/*                     no significance for the user, */
/*                    note that LIMIT = (LENIW-NPTS2)/2. */

/*            WORK  - Double precision */
/*                    Vector of dimension at least LENW */
/*                    on return */
/*                    WORK(1), ..., WORK(LAST) contain the left */
/*                     end points of the subintervals in the */
/*                     partition of (A,B), */
/*                    WORK(LIMIT+1), ..., WORK(LIMIT+LAST) contain */
/*                     the right end points, */
/*                    WORK(LIMIT*2+1), ..., WORK(LIMIT*2+LAST) contain */
/*                     the integral approximations over the subintervals, */
/*                    WORK(LIMIT*3+1), ..., WORK(LIMIT*3+LAST) */
/*                     contain the corresponding error estimates, */
/*                    WORK(LIMIT*4+1), ..., WORK(LIMIT*4+NPTS2) */
/*                     contain the integration limits and the */
/*                     break points sorted in an ascending sequence. */
/*                    note that LIMIT = (LENIW-NPTS2)/2. */

/* ***REFERENCES  (NONE) */
/* ***ROUTINES CALLED  DQAGPE, XERMSG */
/* ***REVISION HISTORY  (YYMMDD) */
/*   800101  DATE WRITTEN */
/*   890831  Modified array declarations.  (WRB) */
/*   890831  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ) */
/* ***END PROLOGUE  DQAGP */




/*         CHECK VALIDITY OF LIMIT AND LENW. */

/* ***FIRST EXECUTABLE STATEMENT  DQAGP */
    /* Parameter adjustments */
    --work;
    --iwork;
    --points;

    /* Function Body */
    *ier = 6;
    *neval = 0;
    *last = 0;
    *result = 0.;
    *abserr = 0.;
    if (*leniw < *npts2 * 3 - 2 || *lenw < (*leniw << 1) - *npts2 || *npts2 < 
	    2) {
	goto L10;
    }

/*         PREPARE CALL FOR DQAGPE. */

    limit = (*leniw - *npts2) / 2;
    l1 = limit + 1;
    l2 = limit + l1;
    l3 = limit + l2;
    l4 = limit + l3;

    dqagpe_((D_fp)f, a, b, npts2, &points[1], epsabs, epsrel, &limit, result, 
	    abserr, neval, ier, &work[1], &work[l1], &work[l2], &work[l3], &
	    work[l4], &iwork[1], &iwork[l1], &iwork[l2], last);

/*         CALL ERROR HANDLER IF NECESSARY. */

    lvl = 0;
L10:
    if (*ier == 6) {
	lvl = 1;
    }
    if (*ier != 0) {
	xermsg_("SLATEC", "DQAGP", "ABNORMAL RETURN", ier, &lvl, (ftnlen)6, (
		ftnlen)5, (ftnlen)15);
    }
    return 0;
} /* dqagp_ */

