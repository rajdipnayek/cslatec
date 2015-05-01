/* dqag.f -- translated by f2c (version 12.02.01).
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

/* DECK DQAG */
/* Subroutine */ int dqag_(D_fp f, doublereal *a, doublereal *b, doublereal *
	epsabs, doublereal *epsrel, integer *key, doublereal *result, 
	doublereal *abserr, integer *neval, integer *ier, integer *limit, 
	integer *lenw, integer *last, integer *iwork, doublereal *work)
{
    static integer l1, l2, l3, lvl;
    extern /* Subroutine */ int dqage_(D_fp, doublereal *, doublereal *, 
	    doublereal *, doublereal *, integer *, integer *, doublereal *, 
	    doublereal *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, integer *, integer *), xermsg_(char *,
	     char *, char *, integer *, integer *, ftnlen, ftnlen, ftnlen);

/* ***BEGIN PROLOGUE  DQAG */
/* ***PURPOSE  The routine calculates an approximation result to a given */
/*            definite integral I = integral of F over (A,B), */
/*            hopefully satisfying following claim for accuracy */
/*            ABS(I-RESULT)LE.MAX(EPSABS,EPSREL*ABS(I)). */
/* ***LIBRARY   SLATEC (QUADPACK) */
/* ***CATEGORY  H2A1A1 */
/* ***TYPE      DOUBLE PRECISION (QAG-S, DQAG-D) */
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
/*        Double precision version */

/*            F      - Double precision */
/*                     Function subprogram defining the integrand */
/*                     Function F(X). The actual name for F needs to be */
/*                     Declared E X T E R N A L in the driver program. */

/*            A      - Double precision */
/*                     Lower limit of integration */

/*            B      - Double precision */
/*                     Upper limit of integration */

/*            EPSABS - Double precision */
/*                     Absolute accuracy requested */
/*            EPSREL - Double precision */
/*                     Relative accuracy requested */
/*                     If  EPSABS.LE.0 */
/*                     And EPSREL.LT.MAX(50*REL.MACH.ACC.,0.5D-28), */
/*                     The routine will end with IER = 6. */

/*            KEY    - Integer */
/*                     Key for choice of local integration rule */
/*                     A GAUSS-KRONROD PAIR is used with */
/*                       7 - 15 POINTS If KEY.LT.2, */
/*                      10 - 21 POINTS If KEY = 2, */
/*                      15 - 31 POINTS If KEY = 3, */
/*                      20 - 41 POINTS If KEY = 4, */
/*                      25 - 51 POINTS If KEY = 5, */
/*                      30 - 61 POINTS If KEY.GT.5. */

/*         ON RETURN */
/*            RESULT - Double precision */
/*                     Approximation to the integral */

/*            ABSERR - Double precision */
/*                     Estimate of the modulus of the absolute error, */
/*                     Which should EQUAL or EXCEED ABS(I-RESULT) */

/*            NEVAL  - Integer */
/*                     Number of integrand evaluations */

/*            IER    - Integer */
/*                     IER = 0 Normal and reliable termination of the */
/*                             routine. It is assumed that the requested */
/*                             accuracy has been achieved. */
/*                     IER.GT.0 Abnormal termination of the routine */
/*                             The estimates for RESULT and ERROR are */
/*                             Less reliable. It is assumed that the */
/*                             requested accuracy has not been achieved. */
/*                      ERROR MESSAGES */
/*                     IER = 1 Maximum number of subdivisions allowed */
/*                             has been achieved. One can allow more */
/*                             subdivisions by increasing the value of */
/*                             LIMIT (and taking the according dimension */
/*                             adjustments into account). HOWEVER, If */
/*                             this yield no improvement it is advised */
/*                             to analyze the integrand in order to */
/*                             determine the integration difficulties. */
/*                             If the position of a local difficulty can */
/*                             be determined (I.E. SINGULARITY, */
/*                             DISCONTINUITY WITHIN THE INTERVAL) One */
/*                             will probably gain from splitting up the */
/*                             interval at this point and calling the */
/*                             INTEGRATOR on the SUBRANGES. If possible, */
/*                             AN APPROPRIATE SPECIAL-PURPOSE INTEGRATOR */
/*                             should be used which is designed for */
/*                             handling the type of difficulty involved. */
/*                         = 2 The occurrence of roundoff error is */
/*                             detected, which prevents the requested */
/*                             tolerance from being achieved. */
/*                         = 3 Extremely bad integrand behaviour occurs */
/*                             at some points of the integration */
/*                             interval. */
/*                         = 6 The input is invalid, because */
/*                             (EPSABS.LE.0 AND */
/*                              EPSREL.LT.MAX(50*REL.MACH.ACC.,0.5D-28)) */
/*                             OR LIMIT.LT.1 OR LENW.LT.LIMIT*4. */
/*                             RESULT, ABSERR, NEVAL, LAST are set */
/*                             to zero. */
/*                             EXCEPT when LENW is invalid, IWORK(1), */
/*                             WORK(LIMIT*2+1) and WORK(LIMIT*3+1) are */
/*                             set to zero, WORK(1) is set to A and */
/*                             WORK(LIMIT+1) to B. */

/*         DIMENSIONING PARAMETERS */
/*            LIMIT - Integer */
/*                    Dimensioning parameter for IWORK */
/*                    Limit determines the maximum number of subintervals */
/*                    in the partition of the given integration interval */
/*                    (A,B), LIMIT.GE.1. */
/*                    If LIMIT.LT.1, the routine will end with IER = 6. */

/*            LENW  - Integer */
/*                    Dimensioning parameter for work */
/*                    LENW must be at least LIMIT*4. */
/*                    IF LENW.LT.LIMIT*4, the routine will end with */
/*                    IER = 6. */

/*            LAST  - Integer */
/*                    On return, LAST equals the number of subintervals */
/*                    produced in the subdivision process, which */
/*                    determines the number of significant elements */
/*                    actually in the WORK ARRAYS. */

/*         WORK ARRAYS */
/*            IWORK - Integer */
/*                    Vector of dimension at least limit, the first K */
/*                    elements of which contain pointers to the error */
/*                    estimates over the subintervals, such that */
/*                    WORK(LIMIT*3+IWORK(1)),... , WORK(LIMIT*3+IWORK(K)) */
/*                    form a decreasing sequence with K = LAST If */
/*                    LAST.LE.(LIMIT/2+2), and K = LIMIT+1-LAST otherwise */

/*            WORK  - Double precision */
/*                    Vector of dimension at least LENW */
/*                    on return */
/*                    WORK(1), ..., WORK(LAST) contain the left end */
/*                    points of the subintervals in the partition of */
/*                     (A,B), */
/*                    WORK(LIMIT+1), ..., WORK(LIMIT+LAST) contain the */
/*                     right end points, */
/*                    WORK(LIMIT*2+1), ..., WORK(LIMIT*2+LAST) contain */
/*                     the integral approximations over the subintervals, */
/*                    WORK(LIMIT*3+1), ..., WORK(LIMIT*3+LAST) contain */
/*                     the error estimates. */

/* ***REFERENCES  (NONE) */
/* ***ROUTINES CALLED  DQAGE, XERMSG */
/* ***REVISION HISTORY  (YYMMDD) */
/*   800101  DATE WRITTEN */
/*   890831  Modified array declarations.  (WRB) */
/*   890831  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ) */
/* ***END PROLOGUE  DQAG */


/* ***FIRST EXECUTABLE STATEMENT  DQAG */
    /* Parameter adjustments */
    --work;
    --iwork;

    /* Function Body */
    *ier = 6;
    *neval = 0;
    *last = 0;
    *result = 0.;
    *abserr = 0.;
    if (*limit >= 1 && *lenw >= *limit << 2) {

/*        PREPARE CALL FOR DQAGE. */

	l1 = *limit + 1;
	l2 = *limit + l1;
	l3 = *limit + l2;

	dqage_((D_fp)f, a, b, epsabs, epsrel, key, limit, result, abserr, 
		neval, ier, &work[1], &work[l1], &work[l2], &work[l3], &iwork[
		1], last);

/*        CALL ERROR HANDLER IF NECESSARY. */

	lvl = 0;
    }

    if (*ier == 6) {
	lvl = 1;
    }
    if (*ier != 0) {
	xermsg_("SLATEC", "DQAG", "ABNORMAL RETURN", ier, &lvl, (ftnlen)6, (
		ftnlen)4, (ftnlen)15);
    }
    return 0;
} /* dqag_ */

