/* dqawo.f -- translated by f2c (version 12.02.01).
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

/* DECK DQAWO */
/* Subroutine */ int dqawo_(D_fp f, doublereal *a, doublereal *b, doublereal *
	omega, integer *integr, doublereal *epsabs, doublereal *epsrel, 
	doublereal *result, doublereal *abserr, integer *neval, integer *ier, 
	integer *leniw, integer *maxp1, integer *lenw, integer *last, integer 
	*iwork, doublereal *work)
{
    static integer l1, l2, l3, l4, lvl, limit;
    extern /* Subroutine */ int dqawoe_(D_fp, doublereal *, doublereal *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    integer *, integer *, doublereal *, doublereal *, integer *, 
	    integer *, integer *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, integer *, integer *, integer *, doublereal *);
    static integer momcom;
    extern /* Subroutine */ int xermsg_(char *, char *, char *, integer *, 
	    integer *, ftnlen, ftnlen, ftnlen);

/* ***BEGIN PROLOGUE  DQAWO */
/* ***PURPOSE  Calculate an approximation to a given definite integral */
/*            I= Integral of F(X)*W(X) over (A,B), where */
/*                   W(X) = COS(OMEGA*X) */
/*               or  W(X) = SIN(OMEGA*X), */
/*            hopefully satisfying the following claim for accuracy */
/*                ABS(I-RESULT).LE.MAX(EPSABS,EPSREL*ABS(I)). */
/* ***LIBRARY   SLATEC (QUADPACK) */
/* ***CATEGORY  H2A2A1 */
/* ***TYPE      DOUBLE PRECISION (QAWO-S, DQAWO-D) */
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

/*        Computation of oscillatory integrals */
/*        Standard fortran subroutine */
/*        Double precision version */

/*        PARAMETERS */
/*         ON ENTRY */
/*            F      - Double precision */
/*                     Function subprogram defining the function */
/*                     F(X).  The actual name for F needs to be */
/*                     declared E X T E R N A L in the driver program. */

/*            A      - Double precision */
/*                     Lower limit of integration */

/*            B      - Double precision */
/*                     Upper limit of integration */

/*            OMEGA  - Double precision */
/*                     Parameter in the integrand weight function */

/*            INTEGR - Integer */
/*                     Indicates which of the weight functions is used */
/*                     INTEGR = 1      W(X) = COS(OMEGA*X) */
/*                     INTEGR = 2      W(X) = SIN(OMEGA*X) */
/*                     If INTEGR.NE.1.AND.INTEGR.NE.2, the routine will */
/*                     end with IER = 6. */

/*            EPSABS - Double precision */
/*                     Absolute accuracy requested */
/*            EPSREL - Double precision */
/*                     Relative accuracy requested */
/*                     If EPSABS.LE.0 and */
/*                     EPSREL.LT.MAX(50*REL.MACH.ACC.,0.5D-28), */
/*                     the routine will end with IER = 6. */

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
/*                   - IER.GT.0 Abnormal termination of the routine. */
/*                             The estimates for integral and error are */
/*                             less reliable. It is assumed that the */
/*                             requested accuracy has not been achieved. */
/*            ERROR MESSAGES */
/*                     IER = 1 Maximum number of subdivisions allowed */
/*                             has been achieved (= LENIW/2). One can */
/*                             allow more subdivisions by increasing the */
/*                             value of LENIW (and taking the according */
/*                             dimension adjustments into account). */
/*                             However, if this yields no improvement it */
/*                             is advised to analyze the integrand in */
/*                             order to determine the integration */
/*                             difficulties. If the position of a local */
/*                             difficulty can be determined (e.g. */
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
/*                             The error may be under-estimated. */
/*                         = 3 Extremely bad integrand behaviour occurs */
/*                             at some interior points of the */
/*                             integration interval. */
/*                         = 4 The algorithm does not converge. */
/*                             Roundoff error is detected in the */
/*                             extrapolation table. It is presumed that */
/*                             the requested tolerance cannot be achieved */
/*                             due to roundoff in the extrapolation */
/*                             table, and that the returned result is */
/*                             the best which can be obtained. */
/*                         = 5 The integral is probably divergent, or */
/*                             slowly convergent. It must be noted that */
/*                             divergence can occur with any other value */
/*                             of IER. */
/*                         = 6 The input is invalid, because */
/*                             (EPSABS.LE.0 and */
/*                              EPSREL.LT.MAX(50*REL.MACH.ACC.,0.5D-28)) */
/*                             or (INTEGR.NE.1 AND INTEGR.NE.2), */
/*                             or LENIW.LT.2 OR MAXP1.LT.1 or */
/*                             LENW.LT.LENIW*2+MAXP1*25. */
/*                             RESULT, ABSERR, NEVAL, LAST are set to */
/*                             zero. Except when LENIW, MAXP1 or LENW are */
/*                             invalid, WORK(LIMIT*2+1), WORK(LIMIT*3+1), */
/*                             IWORK(1), IWORK(LIMIT+1) are set to zero, */
/*                             WORK(1) is set to A and WORK(LIMIT+1) to */
/*                             B. */

/*         DIMENSIONING PARAMETERS */
/*            LENIW  - Integer */
/*                     Dimensioning parameter for IWORK. */
/*                     LENIW/2 equals the maximum number of subintervals */
/*                     allowed in the partition of the given integration */
/*                     interval (A,B), LENIW.GE.2. */
/*                     If LENIW.LT.2, the routine will end with IER = 6. */

/*            MAXP1  - Integer */
/*                     Gives an upper bound on the number of Chebyshev */
/*                     moments which can be stored, i.e. for the */
/*                     intervals of lengths ABS(B-A)*2**(-L), */
/*                     L=0,1, ..., MAXP1-2, MAXP1.GE.1 */
/*                     If MAXP1.LT.1, the routine will end with IER = 6. */

/*            LENW   - Integer */
/*                     Dimensioning parameter for WORK */
/*                     LENW must be at least LENIW*2+MAXP1*25. */
/*                     If LENW.LT.(LENIW*2+MAXP1*25), the routine will */
/*                     end with IER = 6. */

/*            LAST   - Integer */
/*                     On return, LAST equals the number of subintervals */
/*                     produced in the subdivision process, which */
/*                     determines the number of significant elements */
/*                     actually in the WORK ARRAYS. */

/*         WORK ARRAYS */
/*            IWORK  - Integer */
/*                     Vector of dimension at least LENIW */
/*                     on return, the first K elements of which contain */
/*                     pointers to the error estimates over the */
/*                     subintervals, such that WORK(LIMIT*3+IWORK(1)), .. */
/*                     WORK(LIMIT*3+IWORK(K)) form a decreasing */
/*                     sequence, with LIMIT = LENW/2 , and K = LAST */
/*                     if LAST.LE.(LIMIT/2+2), and K = LIMIT+1-LAST */
/*                     otherwise. */
/*                     Furthermore, IWORK(LIMIT+1), ..., IWORK(LIMIT+ */
/*                     LAST) indicate the subdivision levels of the */
/*                     subintervals, such that IWORK(LIMIT+I) = L means */
/*                     that the subinterval numbered I is of length */
/*                     ABS(B-A)*2**(1-L). */

/*            WORK   - Double precision */
/*                     Vector of dimension at least LENW */
/*                     On return */
/*                     WORK(1), ..., WORK(LAST) contain the left */
/*                      end points of the subintervals in the */
/*                      partition of (A,B), */
/*                     WORK(LIMIT+1), ..., WORK(LIMIT+LAST) contain */
/*                      the right end points, */
/*                     WORK(LIMIT*2+1), ..., WORK(LIMIT*2+LAST) contain */
/*                      the integral approximations over the */
/*                      subintervals, */
/*                     WORK(LIMIT*3+1), ..., WORK(LIMIT*3+LAST) */
/*                      contain the error estimates. */
/*                     WORK(LIMIT*4+1), ..., WORK(LIMIT*4+MAXP1*25) */
/*                      Provide space for storing the Chebyshev moments. */
/*                     Note that LIMIT = LENW/2. */

/* ***REFERENCES  (NONE) */
/* ***ROUTINES CALLED  DQAWOE, XERMSG */
/* ***REVISION HISTORY  (YYMMDD) */
/*   800101  DATE WRITTEN */
/*   890831  Modified array declarations.  (WRB) */
/*   890831  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ) */
/* ***END PROLOGUE  DQAWO */




/*         CHECK VALIDITY OF LENIW, MAXP1 AND LENW. */

/* ***FIRST EXECUTABLE STATEMENT  DQAWO */
    /* Parameter adjustments */
    --work;
    --iwork;

    /* Function Body */
    *ier = 6;
    *neval = 0;
    *last = 0;
    *result = 0.;
    *abserr = 0.;
    if (*leniw < 2 || *maxp1 < 1 || *lenw < (*leniw << 1) + *maxp1 * 25) {
	goto L10;
    }

/*         PREPARE CALL FOR DQAWOE */

    limit = *leniw / 2;
    l1 = limit + 1;
    l2 = limit + l1;
    l3 = limit + l2;
    l4 = limit + l3;
    dqawoe_((D_fp)f, a, b, omega, integr, epsabs, epsrel, &limit, &c__1, 
	    maxp1, result, abserr, neval, ier, last, &work[1], &work[l1], &
	    work[l2], &work[l3], &iwork[1], &iwork[l1], &momcom, &work[l4]);

/*         CALL ERROR HANDLER IF NECESSARY */

    lvl = 0;
L10:
    if (*ier == 6) {
	lvl = 0;
    }
    if (*ier != 0) {
	xermsg_("SLATEC", "DQAWO", "ABNORMAL RETURN", ier, &lvl, (ftnlen)6, (
		ftnlen)5, (ftnlen)15);
    }
    return 0;
} /* dqawo_ */

