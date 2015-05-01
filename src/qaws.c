/* qaws.f -- translated by f2c (version 12.02.01).
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

/* DECK QAWS */
/* Subroutine */ int qaws_(E_fp f, real *a, real *b, real *alfa, real *beta, 
	integer *integr, real *epsabs, real *epsrel, real *result, real *
	abserr, integer *neval, integer *ier, integer *limit, integer *lenw, 
	integer *last, integer *iwork, real *work)
{
    static integer l1, l2, l3, lvl;
    extern /* Subroutine */ int qawse_(E_fp, real *, real *, real *, real *, 
	    integer *, real *, real *, integer *, real *, real *, integer *, 
	    integer *, real *, real *, real *, real *, integer *, integer *), 
	    xermsg_(char *, char *, char *, integer *, integer *, ftnlen, 
	    ftnlen, ftnlen);

/* ***BEGIN PROLOGUE  QAWS */
/* ***PURPOSE  The routine calculates an approximation result to a given */
/*            definite integral I = Integral of F*W over (A,B), */
/*            (where W shows a singular behaviour at the end points */
/*            see parameter INTEGR). */
/*            Hopefully satisfying following claim for accuracy */
/*            ABS(I-RESULT).LE.MAX(EPSABS,EPSREL*ABS(I)). */
/* ***LIBRARY   SLATEC (QUADPACK) */
/* ***CATEGORY  H2A2A1 */
/* ***TYPE      SINGLE PRECISION (QAWS-S, DQAWS-D) */
/* ***KEYWORDS  ALGEBRAIC-LOGARITHMIC END POINT SINGULARITIES, */
/*             AUTOMATIC INTEGRATOR, CLENSHAW-CURTIS METHOD, */
/*             GLOBALLY ADAPTIVE, QUADPACK, QUADRATURE, SPECIAL-PURPOSE */
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
/*                     Upper limit of integration, B.GT.A */
/*                     If B.LE.A, the routine will end with IER = 6. */

/*            ALFA   - Real */
/*                     Parameter in the integrand function, ALFA.GT.(-1) */
/*                     If ALFA.LE.(-1), the routine will end with */
/*                     IER = 6. */

/*            BETA   - Real */
/*                     Parameter in the integrand function, BETA.GT.(-1) */
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

/*            EPSABS - Real */
/*                     Absolute accuracy requested */
/*            EPSREL - Real */
/*                     Relative accuracy requested */
/*                     If  EPSABS.LE.0 */
/*                     and EPSREL.LT.MAX(50*REL.MACH.ACC.,0.5D-28), */
/*                     the routine will end with IER = 6. */

/*         ON RETURN */
/*            RESULT - Real */
/*                     Approximation to the integral */

/*            ABSERR - Real */
/*                     Estimate of the modulus of the absolute error, */
/*                     Which should equal or exceed ABS(I-RESULT) */

/*            NEVAL  - Integer */
/*                     Number of integrand evaluations */

/*            IER    - Integer */
/*                     IER = 0 Normal and reliable termination of the */
/*                             routine. It is assumed that the requested */
/*                             accuracy has been achieved. */
/*                     IER.GT.0 Abnormal termination of the routine */
/*                             The estimates for the integral and error */
/*                             are less reliable. It is assumed that the */
/*                             requested accuracy has not been achieved. */
/*            ERROR MESSAGES */
/*                     IER = 1 Maximum number of subdivisions allowed */
/*                             has been achieved. One can allow more */
/*                             subdivisions by increasing the value of */
/*                             LIMIT (and taking the according dimension */
/*                             adjustments into account). However, if */
/*                             this yields no improvement it is advised */
/*                             to analyze the integrand, in order to */
/*                             determine the integration difficulties */
/*                             which prevent the requested tolerance from */
/*                             being achieved. In case of a jump */
/*                             discontinuity or a local singularity */
/*                             of algebraico-logarithmic type at one or */
/*                             more interior points of the integration */
/*                             range, one should proceed by splitting up */
/*                             the interval at these points and calling */
/*                             the integrator on the subranges. */
/*                         = 2 The occurrence of roundoff error is */
/*                             detected, which prevents the requested */
/*                             tolerance from being achieved. */
/*                         = 3 Extremely bad integrand behaviour occurs */
/*                             at some points of the integration */
/*                             interval. */
/*                         = 6 The input is invalid, because */
/*                             B.LE.A or ALFA.LE.(-1) or BETA.LE.(-1) or */
/*                             or INTEGR.LT.1 or INTEGR.GT.4 or */
/*                             (EPSABS.LE.0 and */
/*                              EPSREL.LT.MAX(50*REL.MACH.ACC.,0.5D-28)) */
/*                             or LIMIT.LT.2 or LENW.LT.LIMIT*4. */
/*                             RESULT, ABSERR, NEVAL, LAST are set to */
/*                             zero. Except when LENW or LIMIT is invalid */
/*                             IWORK(1), WORK(LIMIT*2+1) and */
/*                             WORK(LIMIT*3+1) are set to zero, WORK(1) */
/*                             is set to A and WORK(LIMIT+1) to B. */

/*         DIMENSIONING PARAMETERS */
/*            LIMIT  - Integer */
/*                     Dimensioning parameter for IWORK */
/*                     LIMIT determines the maximum number of */
/*                     subintervals in the partition of the given */
/*                     integration interval (A,B), LIMIT.GE.2. */
/*                     If LIMIT.LT.2, the routine will end with IER = 6. */

/*            LENW   - Integer */
/*                     Dimensioning parameter for WORK */
/*                     LENW must be at least LIMIT*4. */
/*                     If LENW.LT.LIMIT*4, the routine will end */
/*                     with IER = 6. */

/*            LAST   - Integer */
/*                     On return, LAST equals the number of */
/*                     subintervals produced in the subdivision process, */
/*                     which determines the significant number of */
/*                     elements actually in the WORK ARRAYS. */

/*         WORK ARRAYS */
/*            IWORK  - Integer */
/*                     Vector of dimension LIMIT, the first K */
/*                     elements of which contain pointers */
/*                     to the error estimates over the subintervals, */
/*                     such that WORK(LIMIT*3+IWORK(1)), ..., */
/*                     WORK(LIMIT*3+IWORK(K)) form a decreasing */
/*                     sequence with K = LAST if LAST.LE.(LIMIT/2+2), */
/*                     and K = LIMIT+1-LAST otherwise */

/*            WORK   - Real */
/*                     Vector of dimension LENW */
/*                     On return */
/*                     WORK(1), ..., WORK(LAST) contain the left */
/*                      end points of the subintervals in the */
/*                      partition of (A,B), */
/*                     WORK(LIMIT+1), ..., WORK(LIMIT+LAST) contain */
/*                      the right end points, */
/*                     WORK(LIMIT*2+1), ..., WORK(LIMIT*2+LAST) */
/*                      contain the integral approximations over */
/*                      the subintervals, */
/*                     WORK(LIMIT*3+1), ..., WORK(LIMIT*3+LAST) */
/*                      contain the error estimates. */

/* ***REFERENCES  (NONE) */
/* ***ROUTINES CALLED  QAWSE, XERMSG */
/* ***REVISION HISTORY  (YYMMDD) */
/*   800101  DATE WRITTEN */
/*   890831  Modified array declarations.  (WRB) */
/*   890831  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ) */
/* ***END PROLOGUE  QAWS */




/*         CHECK VALIDITY OF LIMIT AND LENW. */

/* ***FIRST EXECUTABLE STATEMENT  QAWS */
    /* Parameter adjustments */
    --work;
    --iwork;

    /* Function Body */
    *ier = 6;
    *neval = 0;
    *last = 0;
    *result = 0.f;
    *abserr = 0.f;
    if (*limit < 2 || *lenw < *limit << 2) {
	goto L10;
    }

/*         PREPARE CALL FOR QAWSE. */

    l1 = *limit + 1;
    l2 = *limit + l1;
    l3 = *limit + l2;

    qawse_((E_fp)f, a, b, alfa, beta, integr, epsabs, epsrel, limit, result, 
	    abserr, neval, ier, &work[1], &work[l1], &work[l2], &work[l3], &
	    iwork[1], last);

/*         CALL ERROR HANDLER IF NECESSARY. */

    lvl = 0;
L10:
    if (*ier == 6) {
	lvl = 1;
    }
    if (*ier != 0) {
	xermsg_("SLATEC", "QAWS", "ABNORMAL RETURN", ier, &lvl, (ftnlen)6, (
		ftnlen)4, (ftnlen)15);
    }
    return 0;
} /* qaws_ */

