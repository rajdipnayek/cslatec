/* qawc.f -- translated by f2c (version 12.02.01).
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

/* DECK QAWC */
/* Subroutine */ int qawc_(E_fp f, real *a, real *b, real *c__, real *epsabs, 
	real *epsrel, real *result, real *abserr, integer *neval, integer *
	ier, integer *limit, integer *lenw, integer *last, integer *iwork, 
	real *work)
{
    static integer l1, l2, l3, lvl;
    extern /* Subroutine */ int qawce_(E_fp, real *, real *, real *, real *, 
	    real *, integer *, real *, real *, integer *, integer *, real *, 
	    real *, real *, real *, integer *, integer *), xermsg_(char *, 
	    char *, char *, integer *, integer *, ftnlen, ftnlen, ftnlen);

/* ***BEGIN PROLOGUE  QAWC */
/* ***PURPOSE  The routine calculates an approximation result to a */
/*            Cauchy principal value I = INTEGRAL of F*W over (A,B) */
/*            (W(X) = 1/((X-C), C.NE.A, C.NE.B), hopefully satisfying */
/*            following claim for accuracy */
/*            ABS(I-RESULT).LE.MAX(EPSABE,EPSREL*ABS(I)). */
/* ***LIBRARY   SLATEC (QUADPACK) */
/* ***CATEGORY  H2A2A1, J4 */
/* ***TYPE      SINGLE PRECISION (QAWC-S, DQAWC-D) */
/* ***KEYWORDS  AUTOMATIC INTEGRATOR, CAUCHY PRINCIPAL VALUE, */
/*             CLENSHAW-CURTIS METHOD, GLOBALLY ADAPTIVE, QUADPACK, */
/*             QUADRATURE, SPECIAL-PURPOSE */
/* ***AUTHOR  Piessens, Robert */
/*             Applied Mathematics and Programming Division */
/*             K. U. Leuven */
/*           de Doncker, Elise */
/*             Applied Mathematics and Programming Division */
/*             K. U. Leuven */
/* ***DESCRIPTION */

/*        Computation of a Cauchy principal value */
/*        Standard fortran subroutine */
/*        Real version */


/*        PARAMETERS */
/*         ON ENTRY */
/*            F      - Real */
/*                     Function subprogram defining the integrand */
/*                     Function F(X). The actual name for F needs to be */
/*                     declared E X T E R N A L in the driver program. */

/*            A      - Real */
/*                     Under limit of integration */

/*            B      - Real */
/*                     Upper limit of integration */

/*            C      - Parameter in the weight function, C.NE.A, C.NE.B. */
/*                     If C = A or C = B, the routine will end with */
/*                     IER = 6 . */

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
/*                     Estimate or the modulus of the absolute error, */
/*                     Which should equal or exceed ABS(I-RESULT) */

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
/*                             divisions by increasing the value of LIMIT */
/*                             (and taking the according dimension */
/*                             adjustments into account). However, if */
/*                             this yields no improvement it is advised */
/*                             to analyze the integrand in order to */
/*                             determine the integration difficulties. */
/*                             If the position of a local difficulty */
/*                             can be determined (e.g. SINGULARITY, */
/*                             DISCONTINUITY within the interval) one */
/*                             will probably gain from splitting up the */
/*                             interval at this point and calling */
/*                             appropriate integrators on the subranges. */
/*                         = 2 The occurrence of roundoff error is detec- */
/*                             ted, which prevents the requested */
/*                             tolerance from being achieved. */
/*                         = 3 Extremely bad integrand behaviour occurs */
/*                             at some points of the integration */
/*                             interval. */
/*                         = 6 The input is invalid, because */
/*                             C = A or C = B or */
/*                             (EPSABS.LE.0 and */
/*                              EPSREL.LT.MAX(50*REL.MACH.ACC.,0.5D-28)) */
/*                             or LIMIT.LT.1 or LENW.LT.LIMIT*4. */
/*                             RESULT, ABSERR, NEVAL, LAST are set to */
/*                             zero.  Except when LENW or LIMIT is */
/*                             invalid, IWORK(1), WORK(LIMIT*2+1) and */
/*                             WORK(LIMIT*3+1) are set to zero, WORK(1) */
/*                             is set to A and WORK(LIMIT+1) to B. */

/*         DIMENSIONING PARAMETERS */
/*            LIMIT - Integer */
/*                    Dimensioning parameter for IWORK */
/*                    LIMIT determines the maximum number of subintervals */
/*                    in the partition of the given integration interval */
/*                    (A,B), LIMIT.GE.1. */
/*                    If LIMIT.LT.1, the routine will end with IER = 6. */

/*           LENW   - Integer */
/*                    Dimensioning parameter for WORK */
/*                    LENW must be at least LIMIT*4. */
/*                    If LENW.LT.LIMIT*4, the routine will end with */
/*                    IER = 6. */

/*            LAST  - Integer */
/*                    On return, LAST equals the number of subintervals */
/*                    produced in the subdivision process, which */
/*                    determines the number of significant elements */
/*                    actually in the WORK ARRAYS. */

/*         WORK ARRAYS */
/*            IWORK - Integer */
/*                    Vector of dimension at least LIMIT, the first K */
/*                    elements of which contain pointers */
/*                    to the error estimates over the subintervals, */
/*                    such that WORK(LIMIT*3+IWORK(1)), ... , */
/*                    WORK(LIMIT*3+IWORK(K)) form a decreasing */
/*                    sequence, with K = LAST if LAST.LE.(LIMIT/2+2), */
/*                    and K = LIMIT+1-LAST otherwise */

/*            WORK  - Real */
/*                    Vector of dimension at least LENW */
/*                    On return */
/*                    WORK(1), ..., WORK(LAST) contain the left */
/*                     end points of the subintervals in the */
/*                     partition of (A,B), */
/*                    WORK(LIMIT+1), ..., WORK(LIMIT+LAST) contain */
/*                     the right end points, */
/*                    WORK(LIMIT*2+1), ..., WORK(LIMIT*2+LAST) contain */
/*                     the integral approximations over the subintervals, */
/*                    WORK(LIMIT*3+1), ..., WORK(LIMIT*3+LAST) */
/*                     contain the error estimates. */

/* ***REFERENCES  (NONE) */
/* ***ROUTINES CALLED  QAWCE, XERMSG */
/* ***REVISION HISTORY  (YYMMDD) */
/*   800101  DATE WRITTEN */
/*   890831  Modified array declarations.  (WRB) */
/*   890831  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ) */
/* ***END PROLOGUE  QAWC */




/*         CHECK VALIDITY OF LIMIT AND LENW. */

/* ***FIRST EXECUTABLE STATEMENT  QAWC */
    /* Parameter adjustments */
    --work;
    --iwork;

    /* Function Body */
    *ier = 6;
    *neval = 0;
    *last = 0;
    *result = 0.f;
    *abserr = 0.f;
    if (*limit < 1 || *lenw < *limit << 2) {
	goto L10;
    }

/*         PREPARE CALL FOR QAWCE. */

    l1 = *limit + 1;
    l2 = *limit + l1;
    l3 = *limit + l2;
    qawce_((E_fp)f, a, b, c__, epsabs, epsrel, limit, result, abserr, neval, 
	    ier, &work[1], &work[l1], &work[l2], &work[l3], &iwork[1], last);

/*         CALL ERROR HANDLER IF NECESSARY. */

    lvl = 0;
L10:
    if (*ier == 6) {
	lvl = 1;
    }
    if (*ier != 0) {
	xermsg_("SLATEC", "QAWC", "ABNORMAL RETURN", ier, &lvl, (ftnlen)6, (
		ftnlen)4, (ftnlen)15);
    }
    return 0;
} /* qawc_ */

