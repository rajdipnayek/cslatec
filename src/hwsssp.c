/* hwsssp.f -- translated by f2c (version 12.02.01).
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

/* DECK HWSSSP */
/* Subroutine */ int hwsssp_(real *ts, real *tf, integer *m, integer *mbdcnd, 
	real *bdts, real *bdtf, real *ps, real *pf, integer *n, integer *
	nbdcnd, real *bdps, real *bdpf, real *elmbda, real *f, integer *idimf,
	 real *pertrb, integer *ierror, real *w)
{
    /* System generated locals */
    integer f_dim1, f_offset;

    /* Local variables */
    static real pi, dum, tpi;
    extern /* Subroutine */ int hwsss1_(real *, real *, integer *, integer *, 
	    real *, real *, real *, real *, integer *, integer *, real *, 
	    real *, real *, real *, integer *, real *, real *, real *, real *,
	     real *, real *, real *, real *);
    extern doublereal pimach_(real *);

/* ***BEGIN PROLOGUE  HWSSSP */
/* ***PURPOSE  Solve a finite difference approximation to the Helmholtz */
/*            equation in spherical coordinates and on the surface of the */
/*            unit sphere (radius of 1). */
/* ***LIBRARY   SLATEC (FISHPACK) */
/* ***CATEGORY  I2B1A1A */
/* ***TYPE      SINGLE PRECISION (HWSSSP-S) */
/* ***KEYWORDS  ELLIPTIC, FISHPACK, HELMHOLTZ, PDE, SPHERICAL */
/* ***AUTHOR  Adams, J., (NCAR) */
/*           Swarztrauber, P. N., (NCAR) */
/*           Sweet, R., (NCAR) */
/* ***DESCRIPTION */

/*     Subroutine HWSSSP solves a finite difference approximation to the */
/*     Helmholtz equation in spherical coordinates and on the surface of */
/*     the unit sphere (radius of 1): */

/*          (1/SIN(THETA))(d/dTHETA)(SIN(THETA)(dU/dTHETA)) */

/*             + (1/SIN(THETA)**2)(d/dPHI)(dU/dPHI) */

/*             + LAMBDA*U = F(THETA,PHI) */

/*     Where THETA is colatitude and PHI is longitude. */

/*     * * * * * * * *    Parameter Description     * * * * * * * * * * */

/*             * * * * * *   On Input    * * * * * * */

/*     TS,TF */
/*       The range of THETA (colatitude), i.e., TS .LE. THETA .LE. TF. */
/*       TS must be less than TF.  TS and TF are in radians.  A TS of */
/*       zero corresponds to the north pole and a TF of PI corresponds to */
/*       the south pole. */

/*     * * * * * * * * * * * * * * IMPORTANT * * * * * * * * * * * * * * */

/*     If TF is equal to PI then it must be computed using the statement */
/*     TF = PIMACH(DUM). This insures that TF in the users program is */
/*     equal to PI in this program which permits several tests of the */
/*     input parameters that otherwise would not be possible. */


/*     M */
/*       The number of panels into which the interval (TS,TF) is */
/*       subdivided.  Hence, there will be M+1 grid points in the */
/*       THETA-direction given by THETA(I) = (I-1)DTHETA+TS for */
/*       I = 1,2,...,M+1, where DTHETA = (TF-TS)/M is the panel width. */
/*       M must be greater than 5. */

/*     MBDCND */
/*       Indicates the type of boundary condition at THETA = TS and */
/*       THETA = TF. */

/*       = 1  If the solution is specified at THETA = TS and THETA = TF. */
/*       = 2  If the solution is specified at THETA = TS and the */
/*            derivative of the solution with respect to THETA is */
/*            specified at THETA = TF (see note 2 below). */
/*       = 3  If the derivative of the solution with respect to THETA is */
/*            specified at THETA = TS and THETA = TF (see notes 1,2 */
/*            below). */
/*       = 4  If the derivative of the solution with respect to THETA is */
/*            specified at THETA = TS (see note 1 below) and the */
/*            solution is specified at THETA = TF. */
/*       = 5  If the solution is unspecified at THETA = TS = 0 and the */
/*            solution is specified at THETA = TF. */
/*       = 6  If the solution is unspecified at THETA = TS = 0 and the */
/*            derivative of the solution with respect to THETA is */
/*            specified at THETA = TF (see note 2 below). */
/*       = 7  If the solution is specified at THETA = TS and the */
/*            solution is unspecified at THETA = TF = PI. */
/*       = 8  If the derivative of the solution with respect to THETA is */
/*            specified at THETA = TS (see note 1 below) and the */
/*            solution is unspecified at THETA = TF = PI. */
/*       = 9  If the solution is unspecified at THETA = TS = 0 and */
/*            THETA = TF = PI. */

/*       NOTES:  1.  If TS = 0, do not use MBDCND = 3,4, or 8, but */
/*                   instead use MBDCND = 5,6, or 9  . */
/*               2.  If TF = PI, do not use MBDCND = 2,3, or 6, but */
/*                   instead use MBDCND = 7,8, or 9  . */

/*     BDTS */
/*       A one-dimensional array of length N+1 that specifies the values */
/*       of the derivative of the solution with respect to THETA at */
/*       THETA = TS.  When MBDCND = 3,4, or 8, */

/*            BDTS(J) = (d/dTHETA)U(TS,PHI(J)), J = 1,2,...,N+1  . */

/*       When MBDCND has any other value, BDTS is a dummy variable. */

/*     BDTF */
/*       A one-dimensional array of length N+1 that specifies the values */
/*       of the derivative of the solution with respect to THETA at */
/*       THETA = TF.  When MBDCND = 2,3, or 6, */

/*            BDTF(J) = (d/dTHETA)U(TF,PHI(J)), J = 1,2,...,N+1  . */

/*       When MBDCND has any other value, BDTF is a dummy variable. */

/*     PS,PF */
/*       The range of PHI (longitude), i.e., PS .LE. PHI .LE. PF.  PS */
/*       must be less than PF.  PS and PF are in radians.  If PS = 0 and */
/*       PF = 2*PI, periodic boundary conditions are usually prescribed. */

/*     * * * * * * * * * * * * * * IMPORTANT * * * * * * * * * * * * * * */

/*     If PF is equal to 2*PI then it must be computed using the */
/*     statement PF = 2.*PIMACH(DUM). This insures that PF in the users */
/*     program is equal to 2*PI in this program which permits tests of */
/*     the input parameters that otherwise would not be possible. */


/*     N */
/*       The number of panels into which the interval (PS,PF) is */
/*       subdivided.  Hence, there will be N+1 grid points in the */
/*       PHI-direction given by PHI(J) = (J-1)DPHI+PS  for */
/*       J = 1,2,...,N+1, where DPHI = (PF-PS)/N is the panel width. */
/*       N must be greater than 4. */

/*     NBDCND */
/*       Indicates the type of boundary condition at PHI = PS and */
/*       PHI = PF. */

/*       = 0  If the solution is periodic in PHI, i.e., */
/*            U(I,J) = U(I,N+J). */
/*       = 1  If the solution is specified at PHI = PS and PHI = PF */
/*            (see note below). */
/*       = 2  If the solution is specified at PHI = PS (see note below) */
/*            and the derivative of the solution with respect to PHI is */
/*            specified at PHI = PF. */
/*       = 3  If the derivative of the solution with respect to PHI is */
/*            specified at PHI = PS and PHI = PF. */
/*       = 4  If the derivative of the solution with respect to PHI is */
/*            specified at PS and the solution is specified at PHI = PF */
/*            (see note below). */

/*       NOTE:  NBDCND = 1,2, or 4 cannot be used with */
/*              MBDCND = 5,6,7,8, or 9 (the former indicates that the */
/*                       solution is specified at a pole, the latter */
/*                       indicates that the solution is unspecified). */
/*                       Use instead */
/*              MBDCND = 1 or 2  . */

/*     BDPS */
/*       A one-dimensional array of length M+1 that specifies the values */
/*       of the derivative of the solution with respect to PHI at */
/*       PHI = PS.  When NBDCND = 3 or 4, */

/*            BDPS(I) = (d/dPHI)U(THETA(I),PS), I = 1,2,...,M+1  . */

/*       When NBDCND has any other value, BDPS is a dummy variable. */

/*     BDPF */
/*       A one-dimensional array of length M+1 that specifies the values */
/*       of the derivative of the solution with respect to PHI at */
/*       PHI = PF.  When NBDCND = 2 or 3, */

/*            BDPF(I) = (d/dPHI)U(THETA(I),PF), I = 1,2,...,M+1  . */

/*       When NBDCND has any other value, BDPF is a dummy variable. */

/*     ELMBDA */
/*       The constant LAMBDA in the Helmholtz equation.  If */
/*       LAMBDA .GT. 0, a solution may not exist.  However, HWSSSP will */
/*       attempt to find a solution. */

/*     F */
/*       A two-dimensional array that specifies the value of the right */
/*       side of the Helmholtz equation and boundary values (if any). */
/*       For I = 2,3,...,M  and  J = 2,3,...,N */

/*            F(I,J) = F(THETA(I),PHI(J)). */

/*       On the boundaries F is defined by */

/*            MBDCND   F(1,J)            F(M+1,J) */
/*            ------   ------------      ------------ */

/*              1      U(TS,PHI(J))      U(TF,PHI(J)) */
/*              2      U(TS,PHI(J))      F(TF,PHI(J)) */
/*              3      F(TS,PHI(J))      F(TF,PHI(J)) */
/*              4      F(TS,PHI(J))      U(TF,PHI(J)) */
/*              5      F(0,PS)           U(TF,PHI(J))   J = 1,2,...,N+1 */
/*              6      F(0,PS)           F(TF,PHI(J)) */
/*              7      U(TS,PHI(J))      F(PI,PS) */
/*              8      F(TS,PHI(J))      F(PI,PS) */
/*              9      F(0,PS)           F(PI,PS) */

/*            NBDCND   F(I,1)            F(I,N+1) */
/*            ------   --------------    -------------- */

/*              0      F(THETA(I),PS)    F(THETA(I),PS) */
/*              1      U(THETA(I),PS)    U(THETA(I),PF) */
/*              2      U(THETA(I),PS)    F(THETA(I),PF)   I = 1,2,...,M+1 */
/*              3      F(THETA(I),PS)    F(THETA(I),PF) */
/*              4      F(THETA(I),PS)    U(THETA(I),PF) */

/*       F must be dimensioned at least (M+1)*(N+1). */

/*      *NOTE* */

/*       If the table calls for both the solution U and the right side F */
/*       at a corner then the solution must be specified. */


/*     IDIMF */
/*       The row (or first) dimension of the array F as it appears in the */
/*       program calling HWSSSP.  This parameter is used to specify the */
/*       variable dimension of F.  IDIMF must be at least M+1  . */

/*     W */
/*       A one-dimensional array that must be provided by the user for */
/*       work space. W may require up to 4*(N+1)+(16+INT(log2(N+1)))(M+1) */
/*       locations. The actual number of locations used is computed by */
/*       HWSSSP and is output in location W(1). INT( ) denotes the */
/*       FORTRAN integer function. */


/*     * * * * * * * * * *     On Output     * * * * * * * * * * */

/*     F */
/*       Contains the solution U(I,J) of the finite difference */
/*       approximation for the grid point (THETA(I),PHI(J)), */
/*       I = 1,2,...,M+1,   J = 1,2,...,N+1  . */

/*     PERTRB */
/*       If one specifies a combination of periodic, derivative or */
/*       unspecified boundary conditions for a Poisson equation */
/*       (LAMBDA = 0), a solution may not exist.  PERTRB is a constant, */
/*       calculated and subtracted from F, which ensures that a solution */
/*       exists.  HWSSSP then computes this solution, which is a least */
/*       squares solution to the original approximation.  This solution */
/*       is not unique and is unnormalized. The value of PERTRB should */
/*       be small compared to the right side F. Otherwise , a solution */
/*       is obtained to an essentially different problem. This comparison */
/*       should always be made to insure that a meaningful solution has */
/*       been obtained. */

/*     IERROR */
/*       An error flag that indicates invalid input parameters.  Except */
/*       for numbers 0 and 8, a solution is not attempted. */

/*       = 0  No error */
/*       = 1  TS.LT.0 or TF.GT.PI */
/*       = 2  TS.GE.TF */
/*       = 3  MBDCND.LT.1 or MBDCND.GT.9 */
/*       = 4  PS.LT.0 or PS.GT.PI+PI */
/*       = 5  PS.GE.PF */
/*       = 6  N.LT.5 */
/*       = 7  M.LT.5 */
/*       = 8  NBDCND.LT.0 or NBDCND.GT.4 */
/*       = 9  ELMBDA.GT.0 */
/*       = 10 IDIMF.LT.M+1 */
/*       = 11 NBDCND equals 1,2 or 4 and MBDCND.GE.5 */
/*       = 12 TS.EQ.0 and MBDCND equals 3,4 or 8 */
/*       = 13 TF.EQ.PI and MBDCND equals 2,3 or 6 */
/*       = 14 MBDCND equals 5,6 or 9 and TS.NE.0 */
/*       = 15 MBDCND.GE.7 and TF.NE.PI */

/*       Since this is the only means of indicating a possibly incorrect */
/*       call to HWSSSP, the user should test IERROR after a call. */

/*     W */
/*       Contains intermediate values that must not be destroyed if */
/*       HWSSSP will be called again with INTL = 1. W(1) contains the */
/*       required length of W . */

/* *Long Description: */

/*     * * * * * * *   Program Specifications    * * * * * * * * * * * * */

/*     Dimension of   BDTS(N+1),BDTF(N+1),BDPS(M+1),BDPF(M+1), */
/*     Arguments      F(IDIMF,N+1),W(see argument list) */

/*     Latest         January 1978 */
/*     Revision */


/*     Subprograms    HWSSSP,HWSSS1,GENBUN,POISD2,POISN2,POISP2,COSGEN,ME */
/*     Required       TRIX,TRI3,PIMACH */

/*     Special        NONE */
/*     Conditions */

/*     Common         NONE */
/*     Blocks */

/*     I/O            NONE */

/*     Precision      Single */

/*     Specialist     Paul Swarztrauber */

/*     Language       FORTRAN */

/*     History        Version 1 - September 1973 */
/*                    Version 2 - April     1976 */
/*                    Version 3 - January   1978 */

/*     Algorithm      The routine defines the finite difference */
/*                    equations, incorporates boundary data, and adjusts */
/*                    the right side of singular systems and then calls */
/*                    GENBUN to solve the system. */

/*     Space */
/*     Required       CONTROL DATA 7600 */

/*     Timing and        The execution time T on the NCAR Control Data */
/*     Accuracy       7600 for subroutine HWSSSP is roughly proportional */
/*                    to M*N*log2(N), but also depends on the input */
/*                    parameters NBDCND and MBDCND.  Some typical values */
/*                    are listed in the table below. */
/*                       The solution process employed results in a loss */
/*                    of no more than three significant digits for N and */
/*                    M as large as 64.  More detailed information about */
/*                    accuracy can be found in the documentation for */
/*                    subroutine GENBUN which is the routine that */
/*                    solves the finite difference equations. */


/*                       M(=N)    MBDCND    NBDCND    T(MSECS) */
/*                       -----    ------    ------    -------- */

/*                        32        0         0          31 */
/*                        32        1         1          23 */
/*                        32        3         3          36 */
/*                        64        0         0         128 */
/*                        64        1         1          96 */
/*                        64        3         3         142 */

/*     Portability    American National Standards Institute FORTRAN. */
/*                    The machine dependent constant PI is defined in */
/*                    function PIMACH. */

/*     Required       SIN,COS */
/*     Resident */
/*     Routines */

/*     References     P. N. Swarztrauber,'The Direct Solution Of The */
/*                    Discrete Poisson Equation On The Surface Of a */
/*                    Sphere, SIAM J. Numer. Anal.,15(1974), pp 212-215 */

/*                    Swarztrauber,P. and R. Sweet, 'Efficient FORTRAN */
/*                    Subprograms for The Solution of Elliptic Equations' */
/*                    NCAR TN/IA-109, July, 1975, 138 pp. */

/*     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/* ***REFERENCES  P. N. Swarztrauber and R. Sweet, Efficient Fortran */
/*                 subprograms for the solution of elliptic equations, */
/*                 NCAR TN/IA-109, July 1975, 138 pp. */
/*               P. N. Swarztrauber, The direct solution of the discrete */
/*                 Poisson equation on the surface of a sphere, SIAM */
/*                 Journal on Numerical Analysis 15 (1974), pp. 212-215. */
/* ***ROUTINES CALLED  HWSSS1, PIMACH */
/* ***REVISION HISTORY  (YYMMDD) */
/*   801001  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   891009  Removed unreferenced variable.  (WRB) */
/*   891009  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   920501  Reformatted the REFERENCES section.  (WRB) */
/* ***END PROLOGUE  HWSSSP */

/* ***FIRST EXECUTABLE STATEMENT  HWSSSP */
    /* Parameter adjustments */
    --bdts;
    --bdtf;
    --bdps;
    --bdpf;
    f_dim1 = *idimf;
    f_offset = 1 + f_dim1;
    f -= f_offset;
    --w;

    /* Function Body */
    pi = pimach_(&dum);
    tpi = pi * 2.f;
    *ierror = 0;
    if (*ts < 0.f || *tf > pi) {
	*ierror = 1;
    }
    if (*ts >= *tf) {
	*ierror = 2;
    }
    if (*mbdcnd < 1 || *mbdcnd > 9) {
	*ierror = 3;
    }
    if (*ps < 0.f || *pf > tpi) {
	*ierror = 4;
    }
    if (*ps >= *pf) {
	*ierror = 5;
    }
    if (*n < 5) {
	*ierror = 6;
    }
    if (*m < 5) {
	*ierror = 7;
    }
    if (*nbdcnd < 0 || *nbdcnd > 4) {
	*ierror = 8;
    }
    if (*elmbda > 0.f) {
	*ierror = 9;
    }
    if (*idimf < *m + 1) {
	*ierror = 10;
    }
    if ((*nbdcnd == 1 || *nbdcnd == 2 || *nbdcnd == 4) && *mbdcnd >= 5) {
	*ierror = 11;
    }
    if (*ts == 0.f && (*mbdcnd == 3 || *mbdcnd == 4 || *mbdcnd == 8)) {
	*ierror = 12;
    }
    if (*tf == pi && (*mbdcnd == 2 || *mbdcnd == 3 || *mbdcnd == 6)) {
	*ierror = 13;
    }
    if ((*mbdcnd == 5 || *mbdcnd == 6 || *mbdcnd == 9) && *ts != 0.f) {
	*ierror = 14;
    }
    if (*mbdcnd >= 7 && *tf != pi) {
	*ierror = 15;
    }
    if (*ierror != 0 && *ierror != 9) {
	return 0;
    }
    hwsss1_(ts, tf, m, mbdcnd, &bdts[1], &bdtf[1], ps, pf, n, nbdcnd, &bdps[1]
	    , &bdpf[1], elmbda, &f[f_offset], idimf, pertrb, &w[1], &w[*m + 2]
	    , &w[(*m << 1) + 3], &w[*m * 3 + 4], &w[(*m << 2) + 5], &w[*m * 5 
	    + 6], &w[*m * 6 + 7]);
    w[1] = w[*m * 6 + 7] + (*m + 1) * 6;
    return 0;
} /* hwsssp_ */

