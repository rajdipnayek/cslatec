/* hstcsp.f -- translated by f2c (version 12.02.01).
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

/* DECK HSTCSP */
/* Subroutine */ int hstcsp_(integer *intl, real *a, real *b, integer *m, 
	integer *mbdcnd, real *bda, real *bdb, real *c__, real *d__, integer *
	n, integer *nbdcnd, real *bdc, real *bdd, real *elmbda, real *f, 
	integer *idimf, real *pertrb, integer *ierror, real *w)
{
    /* System generated locals */
    integer f_dim1, f_offset;

    /* Local variables */
    static real pi, dum;
    static integer iwbm, iwcm, iwan, iwbn, iwcn, ierr1, iwwrk, iwrsq;
    extern /* Subroutine */ int hstcs1_(integer *, real *, real *, integer *, 
	    integer *, real *, real *, real *, real *, integer *, integer *, 
	    real *, real *, real *, real *, integer *, real *, integer *, 
	    real *, real *, real *, real *, real *, real *, real *, real *, 
	    real *);
    extern doublereal pimach_(real *);
    static integer iwsnth;

/* ***BEGIN PROLOGUE  HSTCSP */
/* ***PURPOSE  Solve the standard five-point finite difference */
/*            approximation on a staggered grid to the modified Helmholtz */
/*            equation in spherical coordinates assuming axisymmetry */
/*            (no dependence on longitude). */
/* ***LIBRARY   SLATEC (FISHPACK) */
/* ***CATEGORY  I2B1A1A */
/* ***TYPE      SINGLE PRECISION (HSTCSP-S) */
/* ***KEYWORDS  ELLIPTIC, FISHPACK, HELMHOLTZ, PDE, SPHERICAL */
/* ***AUTHOR  Adams, J., (NCAR) */
/*           Swarztrauber, P. N., (NCAR) */
/*           Sweet, R., (NCAR) */
/* ***DESCRIPTION */

/*     HSTCSP solves the standard five-point finite difference */
/*     approximation on a staggered grid to the modified Helmholtz */
/*     equation spherical coordinates assuming axisymmetry (no dependence */
/*     on longitude). */

/*                  (1/R**2)(d/dR)(R**2(dU/dR)) + */

/*       1/(R**2*SIN(THETA))(d/dTHETA)(SIN(THETA)(dU/dTHETA)) + */

/*            (LAMBDA/(R*SIN(THETA))**2)U  =  F(THETA,R) */

/*     where THETA is colatitude and R is the radial coordinate. */
/*     This two-dimensional modified Helmholtz equation results from */
/*     the Fourier transform of the three-dimensional Poisson equation. */

/*    * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */


/*    * * * * * * * *    Parameter Description     * * * * * * * * * * */


/*            * * * * * *   On Input    * * * * * * */

/*   INTL */
/*     = 0  On initial entry to HSTCSP or if any of the arguments */
/*          C, D, N, or NBDCND are changed from a previous call. */

/*     = 1  If C, D, N, and NBDCND are all unchanged from previous */
/*          call to HSTCSP. */

/*     NOTE:  A call with INTL = 0 takes approximately 1.5 times as much */
/*            time as a call with INTL = 1.  Once a call with INTL = 0 */
/*            has been made then subsequent solutions corresponding to */
/*            different F, BDA, BDB, BDC, and BDD can be obtained */
/*            faster with INTL = 1 since initialization is not repeated. */

/*   A,B */
/*     The range of THETA (colatitude), i.e. A .LE. THETA .LE. B.  A */
/*     must be less than B and A must be non-negative.  A and B are in */
/*     radians.  A = 0 corresponds to the north pole and B = PI */
/*     corresponds to the south pole. */

/*                  * * *  IMPORTANT  * * * */

/*     If B is equal to PI, then B must be computed using the statement */

/*     B = PIMACH(DUM) */

/*     This insures that B in the user's program is equal to PI in this */
/*     program which permits several tests of the input parameters that */
/*     otherwise would not be possible. */

/*                  * * * * * * * * * * * * */

/*   M */
/*     The number of grid points in the interval (A,B).  The grid points */
/*     in the THETA-direction are given by THETA(I) = A + (I-0.5)DTHETA */
/*     for I=1,2,...,M where DTHETA =(B-A)/M.  M must be greater than 4. */

/*   MBDCND */
/*     Indicates the type of boundary conditions at THETA = A and */
/*     THETA = B. */

/*     = 1  If the solution is specified at THETA = A and THETA = B. */
/*          (See notes 1, 2 below) */

/*     = 2  If the solution is specified at THETA = A and the derivative */
/*          of the solution with respect to THETA is specified at */
/*          THETA = B (See notes 1, 2 below). */

/*     = 3  If the derivative of the solution with respect to THETA is */
/*          specified at THETA = A (See notes 1, 2 below) and THETA = B. */

/*     = 4  If the derivative of the solution with respect to THETA is */
/*          specified at THETA = A (See notes 1, 2 below) and the */
/*          solution is specified at THETA = B. */

/*     = 5  If the solution is unspecified at THETA = A = 0 and the */
/*          solution is specified at THETA = B. (See note 2 below) */

/*     = 6  If the solution is unspecified at THETA = A = 0 and the */
/*          derivative of the solution with respect to THETA is */
/*          specified at THETA = B (See note 2 below). */

/*     = 7  If the solution is specified at THETA = A and the */
/*          solution is unspecified at THETA = B = PI. */

/*     = 8  If the derivative of the solution with respect to */
/*          THETA is specified at THETA = A (See note 1 below) */
/*          and the solution is unspecified at THETA = B = PI. */

/*     = 9  If the solution is unspecified at THETA = A = 0 and */
/*          THETA = B = PI. */

/*     NOTES:  1.  If A = 0, do not use MBDCND = 1,2,3,4,7 or 8, */
/*                 but instead use MBDCND = 5, 6, or 9. */

/*             2.  if B = PI, do not use MBDCND = 1,2,3,4,5 or 6, */
/*                 but instead use MBDCND = 7, 8, or 9. */

/*             When A = 0  and/or B = PI the only meaningful boundary */
/*             condition is dU/dTHETA = 0.  (See D. Greenspan, 'Numerical */
/*             Analysis of Elliptic Boundary Value Problems,' Harper and */
/*             Row, 1965, Chapter 5.) */

/*   BDA */
/*     A one-dimensional array of length N that specifies the boundary */
/*     values (if any) of the solution at THETA = A.  When */
/*     MBDCND = 1, 2, or 7, */

/*              BDA(J) = U(A,R(J)) ,              J=1,2,...,N. */

/*     When MBDCND = 3, 4, or 8, */

/*              BDA(J) = (d/dTHETA)U(A,R(J)) ,    J=1,2,...,N. */

/*     When MBDCND has any other value, BDA is a dummy variable. */

/*   BDB */
/*     A one-dimensional array of length N that specifies the boundary */
/*     values of the solution at THETA = B.  When MBDCND = 1, 4, or 5, */

/*              BDB(J) = U(B,R(J)) ,              J=1,2,...,N. */

/*     When MBDCND = 2,3, or 6, */

/*              BDB(J) = (d/dTHETA)U(B,R(J)) ,    J=1,2,...,N. */

/*     When MBDCND has any other value, BDB is a dummy variable. */

/*   C,D */
/*     The range of R , i.e. C .LE. R .LE. D. */
/*     C must be less than D.  C must be non-negative. */

/*   N */
/*     The number of unknowns in the interval (C,D).  The unknowns in */
/*     the R-direction are given by R(J) = C + (J-0.5)DR, */
/*     J=1,2,...,N, where DR = (D-C)/N.  N must be greater than 4. */

/*   NBDCND */
/*     Indicates the type of boundary conditions at R = C */
/*     and R = D. */

/*     = 1  If the solution is specified at R = C and R = D. */

/*     = 2  If the solution is specified at R = C and the derivative */
/*          of the solution with respect to R is specified at */
/*          R = D. (See note 1 below) */

/*     = 3  If the derivative of the solution with respect to R is */
/*          specified at R = C and R = D. */

/*     = 4  If the derivative of the solution with respect to R is */
/*          specified at R = C and the solution is specified at */
/*          R = D. */

/*     = 5  If the solution is unspecified at R = C = 0 (See note 2 */
/*          below) and the solution is specified at R = D. */

/*     = 6  If the solution is unspecified at R = C = 0 (See note 2 */
/*          below) and the derivative of the solution with respect to R */
/*          is specified at R = D. */

/*     NOTE 1:  If C = 0 and MBDCND = 3,6,8 or 9, the system of equations */
/*              to be solved is singular.  The unique solution is */
/*              determined by extrapolation to the specification of */
/*              U(THETA(1),C).  But in these cases the right side of the */
/*              system will be perturbed by the constant PERTRB. */

/*     NOTE 2:  NBDCND = 5 or 6 cannot be used with MBDCND = 1, 2, 4, 5, */
/*              or 7 (the former indicates that the solution is */
/*              unspecified at R = 0; the latter indicates that the */
/*              solution is specified).  Use instead NBDCND = 1 or 2. */

/*   BDC */
/*     A one dimensional array of length M that specifies the boundary */
/*     values of the solution at R = C.   When NBDCND = 1 or 2, */

/*              BDC(I) = U(THETA(I),C) ,              I=1,2,...,M. */

/*     When NBDCND = 3 or 4, */

/*              BDC(I) = (d/dR)U(THETA(I),C),         I=1,2,...,M. */

/*     When NBDCND has any other value, BDC is a dummy variable. */

/*   BDD */
/*     A one-dimensional array of length M that specifies the boundary */
/*     values of the solution at R = D.  When NBDCND = 1 or 4, */

/*              BDD(I) = U(THETA(I),D) ,              I=1,2,...,M. */

/*     When NBDCND = 2 or 3, */

/*              BDD(I) = (d/dR)U(THETA(I),D) ,        I=1,2,...,M. */

/*     When NBDCND has any other value, BDD is a dummy variable. */

/*   ELMBDA */
/*     The constant LAMBDA in the modified Helmholtz equation.  If */
/*     LAMBDA is greater than 0, a solution may not exist.  However, */
/*     HSTCSP will attempt to find a solution. */

/*   F */
/*     A two-dimensional array that specifies the values of the right */
/*     side of the modified Helmholtz equation.  For I=1,2,...,M and */
/*     J=1,2,...,N */

/*              F(I,J) = F(THETA(I),R(J)) . */

/*     F must be dimensioned at least M X N. */

/*   IDIMF */
/*     The row (or first) dimension of the array F as it appears in the */
/*     program calling HSTCSP.  This parameter is used to specify the */
/*     variable dimension of F.  IDIMF must be at least M. */

/*   W */
/*     A one-dimensional array that must be provided by the user for */
/*     work space.  With K = INT(log2(N))+1 and L = 2**(K+1), W may */
/*     require up to (K-2)*L+K+MAX(2N,6M)+4(N+M)+5 locations.  The */
/*     actual number of locations used is computed by HSTCSP and is */
/*     returned in the location W(1). */


/*            * * * * * *   On Output   * * * * * * */

/*   F */
/*     Contains the solution U(I,J) of the finite difference */
/*     approximation for the grid point (THETA(I),R(J)) for */
/*     I=1,2,...,M, J=1,2,...,N. */

/*   PERTRB */
/*     If a combination of periodic, derivative, or unspecified */
/*     boundary conditions is specified for a Poisson equation */
/*     (LAMBDA = 0), a solution may not exist.  PERTRB is a con- */
/*     stant, calculated and subtracted from F, which ensures */
/*     that a solution exists.  HSTCSP then computes this */
/*     solution, which is a least squares solution to the */
/*     original approximation.  This solution plus any constant is also */
/*     a solution; hence, the solution is not unique.  The value of */
/*     PERTRB should be small compared to the right side F. */
/*     Otherwise, a solution is obtained to an essentially different */
/*     problem.  This comparison should always be made to insure that */
/*     a meaningful solution has been obtained. */

/*   IERROR */
/*     An error flag that indicates invalid input parameters. */
/*     Except for numbers 0 and 10, a solution is not attempted. */

/*     =  0  No error */

/*     =  1  A .LT. 0 or B .GT. PI */

/*     =  2  A .GE. B */

/*     =  3  MBDCND .LT. 1 or MBDCND .GT. 9 */

/*     =  4  C .LT. 0 */

/*     =  5  C .GE. D */

/*     =  6  NBDCND .LT. 1 or NBDCND .GT. 6 */

/*     =  7  N .LT. 5 */

/*     =  8  NBDCND = 5 or 6 and MBDCND = 1, 2, 4, 5, or 7 */

/*     =  9  C .GT. 0 and NBDCND .GE. 5 */

/*     = 10  ELMBDA .GT. 0 */

/*     = 11  IDIMF .LT. M */

/*     = 12  M .LT. 5 */

/*     = 13  A = 0 and MBDCND =1,2,3,4,7 or 8 */

/*     = 14  B = PI and MBDCND .LE. 6 */

/*     = 15  A .GT. 0 and MBDCND = 5, 6, or 9 */

/*     = 16  B .LT. PI and MBDCND .GE. 7 */

/*     = 17  LAMBDA .NE. 0 and NBDCND .GE. 5 */

/*     Since this is the only means of indicating a possibly */
/*     incorrect call to HSTCSP, the user should test IERROR after */
/*     the call. */

/*   W */
/*     W(1) contains the required length of W.  Also  W contains */
/*     intermediate values that must not be destroyed if HSTCSP */
/*     will be called again with INTL = 1. */

/* *Long Description: */

/*    * * * * * * *   Program Specifications    * * * * * * * * * * * * */

/*    Dimension of   BDA(N),BDB(N),BDC(M),BDD(M),F(IDIMF,N), */
/*    Arguments      W(See argument list) */

/*    Latest         June 1979 */
/*    Revision */

/*    Subprograms    HSTCSP,HSTCS1,BLKTRI,BLKTR1,INDXA,INDXB,INDXC, */
/*    Required       PROD,PRODP,CPROD,CPRODP,PPADD,PSGF,BSRH,PPSGF, */
/*                   PPSPF,COMPB,TEVLS,R1MACH */

/*    Special        NONE */
/*    Conditions */

/*    Common         CBLKT */
/*    Blocks */

/*    I/O            NONE */

/*    Precision      Single */

/*    Specialist     Roland Sweet */

/*    Language       FORTRAN */

/*    History        Written by Roland Sweet at NCAR in May, 1977 */

/*    Algorithm      This subroutine defines the finite-difference */
/*                   equations, incorporates boundary data, adjusts the */
/*                   right side when the system is singular and calls */
/*                   BLKTRI which solves the linear system of equations. */

/*    Space          5269(decimal) = 12225(octal) locations on the */
/*    Required       NCAR Control Data 7600 */

/*    Timing and        The execution time T on the NCAR Control Data */
/*    Accuracy       7600 for subroutine HSTCSP is roughly proportional */
/*                   to M*N*log2(N), but depends on the input parameter */
/*                   INTL.  Some values are listed in the table below. */
/*                      The solution process employed results in a loss */
/*                   of no more than FOUR significant digits for N and M */
/*                   as large as 64.  More detailed information about */
/*                   accuracy can be found in the documentation for */
/*                   subroutine BLKTRI which is the routine that */
/*                   actually solves the finite difference equations. */


/*                      M(=N)     INTL      MBDCND(=NBDCND)     T(MSECS) */
/*                      -----     ----      ---------------     -------- */

/*                       32        0              1-6             132 */
/*                       32        1              1-6              88 */
/*                       64        0              1-6             546 */
/*                       64        1              1-6             380 */

/*    Portability    American National Standards Institute Fortran. */
/*                   The machine accuracy is set using function R1MACH. */

/*    Required       COS,SIN,ABS,SQRT */
/*    Resident */
/*    Routines */

/*    Reference      Swarztrauber, P.N., 'A Direct Method For The */
/*                   Discrete Solution Of Separable Elliptic Equations,' */
/*                   SIAM J. Numer. Anal. 11(1974), pp. 1136-1150. */

/*    * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/* ***REFERENCES  P. N. Swarztrauber and R. Sweet, Efficient Fortran */
/*                 subprograms for the solution of elliptic equations, */
/*                 NCAR TN/IA-109, July 1975, 138 pp. */
/*               P. N. Swarztrauber, A direct method for the discrete */
/*                 solution of separable elliptic equations, SIAM Journal */
/*                 on Numerical Analysis 11, (1974), pp. 1136-1150. */
/* ***ROUTINES CALLED  HSTCS1, PIMACH */
/* ***REVISION HISTORY  (YYMMDD) */
/*   801001  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890531  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   920501  Reformatted the REFERENCES section.  (WRB) */
/* ***END PROLOGUE  HSTCSP */


/* ***FIRST EXECUTABLE STATEMENT  HSTCSP */
    /* Parameter adjustments */
    --bda;
    --bdb;
    --bdc;
    --bdd;
    f_dim1 = *idimf;
    f_offset = 1 + f_dim1;
    f -= f_offset;
    --w;

    /* Function Body */
    pi = pimach_(&dum);

/*     CHECK FOR INVALID INPUT PARAMETERS */

    *ierror = 0;
    if (*a < 0.f || *b > pi) {
	*ierror = 1;
    }
    if (*a >= *b) {
	*ierror = 2;
    }
    if (*mbdcnd < 1 || *mbdcnd > 9) {
	*ierror = 3;
    }
    if (*c__ < 0.f) {
	*ierror = 4;
    }
    if (*c__ >= *d__) {
	*ierror = 5;
    }
    if (*nbdcnd < 1 || *nbdcnd > 6) {
	*ierror = 6;
    }
    if (*n < 5) {
	*ierror = 7;
    }
    if ((*nbdcnd == 5 || *nbdcnd == 6) && (*mbdcnd == 1 || *mbdcnd == 2 || *
	    mbdcnd == 4 || *mbdcnd == 5 || *mbdcnd == 7)) {
	*ierror = 8;
    }
    if (*c__ > 0.f && *nbdcnd >= 5) {
	*ierror = 9;
    }
    if (*idimf < *m) {
	*ierror = 11;
    }
    if (*m < 5) {
	*ierror = 12;
    }
    if (*a == 0.f && *mbdcnd != 5 && *mbdcnd != 6 && *mbdcnd != 9) {
	*ierror = 13;
    }
    if (*b == pi && *mbdcnd <= 6) {
	*ierror = 14;
    }
    if (*a > 0.f && (*mbdcnd == 5 || *mbdcnd == 6 || *mbdcnd == 9)) {
	*ierror = 15;
    }
    if (*b < pi && *mbdcnd >= 7) {
	*ierror = 16;
    }
    if (*elmbda != 0.f && *nbdcnd >= 5) {
	*ierror = 17;
    }
    if (*ierror != 0) {
	goto L101;
    }
    iwbm = *m + 1;
    iwcm = iwbm + *m;
    iwan = iwcm + *m;
    iwbn = iwan + *n;
    iwcn = iwbn + *n;
    iwsnth = iwcn + *n;
    iwrsq = iwsnth + *m;
    iwwrk = iwrsq + *n;
    ierr1 = 0;
    hstcs1_(intl, a, b, m, mbdcnd, &bda[1], &bdb[1], c__, d__, n, nbdcnd, &
	    bdc[1], &bdd[1], elmbda, &f[f_offset], idimf, pertrb, &ierr1, &w[
	    1], &w[iwbm], &w[iwcm], &w[iwan], &w[iwbn], &w[iwcn], &w[iwsnth], 
	    &w[iwrsq], &w[iwwrk]);
    w[1] = w[iwwrk] + iwwrk - 1;
    *ierror = ierr1;
L101:
    return 0;
} /* hstcsp_ */

