/* hstssp.f -- translated by f2c (version 12.02.01).
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

/* DECK HSTSSP */
/* Subroutine */ int hstssp_(real *a, real *b, integer *m, integer *mbdcnd, 
	real *bda, real *bdb, real *c__, real *d__, integer *n, integer *
	nbdcnd, real *bdc, real *bdd, real *elmbda, real *f, integer *idimf, 
	real *pertrb, integer *ierror, real *w)
{
    /* System generated locals */
    integer f_dim1, f_offset, i__1, i__2;
    real r__1;

    /* Local variables */
    static integer i__, j, k;
    static real a1, a2, a3;
    static integer mb;
    static real pi;
    static integer lp, np, mm1, iwb, iwc;
    static real dum;
    static integer iwr, isw, jsw, iws, ierr1;
    static real dlrsq;
    extern doublereal pimach_(real *);
    static real deltar;
    extern /* Subroutine */ int genbun_(integer *, integer *, integer *, 
	    integer *, real *, real *, real *, integer *, real *, integer *, 
	    real *);
    static real deltht, dlthsq;
    extern /* Subroutine */ int poistg_(integer *, integer *, integer *, 
	    integer *, real *, real *, real *, integer *, real *, integer *, 
	    real *);

/* ***BEGIN PROLOGUE  HSTSSP */
/* ***PURPOSE  Solve the standard five-point finite difference */
/*            approximation on a staggered grid to the Helmholtz */
/*            equation in spherical coordinates and on the surface of */
/*            the unit sphere (radius of 1). */
/* ***LIBRARY   SLATEC (FISHPACK) */
/* ***CATEGORY  I2B1A1A */
/* ***TYPE      SINGLE PRECISION (HSTSSP-S) */
/* ***KEYWORDS  ELLIPTIC, FISHPACK, HELMHOLTZ, PDE, SPHERICAL */
/* ***AUTHOR  Adams, J., (NCAR) */
/*           Swarztrauber, P. N., (NCAR) */
/*           Sweet, R., (NCAR) */
/* ***DESCRIPTION */

/*     HSTSSP solves the standard five-point finite difference */
/*     approximation on a staggered grid to the Helmholtz equation in */
/*     spherical coordinates and on the surface of the unit sphere */
/*     (radius of 1) */

/*             (1/SIN(THETA))(d/dTHETA)(SIN(THETA)(dU/dTHETA)) + */

/*       (1/SIN(THETA)**2)(d/dPHI)(dU/dPHI) + LAMBDA*U = F(THETA,PHI) */

/*     where THETA is colatitude and PHI is longitude. */

/*    * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/*    * * * * * * * *    Parameter Description     * * * * * * * * * * */

/*            * * * * * *   On Input    * * * * * * */

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
/*     for I=1,2,...,M where DTHETA =(B-A)/M.  M must be greater than 2. */

/*   MBDCND */
/*     Indicates the type of boundary conditions at THETA = A and */
/*     THETA = B. */

/*     = 1  If the solution is specified at THETA = A and THETA = B. */
/*          (see note 3 below) */

/*     = 2  If the solution is specified at THETA = A and the derivative */
/*          of the solution with respect to THETA is specified at */
/*          THETA = B (see notes 2 and 3 below). */

/*     = 3  If the derivative of the solution with respect to THETA is */
/*          specified at THETA = A (see notes 1, 2 below) and THETA = B. */

/*     = 4  If the derivative of the solution with respect to THETA is */
/*          specified at THETA = A (see notes 1 and 2 below) and the */
/*          solution is specified at THETA = B. */

/*     = 5  If the solution is unspecified at THETA = A = 0 and the */
/*          solution is specified at THETA = B.  (see note 3 below) */

/*     = 6  If the solution is unspecified at THETA = A = 0 and the */
/*          derivative of the solution with respect to THETA is */
/*          specified at THETA = B (see note 2 below). */

/*     = 7  If the solution is specified at THETA = A and the */
/*          solution is unspecified at THETA = B = PI. (see note 3 below) */

/*     = 8  If the derivative of the solution with respect to */
/*          THETA is specified at THETA = A (see note 1 below) */
/*          and the solution is unspecified at THETA = B = PI. */

/*     = 9  If the solution is unspecified at THETA = A = 0 and */
/*          THETA = B = PI. */

/*     NOTES:  1.  If A = 0, do not use MBDCND = 3, 4, or 8, */
/*                 but instead use MBDCND = 5, 6, or 9. */

/*             2.  If B = PI, do not use MBDCND = 2, 3, or 6, */
/*                 but instead use MBDCND = 7, 8, or 9. */

/*             3.  When the solution is specified at THETA = 0 and/or */
/*                 THETA = PI and the other boundary conditions are */
/*                 combinations of unspecified, normal derivative, or */
/*                 periodicity a singular system results.  The unique */
/*                 solution is determined by extrapolation to the */
/*                 specification of the solution at either THETA = 0 or */
/*                 THETA = PI.  But in these cases the right side of the */
/*                 system will be perturbed by the constant PERTRB. */

/*   BDA */
/*     A one-dimensional array of length N that specifies the boundary */
/*     values (if any) of the solution at THETA = A.  When */
/*     MBDCND = 1, 2, or 7, */

/*              BDA(J) = U(A,PHI(J)) ,              J=1,2,...,N. */

/*     When MBDCND = 3, 4, or 8, */

/*              BDA(J) = (d/dTHETA)U(A,PHI(J)) ,    J=1,2,...,N. */

/*     When MBDCND has any other value, BDA is a dummy variable. */

/*   BDB */
/*     A one-dimensional array of length N that specifies the boundary */
/*     values of the solution at THETA = B.  When MBDCND = 1,4, or 5, */

/*              BDB(J) = U(B,PHI(J)) ,              J=1,2,...,N. */

/*     When MBDCND = 2,3, or 6, */

/*              BDB(J) = (d/dTHETA)U(B,PHI(J)) ,    J=1,2,...,N. */

/*     When MBDCND has any other value, BDB is a dummy variable. */

/*   C,D */
/*     The range of PHI (longitude), i.e. C .LE. PHI .LE. D. */
/*     C must be less than D.  If D-C = 2*PI, periodic boundary */
/*     conditions are usually prescribed. */

/*   N */
/*     The number of unknowns in the interval (C,D).  The unknowns in */
/*     the PHI-direction are given by PHI(J) = C + (J-0.5)DPHI, */
/*     J=1,2,...,N, where DPHI = (D-C)/N.  N must be greater than 2. */

/*   NBDCND */
/*     Indicates the type of boundary conditions at PHI = C */
/*     and PHI = D. */

/*     = 0  If the solution is periodic in PHI, i.e. */
/*          U(I,J) = U(I,N+J). */

/*     = 1  If the solution is specified at PHI = C and PHI = D */
/*          (see note below). */

/*     = 2  If the solution is specified at PHI = C and the derivative */
/*          of the solution with respect to PHI is specified at */
/*          PHI = D (see note below). */

/*     = 3  If the derivative of the solution with respect to PHI is */
/*          specified at PHI = C and PHI = D. */

/*     = 4  If the derivative of the solution with respect to PHI is */
/*          specified at PHI = C and the solution is specified at */
/*          PHI = D (see note below). */

/*     NOTE:  When NBDCND = 1, 2, or 4, do not use MBDCND = 5, 6, 7, 8, */
/*     or 9 (the former indicates that the solution is specified at */
/*     a pole; the latter indicates the solution is unspecified).  Use */
/*     instead MBDCND = 1 or 2. */

/*   BDC */
/*     A one dimensional array of length M that specifies the boundary */
/*     values of the solution at PHI = C.   When NBDCND = 1 or 2, */

/*              BDC(I) = U(THETA(I),C) ,              I=1,2,...,M. */

/*     When NBDCND = 3 or 4, */

/*              BDC(I) = (d/dPHI)U(THETA(I),C),       I=1,2,...,M. */

/*     When NBDCND = 0, BDC is a dummy variable. */

/*   BDD */
/*     A one-dimensional array of length M that specifies the boundary */
/*     values of the solution at PHI = D.  When NBDCND = 1 or 4, */

/*              BDD(I) = U(THETA(I),D) ,              I=1,2,...,M. */

/*     When NBDCND = 2 or 3, */

/*              BDD(I) = (d/dPHI)U(THETA(I),D) ,      I=1,2,...,M. */

/*     When NBDCND = 0, BDD is a dummy variable. */

/*   ELMBDA */
/*     The constant LAMBDA in the Helmholtz equation.  If LAMBDA is */
/*     greater than 0, a solution may not exist.  However, HSTSSP will */
/*     attempt to find a solution. */

/*   F */
/*     A two-dimensional array that specifies the values of the right */
/*     side of the Helmholtz equation.  For I=1,2,...,M and J=1,2,...,N */

/*              F(I,J) = F(THETA(I),PHI(J)) . */

/*     F must be dimensioned at least M X N. */

/*   IDIMF */
/*     The row (or first) dimension of the array F as it appears in the */
/*     program calling HSTSSP.  This parameter is used to specify the */
/*     variable dimension of F.  IDIMF must be at least M. */

/*   W */
/*     A one-dimensional array that must be provided by the user for */
/*     work space.  W may require up to 13M + 4N + M*INT(log2(N)) */
/*     locations.  The actual number of locations used is computed by */
/*     HSTSSP and is returned in the location W(1). */


/*            * * * * * *   On Output   * * * * * * */

/*   F */
/*     Contains the solution U(I,J) of the finite difference */
/*     approximation for the grid point (THETA(I),PHI(J)) for */
/*     I=1,2,...,M, J=1,2,...,N. */

/*   PERTRB */
/*     If a combination of periodic, derivative, or unspecified */
/*     boundary conditions is specified for a Poisson equation */
/*     (LAMBDA = 0), a solution may not exist.  PERTRB is a con- */
/*     stant, calculated and subtracted from F, which ensures */
/*     that a solution exists.  HSTSSP then computes this */
/*     solution, which is a least squares solution to the */
/*     original approximation.  This solution plus any constant is also */
/*     a solution; hence, the solution is not unique.  The value of */
/*     PERTRB should be small compared to the right side F. */
/*     Otherwise, a solution is obtained to an essentially different */
/*     problem.  This comparison should always be made to insure that */
/*     a meaningful solution has been obtained. */

/*   IERROR */
/*     An error flag that indicates invalid input parameters. */
/*      Except for numbers 0 and 14, a solution is not attempted. */

/*     =  0  No error */

/*     =  1  A .LT. 0 or B .GT. PI */

/*     =  2  A .GE. B */

/*     =  3  MBDCND .LT. 1 or MBDCND .GT. 9 */

/*     =  4  C .GE. D */

/*     =  5  N .LE. 2 */

/*     =  6  NBDCND .LT. 0 or NBDCND .GT. 4 */

/*     =  7  A .GT. 0 and MBDCND = 5, 6, or 9 */

/*     =  8  A = 0 and MBDCND = 3, 4, or 8 */

/*     =  9  B .LT. PI and MBDCND .GE. 7 */

/*     = 10  B = PI and MBDCND = 2,3, or 6 */

/*     = 11  MBDCND .GE. 5 and NDBCND = 1, 2, or 4 */

/*     = 12  IDIMF .LT. M */

/*     = 13  M .LE. 2 */

/*     = 14  LAMBDA .GT. 0 */

/*     Since this is the only means of indicating a possibly */
/*     incorrect call to HSTSSP, the user should test IERROR after */
/*     the call. */

/*   W */
/*     W(1) contains the required length of W. */

/* *Long Description: */

/*    * * * * * * *   Program Specifications    * * * * * * * * * * * * */

/*    Dimension of   BDA(N),BDB(N),BDC(M),BDD(M),F(IDIMF,N), */
/*    Arguments      W(see argument list) */

/*    Latest         June 1, 1977 */
/*    Revision */

/*    Subprograms    HSTSSP,POISTG,POSTG2,GENBUN,POISD2,POISN2,POISP2, */
/*    Required       COSGEN,MERGE,TRIX,TRI3,PIMACH */

/*    Special        NONE */
/*    Conditions */

/*    Common         NONE */
/*    Blocks */

/*    I/O            NONE */

/*    Precision      Single */

/*    Specialist     Roland Sweet */

/*    Language       FORTRAN */

/*    History        Written by Roland Sweet at NCAR in April, 1977 */

/*    Algorithm      This subroutine defines the finite-difference */
/*                   equations, incorporates boundary data, adjusts the */
/*                   right side when the system is singular and calls */
/*                   either POISTG or GENBUN which solves the linear */
/*                   system of equations. */

/*    Space          8427(decimal) = 20353(octal) locations on the */
/*    Required       NCAR Control Data 7600 */

/*     Timing and        The execution time T on the NCAR Control Data */
/*     Accuracy       7600 for subroutine HSTSSP is roughly proportional */
/*                    to M*N*log2(N).  Some typical values are listed in */
/*                    the table below. */
/*                       The solution process employed results in a loss */
/*                    of no more than four significant digits for N and M */
/*                    as large as 64.  More detailed information about */
/*                    accuracy can be found in the documentation for */
/*                    subroutine POISTG which is the routine that */
/*                    actually solves the finite difference equations. */


/*                       M(=N)    MBDCND    NBDCND    T(MSECS) */
/*                       -----    ------    ------    -------- */

/*                        32       1-9       1-4         56 */
/*                        64       1-9       1-4        230 */

/*    Portability     American National Standards Institute FORTRAN. */
/*                    The machine dependent constant PI is defined in */
/*                    function PIMACH. */

/*    Required       COS */
/*    Resident */
/*    Routines */

/*    Reference      Schumann, U. and R. Sweet,'A Direct Method For */
/*                   The Solution Of Poisson's Equation With Neumann */
/*                   Boundary Conditions On A Staggered Grid Of */
/*                   Arbitrary Size,' J. Comp. Phys. 20(1976), */
/*                   pp. 171-182. */

/*    * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/* ***REFERENCES  U. Schumann and R. Sweet, A direct method for the */
/*                 solution of Poisson's equation with Neumann boundary */
/*                 conditions on a staggered grid of arbitrary size, */
/*                 Journal of Computational Physics 20, (1976), */
/*                 pp. 171-182. */
/* ***ROUTINES CALLED  GENBUN, PIMACH, POISTG */
/* ***REVISION HISTORY  (YYMMDD) */
/*   801001  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890531  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   920501  Reformatted the REFERENCES section.  (WRB) */
/* ***END PROLOGUE  HSTSSP */


/* ***FIRST EXECUTABLE STATEMENT  HSTSSP */
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
    *ierror = 0;
    pi = pimach_(&dum);
    if (*a < 0.f || *b > pi) {
	*ierror = 1;
    }
    if (*a >= *b) {
	*ierror = 2;
    }
    if (*mbdcnd <= 0 || *mbdcnd > 9) {
	*ierror = 3;
    }
    if (*c__ >= *d__) {
	*ierror = 4;
    }
    if (*n <= 2) {
	*ierror = 5;
    }
    if (*nbdcnd < 0 || *nbdcnd >= 5) {
	*ierror = 6;
    }
    if (*a > 0.f && (*mbdcnd == 5 || *mbdcnd == 6 || *mbdcnd == 9)) {
	*ierror = 7;
    }
    if (*a == 0.f && (*mbdcnd == 3 || *mbdcnd == 4 || *mbdcnd == 8)) {
	*ierror = 8;
    }
    if (*b < pi && *mbdcnd >= 7) {
	*ierror = 9;
    }
    if (*b == pi && (*mbdcnd == 2 || *mbdcnd == 3 || *mbdcnd == 6)) {
	*ierror = 10;
    }
    if (*mbdcnd >= 5 && (*nbdcnd == 1 || *nbdcnd == 2 || *nbdcnd == 4)) {
	*ierror = 11;
    }
    if (*idimf < *m) {
	*ierror = 12;
    }
    if (*m <= 2) {
	*ierror = 13;
    }
    if (*ierror != 0) {
	return 0;
    }
    deltar = (*b - *a) / *m;
/* Computing 2nd power */
    r__1 = deltar;
    dlrsq = r__1 * r__1;
    deltht = (*d__ - *c__) / *n;
/* Computing 2nd power */
    r__1 = deltht;
    dlthsq = r__1 * r__1;
    np = *nbdcnd + 1;
    isw = 1;
    jsw = 1;
    mb = *mbdcnd;
    if (*elmbda != 0.f) {
	goto L105;
    }
    switch (*mbdcnd) {
	case 1:  goto L101;
	case 2:  goto L102;
	case 3:  goto L105;
	case 4:  goto L103;
	case 5:  goto L101;
	case 6:  goto L105;
	case 7:  goto L101;
	case 8:  goto L105;
	case 9:  goto L105;
    }
L101:
    if (*a != 0.f || *b != pi) {
	goto L105;
    }
    mb = 9;
    goto L104;
L102:
    if (*a != 0.f) {
	goto L105;
    }
    mb = 6;
    goto L104;
L103:
    if (*b != pi) {
	goto L105;
    }
    mb = 8;
L104:
    jsw = 2;
L105:

/*     DEFINE A,B,C COEFFICIENTS IN W-ARRAY. */

    iwb = *m;
    iwc = iwb + *m;
    iwr = iwc + *m;
    iws = iwr + *m;
    i__1 = *m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	j = iwr + i__;
	w[j] = sin(*a + (i__ - .5f) * deltar);
	w[i__] = sin(*a + (i__ - 1) * deltar) / dlrsq;
/* L106: */
    }
    mm1 = *m - 1;
    i__1 = mm1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	k = iwc + i__;
	w[k] = w[i__ + 1];
	j = iwr + i__;
	k = iwb + i__;
	w[k] = *elmbda * w[j] - (w[i__] + w[i__ + 1]);
/* L107: */
    }
    w[iwr] = sin(*b) / dlrsq;
    w[iwc] = *elmbda * w[iws] - (w[*m] + w[iwr]);
    i__1 = *m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	j = iwr + i__;
	a1 = w[j];
	i__2 = *n;
	for (j = 1; j <= i__2; ++j) {
	    f[i__ + j * f_dim1] = a1 * f[i__ + j * f_dim1];
/* L108: */
	}
/* L109: */
    }

/*     ENTER BOUNDARY DATA FOR THETA-BOUNDARIES. */

    switch (mb) {
	case 1:  goto L110;
	case 2:  goto L110;
	case 3:  goto L112;
	case 4:  goto L112;
	case 5:  goto L114;
	case 6:  goto L114;
	case 7:  goto L110;
	case 8:  goto L112;
	case 9:  goto L114;
    }
L110:
    a1 = w[1] * 2.f;
    w[iwb + 1] -= w[1];
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	f[j * f_dim1 + 1] -= a1 * bda[j];
/* L111: */
    }
    goto L114;
L112:
    a1 = deltar * w[1];
    w[iwb + 1] += w[1];
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	f[j * f_dim1 + 1] += a1 * bda[j];
/* L113: */
    }
L114:
    switch (mb) {
	case 1:  goto L115;
	case 2:  goto L117;
	case 3:  goto L117;
	case 4:  goto L115;
	case 5:  goto L115;
	case 6:  goto L117;
	case 7:  goto L119;
	case 8:  goto L119;
	case 9:  goto L119;
    }
L115:
    a1 = w[iwr] * 2.f;
    w[iwc] -= w[iwr];
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	f[*m + j * f_dim1] -= a1 * bdb[j];
/* L116: */
    }
    goto L119;
L117:
    a1 = deltar * w[iwr];
    w[iwc] += w[iwr];
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	f[*m + j * f_dim1] -= a1 * bdb[j];
/* L118: */
    }

/*     ENTER BOUNDARY DATA FOR PHI-BOUNDARIES. */

L119:
    a1 = 2.f / dlthsq;
    switch (np) {
	case 1:  goto L129;
	case 2:  goto L120;
	case 3:  goto L120;
	case 4:  goto L122;
	case 5:  goto L122;
    }
L120:
    i__1 = *m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	j = iwr + i__;
	f[i__ + f_dim1] -= a1 * bdc[i__] / w[j];
/* L121: */
    }
    goto L124;
L122:
    a1 = 1.f / deltht;
    i__1 = *m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	j = iwr + i__;
	f[i__ + f_dim1] += a1 * bdc[i__] / w[j];
/* L123: */
    }
L124:
    a1 = 2.f / dlthsq;
    switch (np) {
	case 1:  goto L129;
	case 2:  goto L125;
	case 3:  goto L127;
	case 4:  goto L127;
	case 5:  goto L125;
    }
L125:
    i__1 = *m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	j = iwr + i__;
	f[i__ + *n * f_dim1] -= a1 * bdd[i__] / w[j];
/* L126: */
    }
    goto L129;
L127:
    a1 = 1.f / deltht;
    i__1 = *m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	j = iwr + i__;
	f[i__ + *n * f_dim1] -= a1 * bdd[i__] / w[j];
/* L128: */
    }
L129:

/*     ADJUST RIGHT SIDE OF SINGULAR PROBLEMS TO INSURE EXISTENCE OF A */
/*     SOLUTION. */

    *pertrb = 0.f;
    if (*elmbda < 0.f) {
	goto L139;
    } else if (*elmbda == 0) {
	goto L131;
    } else {
	goto L130;
    }
L130:
    *ierror = 14;
    goto L139;
L131:
    switch (mb) {
	case 1:  goto L139;
	case 2:  goto L139;
	case 3:  goto L132;
	case 4:  goto L139;
	case 5:  goto L139;
	case 6:  goto L132;
	case 7:  goto L139;
	case 8:  goto L132;
	case 9:  goto L132;
    }
L132:
    switch (np) {
	case 1:  goto L133;
	case 2:  goto L139;
	case 3:  goto L139;
	case 4:  goto L133;
	case 5:  goto L139;
    }
L133:
    isw = 2;
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	i__2 = *m;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    *pertrb += f[i__ + j * f_dim1];
/* L134: */
	}
/* L135: */
    }
    a1 = *n * (cos(*a) - cos(*b)) / (sin(deltar * .5f) * 2.f);
    *pertrb /= a1;
    i__1 = *m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	j = iwr + i__;
	a1 = *pertrb * w[j];
	i__2 = *n;
	for (j = 1; j <= i__2; ++j) {
	    f[i__ + j * f_dim1] -= a1;
/* L136: */
	}
/* L137: */
    }
    a2 = 0.f;
    a3 = 0.f;
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	a2 += f[j * f_dim1 + 1];
	a3 += f[*m + j * f_dim1];
/* L138: */
    }
    a2 /= w[iwr + 1];
    a3 /= w[iws];
L139:

/*     MULTIPLY I-TH EQUATION THROUGH BY  R(I)*DELTHT**2 */

    i__1 = *m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	j = iwr + i__;
	a1 = dlthsq * w[j];
	w[i__] = a1 * w[i__];
	j = iwc + i__;
	w[j] = a1 * w[j];
	j = iwb + i__;
	w[j] = a1 * w[j];
	i__2 = *n;
	for (j = 1; j <= i__2; ++j) {
	    f[i__ + j * f_dim1] = a1 * f[i__ + j * f_dim1];
/* L140: */
	}
/* L141: */
    }
    lp = *nbdcnd;
    w[1] = 0.f;
    w[iwr] = 0.f;

/*     CALL POISTG OR GENBUN TO SOLVE THE SYSTEM OF EQUATIONS. */

    if (*nbdcnd == 0) {
	goto L142;
    }
    poistg_(&lp, n, &c__1, m, &w[1], &w[iwb + 1], &w[iwc + 1], idimf, &f[
	    f_offset], &ierr1, &w[iwr + 1]);
    goto L143;
L142:
    genbun_(&lp, n, &c__1, m, &w[1], &w[iwb + 1], &w[iwc + 1], idimf, &f[
	    f_offset], &ierr1, &w[iwr + 1]);
L143:
    w[1] = w[iwr + 1] + *m * 3;
    if (isw != 2 || jsw != 2) {
	goto L150;
    }
    if (mb != 8) {
	goto L145;
    }
    a1 = 0.f;
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	a1 += f[*m + j * f_dim1];
/* L144: */
    }
    a1 = (a1 - dlrsq * a3 / 16.f) / *n;
    if (*nbdcnd == 3) {
	a1 += (bdd[*m] - bdc[*m]) / (*d__ - *c__);
    }
    a1 = bdb[1] - a1;
    goto L147;
L145:
    a1 = 0.f;
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	a1 += f[j * f_dim1 + 1];
/* L146: */
    }
    a1 = (a1 - dlrsq * a2 / 16.f) / *n;
    if (*nbdcnd == 3) {
	a1 += (bdd[1] - bdc[1]) / (*d__ - *c__);
    }
    a1 = bda[1] - a1;
L147:
    i__1 = *m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = *n;
	for (j = 1; j <= i__2; ++j) {
	    f[i__ + j * f_dim1] += a1;
/* L148: */
	}
/* L149: */
    }
L150:
    return 0;
} /* hstssp_ */

