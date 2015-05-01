/* hstcyl.f -- translated by f2c (version 12.02.01).
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

/* DECK HSTCYL */
/* Subroutine */ int hstcyl_(real *a, real *b, integer *m, integer *mbdcnd, 
	real *bda, real *bdb, real *c__, real *d__, integer *n, integer *
	nbdcnd, real *bdc, real *bdd, real *elmbda, real *f, integer *idimf, 
	real *pertrb, integer *ierror, real *w)
{
    /* System generated locals */
    integer f_dim1, f_offset, i__1, i__2;
    real r__1;

    /* Local variables */
    static integer i__, j, k;
    static real a1;
    static integer lp, np, iwb, iwc, iwr, ierr1;
    static real dlrsq, deltar;
    extern /* Subroutine */ int genbun_(integer *, integer *, integer *, 
	    integer *, real *, real *, real *, integer *, real *, integer *, 
	    real *);
    static real deltht, dlthsq;
    extern /* Subroutine */ int poistg_(integer *, integer *, integer *, 
	    integer *, real *, real *, real *, integer *, real *, integer *, 
	    real *);

/* ***BEGIN PROLOGUE  HSTCYL */
/* ***PURPOSE  Solve the standard five-point finite difference */
/*            approximation on a staggered grid to the modified */
/*            Helmholtz equation in cylindrical coordinates. */
/* ***LIBRARY   SLATEC (FISHPACK) */
/* ***CATEGORY  I2B1A1A */
/* ***TYPE      SINGLE PRECISION (HSTCYL-S) */
/* ***KEYWORDS  CYLINDRICAL, ELLIPTIC, FISHPACK, HELMHOLTZ, PDE */
/* ***AUTHOR  Adams, J., (NCAR) */
/*           Swarztrauber, P. N., (NCAR) */
/*           Sweet, R., (NCAR) */
/* ***DESCRIPTION */

/*      HSTCYL solves the standard five-point finite difference */
/*      approximation on a staggered grid to the modified Helmholtz */
/*      equation in cylindrical coordinates */

/*          (1/R)(d/dR)(R(dU/dR)) + (d/dZ)(dU/dZ)C */
/*                      + LAMBDA*(1/R**2)*U = F(R,Z) */

/*      This two-dimensional modified Helmholtz equation results */
/*      from the Fourier transform of a three-dimensional Poisson */
/*      equation. */

/*     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/*     * * * * * * * *    Parameter Description     * * * * * * * * * * */

/*             * * * * * *   On Input    * * * * * * */

/*    A,B */
/*      The range of R, i.e. A .LE. R .LE. B.  A must be less than B and */
/*      A must be non-negative. */

/*    M */
/*      The number of grid points in the interval (A,B).  The grid points */
/*      in the R-direction are given by R(I) = A + (I-0.5)DR for */
/*      I=1,2,...,M where DR =(B-A)/M.  M must be greater than 2. */

/*    MBDCND */
/*      Indicates the type of boundary conditions at R = A and R = B. */

/*      = 1  If the solution is specified at R = A (see note below) and */
/*           R = B. */

/*      = 2  If the solution is specified at R = A (see note below) and */
/*           the derivative of the solution with respect to R is */
/*           specified at R = B. */

/*      = 3  If the derivative of the solution with respect to R is */
/*           specified at R = A (see note below) and R = B. */

/*      = 4  If the derivative of the solution with respect to R is */
/*           specified at R = A (see note below) and the solution is */
/*           specified at R = B. */

/*      = 5  If the solution is unspecified at R = A = 0 and the solution */
/*           is specified at R = B. */

/*      = 6  If the solution is unspecified at R = A = 0 and the */
/*           derivative of the solution with respect to R is specified at */
/*           R = B. */

/*      NOTE:  If A = 0, do not use MBDCND = 1,2,3, or 4, but instead */
/*             use MBDCND = 5 or 6.  The resulting approximation gives */
/*             the only meaningful boundary condition, i.e. dU/dR = 0. */
/*             (see D. Greenspan, 'Introductory Numerical Analysis Of */
/*             Elliptic Boundary Value Problems,' Harper and Row, 1965, */
/*             Chapter 5.) */

/*    BDA */
/*      A one-dimensional array of length N that specifies the boundary */
/*      values (if any) of the solution at R = A.  When MBDCND = 1 or 2, */

/*               BDA(J) = U(A,Z(J)) ,          J=1,2,...,N. */

/*      When MBDCND = 3 or 4, */

/*               BDA(J) = (d/dR)U(A,Z(J)) ,    J=1,2,...,N. */

/*      When MBDCND = 5 or 6, BDA is a dummy variable. */

/*    BDB */
/*      A one-dimensional array of length N that specifies the boundary */
/*      values of the solution at R = B.  When MBDCND = 1,4, or 5, */

/*               BDB(J) = U(B,Z(J)) ,          J=1,2,...,N. */

/*      When MBDCND = 2,3, or 6, */

/*               BDB(J) = (d/dR)U(B,Z(J)) ,    J=1,2,...,N. */

/*    C,D */
/*      The range of Z, i.e. C .LE. Z .LE. D.  C must be less */
/*      than D. */

/*    N */
/*      The number of unknowns in the interval (C,D).  The unknowns in */
/*      the Z-direction are given by Z(J) = C + (J-0.5)DZ, */
/*      J=1,2,...,N, where DZ = (D-C)/N.  N must be greater than 2. */

/*    NBDCND */
/*      Indicates the type of boundary conditions at Z = C */
/*      and Z = D. */

/*      = 0  If the solution is periodic in Z, i.e. */
/*           U(I,J) = U(I,N+J). */

/*      = 1  If the solution is specified at Z = C and Z = D. */

/*      = 2  If the solution is specified at Z = C and the derivative */
/*           of the solution with respect to Z is specified at */
/*           Z = D. */

/*      = 3  If the derivative of the solution with respect to Z is */
/*           specified at Z = C and Z = D. */

/*      = 4  If the derivative of the solution with respect to Z is */
/*           specified at Z = C and the solution is specified at */
/*           Z = D. */

/*    BDC */
/*      A one dimensional array of length M that specifies the boundary */
/*      values of the solution at Z = C.   When NBDCND = 1 or 2, */

/*               BDC(I) = U(R(I),C) ,              I=1,2,...,M. */

/*      When NBDCND = 3 or 4, */

/*               BDC(I) = (d/dZ)U(R(I),C),         I=1,2,...,M. */

/*      When NBDCND = 0, BDC is a dummy variable. */

/*    BDD */
/*      A one-dimensional array of length M that specifies the boundary */
/*      values of the solution at Z = D.  when NBDCND = 1 or 4, */

/*               BDD(I) = U(R(I),D) ,              I=1,2,...,M. */

/*      When NBDCND = 2 or 3, */

/*               BDD(I) = (d/dZ)U(R(I),D) ,        I=1,2,...,M. */

/*      When NBDCND = 0, BDD is a dummy variable. */

/*    ELMBDA */
/*      The constant LAMBDA in the modified Helmholtz equation.  If */
/*      LAMBDA is greater than 0, a solution may not exist.  However, */
/*      HSTCYL will attempt to find a solution.  LAMBDA must be zero */
/*      when MBDCND = 5 or 6. */

/*    F */
/*      A two-dimensional array that specifies the values of the right */
/*      side of the modified Helmholtz equation.  For I=1,2,...,M */
/*      and J=1,2,...,N */

/*               F(I,J) = F(R(I),Z(J)) . */

/*      F must be dimensioned at least M X N. */

/*    IDIMF */
/*      The row (or first) dimension of the array F as it appears in the */
/*      program calling HSTCYL.  This parameter is used to specify the */
/*      variable dimension of F.  IDIMF must be at least M. */

/*    W */
/*      A one-dimensional array that must be provided by the user for */
/*      work space.  W may require up to 13M + 4N + M*INT(log2(N)) */
/*      locations.  The actual number of locations used is computed by */
/*      HSTCYL and is returned in the location W(1). */


/*             * * * * * *   On Output   * * * * * * */

/*    F */
/*      Contains the solution U(I,J) of the finite difference */
/*      approximation for the grid point (R(I),Z(J)) for */
/*      I=1,2,...,M, J=1,2,...,N. */

/*    PERTRB */
/*      If a combination of periodic, derivative, or unspecified */
/*      boundary conditions is specified for a Poisson equation */
/*      (LAMBDA = 0), a solution may not exist.  PERTRB is a con- */
/*      stant, calculated and subtracted from F, which ensures */
/*      that a solution exists.  HSTCYL then computes this */
/*      solution, which is a least squares solution to the */
/*      original approximation.  This solution plus any constant is also */
/*      a solution; hence, the solution is not unique.  The value of */
/*      PERTRB should be small compared to the right side F. */
/*      Otherwise, a solution is obtained to an essentially different */
/*      problem.  This comparison should always be made to insure that */
/*      a meaningful solution has been obtained. */

/*    IERROR */
/*      An error flag that indicates invalid input parameters. */
/*      Except for numbers 0 and 11, a solution is not attempted. */

/*      =  0  No error */

/*      =  1  A .LT. 0 */

/*      =  2  A .GE. B */

/*      =  3  MBDCND .LT. 1 or MBDCND .GT. 6 */

/*      =  4  C .GE. D */

/*      =  5  N .LE. 2 */

/*      =  6  NBDCND .LT. 0 or NBDCND .GT. 4 */

/*      =  7  A = 0 and MBDCND = 1,2,3, or 4 */

/*      =  8  A .GT. 0 and MBDCND .GE. 5 */

/*      =  9  M .LE. 2 */

/*      = 10  IDIMF .LT. M */

/*      = 11  LAMBDA .GT. 0 */

/*      = 12  A=0, MBDCND .GE. 5, ELMBDA .NE. 0 */

/*      Since this is the only means of indicating a possibly */
/*      incorrect call to HSTCYL, the user should test IERROR after */
/*      the call. */

/*    W */
/*      W(1) contains the required length of W. */

/* *Long Description: */

/*     * * * * * * *   Program Specifications    * * * * * * * * * * * * */

/*     Dimension OF   BDA(N),BDB(N),BDC(M),BDD(M),F(IDIMF,N), */
/*     Arguments      W(see argument list) */

/*     Latest         June 1, 1977 */
/*     Revision */

/*     Subprograms    HSTCYL,POISTG,POSTG2,GENBUN,POISD2,POISN2,POISP2, */
/*     Required       COSGEN,MERGE,TRIX,TRI3,PIMACH */

/*     Special        NONE */
/*     Conditions */

/*     Common         NONE */
/*     Blocks */

/*     I/O            NONE */

/*     Precision      Single */

/*     Specialist     Roland Sweet */

/*     Language       FORTRAN */

/*     History        Written by Roland Sweet at NCAR in March, 1977 */

/*     Algorithm      This subroutine defines the finite-difference */
/*                    equations, incorporates boundary data, adjusts the */
/*                    right side when the system is singular and calls */
/*                    either POISTG or GENBUN which solves the linear */
/*                    system of equations. */

/*     Space          8228(decimal) = 20044(octal) locations on the */
/*     Required       NCAR Control Data 7600 */

/*     Timing and        The execution time T on the NCAR Control Data */
/*     Accuracy       7600 for subroutine HSTCYL is roughly proportional */
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

/*                        32       1-6       1-4         56 */
/*                        64       1-6       1-4        230 */

/*     Portability    American National Standards Institute Fortran. */
/*                    The machine dependent constant PI is defined in */
/*                    function PIMACH. */

/*     Required       COS */
/*     Resident */
/*     Routines */

/*     Reference      Schumann, U. and R. Sweet,'A Direct Method For */
/*                    The Solution of Poisson's Equation With Neumann */
/*                    Boundary Conditions On A Staggered Grid Of */
/*                    Arbitrary Size,' J. Comp. Phys. 20(1976), */
/*                    pp. 171-182. */

/*     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/* ***REFERENCES  U. Schumann and R. Sweet, A direct method for the */
/*                 solution of Poisson's equation with Neumann boundary */
/*                 conditions on a staggered grid of arbitrary size, */
/*                 Journal of Computational Physics 20, (1976), */
/*                 pp. 171-182. */
/* ***ROUTINES CALLED  GENBUN, POISTG */
/* ***REVISION HISTORY  (YYMMDD) */
/*   801001  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890531  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   920501  Reformatted the REFERENCES section.  (WRB) */
/* ***END PROLOGUE  HSTCYL */


/* ***FIRST EXECUTABLE STATEMENT  HSTCYL */
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
    if (*a < 0.f) {
	*ierror = 1;
    }
    if (*a >= *b) {
	*ierror = 2;
    }
    if (*mbdcnd <= 0 || *mbdcnd >= 7) {
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
    if (*a == 0.f && *mbdcnd != 5 && *mbdcnd != 6) {
	*ierror = 7;
    }
    if (*a > 0.f && *mbdcnd >= 5) {
	*ierror = 8;
    }
    if (*idimf < *m) {
	*ierror = 10;
    }
    if (*m <= 2) {
	*ierror = 9;
    }
    if (*a == 0.f && *mbdcnd >= 5 && *elmbda != 0.f) {
	*ierror = 12;
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

/*     DEFINE A,B,C COEFFICIENTS IN W-ARRAY. */

    iwb = *m;
    iwc = iwb + *m;
    iwr = iwc + *m;
    i__1 = *m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	j = iwr + i__;
	w[j] = *a + (i__ - .5f) * deltar;
	w[i__] = (*a + (i__ - 1) * deltar) / (dlrsq * w[j]);
	k = iwc + i__;
	w[k] = (*a + i__ * deltar) / (dlrsq * w[j]);
	k = iwb + i__;
/* Computing 2nd power */
	r__1 = w[j];
	w[k] = *elmbda / (r__1 * r__1) - 2.f / dlrsq;
/* L101: */
    }

/*     ENTER BOUNDARY DATA FOR R-BOUNDARIES. */

    switch (*mbdcnd) {
	case 1:  goto L102;
	case 2:  goto L102;
	case 3:  goto L104;
	case 4:  goto L104;
	case 5:  goto L106;
	case 6:  goto L106;
    }
L102:
    a1 = w[1] * 2.f;
    w[iwb + 1] -= w[1];
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	f[j * f_dim1 + 1] -= a1 * bda[j];
/* L103: */
    }
    goto L106;
L104:
    a1 = deltar * w[1];
    w[iwb + 1] += w[1];
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	f[j * f_dim1 + 1] += a1 * bda[j];
/* L105: */
    }
L106:
    switch (*mbdcnd) {
	case 1:  goto L107;
	case 2:  goto L109;
	case 3:  goto L109;
	case 4:  goto L107;
	case 5:  goto L107;
	case 6:  goto L109;
    }
L107:
    w[iwc] -= w[iwr];
    a1 = w[iwr] * 2.f;
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	f[*m + j * f_dim1] -= a1 * bdb[j];
/* L108: */
    }
    goto L111;
L109:
    w[iwc] += w[iwr];
    a1 = deltar * w[iwr];
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	f[*m + j * f_dim1] -= a1 * bdb[j];
/* L110: */
    }

/*     ENTER BOUNDARY DATA FOR THETA-BOUNDARIES. */

L111:
    a1 = 2.f / dlthsq;
    switch (np) {
	case 1:  goto L121;
	case 2:  goto L112;
	case 3:  goto L112;
	case 4:  goto L114;
	case 5:  goto L114;
    }
L112:
    i__1 = *m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	f[i__ + f_dim1] -= a1 * bdc[i__];
/* L113: */
    }
    goto L116;
L114:
    a1 = 1.f / deltht;
    i__1 = *m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	f[i__ + f_dim1] += a1 * bdc[i__];
/* L115: */
    }
L116:
    a1 = 2.f / dlthsq;
    switch (np) {
	case 1:  goto L121;
	case 2:  goto L117;
	case 3:  goto L119;
	case 4:  goto L119;
	case 5:  goto L117;
    }
L117:
    i__1 = *m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	f[i__ + *n * f_dim1] -= a1 * bdd[i__];
/* L118: */
    }
    goto L121;
L119:
    a1 = 1.f / deltht;
    i__1 = *m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	f[i__ + *n * f_dim1] -= a1 * bdd[i__];
/* L120: */
    }
L121:

/*     ADJUST RIGHT SIDE OF SINGULAR PROBLEMS TO INSURE EXISTENCE OF A */
/*     SOLUTION. */

    *pertrb = 0.f;
    if (*elmbda < 0.f) {
	goto L130;
    } else if (*elmbda == 0) {
	goto L123;
    } else {
	goto L122;
    }
L122:
    *ierror = 11;
    goto L130;
L123:
    switch (*mbdcnd) {
	case 1:  goto L130;
	case 2:  goto L130;
	case 3:  goto L124;
	case 4:  goto L130;
	case 5:  goto L130;
	case 6:  goto L124;
    }
L124:
    switch (np) {
	case 1:  goto L125;
	case 2:  goto L130;
	case 3:  goto L130;
	case 4:  goto L125;
	case 5:  goto L130;
    }
L125:
    i__1 = *m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	a1 = 0.f;
	i__2 = *n;
	for (j = 1; j <= i__2; ++j) {
	    a1 += f[i__ + j * f_dim1];
/* L126: */
	}
	j = iwr + i__;
	*pertrb += a1 * w[j];
/* L127: */
    }
    *pertrb /= *m * *n * .5f * (*a + *b);
    i__1 = *m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = *n;
	for (j = 1; j <= i__2; ++j) {
	    f[i__ + j * f_dim1] -= *pertrb;
/* L128: */
	}
/* L129: */
    }
L130:

/*     MULTIPLY I-TH EQUATION THROUGH BY  DELTHT**2 */

    i__1 = *m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	w[i__] *= dlthsq;
	j = iwc + i__;
	w[j] *= dlthsq;
	j = iwb + i__;
	w[j] *= dlthsq;
	i__2 = *n;
	for (j = 1; j <= i__2; ++j) {
	    f[i__ + j * f_dim1] *= dlthsq;
/* L131: */
	}
/* L132: */
    }
    lp = *nbdcnd;
    w[1] = 0.f;
    w[iwr] = 0.f;

/*     CALL GENBUN TO SOLVE THE SYSTEM OF EQUATIONS. */

    if (*nbdcnd == 0) {
	goto L133;
    }
    poistg_(&lp, n, &c__1, m, &w[1], &w[iwb + 1], &w[iwc + 1], idimf, &f[
	    f_offset], &ierr1, &w[iwr + 1]);
    goto L134;
L133:
    genbun_(&lp, n, &c__1, m, &w[1], &w[iwb + 1], &w[iwc + 1], idimf, &f[
	    f_offset], &ierr1, &w[iwr + 1]);
L134:
    w[1] = w[iwr + 1] + *m * 3;
    return 0;
} /* hstcyl_ */

