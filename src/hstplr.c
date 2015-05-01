/* hstplr.f -- translated by f2c (version 12.02.01).
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

/* DECK HSTPLR */
/* Subroutine */ int hstplr_(real *a, real *b, integer *m, integer *mbdcnd, 
	real *bda, real *bdb, real *c__, real *d__, integer *n, integer *
	nbdcnd, real *bdc, real *bdd, real *elmbda, real *f, integer *idimf, 
	real *pertrb, integer *ierror, real *w)
{
    /* System generated locals */
    integer f_dim1, f_offset, i__1, i__2;
    real r__1;

    /* Local variables */
    static integer i__, j, k;
    static real a1, a2;
    static integer mb, lp, np, iwb, iwc, iwr, isw, ierr1;
    static real dlrsq, deltar;
    extern /* Subroutine */ int genbun_(integer *, integer *, integer *, 
	    integer *, real *, real *, real *, integer *, real *, integer *, 
	    real *);
    static real deltht, dlthsq;
    extern /* Subroutine */ int poistg_(integer *, integer *, integer *, 
	    integer *, real *, real *, real *, integer *, real *, integer *, 
	    real *);

/* ***BEGIN PROLOGUE  HSTPLR */
/* ***PURPOSE  Solve the standard five-point finite difference */
/*            approximation on a staggered grid to the Helmholtz equation */
/*            in polar coordinates. */
/* ***LIBRARY   SLATEC (FISHPACK) */
/* ***CATEGORY  I2B1A1A */
/* ***TYPE      SINGLE PRECISION (HSTPLR-S) */
/* ***KEYWORDS  ELLIPTIC, FISHPACK, HELMHOLTZ, PDE, POLAR */
/* ***AUTHOR  Adams, J., (NCAR) */
/*           Swarztrauber, P. N., (NCAR) */
/*           Sweet, R., (NCAR) */
/* ***DESCRIPTION */

/*      HSTPLR solves the standard five-point finite difference */
/*      approximation on a staggered grid to the Helmholtz equation in */
/*      polar coordinates */

/*      (1/R)(d/DR)(R(dU/DR)) + (1/R**2)(d/dTHETA)(dU/dTHETA) */

/*                      + LAMBDA*U = F(R,THETA) */

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

/*      = 1  If the solution is specified at R = A and R = B. */

/*      = 2  If the solution is specified at R = A and the derivative */
/*           of the solution with respect to R is specified at R = B. */
/*           (see note 1 below) */

/*      = 3  If the derivative of the solution with respect to R is */
/*           specified at R = A (see note 2 below) and R = B. */

/*      = 4  If the derivative of the solution with respect to R is */
/*           specified at R = A (see note 2 below) and the solution is */
/*           specified at R = B. */

/*      = 5  If the solution is unspecified at R = A = 0 and the solution */
/*           is specified at R = B. */

/*      = 6  If the solution is unspecified at R = A = 0 and the */
/*           derivative of the solution with respect to R is specified at */
/*           R = B. */

/*      NOTE 1:  If A = 0, MBDCND = 2, and NBDCND = 0 or 3, the system of */
/*               equations to be solved is singular.  The unique solution */
/*               is determined by extrapolation to the specification of */
/*               U(0,THETA(1)).  But in this case the right side of the */
/*               system will be perturbed by the constant PERTRB. */

/*      NOTE 2:  If A = 0, do not use MBDCND = 3 or 4, but instead use */
/*               MBDCND = 1,2,5, or 6. */

/*    BDA */
/*      A one-dimensional array of length N that specifies the boundary */
/*      values (if any) of the solution at R = A.  When MBDCND = 1 or 2, */

/*               BDA(J) = U(A,THETA(J)) ,          J=1,2,...,N. */

/*      When MBDCND = 3 or 4, */

/*               BDA(J) = (d/dR)U(A,THETA(J)) ,    J=1,2,...,N. */

/*      When MBDCND = 5 or 6, BDA is a dummy variable. */

/*    BDB */
/*      A one-dimensional array of length N that specifies the boundary */
/*      values of the solution at R = B.  When MBDCND = 1,4, or 5, */

/*               BDB(J) = U(B,THETA(J)) ,          J=1,2,...,N. */

/*      When MBDCND = 2,3, or 6, */

/*               BDB(J) = (d/dR)U(B,THETA(J)) ,    J=1,2,...,N. */

/*    C,D */
/*      The range of THETA, i.e. C .LE. THETA .LE. D.  C must be less */
/*      than D. */

/*    N */
/*      The number of unknowns in the interval (C,D).  The unknowns in */
/*      the THETA-direction are given by THETA(J) = C + (J-0.5)DT, */
/*      J=1,2,...,N, where DT = (D-C)/N.  N must be greater than 2. */

/*    NBDCND */
/*      Indicates the type of boundary conditions at THETA = C */
/*      and THETA = D. */

/*      = 0  If the solution is periodic in THETA, i.e. */
/*           U(I,J) = U(I,N+J). */

/*      = 1  If the solution is specified at THETA = C and THETA = D */
/*           (see note below). */

/*      = 2  If the solution is specified at THETA = C and the derivative */
/*           of the solution with respect to THETA is specified at */
/*           THETA = D (see note below). */

/*      = 3  If the derivative of the solution with respect to THETA is */
/*           specified at THETA = C and THETA = D. */

/*      = 4  If the derivative of the solution with respect to THETA is */
/*           specified at THETA = C and the solution is specified at */
/*           THETA = d (see note below). */

/*      NOTE:  When NBDCND = 1, 2, or 4, do not use MBDCND = 5 or 6 (the */
/*      former indicates that the solution is specified at R =  0; the */
/*      latter indicates the solution is unspecified at R = 0).  Use */
/*      instead MBDCND = 1 or 2. */

/*    BDC */
/*      A one dimensional array of length M that specifies the boundary */
/*      values of the solution at THETA = C.   When NBDCND = 1 or 2, */

/*               BDC(I) = U(R(I),C) ,              I=1,2,...,M. */

/*      When NBDCND = 3 or 4, */

/*               BDC(I) = (d/dTHETA)U(R(I),C),     I=1,2,...,M. */

/*      When NBDCND = 0, BDC is a dummy variable. */

/*    BDD */
/*      A one-dimensional array of length M that specifies the boundary */
/*      values of the solution at THETA = D.  When NBDCND = 1 or 4, */

/*               BDD(I) = U(R(I),D) ,              I=1,2,...,M. */

/*      When NBDCND = 2 or 3, */

/*               BDD(I) = (d/dTHETA)U(R(I),D) ,    I=1,2,...,M. */

/*      When NBDCND = 0, BDD is a dummy variable. */

/*    ELMBDA */
/*      The constant LAMBDA in the Helmholtz equation.  If LAMBDA is */
/*      greater than 0, a solution may not exist.  However, HSTPLR will */
/*      attempt to find a solution. */

/*    F */
/*      A two-dimensional array that specifies the values of the right */
/*      side of the Helmholtz equation.  For I=1,2,...,M and J=1,2,...,N */

/*               F(I,J) = F(R(I),THETA(J)) . */

/*      F must be dimensioned at least M X N. */

/*    IDIMF */
/*      The row (or first) dimension of the array F as it appears in the */
/*      program calling HSTPLR.  This parameter is used to specify the */
/*      variable dimension of F.  IDIMF must be at least M. */

/*    W */
/*      A one-dimensional array that must be provided by the user for */
/*      work space.  W may require up to 13M + 4N + M*INT(log2(N)) */
/*      locations.  The actual number of locations used is computed by */
/*      HSTPLR and is returned in the location W(1). */


/*             * * * * * *   On Output   * * * * * * */

/*    F */
/*      Contains the solution U(I,J) of the finite difference */
/*      approximation for the grid point (R(I),THETA(J)) for */
/*      I=1,2,...,M, J=1,2,...,N. */

/*    PERTRB */
/*      If a combination of periodic, derivative, or unspecified */
/*      boundary conditions is specified for a Poisson equation */
/*      (LAMBDA = 0), a solution may not exist.  PERTRB is a con- */
/*      stant, calculated and subtracted from F, which ensures */
/*      that a solution exists.  HSTPLR then computes this */
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

/*      =  7  A = 0 and MBDCND = 3 or 4 */

/*      =  8  A .GT. 0 and MBDCND .GE. 5 */

/*      =  9  MBDCND .GE. 5 and NBDCND .NE. 0 or 3 */

/*      = 10  IDIMF .LT. M */

/*      = 11  LAMBDA .GT. 0 */

/*      = 12  M .LE. 2 */

/*      Since this is the only means of indicating a possibly */
/*      incorrect call to HSTPLR, the user should test IERROR after */
/*      the call. */

/*    W */
/*      W(1) contains the required length of W. */

/* *Long Description: */

/*     * * * * * * *   Program Specifications    * * * * * * * * * * * * */

/*     Dimension of   BDA(N),BDB(N),BDC(M),BDD(M),F(IDIMF,N), */
/*     Arguments      W(see ARGUMENT LIST) */

/*     Latest         June 1, 1977 */
/*     Revision */

/*     Subprograms    HSTPLR,POISTG,POSTG2,GENBUN,POISD2,POISN2,POISP2, */
/*     Required       COSGEN,MERGE,TRIX,TRI3,PIMACH */

/*     Special        NONE */
/*     Conditions */

/*     Common         NONE */
/*     Blocks */

/*     I/O            NONE */

/*     Precision      Single */

/*     Specialist     Roland Sweet */

/*     Language       FORTRAN */

/*     History        Written by Roland Sweet at NCAR in February, 1977 */

/*     Algorithm      This subroutine defines the finite-difference */
/*                    equations, incorporates boundary data, adjusts the */
/*                    right side when the system is singular and calls */
/*                    either POISTG or GENBUN which solves the linear */
/*                    system of equations. */

/*     Space          8265(decimal) = 20111(octal) LOCATIONS ON THE */
/*     Required       NCAR Control Data 7600 */

/*     Timing and        The execution time T on the NCAR Control Data */
/*     Accuracy       7600 for subroutine HSTPLR is roughly proportional */
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
/*                    The Solution Of Poisson's Equation With Neumann */
/*                    Boundary Conditions On A Staggered Grid of */
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
/* ***END PROLOGUE  HSTPLR */


/* ***FIRST EXECUTABLE STATEMENT  HSTPLR */
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
    if (*a == 0.f && (*mbdcnd == 3 || *mbdcnd == 4)) {
	*ierror = 7;
    }
    if (*a > 0.f && *mbdcnd >= 5) {
	*ierror = 8;
    }
    if (*mbdcnd >= 5 && *nbdcnd != 0 && *nbdcnd != 3) {
	*ierror = 9;
    }
    if (*idimf < *m) {
	*ierror = 10;
    }
    if (*m <= 2) {
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
    isw = 1;
    mb = *mbdcnd;
    if (*a == 0.f && *mbdcnd == 2) {
	mb = 6;
    }

/*     DEFINE A,B,C COEFFICIENTS IN W-ARRAY. */

    iwb = *m;
    iwc = iwb + *m;
    iwr = iwc + *m;
    i__1 = *m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	j = iwr + i__;
	w[j] = *a + (i__ - .5f) * deltar;
	w[i__] = (*a + (i__ - 1) * deltar) / dlrsq;
	k = iwc + i__;
	w[k] = (*a + i__ * deltar) / dlrsq;
	k = iwb + i__;
	w[k] = (*elmbda - 2.f / dlrsq) * w[j];
/* L101: */
    }
    i__1 = *m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	j = iwr + i__;
	a1 = w[j];
	i__2 = *n;
	for (j = 1; j <= i__2; ++j) {
	    f[i__ + j * f_dim1] = a1 * f[i__ + j * f_dim1];
/* L102: */
	}
/* L103: */
    }

/*     ENTER BOUNDARY DATA FOR R-BOUNDARIES. */

    switch (mb) {
	case 1:  goto L104;
	case 2:  goto L104;
	case 3:  goto L106;
	case 4:  goto L106;
	case 5:  goto L108;
	case 6:  goto L108;
    }
L104:
    a1 = w[1] * 2.f;
    w[iwb + 1] -= w[1];
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	f[j * f_dim1 + 1] -= a1 * bda[j];
/* L105: */
    }
    goto L108;
L106:
    a1 = deltar * w[1];
    w[iwb + 1] += w[1];
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	f[j * f_dim1 + 1] += a1 * bda[j];
/* L107: */
    }
L108:
    switch (mb) {
	case 1:  goto L109;
	case 2:  goto L111;
	case 3:  goto L111;
	case 4:  goto L109;
	case 5:  goto L109;
	case 6:  goto L111;
    }
L109:
    a1 = w[iwr] * 2.f;
    w[iwc] -= w[iwr];
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	f[*m + j * f_dim1] -= a1 * bdb[j];
/* L110: */
    }
    goto L113;
L111:
    a1 = deltar * w[iwr];
    w[iwc] += w[iwr];
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	f[*m + j * f_dim1] -= a1 * bdb[j];
/* L112: */
    }

/*     ENTER BOUNDARY DATA FOR THETA-BOUNDARIES. */

L113:
    a1 = 2.f / dlthsq;
    switch (np) {
	case 1:  goto L123;
	case 2:  goto L114;
	case 3:  goto L114;
	case 4:  goto L116;
	case 5:  goto L116;
    }
L114:
    i__1 = *m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	j = iwr + i__;
	f[i__ + f_dim1] -= a1 * bdc[i__] / w[j];
/* L115: */
    }
    goto L118;
L116:
    a1 = 1.f / deltht;
    i__1 = *m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	j = iwr + i__;
	f[i__ + f_dim1] += a1 * bdc[i__] / w[j];
/* L117: */
    }
L118:
    a1 = 2.f / dlthsq;
    switch (np) {
	case 1:  goto L123;
	case 2:  goto L119;
	case 3:  goto L121;
	case 4:  goto L121;
	case 5:  goto L119;
    }
L119:
    i__1 = *m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	j = iwr + i__;
	f[i__ + *n * f_dim1] -= a1 * bdd[i__] / w[j];
/* L120: */
    }
    goto L123;
L121:
    a1 = 1.f / deltht;
    i__1 = *m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	j = iwr + i__;
	f[i__ + *n * f_dim1] -= a1 * bdd[i__] / w[j];
/* L122: */
    }
L123:

/*     ADJUST RIGHT SIDE OF SINGULAR PROBLEMS TO INSURE EXISTENCE OF A */
/*     SOLUTION. */

    *pertrb = 0.f;
    if (*elmbda < 0.f) {
	goto L133;
    } else if (*elmbda == 0) {
	goto L125;
    } else {
	goto L124;
    }
L124:
    *ierror = 11;
    goto L133;
L125:
    switch (mb) {
	case 1:  goto L133;
	case 2:  goto L133;
	case 3:  goto L126;
	case 4:  goto L133;
	case 5:  goto L133;
	case 6:  goto L126;
    }
L126:
    switch (np) {
	case 1:  goto L127;
	case 2:  goto L133;
	case 3:  goto L133;
	case 4:  goto L127;
	case 5:  goto L133;
    }
L127:
    isw = 2;
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	i__2 = *m;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    *pertrb += f[i__ + j * f_dim1];
/* L128: */
	}
/* L129: */
    }
    *pertrb /= *m * *n * .5f * (*a + *b);
    i__1 = *m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	j = iwr + i__;
	a1 = *pertrb * w[j];
	i__2 = *n;
	for (j = 1; j <= i__2; ++j) {
	    f[i__ + j * f_dim1] -= a1;
/* L130: */
	}
/* L131: */
    }
    a2 = 0.f;
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	a2 += f[j * f_dim1 + 1];
/* L132: */
    }
    a2 /= w[iwr + 1];
L133:

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
/* L134: */
	}
/* L135: */
    }
    lp = *nbdcnd;
    w[1] = 0.f;
    w[iwr] = 0.f;

/*     CALL POISTG OR GENBUN TO SOLVE THE SYSTEM OF EQUATIONS. */

    if (lp == 0) {
	goto L136;
    }
    poistg_(&lp, n, &c__1, m, &w[1], &w[iwb + 1], &w[iwc + 1], idimf, &f[
	    f_offset], &ierr1, &w[iwr + 1]);
    goto L137;
L136:
    genbun_(&lp, n, &c__1, m, &w[1], &w[iwb + 1], &w[iwc + 1], idimf, &f[
	    f_offset], &ierr1, &w[iwr + 1]);
L137:
    w[1] = w[iwr + 1] + *m * 3;
    if (*a != 0.f || *mbdcnd != 2 || isw != 2) {
	goto L141;
    }
    a1 = 0.f;
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	a1 += f[j * f_dim1 + 1];
/* L138: */
    }
    a1 = (a1 - dlrsq * a2 / 16.f) / *n;
    if (*nbdcnd == 3) {
	a1 += (bdd[1] - bdc[1]) / (*d__ - *c__);
    }
    a1 = bda[1] - a1;
    i__1 = *m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = *n;
	for (j = 1; j <= i__2; ++j) {
	    f[i__ + j * f_dim1] += a1;
/* L139: */
	}
/* L140: */
    }
L141:
    return 0;
} /* hstplr_ */

