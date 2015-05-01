/* hstcrt.f -- translated by f2c (version 12.02.01).
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

/* DECK HSTCRT */
/* Subroutine */ int hstcrt_(real *a, real *b, integer *m, integer *mbdcnd, 
	real *bda, real *bdb, real *c__, real *d__, integer *n, integer *
	nbdcnd, real *bdc, real *bdd, real *elmbda, real *f, integer *idimf, 
	real *pertrb, integer *ierror, real *w)
{
    /* System generated locals */
    integer f_dim1, f_offset, i__1, i__2;
    real r__1;

    /* Local variables */
    static integer i__, j;
    static real s;
    static integer mp, np, id2, id3, id4;
    static real st2;
    static integer ierr1;
    extern /* Subroutine */ int genbun_(integer *, integer *, integer *, 
	    integer *, real *, real *, real *, integer *, real *, integer *, 
	    real *);
    static real deltax, deltay;
    static integer mperod, nperod;
    static real delxsq, delysq;
    extern /* Subroutine */ int poistg_(integer *, integer *, integer *, 
	    integer *, real *, real *, real *, integer *, real *, integer *, 
	    real *);
    static real twdelx, twdely, twdysq;

/* ***BEGIN PROLOGUE  HSTCRT */
/* ***PURPOSE  Solve the standard five-point finite difference */
/*            approximation on a staggered grid to the Helmholtz equation */
/*            in Cartesian coordinates. */
/* ***LIBRARY   SLATEC (FISHPACK) */
/* ***CATEGORY  I2B1A1A */
/* ***TYPE      SINGLE PRECISION (HSTCRT-S) */
/* ***KEYWORDS  ELLIPTIC, FISHPACK, HELMHOLTZ, PDE */
/* ***AUTHOR  Adams, J., (NCAR) */
/*           Swarztrauber, P. N., (NCAR) */
/*           Sweet, R., (NCAR) */
/* ***DESCRIPTION */

/*      HSTCRT solves the standard five-point finite difference */
/*      approximation on a staggered grid to the Helmholtz equation in */
/*      Cartesian coordinates */

/*      (d/dX)(dU/dX) + (d/dY)(dU/dY) + LAMBDA*U = F(X,Y) */

/*     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/*     * * * * * * * *    Parameter Description     * * * * * * * * * * */

/*             * * * * * *   On Input    * * * * * * */

/*    A,B */
/*      The range of X, i.e. A .LE. X .LE. B.  A must be less than B. */

/*    M */
/*      The number of grid points in the interval (A,B).  The grid points */
/*      in the X-direction are given by X(I) = A + (I-0.5)dX for */
/*      I=1,2,...,M where dX =(B-A)/M.  M must be greater than 2. */

/*    MBDCND */
/*      Indicates the type of boundary conditions at X = A and X = B. */

/*      = 0  If the solution is periodic in X, */
/*           U(M+I,J) = U(I,J). */

/*      = 1  If the solution is specified at X = A and X = B. */

/*      = 2  If the solution is specified at X = A and the derivative */
/*           of the solution with respect to X is specified at X = B. */

/*      = 3  If the derivative of the solution with respect to X is */
/*           specified at X = A  and X = B. */

/*      = 4  If the derivative of the solution with respect to X is */
/*           specified at X = A  and the solution is specified at X = B. */

/*    BDA */
/*      A one-dimensional array of length N that specifies the boundary */
/*      values (if any) of the solution at X = A.  When MBDCND = 1 or 2, */

/*               BDA(J) = U(A,Y(J)) ,          J=1,2,...,N. */

/*      When MBDCND = 3 or 4, */

/*               BDA(J) = (d/dX)U(A,Y(J)) ,    J=1,2,...,N. */

/*    BDB */
/*      A one-dimensional array of length N that specifies the boundary */
/*      values of the solution at X = B.  When MBDCND = 1 or 4 */

/*               BDB(J) = U(B,Y(J)) ,          J=1,2,...,N. */

/*      When MBDCND = 2 or 3 */

/*               BDB(J) = (d/dX)U(B,Y(J)) ,    J=1,2,...,N. */

/*    C,D */
/*      The range of Y, i.e. C .LE. Y .LE. D.  C must be less */
/*      than D. */

/*    N */
/*      The number of unknowns in the interval (C,D).  The unknowns in */
/*      the Y-direction are given by Y(J) = C + (J-0.5)DY, */
/*      J=1,2,...,N, where DY = (D-C)/N.  N must be greater than 2. */

/*    NBDCND */
/*      Indicates the type of boundary conditions at Y = C */
/*      and Y = D. */

/*      = 0  If the solution is periodic in Y, i.e. */
/*           U(I,J) = U(I,N+J). */

/*      = 1  If the solution is specified at Y = C and Y = D. */

/*      = 2  If the solution is specified at Y = C and the derivative */
/*           of the solution with respect to Y is specified at Y = D. */

/*      = 3  If the derivative of the solution with respect to Y is */
/*           specified at Y = C and Y = D. */

/*      = 4  If the derivative of the solution with respect to Y is */
/*           specified at Y = C and the solution is specified at Y = D. */

/*    BDC */
/*      A one dimensional array of length M that specifies the boundary */
/*      values of the solution at Y = C.   When NBDCND = 1 or 2, */

/*               BDC(I) = U(X(I),C) ,              I=1,2,...,M. */

/*      When NBDCND = 3 or 4, */

/*               BDC(I) = (d/dY)U(X(I),C),     I=1,2,...,M. */

/*      When NBDCND = 0, BDC is a dummy variable. */

/*    BDD */
/*      A one-dimensional array of length M that specifies the boundary */
/*      values of the solution at Y = D.  When NBDCND = 1 or 4, */

/*               BDD(I) = U(X(I),D) ,              I=1,2,...,M. */

/*      When NBDCND = 2 or 3, */

/*               BDD(I) = (d/dY)U(X(I),D) ,    I=1,2,...,M. */

/*      When NBDCND = 0, BDD is a dummy variable. */

/*    ELMBDA */
/*      The constant LAMBDA in the Helmholtz equation.  If LAMBDA is */
/*      greater than 0, a solution may not exist.  However, HSTCRT will */
/*      attempt to find a solution. */

/*    F */
/*      A two-dimensional array that specifies the values of the right */
/*      side of the Helmholtz equation.  For I=1,2,...,M and J=1,2,...,N */

/*               F(I,J) = F(X(I),Y(J)) . */

/*      F must be dimensioned at least M X N. */

/*    IDIMF */
/*      The row (or first) dimension of the array F as it appears in the */
/*      program calling HSTCRT.  This parameter is used to specify the */
/*      variable dimension of F.  IDIMF must be at least M. */

/*    W */
/*      A one-dimensional array that must be provided by the user for */
/*      work space.  W may require up to 13M + 4N + M*INT(log2(N)) */
/*      locations.  The actual number of locations used is computed by */
/*      HSTCRT and is returned in the location W(1). */


/*             * * * * * *   On Output   * * * * * * */

/*    F */
/*      Contains the solution U(I,J) of the finite difference */
/*      approximation for the grid point (X(I),Y(J)) for */
/*      I=1,2,...,M, J=1,2,...,N. */

/*    PERTRB */
/*      If a combination of periodic or derivative boundary conditions is */
/*      specified for a Poisson equation (LAMBDA = 0), a solution may not */
/*      exist.  PERTRB is a constant, calculated and subtracted from F, */
/*      which ensures that a solution exists.  HSTCRT then computes this */
/*      solution, which is a least squares solution to the original */
/*      approximation.  This solution plus any constant is also a */
/*      solution; hence, the solution is not unique.  The value of PERTRB */
/*      should be small compared to the right side F.  Otherwise, a */
/*      solution is obtained to an essentially different problem.  This */
/*      comparison should always be made to insure that a meaningful */
/*      solution has been obtained. */

/*    IERROR */
/*      An error flag that indicates invalid input parameters. */
/*       Except for numbers 0 and  6, a solution is not attempted. */

/*      =  0  No error */

/*      =  1  A .GE. B */

/*      =  2  MBDCND .LT. 0 or MBDCND .GT. 4 */

/*      =  3  C .GE. D */

/*      =  4  N .LE. 2 */

/*      =  5  NBDCND .LT. 0 or NBDCND .GT. 4 */

/*      =  6  LAMBDA .GT. 0 */

/*      =  7  IDIMF .LT. M */

/*      =  8  M .LE. 2 */

/*      Since this is the only means of indicating a possibly */
/*      incorrect call to HSTCRT, the user should test IERROR after */
/*      the call. */

/*    W */
/*      W(1) contains the required length of W. */

/* *Long Description: */

/*     * * * * * * *   Program Specifications    * * * * * * * * * * * * */

/*     Dimension of   BDA(N),BDB(N),BDC(M),BDD(M),F(IDIMF,N), */
/*     Arguments      W(See argument list) */

/*     Latest         June 1, 1977 */
/*     Revision */

/*     Subprograms    HSTCRT,POISTG,POSTG2,GENBUN,POISD2,POISN2,POISP2, */
/*     Required       COSGEN,MERGE,TRIX,TRI3,PIMACH */

/*     Special        NONE */
/*     Conditions */

/*     Common         NONE */
/*     Blocks */

/*     I/O            NONE */

/*     Precision      Single */

/*     Specialist     Roland Sweet */

/*     Language       FORTRAN */

/*     History        Written by Roland Sweet at NCAR in January , 1977 */

/*     Algorithm      This subroutine defines the finite-difference */
/*                    equations, incorporates boundary data, adjusts the */
/*                    right side when the system is singular and calls */
/*                    either POISTG or GENBUN which solves the linear */
/*                    system of equations. */

/*     Space          8131(decimal) = 17703(octal) locations on the */
/*     Required       NCAR Control Data 7600 */

/*     Timing and        The execution time T on the NCAR Control Data */
/*     Accuracy       7600 for subroutine HSTCRT is roughly proportional */
/*                    to M*N*log2(N).  Some typical values are listed in */
/*                    the table below. */
/*                       The solution process employed results in a loss */
/*                    of no more than FOUR significant digits for N and M */
/*                    as large as 64.  More detailed information about */
/*                    accuracy can be found in the documentation for */
/*                    subroutine POISTG which is the routine that */
/*                    actually solves the finite difference equations. */


/*                       M(=N)    MBDCND    NBDCND    T(MSECS) */
/*                       -----    ------    ------    -------- */

/*                        32       1-4       1-4         56 */
/*                        64       1-4       1-4        230 */

/*     Portability    American National Standards Institute Fortran. */
/*                    The machine dependent constant PI is defined in */
/*                    function PIMACH. */

/*     Required       COS */
/*     Resident */
/*     Routines */

/*     Reference      Schumann, U. and R. Sweet,'A Direct Method For */
/*                    The Solution Of Poisson's Equation With Neumann */
/*                    Boundary Conditions On A Staggered Grid Of */
/*                    Arbitrary Size,' J. COMP. PHYS. 20(1976), */
/*                    PP. 171-182. */

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
/* ***END PROLOGUE  HSTCRT */


/* ***FIRST EXECUTABLE STATEMENT  HSTCRT */
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
    if (*a >= *b) {
	*ierror = 1;
    }
    if (*mbdcnd < 0 || *mbdcnd > 4) {
	*ierror = 2;
    }
    if (*c__ >= *d__) {
	*ierror = 3;
    }
    if (*n <= 2) {
	*ierror = 4;
    }
    if (*nbdcnd < 0 || *nbdcnd > 4) {
	*ierror = 5;
    }
    if (*idimf < *m) {
	*ierror = 7;
    }
    if (*m <= 2) {
	*ierror = 8;
    }
    if (*ierror != 0) {
	return 0;
    }
    nperod = *nbdcnd;
    mperod = 0;
    if (*mbdcnd > 0) {
	mperod = 1;
    }
    deltax = (*b - *a) / *m;
    twdelx = 1.f / deltax;
/* Computing 2nd power */
    r__1 = deltax;
    delxsq = 2.f / (r__1 * r__1);
    deltay = (*d__ - *c__) / *n;
    twdely = 1.f / deltay;
/* Computing 2nd power */
    r__1 = deltay;
    delysq = r__1 * r__1;
    twdysq = 2.f / delysq;
    np = *nbdcnd + 1;
    mp = *mbdcnd + 1;

/*     DEFINE THE A,B,C COEFFICIENTS IN W-ARRAY. */

    id2 = *m;
    id3 = id2 + *m;
    id4 = id3 + *m;
/* Computing 2nd power */
    r__1 = deltay / deltax;
    s = r__1 * r__1;
    st2 = s * 2.f;
    i__1 = *m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	w[i__] = s;
	j = id2 + i__;
	w[j] = -st2 + *elmbda * delysq;
	j = id3 + i__;
	w[j] = s;
/* L101: */
    }

/*     ENTER BOUNDARY DATA FOR X-BOUNDARIES. */

    switch (mp) {
	case 1:  goto L111;
	case 2:  goto L102;
	case 3:  goto L102;
	case 4:  goto L104;
	case 5:  goto L104;
    }
L102:
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	f[j * f_dim1 + 1] -= bda[j] * delxsq;
/* L103: */
    }
    w[id2 + 1] -= w[1];
    goto L106;
L104:
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	f[j * f_dim1 + 1] += bda[j] * twdelx;
/* L105: */
    }
    w[id2 + 1] += w[1];
L106:
    switch (mp) {
	case 1:  goto L111;
	case 2:  goto L107;
	case 3:  goto L109;
	case 4:  goto L109;
	case 5:  goto L107;
    }
L107:
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	f[*m + j * f_dim1] -= bdb[j] * delxsq;
/* L108: */
    }
    w[id3] -= w[1];
    goto L111;
L109:
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	f[*m + j * f_dim1] -= bdb[j] * twdelx;
/* L110: */
    }
    w[id3] += w[1];
L111:

/*     ENTER BOUNDARY DATA FOR Y-BOUNDARIES. */

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
	f[i__ + f_dim1] -= bdc[i__] * twdysq;
/* L113: */
    }
    goto L116;
L114:
    i__1 = *m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	f[i__ + f_dim1] += bdc[i__] * twdely;
/* L115: */
    }
L116:
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
	f[i__ + *n * f_dim1] -= bdd[i__] * twdysq;
/* L118: */
    }
    goto L121;
L119:
    i__1 = *m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	f[i__ + *n * f_dim1] -= bdd[i__] * twdely;
/* L120: */
    }
L121:
    i__1 = *m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = *n;
	for (j = 1; j <= i__2; ++j) {
	    f[i__ + j * f_dim1] *= delysq;
/* L122: */
	}
/* L123: */
    }
    if (mperod == 0) {
	goto L124;
    }
    w[1] = 0.f;
    w[id4] = 0.f;
L124:
    *pertrb = 0.f;
    if (*elmbda < 0.f) {
	goto L133;
    } else if (*elmbda == 0) {
	goto L126;
    } else {
	goto L125;
    }
L125:
    *ierror = 6;
    goto L133;
L126:
    switch (mp) {
	case 1:  goto L127;
	case 2:  goto L133;
	case 3:  goto L133;
	case 4:  goto L127;
	case 5:  goto L133;
    }
L127:
    switch (np) {
	case 1:  goto L128;
	case 2:  goto L133;
	case 3:  goto L133;
	case 4:  goto L128;
	case 5:  goto L133;
    }

/*     FOR SINGULAR PROBLEMS MUST ADJUST DATA TO INSURE THAT A SOLUTION */
/*     WILL EXIST. */

L128:
    s = 0.f;
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	i__2 = *m;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    s += f[i__ + j * f_dim1];
/* L129: */
	}
/* L130: */
    }
    *pertrb = s / (*m * *n);
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	i__2 = *m;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    f[i__ + j * f_dim1] -= *pertrb;
/* L131: */
	}
/* L132: */
    }
    *pertrb /= delysq;

/*     SOLVE THE EQUATION. */

L133:
    if (nperod == 0) {
	goto L134;
    }
    poistg_(&nperod, n, &mperod, m, &w[1], &w[id2 + 1], &w[id3 + 1], idimf, &
	    f[f_offset], &ierr1, &w[id4 + 1]);
    goto L135;
L134:
    genbun_(&nperod, n, &mperod, m, &w[1], &w[id2 + 1], &w[id3 + 1], idimf, &
	    f[f_offset], &ierr1, &w[id4 + 1]);
L135:
    w[1] = w[id4 + 1] + *m * 3;
    return 0;
} /* hstcrt_ */

