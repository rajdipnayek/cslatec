/* hwscrt.f -- translated by f2c (version 12.02.01).
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

/* DECK HWSCRT */
/* Subroutine */ int hwscrt_(real *a, real *b, integer *m, integer *mbdcnd, 
	real *bda, real *bdb, real *c__, real *d__, integer *n, integer *
	nbdcnd, real *bdc, real *bdd, real *elmbda, real *f, integer *idimf, 
	real *pertrb, integer *ierror, real *w)
{
    /* System generated locals */
    integer f_dim1, f_offset, i__1, i__2;
    real r__1;

    /* Local variables */
    static integer i__, j;
    static real s, a1, a2, s1;
    static integer mp, np, id2, id3, id4, mp1, np1;
    static real st2;
    static integer msp1, nsp1, munk, nunk, ierr1, mstm1, nstm1, mskip, nskip, 
	    mstop, nstop;
    extern /* Subroutine */ int genbun_(integer *, integer *, integer *, 
	    integer *, real *, real *, real *, integer *, real *, integer *, 
	    real *);
    static real deltax, deltay;
    static integer mperod, nperod;
    static real delxsq, delysq, twdelx, twdely;
    static integer nstart, mstart;

/* ***BEGIN PROLOGUE  HWSCRT */
/* ***PURPOSE  Solves the standard five-point finite difference */
/*            approximation to the Helmholtz equation in Cartesian */
/*            coordinates. */
/* ***LIBRARY   SLATEC (FISHPACK) */
/* ***CATEGORY  I2B1A1A */
/* ***TYPE      SINGLE PRECISION (HWSCRT-S) */
/* ***KEYWORDS  CARTESIAN, ELLIPTIC, FISHPACK, HELMHOLTZ, PDE */
/* ***AUTHOR  Adams, J., (NCAR) */
/*           Swarztrauber, P. N., (NCAR) */
/*           Sweet, R., (NCAR) */
/* ***DESCRIPTION */

/*     Subroutine HWSCRT solves the standard five-point finite */
/*     difference approximation to the Helmholtz equation in Cartesian */
/*     coordinates: */

/*          (d/dX)(dU/dX) + (d/dY)(dU/dY) + LAMBDA*U = F(X,Y). */



/*     * * * * * * * *    Parameter Description     * * * * * * * * * * */

/*             * * * * * *   On Input    * * * * * * */

/*     A,B */
/*       The range of X, i.e., A .LE. X .LE. B.  A must be less than B. */

/*     M */
/*       The number of panels into which the interval (A,B) is */
/*       subdivided.  Hence, there will be M+1 grid points in the */
/*       X-direction given by X(I) = A+(I-1)DX for I = 1,2,...,M+1, */
/*       where DX = (B-A)/M is the panel width. M must be greater than 3. */

/*     MBDCND */
/*       Indicates the type of boundary conditions at X = A and X = B. */

/*       = 0  If the solution is periodic in X, i.e., U(I,J) = U(M+I,J). */
/*       = 1  If the solution is specified at X = A and X = B. */
/*       = 2  If the solution is specified at X = A and the derivative of */
/*            the solution with respect to X is specified at X = B. */
/*       = 3  If the derivative of the solution with respect to X is */
/*            specified at X = A and X = B. */
/*       = 4  If the derivative of the solution with respect to X is */
/*            specified at X = A and the solution is specified at X = B. */

/*     BDA */
/*       A one-dimensional array of length N+1 that specifies the values */
/*       of the derivative of the solution with respect to X at X = A. */
/*       When MBDCND = 3 or 4, */

/*            BDA(J) = (d/dX)U(A,Y(J)), J = 1,2,...,N+1  . */

/*       When MBDCND has any other value, BDA is a dummy variable. */

/*     BDB */
/*       A one-dimensional array of length N+1 that specifies the values */
/*       of the derivative of the solution with respect to X at X = B. */
/*       When MBDCND = 2 or 3, */

/*            BDB(J) = (d/dX)U(B,Y(J)), J = 1,2,...,N+1  . */

/*       When MBDCND has any other value BDB is a dummy variable. */

/*     C,D */
/*       The range of Y, i.e., C .LE. Y .LE. D.  C must be less than D. */

/*     N */
/*       The number of panels into which the interval (C,D) is */
/*       subdivided.  Hence, there will be N+1 grid points in the */
/*       Y-direction given by Y(J) = C+(J-1)DY for J = 1,2,...,N+1, where */
/*       DY = (D-C)/N is the panel width.  N must be greater than 3. */

/*     NBDCND */
/*       Indicates the type of boundary conditions at Y = C and Y = D. */

/*       = 0  If the solution is periodic in Y, i.e., U(I,J) = U(I,N+J). */
/*       = 1  If the solution is specified at Y = C and Y = D. */
/*       = 2  If the solution is specified at Y = C and the derivative of */
/*            the solution with respect to Y is specified at Y = D. */
/*       = 3  If the derivative of the solution with respect to Y is */
/*            specified at Y = C and Y = D. */
/*       = 4  If the derivative of the solution with respect to Y is */
/*            specified at Y = C and the solution is specified at Y = D. */

/*     BDC */
/*       A one-dimensional array of length M+1 that specifies the values */
/*       of the derivative of the solution with respect to Y at Y = C. */
/*       When NBDCND = 3 or 4, */

/*            BDC(I) = (d/dY)U(X(I),C), I = 1,2,...,M+1  . */

/*       When NBDCND has any other value, BDC is a dummy variable. */

/*     BDD */
/*       A one-dimensional array of length M+1 that specifies the values */
/*       of the derivative of the solution with respect to Y at Y = D. */
/*       When NBDCND = 2 or 3, */

/*            BDD(I) = (d/dY)U(X(I),D), I = 1,2,...,M+1  . */

/*       When NBDCND has any other value, BDD is a dummy variable. */

/*     ELMBDA */
/*       The constant LAMBDA in the Helmholtz equation.  If */
/*       LAMBDA .GT. 0, a solution may not exist.  However, HWSCRT will */
/*       attempt to find a solution. */

/*     F */
/*       A two-dimensional array which specifies the values of the right */
/*       side of the Helmholtz equation and boundary values (if any). */
/*       For I = 2,3,...,M and J = 2,3,...,N */

/*            F(I,J) = F(X(I),Y(J)). */

/*       On the boundaries F is defined by */

/*            MBDCND     F(1,J)        F(M+1,J) */
/*            ------     ---------     -------- */

/*              0        F(A,Y(J))     F(A,Y(J)) */
/*              1        U(A,Y(J))     U(B,Y(J)) */
/*              2        U(A,Y(J))     F(B,Y(J))     J = 1,2,...,N+1 */
/*              3        F(A,Y(J))     F(B,Y(J)) */
/*              4        F(A,Y(J))     U(B,Y(J)) */


/*            NBDCND     F(I,1)        F(I,N+1) */
/*            ------     ---------     -------- */

/*              0        F(X(I),C)     F(X(I),C) */
/*              1        U(X(I),C)     U(X(I),D) */
/*              2        U(X(I),C)     F(X(I),D)     I = 1,2,...,M+1 */
/*              3        F(X(I),C)     F(X(I),D) */
/*              4        F(X(I),C)     U(X(I),D) */

/*       F must be dimensioned at least (M+1)*(N+1). */

/*       NOTE: */

/*       If the table calls for both the solution U and the right side F */
/*       at a corner then the solution must be specified. */

/*     IDIMF */
/*       The row (or first) dimension of the array F as it appears in the */
/*       program calling HWSCRT.  This parameter is used to specify the */
/*       variable dimension of F.  IDIMF must be at least M+1  . */

/*     W */
/*       A one-dimensional array that must be provided by the user for */
/*       work space.  W may require up to 4*(N+1) + */
/*       (13 + INT(log2(N+1)))*(M+1) locations.  The actual number of */
/*       locations used is computed by HWSCRT and is returned in location */
/*       W(1). */


/*             * * * * * *   On Output     * * * * * * */

/*     F */
/*       Contains the solution U(I,J) of the finite difference */
/*       approximation for the grid point (X(I),Y(J)), I = 1,2,...,M+1, */
/*       J = 1,2,...,N+1  . */

/*     PERTRB */
/*       If a combination of periodic or derivative boundary conditions */
/*       is specified for a Poisson equation (LAMBDA = 0), a solution may */
/*       not exist.  PERTRB is a constant, calculated and subtracted from */
/*       F, which ensures that a solution exists.  HWSCRT then computes */
/*       this solution, which is a least squares solution to the original */
/*       approximation.  This solution plus any constant is also a */
/*       solution.  Hence, the solution is not unique.  The value of */
/*       PERTRB should be small compared to the right side F.  Otherwise, */
/*       a solution is obtained to an essentially different problem. */
/*       This comparison should always be made to insure that a */
/*       meaningful solution has been obtained. */

/*     IERROR */
/*       An error flag that indicates invalid input parameters.  Except */
/*       for numbers 0 and 6, a solution is not attempted. */

/*       = 0  No error. */
/*       = 1  A .GE. B. */
/*       = 2  MBDCND .LT. 0 or MBDCND .GT. 4  . */
/*       = 3  C .GE. D. */
/*       = 4  N .LE. 3 */
/*       = 5  NBDCND .LT. 0 or NBDCND .GT. 4  . */
/*       = 6  LAMBDA .GT. 0  . */
/*       = 7  IDIMF .LT. M+1  . */
/*       = 8  M .LE. 3 */

/*       Since this is the only means of indicating a possibly incorrect */
/*       call to HWSCRT, the user should test IERROR after the call. */

/*     W */
/*       W(1) contains the required length of W. */

/* *Long Description: */

/*     * * * * * * *   Program Specifications    * * * * * * * * * * * * */


/*     Dimension of   BDA(N+1),BDB(N+1),BDC(M+1),BDD(M+1),F(IDIMF,N+1), */
/*     Arguments      W(see argument list) */

/*     Latest         June 1, 1976 */
/*     Revision */

/*     Subprograms    HWSCRT,GENBUN,POISD2,POISN2,POISP2,COSGEN,MERGE, */
/*     Required       TRIX,TRI3,PIMACH */

/*     Special        NONE */
/*     Conditions */

/*     Common         NONE */
/*     Blocks */

/*     I/O            NONE */

/*     Precision      Single */

/*     Specialist     Roland Sweet */

/*     Language       FORTRAN */

/*     History        Standardized September 1, 1973 */
/*                    Revised April 1, 1976 */

/*     Algorithm      The routine defines the finite difference */
/*                    equations, incorporates boundary data, and adjusts */
/*                    the right side of singular systems and then calls */
/*                    GENBUN to solve the system. */

/*     Space          13110(octal) = 5704(decimal) locations on the NCAR */
/*     Required       Control Data 7600 */

/*     Timing and        The execution time T on the NCAR Control Data */
/*     Accuracy       7600 for subroutine HWSCRT is roughly proportional */
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

/*     Reference      Swarztrauber, P. and R. Sweet, 'Efficient FORTRAN */
/*                    Subprograms for The Solution Of Elliptic Equations' */
/*                    NCAR TN/IA-109, July, 1975, 138 pp. */

/*     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/* ***REFERENCES  P. N. Swarztrauber and R. Sweet, Efficient Fortran */
/*                 subprograms for the solution of elliptic equations, */
/*                 NCAR TN/IA-109, July 1975, 138 pp. */
/* ***ROUTINES CALLED  GENBUN */
/* ***REVISION HISTORY  (YYMMDD) */
/*   801001  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890531  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   920501  Reformatted the REFERENCES section.  (WRB) */
/* ***END PROLOGUE  HWSCRT */


/* ***FIRST EXECUTABLE STATEMENT  HWSCRT */
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
    if (*n <= 3) {
	*ierror = 4;
    }
    if (*nbdcnd < 0 || *nbdcnd > 4) {
	*ierror = 5;
    }
    if (*idimf < *m + 1) {
	*ierror = 7;
    }
    if (*m <= 3) {
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
    twdelx = 2.f / deltax;
/* Computing 2nd power */
    r__1 = deltax;
    delxsq = 1.f / (r__1 * r__1);
    deltay = (*d__ - *c__) / *n;
    twdely = 2.f / deltay;
/* Computing 2nd power */
    r__1 = deltay;
    delysq = 1.f / (r__1 * r__1);
    np = *nbdcnd + 1;
    np1 = *n + 1;
    mp = *mbdcnd + 1;
    mp1 = *m + 1;
    nstart = 1;
    nstop = *n;
    nskip = 1;
    switch (np) {
	case 1:  goto L104;
	case 2:  goto L101;
	case 3:  goto L102;
	case 4:  goto L103;
	case 5:  goto L104;
    }
L101:
    nstart = 2;
    goto L104;
L102:
    nstart = 2;
L103:
    nstop = np1;
    nskip = 2;
L104:
    nunk = nstop - nstart + 1;

/*     ENTER BOUNDARY DATA FOR X-BOUNDARIES. */

    mstart = 1;
    mstop = *m;
    mskip = 1;
    switch (mp) {
	case 1:  goto L117;
	case 2:  goto L105;
	case 3:  goto L106;
	case 4:  goto L109;
	case 5:  goto L110;
    }
L105:
    mstart = 2;
    goto L107;
L106:
    mstart = 2;
    mstop = mp1;
    mskip = 2;
L107:
    i__1 = nstop;
    for (j = nstart; j <= i__1; ++j) {
	f[j * f_dim1 + 2] -= f[j * f_dim1 + 1] * delxsq;
/* L108: */
    }
    goto L112;
L109:
    mstop = mp1;
    mskip = 2;
L110:
    i__1 = nstop;
    for (j = nstart; j <= i__1; ++j) {
	f[j * f_dim1 + 1] += bda[j] * twdelx;
/* L111: */
    }
L112:
    switch (mskip) {
	case 1:  goto L113;
	case 2:  goto L115;
    }
L113:
    i__1 = nstop;
    for (j = nstart; j <= i__1; ++j) {
	f[*m + j * f_dim1] -= f[mp1 + j * f_dim1] * delxsq;
/* L114: */
    }
    goto L117;
L115:
    i__1 = nstop;
    for (j = nstart; j <= i__1; ++j) {
	f[mp1 + j * f_dim1] -= bdb[j] * twdelx;
/* L116: */
    }
L117:
    munk = mstop - mstart + 1;

/*     ENTER BOUNDARY DATA FOR Y-BOUNDARIES. */

    switch (np) {
	case 1:  goto L127;
	case 2:  goto L118;
	case 3:  goto L118;
	case 4:  goto L120;
	case 5:  goto L120;
    }
L118:
    i__1 = mstop;
    for (i__ = mstart; i__ <= i__1; ++i__) {
	f[i__ + (f_dim1 << 1)] -= f[i__ + f_dim1] * delysq;
/* L119: */
    }
    goto L122;
L120:
    i__1 = mstop;
    for (i__ = mstart; i__ <= i__1; ++i__) {
	f[i__ + f_dim1] += bdc[i__] * twdely;
/* L121: */
    }
L122:
    switch (nskip) {
	case 1:  goto L123;
	case 2:  goto L125;
    }
L123:
    i__1 = mstop;
    for (i__ = mstart; i__ <= i__1; ++i__) {
	f[i__ + *n * f_dim1] -= f[i__ + np1 * f_dim1] * delysq;
/* L124: */
    }
    goto L127;
L125:
    i__1 = mstop;
    for (i__ = mstart; i__ <= i__1; ++i__) {
	f[i__ + np1 * f_dim1] -= bdd[i__] * twdely;
/* L126: */
    }

/*    MULTIPLY RIGHT SIDE BY DELTAY**2. */

L127:
    delysq = deltay * deltay;
    i__1 = mstop;
    for (i__ = mstart; i__ <= i__1; ++i__) {
	i__2 = nstop;
	for (j = nstart; j <= i__2; ++j) {
	    f[i__ + j * f_dim1] *= delysq;
/* L128: */
	}
/* L129: */
    }

/*     DEFINE THE A,B,C COEFFICIENTS IN W-ARRAY. */

    id2 = munk;
    id3 = id2 + munk;
    id4 = id3 + munk;
    s = delysq * delxsq;
    st2 = s * 2.f;
    i__1 = munk;
    for (i__ = 1; i__ <= i__1; ++i__) {
	w[i__] = s;
	j = id2 + i__;
	w[j] = -st2 + *elmbda * delysq;
	j = id3 + i__;
	w[j] = s;
/* L130: */
    }
    if (mp == 1) {
	goto L131;
    }
    w[1] = 0.f;
    w[id4] = 0.f;
L131:
    switch (mp) {
	case 1:  goto L135;
	case 2:  goto L135;
	case 3:  goto L132;
	case 4:  goto L133;
	case 5:  goto L134;
    }
L132:
    w[id2] = st2;
    goto L135;
L133:
    w[id2] = st2;
L134:
    w[id3 + 1] = st2;
L135:
    *pertrb = 0.f;
    if (*elmbda < 0.f) {
	goto L144;
    } else if (*elmbda == 0) {
	goto L137;
    } else {
	goto L136;
    }
L136:
    *ierror = 6;
    goto L144;
L137:
    if ((*nbdcnd == 0 || *nbdcnd == 3) && (*mbdcnd == 0 || *mbdcnd == 3)) {
	goto L138;
    }
    goto L144;

/*     FOR SINGULAR PROBLEMS MUST ADJUST DATA TO INSURE THAT A SOLUTION */
/*     WILL EXIST. */

L138:
    a1 = 1.f;
    a2 = 1.f;
    if (*nbdcnd == 3) {
	a2 = 2.f;
    }
    if (*mbdcnd == 3) {
	a1 = 2.f;
    }
    s1 = 0.f;
    msp1 = mstart + 1;
    mstm1 = mstop - 1;
    nsp1 = nstart + 1;
    nstm1 = nstop - 1;
    i__1 = nstm1;
    for (j = nsp1; j <= i__1; ++j) {
	s = 0.f;
	i__2 = mstm1;
	for (i__ = msp1; i__ <= i__2; ++i__) {
	    s += f[i__ + j * f_dim1];
/* L139: */
	}
	s1 = s1 + s * a1 + f[mstart + j * f_dim1] + f[mstop + j * f_dim1];
/* L140: */
    }
    s1 = a2 * s1;
    s = 0.f;
    i__1 = mstm1;
    for (i__ = msp1; i__ <= i__1; ++i__) {
	s = s + f[i__ + nstart * f_dim1] + f[i__ + nstop * f_dim1];
/* L141: */
    }
    s1 = s1 + s * a1 + f[mstart + nstart * f_dim1] + f[mstart + nstop * 
	    f_dim1] + f[mstop + nstart * f_dim1] + f[mstop + nstop * f_dim1];
    s = ((nunk - 2) * a2 + 2.f) * ((munk - 2) * a1 + 2.f);
    *pertrb = s1 / s;
    i__1 = nstop;
    for (j = nstart; j <= i__1; ++j) {
	i__2 = mstop;
	for (i__ = mstart; i__ <= i__2; ++i__) {
	    f[i__ + j * f_dim1] -= *pertrb;
/* L142: */
	}
/* L143: */
    }
    *pertrb /= delysq;

/*     SOLVE THE EQUATION. */

L144:
    genbun_(&nperod, &nunk, &mperod, &munk, &w[1], &w[id2 + 1], &w[id3 + 1], 
	    idimf, &f[mstart + nstart * f_dim1], &ierr1, &w[id4 + 1]);
    w[1] = w[id4 + 1] + munk * 3;

/*     FILL IN IDENTICAL VALUES WHEN HAVE PERIODIC BOUNDARY CONDITIONS. */

    if (*nbdcnd != 0) {
	goto L146;
    }
    i__1 = mstop;
    for (i__ = mstart; i__ <= i__1; ++i__) {
	f[i__ + np1 * f_dim1] = f[i__ + f_dim1];
/* L145: */
    }
L146:
    if (*mbdcnd != 0) {
	goto L148;
    }
    i__1 = nstop;
    for (j = nstart; j <= i__1; ++j) {
	f[mp1 + j * f_dim1] = f[j * f_dim1 + 1];
/* L147: */
    }
    if (*nbdcnd == 0) {
	f[mp1 + np1 * f_dim1] = f[np1 * f_dim1 + 1];
    }
L148:
    return 0;
} /* hwscrt_ */

