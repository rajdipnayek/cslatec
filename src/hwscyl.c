/* hwscyl.f -- translated by f2c (version 12.02.01).
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

/* DECK HWSCYL */
/* Subroutine */ int hwscyl_(real *a, real *b, integer *m, integer *mbdcnd, 
	real *bda, real *bdb, real *c__, real *d__, integer *n, integer *
	nbdcnd, real *bdc, real *bdd, real *elmbda, real *f, integer *idimf, 
	real *pertrb, integer *ierror, real *w)
{
    /* System generated locals */
    integer f_dim1, f_offset, i__1, i__2;
    real r__1;

    /* Local variables */
    static integer i__, j, k, l;
    static real r__, s, a1, a2, s1, s2;
    static integer ij, np, id2, id3, id4, id5, id6, mp1, np1, nsp1, munk, 
	    nunk, ierr1, nstm1;
    static real dlrsq;
    static integer mstop, nstop;
    static real dlrby2, deltar;
    extern /* Subroutine */ int genbun_(integer *, integer *, integer *, 
	    integer *, real *, real *, real *, integer *, real *, integer *, 
	    real *);
    static real deltht, dlthsq;
    static integer istart, mstart, nstart;

/* ***BEGIN PROLOGUE  HWSCYL */
/* ***PURPOSE  Solve a standard finite difference approximation */
/*            to the Helmholtz equation in cylindrical coordinates. */
/* ***LIBRARY   SLATEC (FISHPACK) */
/* ***CATEGORY  I2B1A1A */
/* ***TYPE      SINGLE PRECISION (HWSCYL-S) */
/* ***KEYWORDS  CYLINDRICAL, ELLIPTIC, FISHPACK, HELMHOLTZ, PDE */
/* ***AUTHOR  Adams, J., (NCAR) */
/*           Swarztrauber, P. N., (NCAR) */
/*           Sweet, R., (NCAR) */
/* ***DESCRIPTION */

/*     Subroutine HWSCYL solves a finite difference approximation to the */
/*     Helmholtz equation in cylindrical coordinates: */

/*          (1/R)(d/dR)(R(dU/dR)) + (d/dZ)(dU/dZ) */

/*                                + (LAMBDA/R**2)U = F(R,Z) */

/*     This modified Helmholtz equation results from the Fourier */
/*     transform of the three-dimensional Poisson equation. */

/*     * * * * * * * *    Parameter Description     * * * * * * * * * * */

/*             * * * * * *   On Input    * * * * * * */

/*     A,B */
/*       The range of R, i.e., A .LE. R .LE. B.  A must be less than B */
/*       and A must be non-negative. */

/*     M */
/*       The number of panels into which the interval (A,B) is */
/*       subdivided.  Hence, there will be M+1 grid points in the */
/*       R-direction given by R(I) = A+(I-1)DR, for I = 1,2,...,M+1, */
/*       where DR = (B-A)/M is the panel width. M must be greater than 3. */

/*     MBDCND */
/*       Indicates the type of boundary conditions at R = A and R = B. */

/*       = 1  If the solution is specified at R = A and R = B. */
/*       = 2  If the solution is specified at R = A and the derivative of */
/*            the solution with respect to R is specified at R = B. */
/*       = 3  If the derivative of the solution with respect to R is */
/*            specified at R = A (see note below) and R = B. */
/*       = 4  If the derivative of the solution with respect to R is */
/*            specified at R = A (see note below) and the solution is */
/*            specified at R = B. */
/*       = 5  If the solution is unspecified at R = A = 0 and the */
/*            solution is specified at R = B. */
/*       = 6  If the solution is unspecified at R = A = 0 and the */
/*            derivative of the solution with respect to R is specified */
/*            at R = B. */

/*       NOTE:  If A = 0, do not use MBDCND = 3 or 4, but instead use */
/*              MBDCND = 1,2,5, or 6  . */

/*     BDA */
/*       A one-dimensional array of length N+1 that specifies the values */
/*       of the derivative of the solution with respect to R at R = A. */
/*       When MBDCND = 3 or 4, */

/*            BDA(J) = (d/dR)U(A,Z(J)), J = 1,2,...,N+1  . */

/*       When MBDCND has any other value, BDA is a dummy variable. */

/*     BDB */
/*       A one-dimensional array of length N+1 that specifies the values */
/*       of the derivative of the solution with respect to R at R = B. */
/*       When MBDCND = 2,3, or 6, */

/*            BDB(J) = (d/dR)U(B,Z(J)), J = 1,2,...,N+1  . */

/*       When MBDCND has any other value, BDB is a dummy variable. */

/*     C,D */
/*       The range of Z, i.e., C .LE. Z .LE. D.  C must be less than D. */

/*     N */
/*       The number of panels into which the interval (C,D) is */
/*       subdivided.  Hence, there will be N+1 grid points in the */
/*       Z-direction given by Z(J) = C+(J-1)DZ, for J = 1,2,...,N+1, */
/*       where DZ = (D-C)/N is the panel width. N must be greater than 3. */

/*     NBDCND */
/*       Indicates the type of boundary conditions at Z = C and Z = D. */

/*       = 0  If the solution is periodic in Z, i.e., U(I,1) = U(I,N+1). */
/*       = 1  If the solution is specified at Z = C and Z = D. */
/*       = 2  If the solution is specified at Z = C and the derivative of */
/*            the solution with respect to Z is specified at Z = D. */
/*       = 3  If the derivative of the solution with respect to Z is */
/*            specified at Z = C and Z = D. */
/*       = 4  If the derivative of the solution with respect to Z is */
/*            specified at Z = C and the solution is specified at Z = D. */

/*     BDC */
/*       A one-dimensional array of length M+1 that specifies the values */
/*       of the derivative of the solution with respect to Z at Z = C. */
/*       When NBDCND = 3 or 4, */

/*            BDC(I) = (d/dZ)U(R(I),C), I = 1,2,...,M+1  . */

/*       When NBDCND has any other value, BDC is a dummy variable. */

/*     BDD */
/*       A one-dimensional array of length M+1 that specifies the values */
/*       of the derivative of the solution with respect to Z at Z = D. */
/*       When NBDCND = 2 or 3, */

/*            BDD(I) = (d/dZ)U(R(I),D), I = 1,2,...,M+1  . */

/*       When NBDCND has any other value, BDD is a dummy variable. */

/*     ELMBDA */
/*       The constant LAMBDA in the Helmholtz equation.  If */
/*       LAMBDA .GT. 0, a solution may not exist.  However, HWSCYL will */
/*       attempt to find a solution.  LAMBDA must be zero when */
/*       MBDCND = 5 or 6  . */

/*     F */
/*       A two-dimensional array that specifies the values of the right */
/*       side of the Helmholtz equation and boundary data (if any).  For */
/*       I = 2,3,...,M and J = 2,3,...,N */

/*            F(I,J) = F(R(I),Z(J)). */

/*       On the boundaries F is defined by */

/*            MBDCND   F(1,J)            F(M+1,J) */
/*            ------   ---------         --------- */

/*              1      U(A,Z(J))         U(B,Z(J)) */
/*              2      U(A,Z(J))         F(B,Z(J)) */
/*              3      F(A,Z(J))         F(B,Z(J))   J = 1,2,...,N+1 */
/*              4      F(A,Z(J))         U(B,Z(J)) */
/*              5      F(0,Z(J))         U(B,Z(J)) */
/*              6      F(0,Z(J))         F(B,Z(J)) */

/*            NBDCND   F(I,1)            F(I,N+1) */
/*            ------   ---------         --------- */

/*              0      F(R(I),C)         F(R(I),C) */
/*              1      U(R(I),C)         U(R(I),D) */
/*              2      U(R(I),C)         F(R(I),D)   I = 1,2,...,M+1 */
/*              3      F(R(I),C)         F(R(I),D) */
/*              4      F(R(I),C)         U(R(I),D) */

/*       F must be dimensioned at least (M+1)*(N+1). */

/*       NOTE */

/*       If the table calls for both the solution U and the right side F */
/*       at a corner then the solution must be specified. */

/*     IDIMF */
/*       The row (or first) dimension of the array F as it appears in the */
/*       program calling HWSCYL.  This parameter is used to specify the */
/*       variable dimension of F.  IDIMF must be at least M+1  . */

/*     W */
/*       A one-dimensional array that must be provided by the user for */
/*       work space.  W may require up to 4*(N+1) + */
/*       (13 + INT(log2(N+1)))*(M+1) locations.  The actual number of */
/*       locations used is computed by HWSCYL and is returned in location */
/*       W(1). */


/*             * * * * * *   On Output     * * * * * * */

/*     F */
/*       Contains the solution U(I,J) of the finite difference */
/*       approximation for the grid point (R(I),Z(J)), I = 1,2,...,M+1, */
/*       J = 1,2,...,N+1  . */

/*     PERTRB */
/*       If one specifies a combination of periodic, derivative, and */
/*       unspecified boundary conditions for a Poisson equation */
/*       (LAMBDA = 0), a solution may not exist.  PERTRB is a constant, */
/*       calculated and subtracted from F, which ensures that a solution */
/*       exists.  HWSCYL then computes this solution, which is a least */
/*       squares solution to the original approximation.  This solution */
/*       plus any constant is also a solution.  Hence, the solution is */
/*       not unique.  The value of PERTRB should be small compared to the */
/*       right side F.  Otherwise, a solution is obtained to an */
/*       essentially different problem.  This comparison should always */
/*       be made to insure that a meaningful solution has been obtained. */

/*     IERROR */
/*       An error flag which indicates invalid input parameters.  Except */
/*       for numbers 0 and 11, a solution is not attempted. */

/*       =  0  No error. */
/*       =  1  A .LT. 0  . */
/*       =  2  A .GE. B. */
/*       =  3  MBDCND .LT. 1 or MBDCND .GT. 6  . */
/*       =  4  C .GE. D. */
/*       =  5  N .LE. 3 */
/*       =  6  NBDCND .LT. 0 or NBDCND .GT. 4  . */
/*       =  7  A = 0, MBDCND = 3 or 4  . */
/*       =  8  A .GT. 0, MBDCND .GE. 5  . */
/*       =  9  A = 0, LAMBDA .NE. 0, MBDCND .GE. 5  . */
/*       = 10  IDIMF .LT. M+1  . */
/*       = 11  LAMBDA .GT. 0  . */
/*       = 12  M .LE. 3 */

/*       Since this is the only means of indicating a possibly incorrect */
/*       call to HWSCYL, the user should test IERROR after the call. */

/*     W */
/*       W(1) contains the required length of W. */

/* *Long Description: */

/*     * * * * * * *   Program Specifications    * * * * * * * * * * * * */

/*     Dimension of   BDA(N+1),BDB(N+1),BDC(M+1),BDD(M+1),F(IDIMF,N+1), */
/*     Arguments      W(see argument list) */

/*     Latest         June 1, 1976 */
/*     Revision */

/*     Subprograms    HWSCYL,GENBUN,POISD2,POISN2,POISP2,COSGEN,MERGE, */
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

/*     Space          5818(decimal) = 13272(octal) locations on the NCAR */
/*     Required       Control Data 7600 */

/*     Timing and        The execution time T on the NCAR Control Data */
/*     Accuracy       7600 for subroutine HWSCYL is roughly proportional */
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

/*                        32        1         0          31 */
/*                        32        1         1          23 */
/*                        32        3         3          36 */
/*                        64        1         0         128 */
/*                        64        1         1          96 */
/*                        64        3         3         142 */

/*     Portability    American National Standards Institute FORTRAN. */
/*                    The machine dependent constant PI is defined in */
/*                    function PIMACH. */

/*     Required       COS */
/*     Resident */
/*     Routines */

/*     Reference      Swarztrauber, P. and R. Sweet, 'Efficient FORTRAN */
/*                    Subprograms for the Solution of Elliptic Equations' */
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
/* ***END PROLOGUE  HWSCYL */


/* ***FIRST EXECUTABLE STATEMENT  HWSCYL */
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
    if (*n <= 3) {
	*ierror = 5;
    }
    if (*nbdcnd <= -1 || *nbdcnd >= 5) {
	*ierror = 6;
    }
    if (*a == 0.f && (*mbdcnd == 3 || *mbdcnd == 4)) {
	*ierror = 7;
    }
    if (*a > 0.f && *mbdcnd >= 5) {
	*ierror = 8;
    }
    if (*a == 0.f && *elmbda != 0.f && *mbdcnd >= 5) {
	*ierror = 9;
    }
    if (*idimf < *m + 1) {
	*ierror = 10;
    }
    if (*m <= 3) {
	*ierror = 12;
    }
    if (*ierror != 0) {
	return 0;
    }
    mp1 = *m + 1;
    deltar = (*b - *a) / *m;
    dlrby2 = deltar / 2.f;
/* Computing 2nd power */
    r__1 = deltar;
    dlrsq = r__1 * r__1;
    np1 = *n + 1;
    deltht = (*d__ - *c__) / *n;
/* Computing 2nd power */
    r__1 = deltht;
    dlthsq = r__1 * r__1;
    np = *nbdcnd + 1;

/*     DEFINE RANGE OF INDICES I AND J FOR UNKNOWNS U(I,J). */

    mstart = 2;
    mstop = *m;
    switch (*mbdcnd) {
	case 1:  goto L104;
	case 2:  goto L103;
	case 3:  goto L102;
	case 4:  goto L101;
	case 5:  goto L101;
	case 6:  goto L102;
    }
L101:
    mstart = 1;
    goto L104;
L102:
    mstart = 1;
L103:
    mstop = mp1;
L104:
    munk = mstop - mstart + 1;
    nstart = 1;
    nstop = *n;
    switch (np) {
	case 1:  goto L108;
	case 2:  goto L105;
	case 3:  goto L106;
	case 4:  goto L107;
	case 5:  goto L108;
    }
L105:
    nstart = 2;
    goto L108;
L106:
    nstart = 2;
L107:
    nstop = np1;
L108:
    nunk = nstop - nstart + 1;

/*     DEFINE A,B,C COEFFICIENTS IN W-ARRAY. */

    id2 = munk;
    id3 = id2 + munk;
    id4 = id3 + munk;
    id5 = id4 + munk;
    id6 = id5 + munk;
    istart = 1;
    a1 = 2.f / dlrsq;
    ij = 0;
    if (*mbdcnd == 3 || *mbdcnd == 4) {
	ij = 1;
    }
    if (*mbdcnd <= 4) {
	goto L109;
    }
    w[1] = 0.f;
    w[id2 + 1] = a1 * -2.f;
    w[id3 + 1] = a1 * 2.f;
    istart = 2;
    ij = 1;
L109:
    i__1 = munk;
    for (i__ = istart; i__ <= i__1; ++i__) {
	r__ = *a + (i__ - ij) * deltar;
	j = id5 + i__;
	w[j] = r__;
	j = id6 + i__;
/* Computing 2nd power */
	r__1 = r__;
	w[j] = 1.f / (r__1 * r__1);
	w[i__] = (r__ - dlrby2) / (r__ * dlrsq);
	j = id3 + i__;
	w[j] = (r__ + dlrby2) / (r__ * dlrsq);
	k = id6 + i__;
	j = id2 + i__;
	w[j] = -a1 + *elmbda * w[k];
/* L110: */
    }
    switch (*mbdcnd) {
	case 1:  goto L114;
	case 2:  goto L111;
	case 3:  goto L112;
	case 4:  goto L113;
	case 5:  goto L114;
	case 6:  goto L112;
    }
L111:
    w[id2] = a1;
    goto L114;
L112:
    w[id2] = a1;
L113:
    w[id3 + 1] = a1 * istart;
L114:

/*     ENTER BOUNDARY DATA FOR R-BOUNDARIES. */

    switch (*mbdcnd) {
	case 1:  goto L115;
	case 2:  goto L115;
	case 3:  goto L117;
	case 4:  goto L117;
	case 5:  goto L119;
	case 6:  goto L119;
    }
L115:
    a1 = w[1];
    i__1 = nstop;
    for (j = nstart; j <= i__1; ++j) {
	f[j * f_dim1 + 2] -= a1 * f[j * f_dim1 + 1];
/* L116: */
    }
    goto L119;
L117:
    a1 = deltar * 2.f * w[1];
    i__1 = nstop;
    for (j = nstart; j <= i__1; ++j) {
	f[j * f_dim1 + 1] += a1 * bda[j];
/* L118: */
    }
L119:
    switch (*mbdcnd) {
	case 1:  goto L120;
	case 2:  goto L122;
	case 3:  goto L122;
	case 4:  goto L120;
	case 5:  goto L120;
	case 6:  goto L122;
    }
L120:
    a1 = w[id4];
    i__1 = nstop;
    for (j = nstart; j <= i__1; ++j) {
	f[*m + j * f_dim1] -= a1 * f[mp1 + j * f_dim1];
/* L121: */
    }
    goto L124;
L122:
    a1 = deltar * 2.f * w[id4];
    i__1 = nstop;
    for (j = nstart; j <= i__1; ++j) {
	f[mp1 + j * f_dim1] -= a1 * bdb[j];
/* L123: */
    }

/*     ENTER BOUNDARY DATA FOR Z-BOUNDARIES. */

L124:
    a1 = 1.f / dlthsq;
    l = id5 - mstart + 1;
    switch (np) {
	case 1:  goto L134;
	case 2:  goto L125;
	case 3:  goto L125;
	case 4:  goto L127;
	case 5:  goto L127;
    }
L125:
    i__1 = mstop;
    for (i__ = mstart; i__ <= i__1; ++i__) {
	f[i__ + (f_dim1 << 1)] -= a1 * f[i__ + f_dim1];
/* L126: */
    }
    goto L129;
L127:
    a1 = 2.f / deltht;
    i__1 = mstop;
    for (i__ = mstart; i__ <= i__1; ++i__) {
	f[i__ + f_dim1] += a1 * bdc[i__];
/* L128: */
    }
L129:
    a1 = 1.f / dlthsq;
    switch (np) {
	case 1:  goto L134;
	case 2:  goto L130;
	case 3:  goto L132;
	case 4:  goto L132;
	case 5:  goto L130;
    }
L130:
    i__1 = mstop;
    for (i__ = mstart; i__ <= i__1; ++i__) {
	f[i__ + *n * f_dim1] -= a1 * f[i__ + np1 * f_dim1];
/* L131: */
    }
    goto L134;
L132:
    a1 = 2.f / deltht;
    i__1 = mstop;
    for (i__ = mstart; i__ <= i__1; ++i__) {
	f[i__ + np1 * f_dim1] -= a1 * bdd[i__];
/* L133: */
    }
L134:

/*     ADJUST RIGHT SIDE OF SINGULAR PROBLEMS TO INSURE EXISTENCE OF A */
/*     SOLUTION. */

    *pertrb = 0.f;
    if (*elmbda < 0.f) {
	goto L146;
    } else if (*elmbda == 0) {
	goto L136;
    } else {
	goto L135;
    }
L135:
    *ierror = 11;
    goto L146;
L136:
    w[id5 + 1] = (w[id5 + 2] - dlrby2) * .5f;
    switch (*mbdcnd) {
	case 1:  goto L146;
	case 2:  goto L146;
	case 3:  goto L138;
	case 4:  goto L146;
	case 5:  goto L146;
	case 6:  goto L137;
    }
L137:
    w[id5 + 1] *= .5f;
L138:
    switch (np) {
	case 1:  goto L140;
	case 2:  goto L146;
	case 3:  goto L146;
	case 4:  goto L139;
	case 5:  goto L146;
    }
L139:
    a2 = 2.f;
    goto L141;
L140:
    a2 = 1.f;
L141:
    k = id5 + munk;
    w[k] = (w[k - 1] + dlrby2) * .5f;
    s = 0.f;
    i__1 = mstop;
    for (i__ = mstart; i__ <= i__1; ++i__) {
	s1 = 0.f;
	nsp1 = nstart + 1;
	nstm1 = nstop - 1;
	i__2 = nstm1;
	for (j = nsp1; j <= i__2; ++j) {
	    s1 += f[i__ + j * f_dim1];
/* L142: */
	}
	k = i__ + l;
	s += (a2 * s1 + f[i__ + nstart * f_dim1] + f[i__ + nstop * f_dim1]) * 
		w[k];
/* L143: */
    }
    s2 = *m * *a + ((*m - 1) * (*m + 1) + .75f) * dlrby2;
    if (*mbdcnd == 3) {
	s2 += dlrby2 * .25f;
    }
    s1 = (a2 * (nunk - 2) + 2.f) * s2;
    *pertrb = s / s1;
    i__1 = mstop;
    for (i__ = mstart; i__ <= i__1; ++i__) {
	i__2 = nstop;
	for (j = nstart; j <= i__2; ++j) {
	    f[i__ + j * f_dim1] -= *pertrb;
/* L144: */
	}
/* L145: */
    }
L146:

/*     MULTIPLY I-TH EQUATION THROUGH BY DELTHT**2 TO PUT EQUATION INTO */
/*     CORRECT FORM FOR SUBROUTINE GENBUN. */

    i__1 = mstop;
    for (i__ = mstart; i__ <= i__1; ++i__) {
	k = i__ - mstart + 1;
	w[k] *= dlthsq;
	j = id2 + k;
	w[j] *= dlthsq;
	j = id3 + k;
	w[j] *= dlthsq;
	i__2 = nstop;
	for (j = nstart; j <= i__2; ++j) {
	    f[i__ + j * f_dim1] *= dlthsq;
/* L147: */
	}
/* L148: */
    }
    w[1] = 0.f;
    w[id4] = 0.f;

/*     CALL GENBUN TO SOLVE THE SYSTEM OF EQUATIONS. */

    genbun_(nbdcnd, &nunk, &c__1, &munk, &w[1], &w[id2 + 1], &w[id3 + 1], 
	    idimf, &f[mstart + nstart * f_dim1], &ierr1, &w[id4 + 1]);
    w[1] = w[id4 + 1] + munk * 3;
    if (*nbdcnd != 0) {
	goto L150;
    }
    i__1 = mstop;
    for (i__ = mstart; i__ <= i__1; ++i__) {
	f[i__ + np1 * f_dim1] = f[i__ + f_dim1];
/* L149: */
    }
L150:
    return 0;
} /* hwscyl_ */

