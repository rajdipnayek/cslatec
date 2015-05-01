/* hwsplr.f -- translated by f2c (version 12.02.01).
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

/* DECK HWSPLR */
/* Subroutine */ int hwsplr_(real *a, real *b, integer *m, integer *mbdcnd, 
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
    static integer ij, ip, lp, np, id2, id3, id4, id5, id6, mp1, np1, munk, 
	    nunk, ierr1;
    static real dlrsq, ypole;
    static integer mstop, nstop;
    static real dlrby2, deltar;
    extern /* Subroutine */ int genbun_(integer *, integer *, integer *, 
	    integer *, real *, real *, real *, integer *, real *, integer *, 
	    real *);
    static real deltht, dlthsq;
    static integer mstart, nstart, iwstor;

/* ***BEGIN PROLOGUE  HWSPLR */
/* ***PURPOSE  Solve a finite difference approximation to the Helmholtz */
/*            equation in polar coordinates. */
/* ***LIBRARY   SLATEC (FISHPACK) */
/* ***CATEGORY  I2B1A1A */
/* ***TYPE      SINGLE PRECISION (HWSPLR-S) */
/* ***KEYWORDS  ELLIPTIC, FISHPACK, HELMHOLTZ, PDE, POLAR */
/* ***AUTHOR  Adams, J., (NCAR) */
/*           Swarztrauber, P. N., (NCAR) */
/*           Sweet, R., (NCAR) */
/* ***DESCRIPTION */

/*     Subroutine HWSPLR solves a finite difference approximation to the */
/*     Helmholtz equation in polar coordinates: */

/*          (1/R)(d/dR)(R(dU/dR)) + (1/R**2)(d/dTHETA)(dU/dTHETA) */

/*                                + LAMBDA*U = F(R,THETA). */




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
/*       Indicates the type of boundary condition at R = A and R = B. */

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

/*            BDA(J) = (d/dR)U(A,THETA(J)), J = 1,2,...,N+1  . */

/*       When MBDCND has any other value, BDA is a dummy variable. */

/*     BDB */
/*       A one-dimensional array of length N+1 that specifies the values */
/*       of the derivative of the solution with respect to R at R = B. */
/*       When MBDCND = 2,3, or 6, */

/*            BDB(J) = (d/dR)U(B,THETA(J)), J = 1,2,...,N+1  . */

/*       When MBDCND has any other value, BDB is a dummy variable. */

/*     C,D */
/*       The range of THETA, i.e., C .LE. THETA .LE. D.  C must be less */
/*       than D. */

/*     N */
/*       The number of panels into which the interval (C,D) is */
/*       subdivided.  Hence, there will be N+1 grid points in the */
/*       THETA-direction given by THETA(J) = C+(J-1)DTHETA for */
/*       J = 1,2,...,N+1, where DTHETA = (D-C)/N is the panel width.  N */
/*       must be greater than 3. */

/*     NBDCND */
/*       Indicates the type of boundary conditions at THETA = C and */
/*       at THETA = D. */

/*       = 0  If the solution is periodic in THETA, i.e., */
/*            U(I,J) = U(I,N+J). */
/*       = 1  If the solution is specified at THETA = C and THETA = D */
/*            (see note below). */
/*       = 2  If the solution is specified at THETA = C and the */
/*            derivative of the solution with respect to THETA is */
/*            specified at THETA = D (see note below). */
/*       = 4  If the derivative of the solution with respect to THETA is */
/*            specified at THETA = C and the solution is specified at */
/*            THETA = D (see note below). */

/*       NOTE:  When NBDCND = 1,2, or 4, do not use MBDCND = 5 or 6 */
/*              (the former indicates that the solution is specified at */
/*              R = 0, the latter indicates the solution is unspecified */
/*              at R = 0).  Use instead MBDCND = 1 or 2  . */

/*     BDC */
/*       A one-dimensional array of length M+1 that specifies the values */
/*       of the derivative of the solution with respect to THETA at */
/*       THETA = C.  When NBDCND = 3 or 4, */

/*            BDC(I) = (d/dTHETA)U(R(I),C), I = 1,2,...,M+1  . */

/*       When NBDCND has any other value, BDC is a dummy variable. */

/*     BDD */
/*       A one-dimensional array of length M+1 that specifies the values */
/*       of the derivative of the solution with respect to THETA at */
/*       THETA = D.  When NBDCND = 2 or 3, */

/*            BDD(I) = (d/dTHETA)U(R(I),D), I = 1,2,...,M+1  . */

/*       When NBDCND has any other value, BDD is a dummy variable. */

/*     ELMBDA */
/*       The constant LAMBDA in the Helmholtz equation.  If */
/*       LAMBDA .LT. 0, a solution may not exist.  However, HWSPLR will */
/*       attempt to find a solution. */

/*     F */
/*       A two-dimensional array that specifies the values of the right */
/*       side of the Helmholtz equation and boundary values (if any). */
/*       For I = 2,3,...,M and J = 2,3,...,N */

/*            F(I,J) = F(R(I),THETA(J)). */

/*       On the boundaries F is defined by */

/*            MBDCND   F(1,J)            F(M+1,J) */
/*            ------   -------------     ------------- */

/*              1      U(A,THETA(J))     U(B,THETA(J)) */
/*              2      U(A,THETA(J))     F(B,THETA(J)) */
/*              3      F(A,THETA(J))     F(B,THETA(J)) */
/*              4      F(A,THETA(J))     U(B,THETA(J))   J = 1,2,...,N+1 */
/*              5      F(0,0)            U(B,THETA(J)) */
/*              6      F(0,0)            F(B,THETA(J)) */

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
/*       program calling HWSPLR.  This parameter is used to specify the */
/*       variable dimension of F.  IDIMF must be at least M+1  . */

/*     W */
/*       A one-dimensional array that must be provided by the user for */
/*       work space.  W may require up to 4*(N+1) + */
/*       (13 + INT(log2(N+1)))*(M+1) locations.  The actual number of */
/*       locations used is computed by HWSPLR and is returned in location */
/*       W(1). */


/*             * * * * * *   On Output     * * * * * * */

/*     F */
/*       Contains the solution U(I,J) of the finite difference */
/*       approximation for the grid point (R(I),THETA(J)), */
/*       I = 1,2,...,M+1, J = 1,2,...,N+1  . */

/*     PERTRB */
/*       If a combination of periodic, derivative, or unspecified */
/*       boundary conditions is specified for a Poisson equation */
/*       (LAMBDA = 0), a solution may not exist.  PERTRB is a constant, */
/*       calculated and subtracted from F, which ensures that a solution */
/*       exists.  HWSPLR then computes this solution, which is a least */
/*       squares solution to the original approximation.  This solution */
/*       plus any constant is also a solution.  Hence, the solution is */
/*       not unique.  PERTRB should be small compared to the right side. */
/*       Otherwise, a solution is obtained to an essentially different */
/*       problem.  This comparison should always be made to insure that a */
/*       meaningful solution has been obtained. */

/*     IERROR */
/*       An error flag that indicates invalid input parameters.  Except */
/*       for numbers 0 and 11, a solution is not attempted. */

/*       =  0  No error. */
/*       =  1  A .LT. 0  . */
/*       =  2  A .GE. B. */
/*       =  3  MBDCND .LT. 1 or MBDCND .GT. 6  . */
/*       =  4  C .GE. D. */
/*       =  5  N .LE. 3 */
/*       =  6  NBDCND .LT. 0 or .GT. 4  . */
/*       =  7  A = 0, MBDCND = 3 or 4  . */
/*       =  8  A .GT. 0, MBDCND .GE. 5  . */
/*       =  9  MBDCND .GE. 5, NBDCND .NE. 0 and NBDCND .NE. 3  . */
/*       = 10  IDIMF .LT. M+1  . */
/*       = 11  LAMBDA .GT. 0  . */
/*       = 12  M .LE. 3 */

/*       Since this is the only means of indicating a possibly incorrect */
/*       call to HWSPLR, the user should test IERROR after the call. */

/*     W */
/*       W(1) contains the required length of W. */

/* *Long Description: */

/*     * * * * * * *   Program Specifications    * * * * * * * * * * * * */

/*     Dimension of   BDA(N+1),BDB(N+1),BDC(M+1),BDD(M+1),F(IDIMF,N+1), */
/*     Arguments      W(see argument list) */

/*     Latest         June 1, 1976 */
/*     Revision */

/*     Subprograms    HWSPLR,GENBUN,POISD2,POISN2,POISP2,COSGEN,MERGE, */
/*     Required       TRIX,TRI3,PIMACH */

/*     Special        None */
/*     Conditions */

/*     Common         NONE */
/*     Blocks */

/*     I/O */

/*     Precision      Single */

/*     Specialist     Roland Sweet */

/*     Language       FORTRAN */

/*     History        Standardized April 1, 1973 */
/*                    Revised January 1, 1976 */

/*     Algorithm      The routine defines the finite difference */
/*                    equations, incorporates boundary data, and adjusts */
/*                    the right side of singular systems and then calls */
/*                    GENBUN to solve the system. */

/*     Space          13430(octal) = 5912(decimal)  locations on the NCAR */
/*     Required       Control Data 7600 */

/*     Timing and        The execution time T on the NCAR Control Data */
/*     Accuracy       7600 for subroutine HWSPLR is roughly proportional */
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
/*                    Subprograms For The Solution Of Elliptic Equations' */
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
/* ***END PROLOGUE  HWSPLR */


/* ***FIRST EXECUTABLE STATEMENT  HWSPLR */
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
    if (*mbdcnd >= 5 && *nbdcnd != 0 && *nbdcnd != 3) {
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
    mstop = mp1;
    switch (*mbdcnd) {
	case 1:  goto L101;
	case 2:  goto L105;
	case 3:  goto L102;
	case 4:  goto L103;
	case 5:  goto L104;
	case 6:  goto L105;
    }
L101:
    mstop = *m;
    goto L105;
L102:
    mstart = 1;
    goto L105;
L103:
    mstart = 1;
L104:
    mstop = *m;
L105:
    munk = mstop - mstart + 1;
    nstart = 1;
    nstop = *n;
    switch (np) {
	case 1:  goto L109;
	case 2:  goto L106;
	case 3:  goto L107;
	case 4:  goto L108;
	case 5:  goto L109;
    }
L106:
    nstart = 2;
    goto L109;
L107:
    nstart = 2;
L108:
    nstop = np1;
L109:
    nunk = nstop - nstart + 1;

/*     DEFINE A,B,C COEFFICIENTS IN W-ARRAY. */

    id2 = munk;
    id3 = id2 + munk;
    id4 = id3 + munk;
    id5 = id4 + munk;
    id6 = id5 + munk;
    a1 = 2.f / dlrsq;
    ij = 0;
    if (*mbdcnd == 3 || *mbdcnd == 4) {
	ij = 1;
    }
    i__1 = munk;
    for (i__ = 1; i__ <= i__1; ++i__) {
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
	j = id2 + i__;
	w[j] = -a1 + *elmbda;
/* L110: */
    }
    switch (*mbdcnd) {
	case 1:  goto L114;
	case 2:  goto L111;
	case 3:  goto L112;
	case 4:  goto L113;
	case 5:  goto L114;
	case 6:  goto L111;
    }
L111:
    w[id2] = a1;
    goto L114;
L112:
    w[id2] = a1;
L113:
    w[id3 + 1] = a1;
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

/*     ENTER BOUNDARY DATA FOR THETA-BOUNDARIES. */

L124:
    a1 = 1.f / dlthsq;
    l = id5 - mstart + 1;
    lp = id6 - mstart + 1;
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
	j = i__ + lp;
	f[i__ + (f_dim1 << 1)] -= a1 * w[j] * f[i__ + f_dim1];
/* L126: */
    }
    goto L129;
L127:
    a1 = 2.f / deltht;
    i__1 = mstop;
    for (i__ = mstart; i__ <= i__1; ++i__) {
	j = i__ + lp;
	f[i__ + f_dim1] += a1 * w[j] * bdc[i__];
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
	j = i__ + lp;
	f[i__ + *n * f_dim1] -= a1 * w[j] * f[i__ + np1 * f_dim1];
/* L131: */
    }
    goto L134;
L132:
    a1 = 2.f / deltht;
    i__1 = mstop;
    for (i__ = mstart; i__ <= i__1; ++i__) {
	j = i__ + lp;
	f[i__ + np1 * f_dim1] -= a1 * w[j] * bdd[i__];
/* L133: */
    }
L134:

/*     ADJUST RIGHT SIDE OF EQUATION FOR UNKNOWN AT POLE WHEN HAVE */
/*     DERIVATIVE SPECIFIED BOUNDARY CONDITIONS. */

    if (*mbdcnd >= 5 && *nbdcnd == 3) {
	f[f_dim1 + 1] -= (bdd[2] - bdc[2]) * 4.f / (*n * deltht * dlrsq);
    }

/*     ADJUST RIGHT SIDE OF SINGULAR PROBLEMS TO INSURE EXISTENCE OF A */
/*     SOLUTION. */

    *pertrb = 0.f;
    if (*elmbda < 0.f) {
	goto L144;
    } else if (*elmbda == 0) {
	goto L136;
    } else {
	goto L135;
    }
L135:
    *ierror = 11;
    goto L144;
L136:
    if (*nbdcnd != 0 && *nbdcnd != 3) {
	goto L144;
    }
    s2 = 0.f;
    switch (*mbdcnd) {
	case 1:  goto L144;
	case 2:  goto L144;
	case 3:  goto L137;
	case 4:  goto L144;
	case 5:  goto L144;
	case 6:  goto L138;
    }
L137:
    w[id5 + 1] = (w[id5 + 2] - dlrby2) * .5f;
    s2 = deltar * .25f;
L138:
    a2 = 2.f;
    if (*nbdcnd == 0) {
	a2 = 1.f;
    }
    j = id5 + munk;
    w[j] = (w[j - 1] + dlrby2) * .5f;
    s = 0.f;
    i__1 = mstop;
    for (i__ = mstart; i__ <= i__1; ++i__) {
	s1 = 0.f;
	ij = nstart + 1;
	k = nstop - 1;
	i__2 = k;
	for (j = ij; j <= i__2; ++j) {
	    s1 += f[i__ + j * f_dim1];
/* L139: */
	}
	j = i__ + l;
	s += (a2 * s1 + f[i__ + nstart * f_dim1] + f[i__ + nstop * f_dim1]) * 
		w[j];
/* L140: */
    }
    s2 = *m * *a + deltar * ((*m - 1) * (*m + 1) * .5f + .25f) + s2;
    s1 = (a2 * (nunk - 2) + 2.f) * s2;
    if (*mbdcnd == 3) {
	goto L141;
    }
    s2 = *n * a2 * deltar / 8.f;
    s += f[f_dim1 + 1] * s2;
    s1 += s2;
L141:
    *pertrb = s / s1;
    i__1 = mstop;
    for (i__ = mstart; i__ <= i__1; ++i__) {
	i__2 = nstop;
	for (j = nstart; j <= i__2; ++j) {
	    f[i__ + j * f_dim1] -= *pertrb;
/* L142: */
	}
/* L143: */
    }
L144:

/*     MULTIPLY I-TH EQUATION THROUGH BY (R(I)*DELTHT)**2. */

    i__1 = mstop;
    for (i__ = mstart; i__ <= i__1; ++i__) {
	k = i__ - mstart + 1;
	j = i__ + lp;
	a1 = dlthsq / w[j];
	w[k] = a1 * w[k];
	j = id2 + k;
	w[j] = a1 * w[j];
	j = id3 + k;
	w[j] = a1 * w[j];
	i__2 = nstop;
	for (j = nstart; j <= i__2; ++j) {
	    f[i__ + j * f_dim1] = a1 * f[i__ + j * f_dim1];
/* L145: */
	}
/* L146: */
    }
    w[1] = 0.f;
    w[id4] = 0.f;

/*     CALL GENBUN TO SOLVE THE SYSTEM OF EQUATIONS. */

    genbun_(nbdcnd, &nunk, &c__1, &munk, &w[1], &w[id2 + 1], &w[id3 + 1], 
	    idimf, &f[mstart + nstart * f_dim1], &ierr1, &w[id4 + 1]);
    iwstor = w[id4 + 1] + munk * 3;
    switch (*mbdcnd) {
	case 1:  goto L157;
	case 2:  goto L157;
	case 3:  goto L157;
	case 4:  goto L157;
	case 5:  goto L148;
	case 6:  goto L147;
    }

/*     ADJUST THE SOLUTION AS NECESSARY FOR THE PROBLEMS WHERE A = 0. */

L147:
    if (*elmbda != 0.f) {
	goto L148;
    }
    ypole = 0.f;
    goto L155;
L148:
    j = id5 + munk;
    w[j] = w[id2] / w[id3];
    i__1 = munk;
    for (ip = 3; ip <= i__1; ++ip) {
	i__ = munk - ip + 2;
	j = id5 + i__;
	lp = id2 + i__;
	k = id3 + i__;
	w[j] = w[i__] / (w[lp] - w[k] * w[j + 1]);
/* L149: */
    }
    w[id5 + 1] = dlthsq * -.5f / (w[id2 + 1] - w[id3 + 1] * w[id5 + 2]);
    i__1 = munk;
    for (i__ = 2; i__ <= i__1; ++i__) {
	j = id5 + i__;
	w[j] = -w[j] * w[j - 1];
/* L150: */
    }
    s = 0.f;
    i__1 = nstop;
    for (j = nstart; j <= i__1; ++j) {
	s += f[j * f_dim1 + 2];
/* L151: */
    }
    a2 = (real) nunk;
    if (*nbdcnd == 0) {
	goto L152;
    }
    s -= (f[nstart * f_dim1 + 2] + f[nstop * f_dim1 + 2]) * .5f;
    a2 += -1.f;
L152:
    ypole = (dlrsq * .25f * f[f_dim1 + 1] - s / a2) / (w[id5 + 1] - 1.f + *
	    elmbda * dlrsq * .25f);
    i__1 = mstop;
    for (i__ = mstart; i__ <= i__1; ++i__) {
	k = l + i__;
	i__2 = nstop;
	for (j = nstart; j <= i__2; ++j) {
	    f[i__ + j * f_dim1] += ypole * w[k];
/* L153: */
	}
/* L154: */
    }
L155:
    i__1 = np1;
    for (j = 1; j <= i__1; ++j) {
	f[j * f_dim1 + 1] = ypole;
/* L156: */
    }
L157:
    if (*nbdcnd != 0) {
	goto L159;
    }
    i__1 = mstop;
    for (i__ = mstart; i__ <= i__1; ++i__) {
	f[i__ + np1 * f_dim1] = f[i__ + f_dim1];
/* L158: */
    }
L159:
    w[1] = (real) iwstor;
    return 0;
} /* hwsplr_ */

