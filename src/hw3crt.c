/* hw3crt.f -- translated by f2c (version 12.02.01).
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

/* DECK HW3CRT */
/* Subroutine */ int hw3crt_(real *xs, real *xf, integer *l, integer *lbdcnd, 
	real *bdxs, real *bdxf, real *ys, real *yf, integer *m, integer *
	mbdcnd, real *bdys, real *bdyf, real *zs, real *zf, integer *n, 
	integer *nbdcnd, real *bdzs, real *bdzf, real *elmbda, integer *ldimf,
	 integer *mdimf, real *f, real *pertrb, integer *ierror, real *w)
{
    /* System generated locals */
    integer bdxs_dim1, bdxs_offset, bdxf_dim1, bdxf_offset, bdys_dim1, 
	    bdys_offset, bdyf_dim1, bdyf_offset, bdzs_dim1, bdzs_offset, 
	    bdzf_dim1, bdzf_offset, f_dim1, f_dim2, f_offset, i__1, i__2, 
	    i__3;
    real r__1;

    /* Local variables */
    static integer i__, j, k;
    static real s, c1, c2, c3, s1, s2;
    static integer ir;
    static real dx, dy;
    static integer mp;
    static real dz;
    static integer np, lp, lp1, mp1, np1, iwb, iwc;
    static real xlp, ylp, zlp;
    static integer iww, lunk, munk, nunk, lstop, mstop, nstop;
    extern /* Subroutine */ int pois3d_(integer *, integer *, real *, integer 
	    *, integer *, real *, integer *, integer *, real *, real *, real *
	    , integer *, integer *, real *, integer *, real *);
    static integer lstpm1, mstpm1, nstpm1, nperod, lstart, mstart, nstart;
    static real twbydx, twbydy, twbydz;

/* ***BEGIN PROLOGUE  HW3CRT */
/* ***PURPOSE  Solve the standard seven-point finite difference */
/*            approximation to the Helmholtz equation in Cartesian */
/*            coordinates. */
/* ***LIBRARY   SLATEC (FISHPACK) */
/* ***CATEGORY  I2B1A1A */
/* ***TYPE      SINGLE PRECISION (HW3CRT-S) */
/* ***KEYWORDS  CARTESIAN, ELLIPTIC, FISHPACK, HELMHOLTZ, PDE */
/* ***AUTHOR  Adams, J., (NCAR) */
/*           Swarztrauber, P. N., (NCAR) */
/*           Sweet, R., (NCAR) */
/* ***DESCRIPTION */

/*     Subroutine HW3CRT solves the standard seven-point finite */
/*     difference approximation to the Helmholtz equation in Cartesian */
/*     coordinates: */

/*         (d/dX)(dU/dX) + (d/dY)(dU/dY) + (d/dZ)(dU/dZ) */

/*                    + LAMBDA*U = F(X,Y,Z) . */

/*    * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */


/*    * * * * * * * *    Parameter Description     * * * * * * * * * * */


/*            * * * * * *   On Input    * * * * * * */

/*     XS,XF */
/*        The range of X, i.e. XS .LE. X .LE. XF . */
/*        XS must be less than XF. */

/*     L */
/*        The number of panels into which the interval (XS,XF) is */
/*        subdivided.  Hence, there will be L+1 grid points in the */
/*        X-direction given by X(I) = XS+(I-1)DX for I=1,2,...,L+1, */
/*        where DX = (XF-XS)/L is the panel width.  L must be at */
/*        least 5 . */

/*     LBDCND */
/*        Indicates the type of boundary conditions at X = XS and X = XF. */

/*        = 0  If the solution is periodic in X, i.e. */
/*             U(L+I,J,K) = U(I,J,K). */
/*        = 1  If the solution is specified at X = XS and X = XF. */
/*        = 2  If the solution is specified at X = XS and the derivative */
/*             of the solution with respect to X is specified at X = XF. */
/*        = 3  If the derivative of the solution with respect to X is */
/*             specified at X = XS and X = XF. */
/*        = 4  If the derivative of the solution with respect to X is */
/*             specified at X = XS and the solution is specified at X=XF. */

/*     BDXS */
/*        A two-dimensional array that specifies the values of the */
/*        derivative of the solution with respect to X at X = XS. */
/*        when LBDCND = 3 or 4, */

/*             BDXS(J,K) = (d/dX)U(XS,Y(J),Z(K)), J=1,2,...,M+1, */
/*                                                K=1,2,...,N+1. */

/*        When LBDCND has any other value, BDXS is a dummy variable. */
/*        BDXS must be dimensioned at least (M+1)*(N+1). */

/*     BDXF */
/*        A two-dimensional array that specifies the values of the */
/*        derivative of the solution with respect to X at X = XF. */
/*        When LBDCND = 2 or 3, */

/*             BDXF(J,K) = (d/dX)U(XF,Y(J),Z(K)), J=1,2,...,M+1, */
/*                                                K=1,2,...,N+1. */

/*        When LBDCND has any other value, BDXF is a dummy variable. */
/*        BDXF must be dimensioned at least (M+1)*(N+1). */

/*     YS,YF */
/*        The range of Y, i.e. YS .LE. Y .LE. YF. */
/*        YS must be less than YF. */

/*     M */
/*        The number of panels into which the interval (YS,YF) is */
/*        subdivided.  Hence, there will be M+1 grid points in the */
/*        Y-direction given by Y(J) = YS+(J-1)DY for J=1,2,...,M+1, */
/*        where DY = (YF-YS)/M is the panel width.  M must be at */
/*        least 5 . */

/*     MBDCND */
/*        Indicates the type of boundary conditions at Y = YS and Y = YF. */

/*        = 0  If the solution is periodic in Y, i.e. */
/*             U(I,M+J,K) = U(I,J,K). */
/*        = 1  If the solution is specified at Y = YS and Y = YF. */
/*        = 2  If the solution is specified at Y = YS and the derivative */
/*             of the solution with respect to Y is specified at Y = YF. */
/*        = 3  If the derivative of the solution with respect to Y is */
/*             specified at Y = YS and Y = YF. */
/*        = 4  If the derivative of the solution with respect to Y is */
/*             specified at Y = YS and the solution is specified at Y=YF. */

/*     BDYS */
/*        A two-dimensional array that specifies the values of the */
/*        derivative of the solution with respect to Y at Y = YS. */
/*        When MBDCND = 3 or 4, */

/*             BDYS(I,K) = (d/dY)U(X(I),YS,Z(K)), I=1,2,...,L+1, */
/*                                                K=1,2,...,N+1. */

/*        When MBDCND has any other value, BDYS is a dummy variable. */
/*        BDYS must be dimensioned at least (L+1)*(N+1). */

/*     BDYF */
/*        A two-dimensional array that specifies the values of the */
/*        derivative of the solution with respect to Y at Y = YF. */
/*        When MBDCND = 2 or 3, */

/*             BDYF(I,K) = (d/dY)U(X(I),YF,Z(K)), I=1,2,...,L+1, */
/*                                                K=1,2,...,N+1. */

/*        When MBDCND has any other value, BDYF is a dummy variable. */
/*        BDYF must be dimensioned at least (L+1)*(N+1). */

/*     ZS,ZF */
/*        The range of Z, i.e. ZS .LE. Z .LE. ZF. */
/*        ZS must be less than ZF. */

/*     N */
/*        The number of panels into which the interval (ZS,ZF) is */
/*        subdivided.  Hence, there will be N+1 grid points in the */
/*        Z-direction given by Z(K) = ZS+(K-1)DZ for K=1,2,...,N+1, */
/*        where DZ = (ZF-ZS)/N is the panel width.  N must be at least 5. */

/*     NBDCND */
/*        Indicates the type of boundary conditions at Z = ZS and Z = ZF. */

/*        = 0  If the solution is periodic in Z, i.e. */
/*             U(I,J,N+K) = U(I,J,K). */
/*        = 1  If the solution is specified at Z = ZS and Z = ZF. */
/*        = 2  If the solution is specified at Z = ZS and the derivative */
/*             of the solution with respect to Z is specified at Z = ZF. */
/*        = 3  If the derivative of the solution with respect to Z is */
/*             specified at Z = ZS and Z = ZF. */
/*        = 4  If the derivative of the solution with respect to Z is */
/*             specified at Z = ZS and the solution is specified at Z=ZF. */

/*     BDZS */
/*        A two-dimensional array that specifies the values of the */
/*        derivative of the solution with respect to Z at Z = ZS. */
/*        When NBDCND = 3 or 4, */

/*             BDZS(I,J) = (d/dZ)U(X(I),Y(J),ZS), I=1,2,...,L+1, */
/*                                                J=1,2,...,M+1. */

/*        When NBDCND has any other value, BDZS is a dummy variable. */
/*        BDZS must be dimensioned at least (L+1)*(M+1). */

/*     BDZF */
/*        A two-dimensional array that specifies the values of the */
/*        derivative of the solution with respect to Z at Z = ZF. */
/*        When NBDCND = 2 or 3, */

/*             BDZF(I,J) = (d/dZ)U(X(I),Y(J),ZF), I=1,2,...,L+1, */
/*                                                J=1,2,...,M+1. */

/*        When NBDCND has any other value, BDZF is a dummy variable. */
/*        BDZF must be dimensioned at least (L+1)*(M+1). */

/*     ELMBDA */
/*        The constant LAMBDA in the Helmholtz equation. If */
/*        LAMBDA .GT. 0, a solution may not exist.  However, HW3CRT will */
/*        attempt to find a solution. */

/*     F */
/*        A three-dimensional array that specifies the values of the */
/*        right side of the Helmholtz equation and boundary values (if */
/*        any).  For I=2,3,...,L, J=2,3,...,M, and K=2,3,...,N */

/*                   F(I,J,K) = F(X(I),Y(J),Z(K)). */

/*        On the boundaries F is defined by */

/*        LBDCND      F(1,J,K)         F(L+1,J,K) */
/*        ------   ---------------   --------------- */

/*          0      F(XS,Y(J),Z(K))   F(XS,Y(J),Z(K)) */
/*          1      U(XS,Y(J),Z(K))   U(XF,Y(J),Z(K)) */
/*          2      U(XS,Y(J),Z(K))   F(XF,Y(J),Z(K))   J=1,2,...,M+1 */
/*          3      F(XS,Y(J),Z(K))   F(XF,Y(J),Z(K))   K=1,2,...,N+1 */
/*          4      F(XS,Y(J),Z(K))   U(XF,Y(J),Z(K)) */

/*        MBDCND      F(I,1,K)         F(I,M+1,K) */
/*        ------   ---------------   --------------- */

/*          0      F(X(I),YS,Z(K))   F(X(I),YS,Z(K)) */
/*          1      U(X(I),YS,Z(K))   U(X(I),YF,Z(K)) */
/*          2      U(X(I),YS,Z(K))   F(X(I),YF,Z(K))   I=1,2,...,L+1 */
/*          3      F(X(I),YS,Z(K))   F(X(I),YF,Z(K))   K=1,2,...,N+1 */
/*          4      F(X(I),YS,Z(K))   U(X(I),YF,Z(K)) */

/*        NBDCND      F(I,J,1)         F(I,J,N+1) */
/*        ------   ---------------   --------------- */

/*          0      F(X(I),Y(J),ZS)   F(X(I),Y(J),ZS) */
/*          1      U(X(I),Y(J),ZS)   U(X(I),Y(J),ZF) */
/*          2      U(X(I),Y(J),ZS)   F(X(I),Y(J),ZF)   I=1,2,...,L+1 */
/*          3      F(X(I),Y(J),ZS)   F(X(I),Y(J),ZF)   J=1,2,...,M+1 */
/*          4      F(X(I),Y(J),ZS)   U(X(I),Y(J),ZF) */

/*        F must be dimensioned at least (L+1)*(M+1)*(N+1). */

/*        NOTE: */

/*        If the table calls for both the solution U and the right side F */
/*        on a boundary, then the solution must be specified. */

/*     LDIMF */
/*        The row (or first) dimension of the arrays F,BDYS,BDYF,BDZS, */
/*        and BDZF as it appears in the program calling HW3CRT. this */
/*        parameter is used to specify the variable dimension of these */
/*        arrays.  LDIMF must be at least L+1. */

/*     MDIMF */
/*        The column (or second) dimension of the array F and the row (or */
/*        first) dimension of the arrays BDXS and BDXF as it appears in */
/*        the program calling HW3CRT.  This parameter is used to specify */
/*        the variable dimension of these arrays. */
/*        MDIMF must be at least M+1. */

/*     W */
/*        A one-dimensional array that must be provided by the user for */
/*        work space.  The length of W must be at least 30 + L + M + 5*N */
/*        + MAX(L,M,N) + 7*(INT((L+1)/2) + INT((M+1)/2)) */


/*            * * * * * *   On Output   * * * * * * */

/*     F */
/*        Contains the solution U(I,J,K) of the finite difference */
/*        approximation for the grid point (X(I),Y(J),Z(K)) for */
/*        I=1,2,...,L+1, J=1,2,...,M+1, and K=1,2,...,N+1. */

/*     PERTRB */
/*        If a combination of periodic or derivative boundary conditions */
/*        is specified for a Poisson equation (LAMBDA = 0), a solution */
/*        may not exist.  PERTRB is a constant, calculated and subtracted */
/*        from F, which ensures that a solution exists.  PWSCRT then */
/*        computes this solution, which is a least squares solution to */
/*        the original approximation.  This solution is not unique and is */
/*        unnormalized.  The value of PERTRB should be small compared to */
/*        the right side F.  Otherwise, a solution is obtained to an */
/*        essentially different problem.  This comparison should always */
/*        be made to insure that a meaningful solution has been obtained. */

/*     IERROR */
/*        An error flag that indicates invalid input parameters.  Except */
/*        for numbers 0 and 12, a solution is not attempted. */

/*        =  0  No error */
/*        =  1  XS .GE. XF */
/*        =  2  L .LT. 5 */
/*        =  3  LBDCND .LT. 0 .OR. LBDCND .GT. 4 */
/*        =  4  YS .GE. YF */
/*        =  5  M .LT. 5 */
/*        =  6  MBDCND .LT. 0 .OR. MBDCND .GT. 4 */
/*        =  7  ZS .GE. ZF */
/*        =  8  N .LT. 5 */
/*        =  9  NBDCND .LT. 0 .OR. NBDCND .GT. 4 */
/*        = 10  LDIMF .LT. L+1 */
/*        = 11  MDIMF .LT. M+1 */
/*        = 12  LAMBDA .GT. 0 */

/*        Since this is the only means of indicating a possibly incorrect */
/*        call to HW3CRT, the user should test IERROR after the call. */

/* *Long Description: */

/*    * * * * * * *   Program Specifications    * * * * * * * * * * * * */

/*     Dimension of   BDXS(MDIMF,N+1),BDXF(MDIMF,N+1),BDYS(LDIMF,N+1), */
/*     Arguments      BDYF(LDIMF,N+1),BDZS(LDIMF,M+1),BDZF(LDIMF,M+1), */
/*                    F(LDIMF,MDIMF,N+1),W(see argument list) */

/*     Latest         December 1, 1978 */
/*     Revision */

/*     Subprograms    HW3CRT,POIS3D,POS3D1,TRIDQ,RFFTI,RFFTF,RFFTF1, */
/*     Required       RFFTB,RFFTB1,COSTI,COST,SINTI,SINT,COSQI,COSQF, */
/*                    COSQF1,COSQB,COSQB1,SINQI,SINQF,SINQB,CFFTI, */
/*                    CFFTI1,CFFTB,CFFTB1,PASSB2,PASSB3,PASSB4,PASSB, */
/*                    CFFTF,CFFTF1,PASSF1,PASSF2,PASSF3,PASSF4,PASSF, */
/*                    PIMACH */

/*     Special        NONE */
/*     Conditions */

/*     Common         NONE */
/*     Blocks */

/*     I/O            NONE */

/*     Precision      Single */

/*     Specialist     Roland Sweet */

/*     Language       FORTRAN */

/*     History        Written by Roland Sweet at NCAR in July 1977 */

/*     Algorithm      This subroutine defines the finite difference */
/*                    equations, incorporates boundary data, and */
/*                    adjusts the right side of singular systems and */
/*                    then calls POIS3D to solve the system. */

/*     Space          7862(decimal) = 17300(octal) locations on the */
/*     Required       NCAR Control Data 7600 */

/*     Timing and        The execution time T on the NCAR Control Data */
/*     Accuracy       7600 for subroutine HW3CRT is roughly proportional */
/*                    to L*M*N*(log2(L)+log2(M)+5), but also depends on */
/*                    input parameters LBDCND and MBDCND.  Some typical */
/*                    values are listed in the table below. */
/*                       The solution process employed results in a loss */
/*                    of no more than three significant digits for L,M */
/*                    and N as large as 32.  More detailed information */
/*                    about accuracy can be found in the documentation */
/*                    for subroutine POIS3D which is the routine that */
/*                    actually solves the finite difference equations. */


/*                       L(=M=N)     LBDCND(=MBDCND=NBDCND)      T(MSECS) */
/*                       -------     ----------------------      -------- */

/*                         16                  0                    300 */
/*                         16                  1                    302 */
/*                         16                  3                    348 */
/*                         32                  0                   1925 */
/*                         32                  1                   1929 */
/*                         32                  3                   2109 */

/*     Portability    American National Standards Institute FORTRAN. */
/*                    The machine dependent constant PI is defined in */
/*                    function PIMACH. */

/*     Required       COS,SIN,ATAN */
/*     Resident */
/*     Routines */

/*     Reference      NONE */

/*    * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/* ***REFERENCES  (NONE) */
/* ***ROUTINES CALLED  POIS3D */
/* ***REVISION HISTORY  (YYMMDD) */
/*   801001  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890531  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/* ***END PROLOGUE  HW3CRT */


/* ***FIRST EXECUTABLE STATEMENT  HW3CRT */
    /* Parameter adjustments */
    bdzf_dim1 = *ldimf;
    bdzf_offset = 1 + bdzf_dim1;
    bdzf -= bdzf_offset;
    bdzs_dim1 = *ldimf;
    bdzs_offset = 1 + bdzs_dim1;
    bdzs -= bdzs_offset;
    bdyf_dim1 = *ldimf;
    bdyf_offset = 1 + bdyf_dim1;
    bdyf -= bdyf_offset;
    bdys_dim1 = *ldimf;
    bdys_offset = 1 + bdys_dim1;
    bdys -= bdys_offset;
    f_dim1 = *ldimf;
    f_dim2 = *mdimf;
    f_offset = 1 + f_dim1 * (1 + f_dim2);
    f -= f_offset;
    bdxf_dim1 = *mdimf;
    bdxf_offset = 1 + bdxf_dim1;
    bdxf -= bdxf_offset;
    bdxs_dim1 = *mdimf;
    bdxs_offset = 1 + bdxs_dim1;
    bdxs -= bdxs_offset;
    --w;

    /* Function Body */
    *ierror = 0;
    if (*xf <= *xs) {
	*ierror = 1;
    }
    if (*l < 5) {
	*ierror = 2;
    }
    if (*lbdcnd < 0 || *lbdcnd > 4) {
	*ierror = 3;
    }
    if (*yf <= *ys) {
	*ierror = 4;
    }
    if (*m < 5) {
	*ierror = 5;
    }
    if (*mbdcnd < 0 || *mbdcnd > 4) {
	*ierror = 6;
    }
    if (*zf <= *zs) {
	*ierror = 7;
    }
    if (*n < 5) {
	*ierror = 8;
    }
    if (*nbdcnd < 0 || *nbdcnd > 4) {
	*ierror = 9;
    }
    if (*ldimf < *l + 1) {
	*ierror = 10;
    }
    if (*mdimf < *m + 1) {
	*ierror = 11;
    }
    if (*ierror != 0) {
	goto L188;
    }
    dy = (*yf - *ys) / *m;
    twbydy = 2.f / dy;
/* Computing 2nd power */
    r__1 = dy;
    c2 = 1.f / (r__1 * r__1);
    mstart = 1;
    mstop = *m;
    mp1 = *m + 1;
    mp = *mbdcnd + 1;
    switch (mp) {
	case 1:  goto L104;
	case 2:  goto L101;
	case 3:  goto L101;
	case 4:  goto L102;
	case 5:  goto L102;
    }
L101:
    mstart = 2;
L102:
    switch (mp) {
	case 1:  goto L104;
	case 2:  goto L104;
	case 3:  goto L103;
	case 4:  goto L103;
	case 5:  goto L104;
    }
L103:
    mstop = mp1;
L104:
    munk = mstop - mstart + 1;
    dz = (*zf - *zs) / *n;
    twbydz = 2.f / dz;
    np = *nbdcnd + 1;
/* Computing 2nd power */
    r__1 = dz;
    c3 = 1.f / (r__1 * r__1);
    np1 = *n + 1;
    nstart = 1;
    nstop = *n;
    switch (np) {
	case 1:  goto L108;
	case 2:  goto L105;
	case 3:  goto L105;
	case 4:  goto L106;
	case 5:  goto L106;
    }
L105:
    nstart = 2;
L106:
    switch (np) {
	case 1:  goto L108;
	case 2:  goto L108;
	case 3:  goto L107;
	case 4:  goto L107;
	case 5:  goto L108;
    }
L107:
    nstop = np1;
L108:
    nunk = nstop - nstart + 1;
    lp1 = *l + 1;
    dx = (*xf - *xs) / *l;
/* Computing 2nd power */
    r__1 = dx;
    c1 = 1.f / (r__1 * r__1);
    twbydx = 2.f / dx;
    lp = *lbdcnd + 1;
    lstart = 1;
    lstop = *l;

/*     ENTER BOUNDARY DATA FOR X-BOUNDARIES. */

    switch (lp) {
	case 1:  goto L122;
	case 2:  goto L109;
	case 3:  goto L109;
	case 4:  goto L112;
	case 5:  goto L112;
    }
L109:
    lstart = 2;
    i__1 = mstop;
    for (j = mstart; j <= i__1; ++j) {
	i__2 = nstop;
	for (k = nstart; k <= i__2; ++k) {
	    f[(j + k * f_dim2) * f_dim1 + 2] -= c1 * f[(j + k * f_dim2) * 
		    f_dim1 + 1];
/* L110: */
	}
/* L111: */
    }
    goto L115;
L112:
    i__1 = mstop;
    for (j = mstart; j <= i__1; ++j) {
	i__2 = nstop;
	for (k = nstart; k <= i__2; ++k) {
	    f[(j + k * f_dim2) * f_dim1 + 1] += twbydx * bdxs[j + k * 
		    bdxs_dim1];
/* L113: */
	}
/* L114: */
    }
L115:
    switch (lp) {
	case 1:  goto L122;
	case 2:  goto L116;
	case 3:  goto L119;
	case 4:  goto L119;
	case 5:  goto L116;
    }
L116:
    i__1 = mstop;
    for (j = mstart; j <= i__1; ++j) {
	i__2 = nstop;
	for (k = nstart; k <= i__2; ++k) {
	    f[*l + (j + k * f_dim2) * f_dim1] -= c1 * f[lp1 + (j + k * f_dim2)
		     * f_dim1];
/* L117: */
	}
/* L118: */
    }
    goto L122;
L119:
    lstop = lp1;
    i__1 = mstop;
    for (j = mstart; j <= i__1; ++j) {
	i__2 = nstop;
	for (k = nstart; k <= i__2; ++k) {
	    f[lp1 + (j + k * f_dim2) * f_dim1] -= twbydx * bdxf[j + k * 
		    bdxf_dim1];
/* L120: */
	}
/* L121: */
    }
L122:
    lunk = lstop - lstart + 1;

/*     ENTER BOUNDARY DATA FOR Y-BOUNDARIES. */

    switch (mp) {
	case 1:  goto L136;
	case 2:  goto L123;
	case 3:  goto L123;
	case 4:  goto L126;
	case 5:  goto L126;
    }
L123:
    i__1 = lstop;
    for (i__ = lstart; i__ <= i__1; ++i__) {
	i__2 = nstop;
	for (k = nstart; k <= i__2; ++k) {
	    f[i__ + (k * f_dim2 + 2) * f_dim1] -= c2 * f[i__ + (k * f_dim2 + 
		    1) * f_dim1];
/* L124: */
	}
/* L125: */
    }
    goto L129;
L126:
    i__1 = lstop;
    for (i__ = lstart; i__ <= i__1; ++i__) {
	i__2 = nstop;
	for (k = nstart; k <= i__2; ++k) {
	    f[i__ + (k * f_dim2 + 1) * f_dim1] += twbydy * bdys[i__ + k * 
		    bdys_dim1];
/* L127: */
	}
/* L128: */
    }
L129:
    switch (mp) {
	case 1:  goto L136;
	case 2:  goto L130;
	case 3:  goto L133;
	case 4:  goto L133;
	case 5:  goto L130;
    }
L130:
    i__1 = lstop;
    for (i__ = lstart; i__ <= i__1; ++i__) {
	i__2 = nstop;
	for (k = nstart; k <= i__2; ++k) {
	    f[i__ + (*m + k * f_dim2) * f_dim1] -= c2 * f[i__ + (mp1 + k * 
		    f_dim2) * f_dim1];
/* L131: */
	}
/* L132: */
    }
    goto L136;
L133:
    i__1 = lstop;
    for (i__ = lstart; i__ <= i__1; ++i__) {
	i__2 = nstop;
	for (k = nstart; k <= i__2; ++k) {
	    f[i__ + (mp1 + k * f_dim2) * f_dim1] -= twbydy * bdyf[i__ + k * 
		    bdyf_dim1];
/* L134: */
	}
/* L135: */
    }
L136:

/*     ENTER BOUNDARY DATA FOR Z-BOUNDARIES. */

    switch (np) {
	case 1:  goto L150;
	case 2:  goto L137;
	case 3:  goto L137;
	case 4:  goto L140;
	case 5:  goto L140;
    }
L137:
    i__1 = lstop;
    for (i__ = lstart; i__ <= i__1; ++i__) {
	i__2 = mstop;
	for (j = mstart; j <= i__2; ++j) {
	    f[i__ + (j + (f_dim2 << 1)) * f_dim1] -= c3 * f[i__ + (j + f_dim2)
		     * f_dim1];
/* L138: */
	}
/* L139: */
    }
    goto L143;
L140:
    i__1 = lstop;
    for (i__ = lstart; i__ <= i__1; ++i__) {
	i__2 = mstop;
	for (j = mstart; j <= i__2; ++j) {
	    f[i__ + (j + f_dim2) * f_dim1] += twbydz * bdzs[i__ + j * 
		    bdzs_dim1];
/* L141: */
	}
/* L142: */
    }
L143:
    switch (np) {
	case 1:  goto L150;
	case 2:  goto L144;
	case 3:  goto L147;
	case 4:  goto L147;
	case 5:  goto L144;
    }
L144:
    i__1 = lstop;
    for (i__ = lstart; i__ <= i__1; ++i__) {
	i__2 = mstop;
	for (j = mstart; j <= i__2; ++j) {
	    f[i__ + (j + *n * f_dim2) * f_dim1] -= c3 * f[i__ + (j + np1 * 
		    f_dim2) * f_dim1];
/* L145: */
	}
/* L146: */
    }
    goto L150;
L147:
    i__1 = lstop;
    for (i__ = lstart; i__ <= i__1; ++i__) {
	i__2 = mstop;
	for (j = mstart; j <= i__2; ++j) {
	    f[i__ + (j + np1 * f_dim2) * f_dim1] -= twbydz * bdzf[i__ + j * 
		    bdzf_dim1];
/* L148: */
	}
/* L149: */
    }

/*     DEFINE A,B,C COEFFICIENTS IN W-ARRAY. */

L150:
    iwb = nunk + 1;
    iwc = iwb + nunk;
    iww = iwc + nunk;
    i__1 = nunk;
    for (k = 1; k <= i__1; ++k) {
	i__ = iwc + k - 1;
	w[k] = c3;
	w[i__] = c3;
	i__ = iwb + k - 1;
	w[i__] = c3 * -2.f + *elmbda;
/* L151: */
    }
    switch (np) {
	case 1:  goto L155;
	case 2:  goto L155;
	case 3:  goto L153;
	case 4:  goto L152;
	case 5:  goto L152;
    }
L152:
    w[iwc] = c3 * 2.f;
L153:
    switch (np) {
	case 1:  goto L155;
	case 2:  goto L155;
	case 3:  goto L154;
	case 4:  goto L154;
	case 5:  goto L155;
    }
L154:
    w[iwb - 1] = c3 * 2.f;
L155:
    *pertrb = 0.f;

/*     FOR SINGULAR PROBLEMS ADJUST DATA TO INSURE A SOLUTION WILL EXIST. */

    switch (lp) {
	case 1:  goto L156;
	case 2:  goto L172;
	case 3:  goto L172;
	case 4:  goto L156;
	case 5:  goto L172;
    }
L156:
    switch (mp) {
	case 1:  goto L157;
	case 2:  goto L172;
	case 3:  goto L172;
	case 4:  goto L157;
	case 5:  goto L172;
    }
L157:
    switch (np) {
	case 1:  goto L158;
	case 2:  goto L172;
	case 3:  goto L172;
	case 4:  goto L158;
	case 5:  goto L172;
    }
L158:
    if (*elmbda < 0.f) {
	goto L172;
    } else if (*elmbda == 0) {
	goto L160;
    } else {
	goto L159;
    }
L159:
    *ierror = 12;
    goto L172;
L160:
    mstpm1 = mstop - 1;
    lstpm1 = lstop - 1;
    nstpm1 = nstop - 1;
    xlp = (real) ((lp + 2) / 3);
    ylp = (real) ((mp + 2) / 3);
    zlp = (real) ((np + 2) / 3);
    s1 = 0.f;
    i__1 = nstpm1;
    for (k = 2; k <= i__1; ++k) {
	i__2 = mstpm1;
	for (j = 2; j <= i__2; ++j) {
	    i__3 = lstpm1;
	    for (i__ = 2; i__ <= i__3; ++i__) {
		s1 += f[i__ + (j + k * f_dim2) * f_dim1];
/* L161: */
	    }
	    s1 += (f[(j + k * f_dim2) * f_dim1 + 1] + f[lstop + (j + k * 
		    f_dim2) * f_dim1]) / xlp;
/* L162: */
	}
	s2 = 0.f;
	i__2 = lstpm1;
	for (i__ = 2; i__ <= i__2; ++i__) {
	    s2 = s2 + f[i__ + (k * f_dim2 + 1) * f_dim1] + f[i__ + (mstop + k 
		    * f_dim2) * f_dim1];
/* L163: */
	}
	s2 = (s2 + (f[(k * f_dim2 + 1) * f_dim1 + 1] + f[(mstop + k * f_dim2) 
		* f_dim1 + 1] + f[lstop + (k * f_dim2 + 1) * f_dim1] + f[
		lstop + (mstop + k * f_dim2) * f_dim1]) / xlp) / ylp;
	s1 += s2;
/* L164: */
    }
    s = (f[(f_dim2 + 1) * f_dim1 + 1] + f[lstop + (f_dim2 + 1) * f_dim1] + f[(
	    nstop * f_dim2 + 1) * f_dim1 + 1] + f[lstop + (nstop * f_dim2 + 1)
	     * f_dim1] + f[(mstop + f_dim2) * f_dim1 + 1] + f[lstop + (mstop 
	    + f_dim2) * f_dim1] + f[(mstop + nstop * f_dim2) * f_dim1 + 1] + 
	    f[lstop + (mstop + nstop * f_dim2) * f_dim1]) / (xlp * ylp);
    i__1 = mstpm1;
    for (j = 2; j <= i__1; ++j) {
	i__2 = lstpm1;
	for (i__ = 2; i__ <= i__2; ++i__) {
	    s = s + f[i__ + (j + f_dim2) * f_dim1] + f[i__ + (j + nstop * 
		    f_dim2) * f_dim1];
/* L165: */
	}
/* L166: */
    }
    s2 = 0.f;
    i__1 = lstpm1;
    for (i__ = 2; i__ <= i__1; ++i__) {
	s2 = s2 + f[i__ + (f_dim2 + 1) * f_dim1] + f[i__ + (nstop * f_dim2 + 
		1) * f_dim1] + f[i__ + (mstop + f_dim2) * f_dim1] + f[i__ + (
		mstop + nstop * f_dim2) * f_dim1];
/* L167: */
    }
    s = s2 / ylp + s;
    s2 = 0.f;
    i__1 = mstpm1;
    for (j = 2; j <= i__1; ++j) {
	s2 = s2 + f[(j + f_dim2) * f_dim1 + 1] + f[(j + nstop * f_dim2) * 
		f_dim1 + 1] + f[lstop + (j + f_dim2) * f_dim1] + f[lstop + (j 
		+ nstop * f_dim2) * f_dim1];
/* L168: */
    }
    s = s2 / xlp + s;
    *pertrb = (s / zlp + s1) / ((lunk + 1.f - xlp) * (munk + 1.f - ylp) * (
	    nunk + 1.f - zlp));
    i__1 = lunk;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = munk;
	for (j = 1; j <= i__2; ++j) {
	    i__3 = nunk;
	    for (k = 1; k <= i__3; ++k) {
		f[i__ + (j + k * f_dim2) * f_dim1] -= *pertrb;
/* L169: */
	    }
/* L170: */
	}
/* L171: */
    }
L172:
    nperod = 0;
    if (*nbdcnd == 0) {
	goto L173;
    }
    nperod = 1;
    w[1] = 0.f;
    w[iww - 1] = 0.f;
L173:
    pois3d_(lbdcnd, &lunk, &c1, mbdcnd, &munk, &c2, &nperod, &nunk, &w[1], &w[
	    iwb], &w[iwc], ldimf, mdimf, &f[lstart + (mstart + nstart * 
	    f_dim2) * f_dim1], &ir, &w[iww]);

/*     FILL IN SIDES FOR PERIODIC BOUNDARY CONDITIONS. */

    if (lp != 1) {
	goto L180;
    }
    if (mp != 1) {
	goto L175;
    }
    i__1 = nstop;
    for (k = nstart; k <= i__1; ++k) {
	f[(mp1 + k * f_dim2) * f_dim1 + 1] = f[(k * f_dim2 + 1) * f_dim1 + 1];
/* L174: */
    }
    mstop = mp1;
L175:
    if (np != 1) {
	goto L177;
    }
    i__1 = mstop;
    for (j = mstart; j <= i__1; ++j) {
	f[(j + np1 * f_dim2) * f_dim1 + 1] = f[(j + f_dim2) * f_dim1 + 1];
/* L176: */
    }
    nstop = np1;
L177:
    i__1 = mstop;
    for (j = mstart; j <= i__1; ++j) {
	i__2 = nstop;
	for (k = nstart; k <= i__2; ++k) {
	    f[lp1 + (j + k * f_dim2) * f_dim1] = f[(j + k * f_dim2) * f_dim1 
		    + 1];
/* L178: */
	}
/* L179: */
    }
L180:
    if (mp != 1) {
	goto L185;
    }
    if (np != 1) {
	goto L182;
    }
    i__1 = lstop;
    for (i__ = lstart; i__ <= i__1; ++i__) {
	f[i__ + (np1 * f_dim2 + 1) * f_dim1] = f[i__ + (f_dim2 + 1) * f_dim1];
/* L181: */
    }
    nstop = np1;
L182:
    i__1 = lstop;
    for (i__ = lstart; i__ <= i__1; ++i__) {
	i__2 = nstop;
	for (k = nstart; k <= i__2; ++k) {
	    f[i__ + (mp1 + k * f_dim2) * f_dim1] = f[i__ + (k * f_dim2 + 1) * 
		    f_dim1];
/* L183: */
	}
/* L184: */
    }
L185:
    if (np != 1) {
	goto L188;
    }
    i__1 = lstop;
    for (i__ = lstart; i__ <= i__1; ++i__) {
	i__2 = mstop;
	for (j = mstart; j <= i__2; ++j) {
	    f[i__ + (j + np1 * f_dim2) * f_dim1] = f[i__ + (j + f_dim2) * 
		    f_dim1];
/* L186: */
	}
/* L187: */
    }
L188:
    return 0;
} /* hw3crt_ */

