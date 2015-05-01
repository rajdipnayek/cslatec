/* hwscsp.f -- translated by f2c (version 12.02.01).
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

/* DECK HWSCSP */
/* Subroutine */ int hwscsp_(integer *intl, real *ts, real *tf, integer *m, 
	integer *mbdcnd, real *bdts, real *bdtf, real *rs, real *rf, integer *
	n, integer *nbdcnd, real *bdrs, real *bdrf, real *elmbda, real *f, 
	integer *idimf, real *pertrb, integer *ierror, real *w)
{
    /* System generated locals */
    integer f_dim1, f_offset, i__1, i__2;

    /* Local variables */
    static integer k, l, i1, i2, i3, i4, i5, i6, i7, i8, i9, i10;
    static real pi;
    static integer mp1, np1, nck;
    static real dum;
    extern /* Subroutine */ int hwscs1_(integer *, real *, real *, integer *, 
	    integer *, real *, real *, real *, real *, integer *, integer *, 
	    real *, real *, real *, real *, integer *, real *, real *, real *,
	     real *, real *, real *, real *, real *, real *, real *, real *, 
	    real *);
    extern doublereal pimach_(real *);

/* ***BEGIN PROLOGUE  HWSCSP */
/* ***PURPOSE  Solve a finite difference approximation to the modified */
/*            Helmholtz equation in spherical coordinates assuming */
/*            axisymmetry  (no dependence on longitude). */
/* ***LIBRARY   SLATEC (FISHPACK) */
/* ***CATEGORY  I2B1A1A */
/* ***TYPE      SINGLE PRECISION (HWSCSP-S) */
/* ***KEYWORDS  ELLIPTIC, FISHPACK, HELMHOLTZ, PDE, SPHERICAL */
/* ***AUTHOR  Adams, J., (NCAR) */
/*           Swarztrauber, P. N., (NCAR) */
/*           Sweet, R., (NCAR) */
/* ***DESCRIPTION */

/*     Subroutine HWSCSP solves a finite difference approximation to the */
/*       modified Helmholtz equation in spherical coordinates assuming */
/*       axisymmetry  (no dependence on longitude) */

/*          (1/R**2)(d/dR)((R**2)(d/dR)U) */

/*             + (1/(R**2)SIN(THETA))(d/dTHETA)(SIN(THETA)(d/dTHETA)U) */

/*             + (LAMBDA/(RSIN(THETA))**2)U = F(THETA,R). */

/*     This two dimensional modified Helmholtz equation results from */
/*     the Fourier transform of the three dimensional Poisson equation */

/*     * * * * * * * * * *     On Input     * * * * * * * * * * */

/*     INTL */
/*       = 0  On initial entry to HWSCSP or if any of the arguments */
/*            RS, RF, N, NBDCND are changed from a previous call. */
/*       = 1  If RS, RF, N, NBDCND are all unchanged from previous call */
/*            to HWSCSP. */

/*       NOTE   A call with INTL=0 takes approximately 1.5 times as */
/*              much time as a call with INTL = 1.  Once a call with */
/*              INTL = 0 has been made then subsequent solutions */
/*              corresponding to different F, BDTS, BDTF, BDRS, BDRF can */
/*              be obtained faster with INTL = 1 since initialization is */
/*              not repeated. */

/*     TS,TF */
/*       The range of THETA (colatitude), i.e., TS .LE. THETA .LE. TF. */
/*       TS must be less than TF.  TS and TF are in radians.  A TS of */
/*       zero corresponds to the north pole and a TF of PI corresponds */
/*       to the south pole. */

/*     * * * * * * * * * * * * * * IMPORTANT * * * * * * * * * * * * * * */

/*     If TF is equal to PI then it must be computed using the statement */
/*     TF = PIMACH(DUM). This insures that TF in the users program is */
/*     equal to PI in this program which permits several tests of the */
/*     input parameters that otherwise would not be possible. */

/*     M */
/*       The number of panels into which the interval (TS,TF) is */
/*       subdivided.  Hence, there will be M+1 grid points in the */
/*       THETA-direction given by THETA(K) = (I-1)DTHETA+TS for */
/*       I = 1,2,...,M+1, where DTHETA = (TF-TS)/M is the panel width. */

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
/*            specified at THETA = TS (see note 1 below) and the solution */
/*            is unspecified at THETA = TF = PI. */
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

/*            BDTS(J) = (d/dTHETA)U(TS,R(J)), J = 1,2,...,N+1  . */

/*       When MBDCND has any other value, BDTS is a dummy variable. */

/*     BDTF */
/*       A one-dimensional array of length N+1 that specifies the values */
/*       of the derivative of the solution with respect to THETA at */
/*       THETA = TF.  When MBDCND = 2,3, or 6, */

/*            BDTF(J) = (d/dTHETA)U(TF,R(J)), J = 1,2,...,N+1  . */

/*       When MBDCND has any other value, BDTF is a dummy variable. */

/*     RS,RF */
/*       The range of R, i.e., RS .LE. R .LT. RF.  RS must be less than */
/*       RF.  RS must be non-negative. */

/*       N */
/*       The number of panels into which the interval (RS,RF) is */
/*       subdivided.  Hence, there will be N+1 grid points in the */
/*       R-direction given by R(J) = (J-1)DR+RS for J = 1,2,...,N+1, */
/*       where DR = (RF-RS)/N is the panel width. */
/*       N must be greater than 2 */

/*     NBDCND */
/*       Indicates the type of boundary condition at R = RS and R = RF. */

/*       = 1  If the solution is specified at R = RS and R = RF. */
/*       = 2  If the solution is specified at R = RS and the derivative */
/*            of the solution with respect to R is specified at R = RF. */
/*       = 3  If the derivative of the solution with respect to R is */
/*            specified at R = RS and R = RF. */
/*       = 4  If the derivative of the solution with respect to R is */
/*            specified at RS and the solution is specified at R = RF. */
/*       = 5  If the solution is unspecified at R = RS = 0 (see note */
/*            below) and the solution is specified at R = RF. */
/*       = 6  If the solution is unspecified at R = RS = 0 (see note */
/*            below) and the derivative of the solution with respect to */
/*            R is specified at R = RF. */

/*       NOTE:  NBDCND = 5 or 6 cannot be used with */
/*              MBDCND = 1,2,4,5, or 7 (the former indicates that the */
/*                       solution is unspecified at R = 0, the latter */
/*                       indicates that the solution is specified). */
/*                       Use instead */
/*              NBDCND = 1 or 2  . */

/*     BDRS */
/*       A one-dimensional array of length M+1 that specifies the values */
/*       of the derivative of the solution with respect to R at R = RS. */
/*       When NBDCND = 3 or 4, */

/*            BDRS(I) = (d/dR)U(THETA(I),RS), I = 1,2,...,M+1  . */

/*       When NBDCND has any other value, BDRS is a dummy variable. */

/*     BDRF */
/*       A one-dimensional array of length M+1 that specifies the values */
/*       of the derivative of the solution with respect to R at R = RF. */
/*       When NBDCND = 2,3, or 6, */

/*            BDRF(I) = (d/dR)U(THETA(I),RF), I = 1,2,...,M+1  . */

/*       When NBDCND has any other value, BDRF is a dummy variable. */

/*     ELMBDA */
/*       The constant LAMBDA in the Helmholtz equation.  If */
/*       LAMBDA .GT. 0, a solution may not exist.  However, HWSCSP will */
/*       attempt to find a solution.  If NBDCND = 5 or 6 or */
/*       MBDCND = 5,6,7,8, or 9, ELMBDA must be zero. */

/*     F */
/*       A two-dimensional array that specifies the value of the right */
/*       side of the Helmholtz equation and boundary values (if any). */
/*       for I = 2,3,...,M and J = 2,3,...,N */

/*            F(I,J) = F(THETA(I),R(J)). */

/*       On the boundaries F is defined by */

/*            MBDCND   F(1,J)            F(M+1,J) */
/*            ------   ----------        ---------- */

/*              1      U(TS,R(J))        U(TF,R(J)) */
/*              2      U(TS,R(J))        F(TF,R(J)) */
/*              3      F(TS,R(J))        F(TF,R(J)) */
/*              4      F(TS,R(J))        U(TF,R(J)) */
/*              5      F(0,R(J))         U(TF,R(J))   J = 1,2,...,N+1 */
/*              6      F(0,R(J))         F(TF,R(J)) */
/*              7      U(TS,R(J))        F(PI,R(J)) */
/*              8      F(TS,R(J))        F(PI,R(J)) */
/*              9      F(0,R(J))         F(PI,R(J)) */

/*            NBDCND   F(I,1)            F(I,N+1) */
/*            ------   --------------    -------------- */

/*              1      U(THETA(I),RS)    U(THETA(I),RF) */
/*              2      U(THETA(I),RS)    F(THETA(I),RF) */
/*              3      F(THETA(I),RS)    F(THETA(I),RF) */
/*              4      F(THETA(I),RS)    U(THETA(I),RF)   I = 1,2,...,M+1 */
/*              5      F(TS,0)           U(THETA(I),RF) */
/*              6      F(TS,0)           F(THETA(I),RF) */

/*       F must be dimensioned at least (M+1)*(N+1). */

/*       NOTE */

/*       If the table calls for both the solution U and the right side F */
/*       at a corner then the solution must be specified. */

/*     IDIMF */
/*       The row (or first) dimension of the array F as it appears in the */
/*       program calling HWSCSP.  This parameter is used to specify the */
/*       variable dimension of F.  IDIMF must be at least M+1  . */

/*     W */
/*       A one-dimensional array that must be provided by the user for */
/*       work space. Its length can be computed from the formula below */
/*       which depends on the value of NBDCND. */

/*       If NBDCND=2,4 or 6 define NUNK=N */
/*       If NBDCND=1 or 5   define NUNK=N-1 */
/*       If NBDCND=3        define NUNK=N+1 */

/*       Now set K=INT(log2(NUNK))+1 and L=2**(K+1) then W must be */
/*       dimensioned at least (K-2)*L+K+5*(M+N)+MAX(2*N,6*M)+23 */

/*       **IMPORTANT** For purposes of checking, the required length */
/*                     of W is computed by HWSCSP and stored in W(1) */
/*                     in floating point format. */


/*     * * * * * * * * * *     On Output     * * * * * * * * * * */

/*     F */
/*       Contains the solution U(I,J) of the finite difference */
/*       approximation for the grid point (THETA(I),R(J)), */
/*       I = 1,2,...,M+1,   J = 1,2,...,N+1  . */

/*     PERTRB */
/*       If a combination of periodic or derivative boundary conditions */
/*       is specified for a Poisson equation (LAMBDA = 0), a solution may */
/*       not exist.  PERTRB is a constant, calculated and subtracted from */
/*       F, which ensures that a solution exists.  HWSCSP then computes */
/*       this solution, which is a least squares solution to the original */
/*       approximation. This solution is not unique and is unnormalized. */
/*       The value of PERTRB should be small compared to the right side */
/*       F. Otherwise , a solution is obtained to an essentially */
/*       different problem. This comparison should always be made to */
/*       insure that a meaningful solution has been obtained. */

/*     IERROR */
/*       An error flag that indicates invalid input parameters.  Except */
/*       for numbers 0 and 10, a solution is not attempted. */

/*       = 1  TS.LT.0. or TF.GT.PI */
/*       = 2  TS.GE.TF */
/*       = 3  M.LT.5 */
/*       = 4  MBDCND.LT.1 or MBDCND.GT.9 */
/*       = 5  RS.LT.0 */
/*       = 6  RS.GE.RF */
/*       = 7  N.LT.5 */
/*       = 8  NBDCND.LT.1 or NBDCND.GT.6 */
/*       = 9  ELMBDA.GT.0 */
/*       = 10 IDIMF.LT.M+1 */
/*       = 11 ELMBDA.NE.0 and MBDCND.GE.5 */
/*       = 12 ELMBDA.NE.0 and NBDCND equals 5 or 6 */
/*       = 13 MBDCND equals 5,6 or 9 and TS.NE.0 */
/*       = 14 MBDCND.GE.7 and TF.NE.PI */
/*       = 15 TS.EQ.0 and MBDCND equals 3,4 or 8 */
/*       = 16 TF.EQ.PI and MBDCND equals 2,3 or 6 */
/*       = 17 NBDCND.GE.5 and RS.NE.0 */
/*       = 18 NBDCND.GE.5 and MBDCND equals 1,2,4,5 or 7 */

/*       Since this is the only means of indicating a possibly incorrect */
/*       call to HWSCSP, the user should test IERROR after a call. */

/*     W */
/*       Contains intermediate values that must not be destroyed if */
/*       HWSCSP will be called again with INTL = 1.  W(1) contains the */
/*       number of locations which W must have. */

/* *Long Description: */

/*     * * * * * * *   Program Specifications    * * * * * * * * * * * * */

/*     Dimension of   BDTS(N+1),BDTF(N+1),BDRS(M+1),BDRF(M+1), */
/*     Arguments      F(IDIMF,N+1),W(see argument list) */

/*     Latest         June 1979 */
/*     Revision */

/*     Subprograms    HWSCSP,HWSCS1,BLKTRI,BLKTR1,PROD,PRODP,CPROD,CPRODP */
/*     Required       ,COMBP,PPADD,PSGF,BSRH,PPSGF,PPSPF,TEVLS,INDXA, */
/*                    ,INDXB,INDXC,R1MACH */

/*     Special */
/*     Conditions */

/*     Common         CBLKT */
/*     Blocks */

/*     I/O            NONE */

/*     Precision      Single */

/*     Specialist     Paul N Swarztrauber */

/*     Language       FORTRAN */

/*     History        Version 1 September 1973 */
/*                    Version 2 April     1976 */
/*                    Version 3 June      1979 */

/*     Algorithm      The routine defines the finite difference */
/*                    equations, incorporates boundary data, and adjusts */
/*                    the right side of singular systems and then calls */
/*                    BLKTRI to solve the system. */

/*     Space */
/*     Required */

/*     Portability    American National Standards Institute FORTRAN. */
/*                    The machine accuracy is set using function R1MACH. */

/*     Required       NONE */
/*     Resident */
/*     Routines */

/*     Reference      Swarztrauber,P. and R. Sweet, 'Efficient FORTRAN */
/*                    Subprograms for The Solution Of Elliptic Equations' */
/*                    NCAR TN/IA-109, July, 1975, 138 pp. */

/*     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/* ***REFERENCES  P. N. Swarztrauber and R. Sweet, Efficient Fortran */
/*                 subprograms for the solution of elliptic equations, */
/*                 NCAR TN/IA-109, July 1975, 138 pp. */
/* ***ROUTINES CALLED  HWSCS1, PIMACH */
/* ***REVISION HISTORY  (YYMMDD) */
/*   801001  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890531  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   920501  Reformatted the REFERENCES section.  (WRB) */
/* ***END PROLOGUE  HWSCSP */

/* ***FIRST EXECUTABLE STATEMENT  HWSCSP */
    /* Parameter adjustments */
    --bdts;
    --bdtf;
    --bdrs;
    --bdrf;
    f_dim1 = *idimf;
    f_offset = 1 + f_dim1;
    f -= f_offset;
    --w;

    /* Function Body */
    pi = pimach_(&dum);
    *ierror = 0;
    if (*ts < 0.f || *tf > pi) {
	*ierror = 1;
    }
    if (*ts >= *tf) {
	*ierror = 2;
    }
    if (*m < 5) {
	*ierror = 3;
    }
    if (*mbdcnd < 1 || *mbdcnd > 9) {
	*ierror = 4;
    }
    if (*rs < 0.f) {
	*ierror = 5;
    }
    if (*rs >= *rf) {
	*ierror = 6;
    }
    if (*n < 5) {
	*ierror = 7;
    }
    if (*nbdcnd < 1 || *nbdcnd > 6) {
	*ierror = 8;
    }
    if (*elmbda > 0.f) {
	*ierror = 9;
    }
    if (*idimf < *m + 1) {
	*ierror = 10;
    }
    if (*elmbda != 0.f && *mbdcnd >= 5) {
	*ierror = 11;
    }
    if (*elmbda != 0.f && (*nbdcnd == 5 || *nbdcnd == 6)) {
	*ierror = 12;
    }
    if ((*mbdcnd == 5 || *mbdcnd == 6 || *mbdcnd == 9) && *ts != 0.f) {
	*ierror = 13;
    }
    if (*mbdcnd >= 7 && *tf != pi) {
	*ierror = 14;
    }
    if (*ts == 0.f && (*mbdcnd == 4 || *mbdcnd == 8 || *mbdcnd == 3)) {
	*ierror = 15;
    }
    if (*tf == pi && (*mbdcnd == 2 || *mbdcnd == 3 || *mbdcnd == 6)) {
	*ierror = 16;
    }
    if (*nbdcnd >= 5 && *rs != 0.f) {
	*ierror = 17;
    }
    if (*nbdcnd >= 5 && (*mbdcnd == 1 || *mbdcnd == 2 || *mbdcnd == 5 || *
	    mbdcnd == 7)) {
	*ierror = 18;
    }
    if (*ierror != 0 && *ierror != 9) {
	return 0;
    }
    nck = *n;
    switch (*nbdcnd) {
	case 1:  goto L101;
	case 2:  goto L103;
	case 3:  goto L102;
	case 4:  goto L103;
	case 5:  goto L101;
	case 6:  goto L103;
    }
L101:
    --nck;
    goto L103;
L102:
    ++nck;
L103:
    l = 2;
    k = 1;
L104:
    l += l;
    ++k;
    if (nck - l <= 0) {
	goto L105;
    } else {
	goto L104;
    }
L105:
    l += l;
    np1 = *n + 1;
    mp1 = *m + 1;
/* Computing MAX */
    i__1 = *n << 1, i__2 = *m * 6;
    i1 = (k - 2) * l + k + max(i__1,i__2) + 13;
    i2 = i1 + np1;
    i3 = i2 + np1;
    i4 = i3 + np1;
    i5 = i4 + np1;
    i6 = i5 + np1;
    i7 = i6 + mp1;
    i8 = i7 + mp1;
    i9 = i8 + mp1;
    i10 = i9 + mp1;
    w[1] = (real) (i10 + *m);
    hwscs1_(intl, ts, tf, m, mbdcnd, &bdts[1], &bdtf[1], rs, rf, n, nbdcnd, &
	    bdrs[1], &bdrf[1], elmbda, &f[f_offset], idimf, pertrb, &w[2], &w[
	    i1], &w[i2], &w[i3], &w[i4], &w[i5], &w[i6], &w[i7], &w[i8], &w[
	    i9], &w[i10]);
    return 0;
} /* hwscsp_ */

