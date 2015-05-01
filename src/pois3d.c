/* pois3d.f -- translated by f2c (version 12.02.01).
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

/* DECK POIS3D */
/* Subroutine */ int pois3d_(integer *lperod, integer *l, real *c1, integer *
	mperod, integer *m, real *c2, integer *nperod, integer *n, real *a, 
	real *b, real *c__, integer *ldimf, integer *mdimf, real *f, integer *
	ierror, real *w)
{
    /* System generated locals */
    integer f_dim1, f_dim2, f_offset, i__1, i__2, i__3;

    /* Local variables */
    static integer i__, j, k, nh, lp, mp, np, iwd, iwt, iwx, iwy, nhm1, iwbb, 
	    nodd, nhmk;
    static real save[6];
    static integer nhpk;
    extern /* Subroutine */ int pos3d1_(integer *, integer *, integer *, 
	    integer *, integer *, real *, real *, real *, integer *, integer *
	    , real *, real *, real *, real *, real *, real *, real *, real *, 
	    real *, real *);
    static integer iwyrt;

/* ***BEGIN PROLOGUE  POIS3D */
/* ***PURPOSE  Solve a three-dimensional block tridiagonal linear system */
/*            which arises from a finite difference approximation to a */
/*            three-dimensional Poisson equation using the Fourier */
/*            transform package FFTPAK written by Paul Swarztrauber. */
/* ***LIBRARY   SLATEC (FISHPACK) */
/* ***CATEGORY  I2B4B */
/* ***TYPE      SINGLE PRECISION (POIS3D-S) */
/* ***KEYWORDS  ELLIPTIC PDE, FISHPACK, HELMHOLTZ, POISSON */
/* ***AUTHOR  Adams, J., (NCAR) */
/*           Swarztrauber, P. N., (NCAR) */
/*           Sweet, R., (NCAR) */
/* ***DESCRIPTION */

/*     Subroutine POIS3D solves the linear system of equations */

/*       C1*(X(I-1,J,K)-2.*X(I,J,K)+X(I+1,J,K)) */
/*     + C2*(X(I,J-1,K)-2.*X(I,J,K)+X(I,J+1,K)) */
/*     + A(K)*X(I,J,K-1)+B(K)*X(I,J,K)+C(K)*X(I,J,K+1) = F(I,J,K) */

/*     for  I=1,2,...,L , J=1,2,...,M , and K=1,2,...,N . */

/*     The indices K-1 and K+1 are evaluated modulo N, i.e. */
/*     X(I,J,0) = X(I,J,N) and X(I,J,N+1) = X(I,J,1). The unknowns */
/*     X(0,J,K), X(L+1,J,K), X(I,0,K), and X(I,M+1,K) are assumed to take */
/*     on certain prescribed values described below. */

/*    * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */


/*    * * * * * * * *    Parameter Description     * * * * * * * * * * */


/*            * * * * * *   On Input    * * * * * * */

/*     LPEROD   Indicates the values that X(0,J,K) and X(L+1,J,K) are */
/*              assumed to have. */

/*              = 0  If X(0,J,K) = X(L,J,K) and X(L+1,J,K) = X(1,J,K). */
/*              = 1  If X(0,J,K) = X(L+1,J,K) = 0. */
/*              = 2  If X(0,J,K) = 0  and X(L+1,J,K) = X(L-1,J,K). */
/*              = 3  If X(0,J,K) = X(2,J,K) and X(L+1,J,K) = X(L-1,J,K). */
/*              = 4  If X(0,J,K) = X(2,J,K) and X(L+1,J,K) = 0. */

/*     L        The number of unknowns in the I-direction. L must be at */
/*              least 3. */

/*     C1       The real constant that appears in the above equation. */

/*     MPEROD   Indicates the values that X(I,0,K) and X(I,M+1,K) are */
/*              assumed to have. */

/*              = 0  If X(I,0,K) = X(I,M,K) and X(I,M+1,K) = X(I,1,K). */
/*              = 1  If X(I,0,K) = X(I,M+1,K) = 0. */
/*              = 2  If X(I,0,K) = 0 and X(I,M+1,K) = X(I,M-1,K). */
/*              = 3  If X(I,0,K) = X(I,2,K) and X(I,M+1,K) = X(I,M-1,K). */
/*              = 4  If X(I,0,K) = X(I,2,K) and X(I,M+1,K) = 0. */

/*     M        The number of unknowns in the J-direction. M must be at */
/*              least 3. */

/*     C2       The real constant which appears in the above equation. */

/*     NPEROD   = 0  If A(1) and C(N) are not zero. */
/*              = 1  If A(1) = C(N) = 0. */

/*     N        The number of unknowns in the K-direction. N must be at */
/*              least 3. */


/*     A,B,C    One-dimensional arrays of length N that specify the */
/*              coefficients in the linear equations given above. */

/*              If NPEROD = 0 the array elements must not depend upon the */
/*              index K, but must be constant.  Specifically, the */
/*              subroutine checks the following condition */

/*                          A(K) = C(1) */
/*                          C(K) = C(1) */
/*                          B(K) = B(1) */

/*                  for K=1,2,...,N. */

/*     LDIMF    The row (or first) dimension of the three-dimensional */
/*              array F as it appears in the program calling POIS3D. */
/*              This parameter is used to specify the variable dimension */
/*              of F.  LDIMF must be at least L. */

/*     MDIMF    The column (or second) dimension of the three-dimensional */
/*              array F as it appears in the program calling POIS3D. */
/*              This parameter is used to specify the variable dimension */
/*              of F.  MDIMF must be at least M. */

/*     F        A three-dimensional array that specifies the values of */
/*              the right side of the linear system of equations given */
/*              above.  F must be dimensioned at least L x M x N. */

/*     W        A one-dimensional array that must be provided by the */
/*              user for work space.  The length of W must be at least */
/*              30 + L + M + 2*N + MAX(L,M,N) + */
/*              7*(INT((L+1)/2) + INT((M+1)/2)). */


/*            * * * * * *   On Output   * * * * * * */

/*     F        Contains the solution X. */

/*     IERROR   An error flag that indicates invalid input parameters. */
/*              Except for number zero, a solution is not attempted. */
/*              = 0  No error */
/*              = 1  If LPEROD .LT. 0 or .GT. 4 */
/*              = 2  If L .LT. 3 */
/*              = 3  If MPEROD .LT. 0 or .GT. 4 */
/*              = 4  If M .LT. 3 */
/*              = 5  If NPEROD .LT. 0 or .GT. 1 */
/*              = 6  If N .LT. 3 */
/*              = 7  If LDIMF .LT. L */
/*              = 8  If MDIMF .LT. M */
/*              = 9  If A(K) .NE. C(1) or C(K) .NE. C(1) or B(I) .NE.B(1) */
/*                      for some K=1,2,...,N. */
/*              = 10 If NPEROD = 1 and A(1) .NE. 0 or C(N) .NE. 0 */

/*              Since this is the only means of indicating a possibly */
/*              incorrect call to POIS3D, the user should test IERROR */
/*              after the call. */

/* *Long Description: */

/*    * * * * * * *   Program Specifications    * * * * * * * * * * * * */

/*     Dimension of   A(N),B(N),C(N),F(LDIMF,MDIMF,N), */
/*     Arguments      W(see argument list) */

/*     Latest         December 1, 1978 */
/*     Revision */

/*     Subprograms    POIS3D,POS3D1,TRIDQ,RFFTI,RFFTF,RFFTF1,RFFTB, */
/*     Required       RFFTB1,COSTI,COST,SINTI,SINT,COSQI,COSQF,COSQF1 */
/*                    COSQB,COSQB1,SINQI,SINQF,SINQB,CFFTI,CFFTI1, */
/*                    CFFTB,CFFTB1,PASSB2,PASSB3,PASSB4,PASSB,CFFTF, */
/*                    CFFTF1,PASSF1,PASSF2,PASSF3,PASSF4,PASSF,PIMACH, */

/*     Special        NONE */
/*     Conditions */

/*     Common         NONE */
/*     Blocks */

/*     I/O            NONE */

/*     Precision      Single */

/*     Specialist     Roland Sweet */

/*     Language       FORTRAN */

/*     History        Written by Roland Sweet at NCAR in July 1977 */

/*     Algorithm      This subroutine solves three-dimensional block */
/*                    tridiagonal linear systems arising from finite */
/*                    difference approximations to three-dimensional */
/*                    Poisson equations using the Fourier transform */
/*                    package FFTPAK written by Paul Swarztrauber. */

/*     Space          6561(decimal) = 14641(octal) locations on the */
/*     Required       NCAR Control Data 7600 */

/*     Timing and        The execution time T on the NCAR Control Data */
/*     Accuracy       7600 for subroutine POIS3D is roughly proportional */
/*                    to L*M*N*(log2(L)+log2(M)+5), but also depends on */
/*                    input parameters LPEROD and MPEROD.  Some typical */
/*                    values are listed in the table below when NPEROD=0. */
/*                       To measure the accuracy of the algorithm a */
/*                    uniform random number generator was used to create */
/*                    a solution array X for the system given in the */
/*                    'PURPOSE' with */

/*                       A(K) = C(K) = -0.5*B(K) = 1,       K=1,2,...,N */

/*                    and, when NPEROD = 1 */

/*                       A(1) = C(N) = 0 */
/*                       A(N) = C(1) = 2. */

/*                    The solution X was substituted into the given sys- */
/*                    tem and, using double precision, a right side Y was */
/*                    computed.  Using this array Y subroutine POIS3D was */
/*                    called to produce an approximate solution Z.  Then */
/*                    the relative error, defined as */

/*                    E = MAX(ABS(Z(I,J,K)-X(I,J,K)))/MAX(ABS(X(I,J,K))) */

/*                    where the two maxima are taken over I=1,2,...,L, */
/*                    J=1,2,...,M and K=1,2,...,N, was computed.  The */
/*                    value of E is given in the table below for some */
/*                    typical values of L,M and N. */


/*                       L(=M=N)   LPEROD    MPEROD    T(MSECS)    E */
/*                       ------    ------    ------    --------  ------ */

/*                         16        0         0         272     1.E-13 */
/*                         15        1         1         287     4.E-13 */
/*                         17        3         3         338     2.E-13 */
/*                         32        0         0        1755     2.E-13 */
/*                         31        1         1        1894     2.E-12 */
/*                         33        3         3        2042     7.E-13 */


/*     Portability    American National Standards Institute FORTRAN. */
/*                    The machine dependent constant PI is defined in */
/*                    function PIMACH. */

/*     Required       COS,SIN,ATAN */
/*     Resident */
/*     Routines */

/*     Reference      NONE */

/*    * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/* ***REFERENCES  (NONE) */
/* ***ROUTINES CALLED  POS3D1 */
/* ***REVISION HISTORY  (YYMMDD) */
/*   801001  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890531  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/* ***END PROLOGUE  POIS3D */
/* ***FIRST EXECUTABLE STATEMENT  POIS3D */
    /* Parameter adjustments */
    --a;
    --b;
    --c__;
    f_dim1 = *ldimf;
    f_dim2 = *mdimf;
    f_offset = 1 + f_dim1 * (1 + f_dim2);
    f -= f_offset;
    --w;

    /* Function Body */
    lp = *lperod + 1;
    mp = *mperod + 1;
    np = *nperod + 1;

/*     CHECK FOR INVALID INPUT. */

    *ierror = 0;
    if (lp < 1 || lp > 5) {
	*ierror = 1;
    }
    if (*l < 3) {
	*ierror = 2;
    }
    if (mp < 1 || mp > 5) {
	*ierror = 3;
    }
    if (*m < 3) {
	*ierror = 4;
    }
    if (np < 1 || np > 2) {
	*ierror = 5;
    }
    if (*n < 3) {
	*ierror = 6;
    }
    if (*ldimf < *l) {
	*ierror = 7;
    }
    if (*mdimf < *m) {
	*ierror = 8;
    }
    if (np != 1) {
	goto L103;
    }
    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
	if (a[k] != c__[1]) {
	    goto L102;
	}
	if (c__[k] != c__[1]) {
	    goto L102;
	}
	if (b[k] != b[1]) {
	    goto L102;
	}
/* L101: */
    }
    goto L104;
L102:
    *ierror = 9;
L103:
    if (*nperod == 1 && (a[1] != 0.f || c__[*n] != 0.f)) {
	*ierror = 10;
    }
L104:
    if (*ierror != 0) {
	goto L122;
    }
    iwyrt = *l + 1;
    iwt = iwyrt + *m;
/* Computing MAX */
    i__1 = max(*l,*m);
    iwd = iwt + max(i__1,*n) + 1;
    iwbb = iwd + *n;
    iwx = iwbb + *n;
    iwy = iwx + (*l + 1) / 2 * 7 + 15;
    switch (np) {
	case 1:  goto L105;
	case 2:  goto L114;
    }

/*     REORDER UNKNOWNS WHEN NPEROD = 0. */

L105:
    nh = (*n + 1) / 2;
    nhm1 = nh - 1;
    nodd = 1;
    if (nh << 1 == *n) {
	nodd = 2;
    }
    i__1 = *l;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = *m;
	for (j = 1; j <= i__2; ++j) {
	    i__3 = nhm1;
	    for (k = 1; k <= i__3; ++k) {
		nhpk = nh + k;
		nhmk = nh - k;
		w[k] = f[i__ + (j + nhmk * f_dim2) * f_dim1] - f[i__ + (j + 
			nhpk * f_dim2) * f_dim1];
		w[nhpk] = f[i__ + (j + nhmk * f_dim2) * f_dim1] + f[i__ + (j 
			+ nhpk * f_dim2) * f_dim1];
/* L106: */
	    }
	    w[nh] = f[i__ + (j + nh * f_dim2) * f_dim1] * 2.f;
	    switch (nodd) {
		case 1:  goto L108;
		case 2:  goto L107;
	    }
L107:
	    w[*n] = f[i__ + (j + *n * f_dim2) * f_dim1] * 2.f;
L108:
	    i__3 = *n;
	    for (k = 1; k <= i__3; ++k) {
		f[i__ + (j + k * f_dim2) * f_dim1] = w[k];
/* L109: */
	    }
/* L110: */
	}
/* L111: */
    }
    save[0] = c__[nhm1];
    save[1] = a[nh];
    save[2] = c__[nh];
    save[3] = b[nhm1];
    save[4] = b[*n];
    save[5] = a[*n];
    c__[nhm1] = 0.f;
    a[nh] = 0.f;
    c__[nh] *= 2.f;
    switch (nodd) {
	case 1:  goto L112;
	case 2:  goto L113;
    }
L112:
    b[nhm1] -= a[nh - 1];
    b[*n] += a[*n];
    goto L114;
L113:
    a[*n] = c__[nh];
L114:
    pos3d1_(&lp, l, &mp, m, n, &a[1], &b[1], &c__[1], ldimf, mdimf, &f[
	    f_offset], &w[1], &w[iwyrt], &w[iwt], &w[iwd], &w[iwx], &w[iwy], 
	    c1, c2, &w[iwbb]);
    switch (np) {
	case 1:  goto L115;
	case 2:  goto L122;
    }
L115:
    i__1 = *l;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = *m;
	for (j = 1; j <= i__2; ++j) {
	    i__3 = nhm1;
	    for (k = 1; k <= i__3; ++k) {
		nhmk = nh - k;
		nhpk = nh + k;
		w[nhmk] = (f[i__ + (j + nhpk * f_dim2) * f_dim1] + f[i__ + (j 
			+ k * f_dim2) * f_dim1]) * .5f;
		w[nhpk] = (f[i__ + (j + nhpk * f_dim2) * f_dim1] - f[i__ + (j 
			+ k * f_dim2) * f_dim1]) * .5f;
/* L116: */
	    }
	    w[nh] = f[i__ + (j + nh * f_dim2) * f_dim1] * .5f;
	    switch (nodd) {
		case 1:  goto L118;
		case 2:  goto L117;
	    }
L117:
	    w[*n] = f[i__ + (j + *n * f_dim2) * f_dim1] * .5f;
L118:
	    i__3 = *n;
	    for (k = 1; k <= i__3; ++k) {
		f[i__ + (j + k * f_dim2) * f_dim1] = w[k];
/* L119: */
	    }
/* L120: */
	}
/* L121: */
    }
    c__[nhm1] = save[0];
    a[nh] = save[1];
    c__[nh] = save[2];
    b[nhm1] = save[3];
    b[*n] = save[4];
    a[*n] = save[5];
L122:
    return 0;
} /* pois3d_ */

