/* genbun.f -- translated by f2c (version 12.02.01).
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
static integer c__2 = 2;

/* DECK GENBUN */
/* Subroutine */ int genbun_(integer *nperod, integer *n, integer *mperod, 
	integer *m, real *a, real *b, real *c__, integer *idimy, real *y, 
	integer *ierror, real *w)
{
    /* System generated locals */
    integer y_dim1, y_offset, i__1, i__2;

    /* Local variables */
    static integer i__, j, k;
    static real a1;
    static integer mh, mp, np, mp1, iwd, iwp, mhm1, iwb2, iwb3, nby2, iww1, 
	    iww2, iww3, iwba, iwbb, iwbc, modd, mhmi, mhpi, irev, mskip;
    extern /* Subroutine */ int poisd2_(integer *, integer *, integer *, real 
	    *, real *, real *, real *, integer *, real *, real *, real *, 
	    real *, real *), poisn2_(integer *, integer *, integer *, integer 
	    *, real *, real *, real *, real *, integer *, real *, real *, 
	    real *, real *, real *, real *, real *, real *, real *), poisp2_(
	    integer *, integer *, real *, real *, real *, real *, integer *, 
	    real *, real *, real *, real *, real *, real *, real *, real *, 
	    real *);
    static integer iwtcos, ipstor;

/* ***BEGIN PROLOGUE  GENBUN */
/* ***PURPOSE  Solve by a cyclic reduction algorithm the linear system */
/*            of equations that results from a finite difference */
/*            approximation to certain 2-d elliptic PDE's on a centered */
/*            grid . */
/* ***LIBRARY   SLATEC (FISHPACK) */
/* ***CATEGORY  I2B4B */
/* ***TYPE      SINGLE PRECISION (GENBUN-S, CMGNBN-C) */
/* ***KEYWORDS  ELLIPTIC, FISHPACK, PDE, TRIDIAGONAL */
/* ***AUTHOR  Adams, J., (NCAR) */
/*           Swarztrauber, P. N., (NCAR) */
/*           Sweet, R., (NCAR) */
/* ***DESCRIPTION */

/*     Subroutine GENBUN solves the linear system of equations */

/*          A(I)*X(I-1,J) + B(I)*X(I,J) + C(I)*X(I+1,J) */

/*          + X(I,J-1) - 2.*X(I,J) + X(I,J+1) = Y(I,J) */

/*               for I = 1,2,...,M  and  J = 1,2,...,N. */

/*     The indices I+1 and I-1 are evaluated modulo M, i.e., */
/*     X(0,J) = X(M,J) and X(M+1,J) = X(1,J), and X(I,0) may be equal to */
/*     0, X(I,2), or X(I,N) and X(I,N+1) may be equal to 0, X(I,N-1), or */
/*     X(I,1) depending on an input parameter. */


/*     * * * * * * * *    Parameter Description     * * * * * * * * * * */

/*             * * * * * *   On Input    * * * * * * */

/*     NPEROD */
/*       Indicates the values that X(I,0) and X(I,N+1) are assumed to */
/*       have. */

/*       = 0  If X(I,0) = X(I,N) and X(I,N+1) = X(I,1). */
/*       = 1  If X(I,0) = X(I,N+1) = 0  . */
/*       = 2  If X(I,0) = 0 and X(I,N+1) = X(I,N-1). */
/*       = 3  If X(I,0) = X(I,2) and X(I,N+1) = X(I,N-1). */
/*       = 4  If X(I,0) = X(I,2) and X(I,N+1) = 0. */

/*     N */
/*       The number of unknowns in the J-direction.  N must be greater */
/*       than 2. */

/*     MPEROD */
/*       = 0 if A(1) and C(M) are not zero. */
/*       = 1 if A(1) = C(M) = 0. */

/*     M */
/*       The number of unknowns in the I-direction.  M must be greater */
/*       than 2. */

/*     A,B,C */
/*       One-dimensional arrays of length M that specify the */
/*       coefficients in the linear equations given above.  If MPEROD = 0 */
/*       the array elements must not depend upon the index I, but must be */
/*       constant.  Specifically, the subroutine checks the following */
/*       condition */

/*             A(I) = C(1) */
/*             C(I) = C(1) */
/*             B(I) = B(1) */

/*       for I=1,2,...,M. */

/*     IDIMY */
/*       The row (or first) dimension of the two-dimensional array Y as */
/*       it appears in the program calling GENBUN.  This parameter is */
/*       used to specify the variable dimension of Y.  IDIMY must be at */
/*       least M. */

/*     Y */
/*       A two-dimensional array that specifies the values of the right */
/*       side of the linear system of equations given above.  Y must be */
/*       dimensioned at least M*N. */

/*     W */
/*       A one-dimensional array that must be provided by the user for */
/*       work space.  W may require up to 4*N + (10 + INT(log2(N)))*M */
/*       locations.  The actual number of locations used is computed by */
/*       GENBUN and is returned in location W(1). */


/*             * * * * * *   On Output     * * * * * * */

/*     Y */
/*       Contains the solution X. */

/*     IERROR */
/*       An error flag that indicates invalid input parameters.  Except */
/*       for number zero, a solution is not attempted. */

/*       = 0  No error. */
/*       = 1  M .LE. 2 */
/*       = 2  N .LE. 2 */
/*       = 3  IDIMY .LT. M */
/*       = 4  NPEROD .LT. 0 or NPEROD .GT. 4 */
/*       = 5  MPEROD .LT. 0 or MPEROD .GT. 1 */
/*       = 6  A(I) .NE. C(1) or C(I) .NE. C(1) or B(I) .NE. B(1) for */
/*            some I=1,2,...,M. */
/*       = 7  A(1) .NE. 0 or C(M) .NE. 0 and MPEROD = 1 */

/*     W */
/*       W(1) contains the required length of W. */

/* *Long Description: */

/*     * * * * * * *   Program Specifications    * * * * * * * * * * * * */

/*     Dimension of   A(M),B(M),C(M),Y(IDIMY,N),W(see parameter list) */
/*     Arguments */

/*     Latest         June 1, 1976 */
/*     Revision */

/*     Subprograms    GENBUN,POISD2,POISN2,POISP2,COSGEN,MERGE,TRIX,TRI3, */
/*     Required       PIMACH */

/*     Special        NONE */
/*     Conditions */

/*     Common         NONE */
/*     Blocks */

/*     I/O            NONE */

/*     Precision      Single */

/*     Specialist     Roland Sweet */

/*     Language       FORTRAN */

/*     History        Standardized April 1, 1973 */
/*                    Revised August 20,1973 */
/*                    Revised January 1, 1976 */

/*     Algorithm      The linear system is solved by a cyclic reduction */
/*                    algorithm described in the reference. */

/*     Space          4944(decimal) = 11520(octal) locations on the NCAR */
/*     Required       Control Data 7600. */

/*     Timing and        The execution time T on the NCAR Control Data */
/*     Accuracy       7600 for subroutine GENBUN is roughly proportional */
/*                    to M*N*log2(N), but also depends on the input */
/*                    parameter NPEROD.  Some typical values are listed */
/*                    in the table below.  More comprehensive timing */
/*                    charts may be found in the reference. */
/*                       To measure the accuracy of the algorithm a */
/*                    uniform random number generator was used to create */
/*                    a solution array X for the system given in the */
/*                    'PURPOSE' with */

/*                       A(I) = C(I) = -0.5*B(I) = 1,       I=1,2,...,M */

/*                    and, when MPEROD = 1 */

/*                       A(1) = C(M) = 0 */
/*                       A(M) = C(1) = 2. */

/*                    The solution X was substituted into the given sys- */
/*                    tem and, using double precision, a right side Y was */
/*                    computed.  Using this array Y subroutine GENBUN was */
/*                    called to produce an approximate solution Z.  Then */
/*                    the relative error, defined as */

/*                       E = MAX(ABS(Z(I,J)-X(I,J)))/MAX(ABS(X(I,J))) */

/*                    where the two maxima are taken over all I=1,2,...,M */
/*                    and J=1,2,...,N, was computed.  The value of E is */
/*                    given in the table below for some typical values of */
/*                    M and N. */


/*                       M (=N)    MPEROD    NPEROD    T(MSECS)    E */
/*                       ------    ------    ------    --------  ------ */

/*                         31        0         0          36     6.E-14 */
/*                         31        1         1          21     4.E-13 */
/*                         31        1         3          41     3.E-13 */
/*                         32        0         0          29     9.E-14 */
/*                         32        1         1          32     3.E-13 */
/*                         32        1         3          48     1.E-13 */
/*                         33        0         0          36     9.E-14 */
/*                         33        1         1          30     4.E-13 */
/*                         33        1         3          34     1.E-13 */
/*                         63        0         0         150     1.E-13 */
/*                         63        1         1          91     1.E-12 */
/*                         63        1         3         173     2.E-13 */
/*                         64        0         0         122     1.E-13 */
/*                         64        1         1         128     1.E-12 */
/*                         64        1         3         199     6.E-13 */
/*                         65        0         0         143     2.E-13 */
/*                         65        1         1         120     1.E-12 */
/*                         65        1         3         138     4.E-13 */

/*     Portability    American National Standards Institute Fortran. */
/*                    The machine dependent constant PI is defined in */
/*                    function PIMACH. */

/*     Required       COS */
/*     Resident */
/*     Routines */

/*     Reference      Sweet, R., 'A Cyclic Reduction Algorithm For */
/*                    Solving Block Tridiagonal Systems Of Arbitrary */
/*                    Dimensions,' SIAM J. on Numer. Anal., */
/*                    14(Sept., 1977), PP. 706-720. */

/*     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/* ***REFERENCES  R. Sweet, A cyclic reduction algorithm for solving */
/*                 block tridiagonal systems of arbitrary dimensions, */
/*                 SIAM Journal on Numerical Analysis 14, (September */
/*                 1977), pp. 706-720. */
/* ***ROUTINES CALLED  POISD2, POISN2, POISP2 */
/* ***REVISION HISTORY  (YYMMDD) */
/*   801001  DATE WRITTEN */
/*   861211  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   920501  Reformatted the REFERENCES section.  (WRB) */
/* ***END PROLOGUE  GENBUN */


/* ***FIRST EXECUTABLE STATEMENT  GENBUN */
    /* Parameter adjustments */
    --a;
    --b;
    --c__;
    y_dim1 = *idimy;
    y_offset = 1 + y_dim1;
    y -= y_offset;
    --w;

    /* Function Body */
    *ierror = 0;
    if (*m <= 2) {
	*ierror = 1;
    }
    if (*n <= 2) {
	*ierror = 2;
    }
    if (*idimy < *m) {
	*ierror = 3;
    }
    if (*nperod < 0 || *nperod > 4) {
	*ierror = 4;
    }
    if (*mperod < 0 || *mperod > 1) {
	*ierror = 5;
    }
    if (*mperod == 1) {
	goto L102;
    }
    i__1 = *m;
    for (i__ = 2; i__ <= i__1; ++i__) {
	if (a[i__] != c__[1]) {
	    goto L103;
	}
	if (c__[i__] != c__[1]) {
	    goto L103;
	}
	if (b[i__] != b[1]) {
	    goto L103;
	}
/* L101: */
    }
    goto L104;
L102:
    if (a[1] != 0.f || c__[*m] != 0.f) {
	*ierror = 7;
    }
    goto L104;
L103:
    *ierror = 6;
L104:
    if (*ierror != 0) {
	return 0;
    }
    mp1 = *m + 1;
    iwba = mp1;
    iwbb = iwba + *m;
    iwbc = iwbb + *m;
    iwb2 = iwbc + *m;
    iwb3 = iwb2 + *m;
    iww1 = iwb3 + *m;
    iww2 = iww1 + *m;
    iww3 = iww2 + *m;
    iwd = iww3 + *m;
    iwtcos = iwd + *m;
    iwp = iwtcos + (*n << 2);
    i__1 = *m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	k = iwba + i__ - 1;
	w[k] = -a[i__];
	k = iwbc + i__ - 1;
	w[k] = -c__[i__];
	k = iwbb + i__ - 1;
	w[k] = 2.f - b[i__];
	i__2 = *n;
	for (j = 1; j <= i__2; ++j) {
	    y[i__ + j * y_dim1] = -y[i__ + j * y_dim1];
/* L105: */
	}
/* L106: */
    }
    mp = *mperod + 1;
    np = *nperod + 1;
    switch (mp) {
	case 1:  goto L114;
	case 2:  goto L107;
    }
L107:
    switch (np) {
	case 1:  goto L108;
	case 2:  goto L109;
	case 3:  goto L110;
	case 4:  goto L111;
	case 5:  goto L123;
    }
L108:
    poisp2_(m, n, &w[iwba], &w[iwbb], &w[iwbc], &y[y_offset], idimy, &w[1], &
	    w[iwb2], &w[iwb3], &w[iww1], &w[iww2], &w[iww3], &w[iwd], &w[
	    iwtcos], &w[iwp]);
    goto L112;
L109:
    poisd2_(m, n, &c__1, &w[iwba], &w[iwbb], &w[iwbc], &y[y_offset], idimy, &
	    w[1], &w[iww1], &w[iwd], &w[iwtcos], &w[iwp]);
    goto L112;
L110:
    poisn2_(m, n, &c__1, &c__2, &w[iwba], &w[iwbb], &w[iwbc], &y[y_offset], 
	    idimy, &w[1], &w[iwb2], &w[iwb3], &w[iww1], &w[iww2], &w[iww3], &
	    w[iwd], &w[iwtcos], &w[iwp]);
    goto L112;
L111:
    poisn2_(m, n, &c__1, &c__1, &w[iwba], &w[iwbb], &w[iwbc], &y[y_offset], 
	    idimy, &w[1], &w[iwb2], &w[iwb3], &w[iww1], &w[iww2], &w[iww3], &
	    w[iwd], &w[iwtcos], &w[iwp]);
L112:
    ipstor = w[iww1];
    irev = 2;
    if (*nperod == 4) {
	goto L124;
    }
L113:
    switch (mp) {
	case 1:  goto L127;
	case 2:  goto L133;
    }
L114:

/*     REORDER UNKNOWNS WHEN MP =0 */

    mh = (*m + 1) / 2;
    mhm1 = mh - 1;
    modd = 1;
    if (mh << 1 == *m) {
	modd = 2;
    }
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	i__2 = mhm1;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    mhpi = mh + i__;
	    mhmi = mh - i__;
	    w[i__] = y[mhmi + j * y_dim1] - y[mhpi + j * y_dim1];
	    w[mhpi] = y[mhmi + j * y_dim1] + y[mhpi + j * y_dim1];
/* L115: */
	}
	w[mh] = y[mh + j * y_dim1] * 2.f;
	switch (modd) {
	    case 1:  goto L117;
	    case 2:  goto L116;
	}
L116:
	w[*m] = y[*m + j * y_dim1] * 2.f;
L117:
	i__2 = *m;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    y[i__ + j * y_dim1] = w[i__];
/* L118: */
	}
/* L119: */
    }
    k = iwbc + mhm1 - 1;
    i__ = iwba + mhm1;
    w[k] = 0.f;
    w[i__] = 0.f;
    w[k + 1] *= 2.f;
    switch (modd) {
	case 1:  goto L120;
	case 2:  goto L121;
    }
L120:
    k = iwbb + mhm1 - 1;
    w[k] -= w[i__ - 1];
    w[iwbc - 1] += w[iwbb - 1];
    goto L122;
L121:
    w[iwbb - 1] = w[k + 1];
L122:
    goto L107;

/*     REVERSE COLUMNS WHEN NPEROD = 4. */

L123:
    irev = 1;
    nby2 = *n / 2;
L124:
    i__1 = nby2;
    for (j = 1; j <= i__1; ++j) {
	mskip = *n + 1 - j;
	i__2 = *m;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    a1 = y[i__ + j * y_dim1];
	    y[i__ + j * y_dim1] = y[i__ + mskip * y_dim1];
	    y[i__ + mskip * y_dim1] = a1;
/* L125: */
	}
/* L126: */
    }
    switch (irev) {
	case 1:  goto L110;
	case 2:  goto L113;
    }
L127:
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	i__2 = mhm1;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    mhmi = mh - i__;
	    mhpi = mh + i__;
	    w[mhmi] = (y[mhpi + j * y_dim1] + y[i__ + j * y_dim1]) * .5f;
	    w[mhpi] = (y[mhpi + j * y_dim1] - y[i__ + j * y_dim1]) * .5f;
/* L128: */
	}
	w[mh] = y[mh + j * y_dim1] * .5f;
	switch (modd) {
	    case 1:  goto L130;
	    case 2:  goto L129;
	}
L129:
	w[*m] = y[*m + j * y_dim1] * .5f;
L130:
	i__2 = *m;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    y[i__ + j * y_dim1] = w[i__];
/* L131: */
	}
/* L132: */
    }
L133:

/*     RETURN STORAGE REQUIREMENTS FOR W ARRAY. */

    w[1] = (real) (ipstor + iwp - 1);
    return 0;
} /* genbun_ */

