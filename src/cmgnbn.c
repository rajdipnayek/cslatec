/* cmgnbn.f -- translated by f2c (version 12.02.01).
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

/* DECK CMGNBN */
/* Subroutine */ int cmgnbn_(integer *nperod, integer *n, integer *mperod, 
	integer *m, complex *a, complex *b, complex *c__, integer *idimy, 
	complex *y, integer *ierror, complex *w)
{
    /* System generated locals */
    integer y_dim1, y_offset, i__1, i__2, i__3, i__4, i__5;
    real r__1;
    complex q__1, q__2;

    /* Local variables */
    static integer i__, j, k;
    static complex a1;
    static integer mh, mp, np, iwd, iwp, mhm1, iwb2, iwb3, nby2, iww1, iww2, 
	    iww3, iwba, iwbb, iwbc, modd, mhmi, mhpi, irev, mskip;
    extern /* Subroutine */ int cmposd_(integer *, integer *, integer *, 
	    complex *, complex *, complex *, complex *, integer *, complex *, 
	    complex *, complex *, complex *, complex *), cmposn_(integer *, 
	    integer *, integer *, integer *, complex *, complex *, complex *, 
	    complex *, integer *, complex *, complex *, complex *, complex *, 
	    complex *, complex *, complex *, complex *, complex *), cmposp_(
	    integer *, integer *, complex *, complex *, complex *, complex *, 
	    integer *, complex *, complex *, complex *, complex *, complex *, 
	    complex *, complex *, complex *, complex *);
    static integer iwtcos, ipstor;

/* ***BEGIN PROLOGUE  CMGNBN */
/* ***PURPOSE  Solve a complex block tridiagonal linear system of */
/*            equations by a cyclic reduction algorithm. */
/* ***LIBRARY   SLATEC (FISHPACK) */
/* ***CATEGORY  I2B4B */
/* ***TYPE      COMPLEX (GENBUN-S, CMGNBN-C) */
/* ***KEYWORDS  CYCLIC REDUCTION, ELLIPTIC PDE, FISHPACK, */
/*             TRIDIAGONAL LINEAR SYSTEM */
/* ***AUTHOR  Adams, J., (NCAR) */
/*           Swarztrauber, P. N., (NCAR) */
/*           Sweet, R., (NCAR) */
/* ***DESCRIPTION */

/*     Subroutine CMGNBN solves the complex linear system of equations */

/*          A(I)*X(I-1,J) + B(I)*X(I,J) + C(I)*X(I+1,J) */

/*          + X(I,J-1) - 2.*X(I,J) + X(I,J+1) = Y(I,J) */

/*               For I = 1,2,...,M  and  J = 1,2,...,N. */

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
/*       = 0 If A(1) and C(M) are not zero */
/*       = 1 If A(1) = C(M) = 0 */

/*     M */
/*       The number of unknowns in the I-direction.  N must be greater */
/*       than 2. */

/*     A,B,C */
/*       One-dimensional complex arrays of length M that specify the */
/*       coefficients in the linear equations given above.  If MPEROD = 0 */
/*       the array elements must not depend upon the index I, but must be */
/*       constant.  Specifically, the subroutine checks the following */
/*       condition */

/*             A(I) = C(1) */
/*             C(I) = C(1) */
/*             B(I) = B(1) */

/*       For I=1,2,...,M. */

/*     IDIMY */
/*       The row (or first) dimension of the two-dimensional array Y as */
/*       it appears in the program calling CMGNBN.  This parameter is */
/*       used to specify the variable dimension of Y.  IDIMY must be at */
/*       least M. */

/*     Y */
/*       A two-dimensional complex array that specifies the values of the */
/*       right side of the linear system of equations given above.  Y */
/*       must be dimensioned at least M*N. */

/*     W */
/*       A one-dimensional complex array that must be provided by the */
/*       user for work space.  W may require up to 4*N + */
/*       (10 + INT(log2(N)))*M LOCATIONS.  The actual number of locations */
/*       used is computed by CMGNBN and is returned in location W(1). */


/*             * * * * * *   On Output     * * * * * * */

/*     Y */
/*       Contains the solution X. */

/*     IERROR */
/*       An error flag which indicates invalid input parameters.  Except */
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

/*     Latest         June 1979 */
/*     Revision */

/*     Subprograms    CMGNBN,CMPOSD,CMPOSN,CMPOSP,CMPCSG,CMPMRG, */
/*     Required       CMPTRX,CMPTR3,PIMACH */

/*     Special        None */
/*     Conditions */

/*     Common         None */
/*     Blocks */

/*     I/O            None */

/*     Precision      Single */

/*     Specialist     Roland Sweet */

/*     Language       FORTRAN */

/*     History        Written by Roland Sweet at NCAR in June, 1977 */

/*     Algorithm      The linear system is solved by a cyclic reduction */
/*                    algorithm described in the reference. */

/*     Space          4944(DECIMAL) = 11520(octal) locations on the NCAR */
/*     Required       Control Data 7600 */

/*     Timing and      The execution time T on the NCAR Control Data */
/*     Accuracy       7600 for subroutine CMGNBN is roughly proportional */
/*                    to M*N*log2(N), but also depends on the input */
/*                    parameter NPEROD.  Some typical values are listed */
/*                    in the table below. */
/*                       To measure the accuracy of the algorithm a */
/*                    uniform random number generator was used to create */
/*                    a solution array X for the system given in the */
/*                    'PURPOSE' with */

/*                       A(I) = C(I) = -0.5*B(I) = 1,       I=1,2,...,M */

/*                    and, when MPEROD = 1 */

/*                       A(1) = C(M) = 0 */
/*                       A(M) = C(1) = 2. */

/*                    The solution X was substituted into the given sys- */
/*                    tem and a right side Y was computed.  Using this */
/*                    array Y subroutine CMGNBN was called to produce an */
/*                    approximate solution Z.  Then the relative error, */
/*                    defined as */

/*                       E = MAX(ABS(Z(I,J)-X(I,J)))/MAX(ABS(X(I,J))) */

/*                    where the two maxima are taken over all I=1,2,...,M */
/*                    and J=1,2,...,N, was computed.  The value of E is */
/*                    given in the table below for some typical values of */
/*                    M and N. */


/*                       M (=N)    MPEROD    NPEROD    T(MSECS)    E */
/*                       ------    ------    ------    --------  ------ */

/*                         31        0         0          77     1.E-12 */
/*                         31        1         1          45     4.E-13 */
/*                         31        1         3          91     2.E-12 */
/*                         32        0         0          59     7.E-14 */
/*                         32        1         1          65     5.E-13 */
/*                         32        1         3          97     2.E-13 */
/*                         33        0         0          80     6.E-13 */
/*                         33        1         1          67     5.E-13 */
/*                         33        1         3          76     3.E-12 */
/*                         63        0         0         350     5.E-12 */
/*                         63        1         1         215     6.E-13 */
/*                         63        1         3         412     1.E-11 */
/*                         64        0         0         264     1.E-13 */
/*                         64        1         1         287     3.E-12 */
/*                         64        1         3         421     3.E-13 */
/*                         65        0         0         338     2.E-12 */
/*                         65        1         1         292     5.E-13 */
/*                         65        1         3         329     1.E-11 */

/*     Portability    American National Standards Institute Fortran. */
/*                    The machine dependent constant PI is defined in */
/*                    function PIMACH. */

/*     Required       COS */
/*     Resident */
/*     Routines */

/*     Reference      Sweet, R., 'A Cyclic Reduction Algorithm for */
/*                    Solving Block Tridiagonal Systems Of Arbitrary */
/*                    Dimensions,' SIAM J. on Numer. Anal., */
/*                    14(SEPT., 1977), PP. 706-720. */

/*     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/* ***REFERENCES  R. Sweet, A cyclic reduction algorithm for solving */
/*                 block tridiagonal systems of arbitrary dimensions, */
/*                 SIAM Journal on Numerical Analysis 14, (September */
/*                 1977), pp. 706-720. */
/* ***ROUTINES CALLED  CMPOSD, CMPOSN, CMPOSP */
/* ***REVISION HISTORY  (YYMMDD) */
/*   801001  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890531  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   920501  Reformatted the REFERENCES section.  (WRB) */
/* ***END PROLOGUE  CMGNBN */


/* ***FIRST EXECUTABLE STATEMENT  CMGNBN */
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
	i__2 = i__;
	q__1.r = a[i__2].r - c__[1].r, q__1.i = a[i__2].i - c__[1].i;
	if (c_abs(&q__1) != 0.f) {
	    goto L103;
	}
	i__2 = i__;
	q__1.r = c__[i__2].r - c__[1].r, q__1.i = c__[i__2].i - c__[1].i;
	if (c_abs(&q__1) != 0.f) {
	    goto L103;
	}
	i__2 = i__;
	q__1.r = b[i__2].r - b[1].r, q__1.i = b[i__2].i - b[1].i;
	if (c_abs(&q__1) != 0.f) {
	    goto L103;
	}
/* L101: */
    }
    goto L104;
L102:
    if (c_abs(&a[1]) != 0.f && c_abs(&c__[*m]) != 0.f) {
	*ierror = 7;
    }
    goto L104;
L103:
    *ierror = 6;
L104:
    if (*ierror != 0) {
	return 0;
    }
    iwba = *m + 1;
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
	i__2 = k;
	i__3 = i__;
	q__1.r = -a[i__3].r, q__1.i = -a[i__3].i;
	w[i__2].r = q__1.r, w[i__2].i = q__1.i;
	k = iwbc + i__ - 1;
	i__2 = k;
	i__3 = i__;
	q__1.r = -c__[i__3].r, q__1.i = -c__[i__3].i;
	w[i__2].r = q__1.r, w[i__2].i = q__1.i;
	k = iwbb + i__ - 1;
	i__2 = k;
	i__3 = i__;
	q__1.r = 2.f - b[i__3].r, q__1.i = -b[i__3].i;
	w[i__2].r = q__1.r, w[i__2].i = q__1.i;
	i__2 = *n;
	for (j = 1; j <= i__2; ++j) {
	    i__3 = i__ + j * y_dim1;
	    i__4 = i__ + j * y_dim1;
	    q__1.r = -y[i__4].r, q__1.i = -y[i__4].i;
	    y[i__3].r = q__1.r, y[i__3].i = q__1.i;
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
    cmposp_(m, n, &w[iwba], &w[iwbb], &w[iwbc], &y[y_offset], idimy, &w[1], &
	    w[iwb2], &w[iwb3], &w[iww1], &w[iww2], &w[iww3], &w[iwd], &w[
	    iwtcos], &w[iwp]);
    goto L112;
L109:
    cmposd_(m, n, &c__1, &w[iwba], &w[iwbb], &w[iwbc], &y[y_offset], idimy, &
	    w[1], &w[iww1], &w[iwd], &w[iwtcos], &w[iwp]);
    goto L112;
L110:
    cmposn_(m, n, &c__1, &c__2, &w[iwba], &w[iwbb], &w[iwbc], &y[y_offset], 
	    idimy, &w[1], &w[iwb2], &w[iwb3], &w[iww1], &w[iww2], &w[iww3], &
	    w[iwd], &w[iwtcos], &w[iwp]);
    goto L112;
L111:
    cmposn_(m, n, &c__1, &c__1, &w[iwba], &w[iwbb], &w[iwbc], &y[y_offset], 
	    idimy, &w[1], &w[iwb2], &w[iwb3], &w[iww1], &w[iww2], &w[iww3], &
	    w[iwd], &w[iwtcos], &w[iwp]);
L112:
    i__1 = iww1;
    ipstor = w[i__1].r;
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
	    i__3 = i__;
	    i__4 = mhmi + j * y_dim1;
	    i__5 = mhpi + j * y_dim1;
	    q__1.r = y[i__4].r - y[i__5].r, q__1.i = y[i__4].i - y[i__5].i;
	    w[i__3].r = q__1.r, w[i__3].i = q__1.i;
	    i__3 = mhpi;
	    i__4 = mhmi + j * y_dim1;
	    i__5 = mhpi + j * y_dim1;
	    q__1.r = y[i__4].r + y[i__5].r, q__1.i = y[i__4].i + y[i__5].i;
	    w[i__3].r = q__1.r, w[i__3].i = q__1.i;
/* L115: */
	}
	i__2 = mh;
	i__3 = mh + j * y_dim1;
	q__1.r = y[i__3].r * 2.f, q__1.i = y[i__3].i * 2.f;
	w[i__2].r = q__1.r, w[i__2].i = q__1.i;
	switch (modd) {
	    case 1:  goto L117;
	    case 2:  goto L116;
	}
L116:
	i__2 = *m;
	i__3 = *m + j * y_dim1;
	q__1.r = y[i__3].r * 2.f, q__1.i = y[i__3].i * 2.f;
	w[i__2].r = q__1.r, w[i__2].i = q__1.i;
L117:
	i__2 = *m;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    i__3 = i__ + j * y_dim1;
	    i__4 = i__;
	    y[i__3].r = w[i__4].r, y[i__3].i = w[i__4].i;
/* L118: */
	}
/* L119: */
    }
    k = iwbc + mhm1 - 1;
    i__ = iwba + mhm1;
    i__1 = k;
    w[i__1].r = 0.f, w[i__1].i = 0.f;
    i__1 = i__;
    w[i__1].r = 0.f, w[i__1].i = 0.f;
    i__1 = k + 1;
    i__2 = k + 1;
    q__1.r = w[i__2].r * 2.f, q__1.i = w[i__2].i * 2.f;
    w[i__1].r = q__1.r, w[i__1].i = q__1.i;
    switch (modd) {
	case 1:  goto L120;
	case 2:  goto L121;
    }
L120:
    k = iwbb + mhm1 - 1;
    i__1 = k;
    i__2 = k;
    i__3 = i__ - 1;
    q__1.r = w[i__2].r - w[i__3].r, q__1.i = w[i__2].i - w[i__3].i;
    w[i__1].r = q__1.r, w[i__1].i = q__1.i;
    i__1 = iwbc - 1;
    i__2 = iwbc - 1;
    i__3 = iwbb - 1;
    q__1.r = w[i__2].r + w[i__3].r, q__1.i = w[i__2].i + w[i__3].i;
    w[i__1].r = q__1.r, w[i__1].i = q__1.i;
    goto L122;
L121:
    i__1 = iwbb - 1;
    i__2 = k + 1;
    w[i__1].r = w[i__2].r, w[i__1].i = w[i__2].i;
L122:
    goto L107;

/*     REVERSE COLUMNS WHEN NPEROD = 4 */

L123:
    irev = 1;
    nby2 = *n / 2;
L124:
    i__1 = nby2;
    for (j = 1; j <= i__1; ++j) {
	mskip = *n + 1 - j;
	i__2 = *m;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    i__3 = i__ + j * y_dim1;
	    a1.r = y[i__3].r, a1.i = y[i__3].i;
	    i__3 = i__ + j * y_dim1;
	    i__4 = i__ + mskip * y_dim1;
	    y[i__3].r = y[i__4].r, y[i__3].i = y[i__4].i;
	    i__3 = i__ + mskip * y_dim1;
	    y[i__3].r = a1.r, y[i__3].i = a1.i;
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
	    i__3 = mhmi;
	    i__4 = mhpi + j * y_dim1;
	    i__5 = i__ + j * y_dim1;
	    q__2.r = y[i__4].r + y[i__5].r, q__2.i = y[i__4].i + y[i__5].i;
	    q__1.r = q__2.r * .5f, q__1.i = q__2.i * .5f;
	    w[i__3].r = q__1.r, w[i__3].i = q__1.i;
	    i__3 = mhpi;
	    i__4 = mhpi + j * y_dim1;
	    i__5 = i__ + j * y_dim1;
	    q__2.r = y[i__4].r - y[i__5].r, q__2.i = y[i__4].i - y[i__5].i;
	    q__1.r = q__2.r * .5f, q__1.i = q__2.i * .5f;
	    w[i__3].r = q__1.r, w[i__3].i = q__1.i;
/* L128: */
	}
	i__2 = mh;
	i__3 = mh + j * y_dim1;
	q__1.r = y[i__3].r * .5f, q__1.i = y[i__3].i * .5f;
	w[i__2].r = q__1.r, w[i__2].i = q__1.i;
	switch (modd) {
	    case 1:  goto L130;
	    case 2:  goto L129;
	}
L129:
	i__2 = *m;
	i__3 = *m + j * y_dim1;
	q__1.r = y[i__3].r * .5f, q__1.i = y[i__3].i * .5f;
	w[i__2].r = q__1.r, w[i__2].i = q__1.i;
L130:
	i__2 = *m;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    i__3 = i__ + j * y_dim1;
	    i__4 = i__;
	    y[i__3].r = w[i__4].r, y[i__3].i = w[i__4].i;
/* L131: */
	}
/* L132: */
    }
L133:

/*     RETURN STORAGE REQUIREMENTS FOR W ARRAY. */

    r__1 = (real) (ipstor + iwp - 1);
    q__1.r = r__1, q__1.i = 0.f;
    w[1].r = q__1.r, w[1].i = q__1.i;
    return 0;
} /* cmgnbn_ */

