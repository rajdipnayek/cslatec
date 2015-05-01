/* scov.f -- translated by f2c (version 12.02.01).
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

static logical c_false = FALSE_;
static integer c__1 = 1;
static integer c__2 = 2;

/* DECK SCOV */
/* Subroutine */ int scov_(S_fp fcn, integer *iopt, integer *m, integer *n, 
	real *x, real *fvec, real *r__, integer *ldr, integer *info, real *
	wa1, real *wa2, real *wa3, real *wa4)
{
    /* Initialized data */

    static real zero = 0.f;
    static real one = 1.f;

    /* System generated locals */
    integer r_dim1, r_offset, i__1, i__2, i__3;

    /* Local variables */
    static integer i__, j, k, kp1, nm1, idum;
    static logical sing;
    static real temp;
    static integer nrow, iflag;
    extern /* Subroutine */ int qrfac_(integer *, integer *, real *, integer *
	    , logical *, integer *, integer *, real *, real *, real *);
    static real sigma;
    extern doublereal enorm_(integer *, real *);
    extern /* Subroutine */ int fdjac3_(S_fp, integer *, integer *, real *, 
	    real *, real *, integer *, integer *, real *, real *), xermsg_(
	    char *, char *, char *, integer *, integer *, ftnlen, ftnlen, 
	    ftnlen), rwupdt_(integer *, real *, integer *, real *, real *, 
	    real *, real *, real *);

/* ***BEGIN PROLOGUE  SCOV */
/* ***PURPOSE  Calculate the covariance matrix for a nonlinear data */
/*            fitting problem.  It is intended to be used after a */
/*            successful return from either SNLS1 or SNLS1E. */
/* ***LIBRARY   SLATEC */
/* ***CATEGORY  K1B1 */
/* ***TYPE      SINGLE PRECISION (SCOV-S, DCOV-D) */
/* ***KEYWORDS  COVARIANCE MATRIX, NONLINEAR DATA FITTING, */
/*             NONLINEAR LEAST SQUARES */
/* ***AUTHOR  Hiebert, K. L., (SNLA) */
/* ***DESCRIPTION */

/*  1. Purpose. */

/*     SCOV calculates the covariance matrix for a nonlinear data */
/*     fitting problem.  It is intended to be used after a */
/*     successful return from either SNLS1 or SNLS1E. SCOV */
/*     and SNLS1 (and SNLS1E) have compatible parameters.  The */
/*     required external subroutine, FCN, is the same */
/*     for all three codes, SCOV, SNLS1, and SNLS1E. */

/*  2. Subroutine and Type Statements. */

/*     SUBROUTINE SCOV(FCN,IOPT,M,N,X,FVEC,R,LDR,INFO, */
/*                     WA1,WA2,WA3,WA4) */
/*     INTEGER IOPT,M,N,LDR,INFO */
/*     REAL X(N),FVEC(M),R(LDR,N),WA1(N),WA2(N),WA3(N),WA4(M) */
/*     EXTERNAL FCN */

/*  3. Parameters. */

/*       FCN is the name of the user-supplied subroutine which calculates */
/*         the functions.  If the user wants to supply the Jacobian */
/*         (IOPT=2 or 3), then FCN must be written to calculate the */
/*         Jacobian, as well as the functions.  See the explanation */
/*         of the IOPT argument below.  FCN must be declared in an */
/*         EXTERNAL statement in the calling program and should be */
/*         written as follows. */

/*         SUBROUTINE FCN(IFLAG,M,N,X,FVEC,FJAC,LDFJAC) */
/*         INTEGER IFLAG,LDFJAC,M,N */
/*         REAL X(N),FVEC(M) */
/*         ---------- */
/*         FJAC and LDFJAC may be ignored     , if IOPT=1. */
/*         REAL FJAC(LDFJAC,N)                , if IOPT=2. */
/*         REAL FJAC(N)                       , if IOPT=3. */
/*         ---------- */
/*           IFLAG will never be zero when FCN is called by SCOV. */
/*         RETURN */
/*         ---------- */
/*           If IFLAG=1, calculate the functions at X and return */
/*           this vector in FVEC. */
/*         RETURN */
/*         ---------- */
/*           If IFLAG=2, calculate the full Jacobian at X and return */
/*           this matrix in FJAC.  Note that IFLAG will never be 2 unless */
/*           IOPT=2.  FVEC contains the function values at X and must */
/*           not be altered.  FJAC(I,J) must be set to the derivative */
/*           of FVEC(I) with respect to X(J). */
/*         RETURN */
/*         ---------- */
/*           If IFLAG=3, calculate the LDFJAC-th row of the Jacobian */
/*           and return this vector in FJAC.  Note that IFLAG will */
/*           never be 3 unless IOPT=3.  FJAC(J) must be set to */
/*           the derivative of FVEC(LDFJAC) with respect to X(J). */
/*         RETURN */
/*         ---------- */
/*         END */


/*         The value of IFLAG should not be changed by FCN unless the */
/*         user wants to terminate execution of SCOV.  In this case, set */
/*         IFLAG to a negative integer. */


/*    IOPT is an input variable which specifies how the Jacobian will */
/*         be calculated.  If IOPT=2 or 3, then the user must supply the */
/*         Jacobian, as well as the function values, through the */
/*         subroutine FCN.  If IOPT=2, the user supplies the full */
/*         Jacobian with one call to FCN.  If IOPT=3, the user supplies */
/*         one row of the Jacobian with each call.  (In this manner, */
/*         storage can be saved because the full Jacobian is not stored.) */
/*         If IOPT=1, the code will approximate the Jacobian by forward */
/*         differencing. */

/*       M is a positive integer input variable set to the number of */
/*         functions. */

/*       N is a positive integer input variable set to the number of */
/*         variables.  N must not exceed M. */

/*       X is an array of length N.  On input X must contain the value */
/*         at which the covariance matrix is to be evaluated.  This is */
/*         usually the value for X returned from a successful run of */
/*         SNLS1 (or SNLS1E).  The value of X will not be changed. */

/*    FVEC is an output array of length M which contains the functions */
/*         evaluated at X. */

/*       R is an output array.  For IOPT=1 and 2, R is an M by N array. */
/*         For IOPT=3, R is an N by N array.  On output, if INFO=1, */
/*         the upper N by N submatrix of R contains the covariance */
/*         matrix evaluated at X. */

/*     LDR is a positive integer input variable which specifies */
/*         the leading dimension of the array R.  For IOPT=1 and 2, */
/*         LDR must not be less than M.  For IOPT=3, LDR must not */
/*         be less than N. */

/*    INFO is an integer output variable.  If the user has terminated */
/*         execution, INFO is set to the (negative) value of IFLAG.  See */
/*         description of FCN. Otherwise, INFO is set as follows. */

/*         INFO = 0 Improper input parameters (M.LE.0 or N.LE.0). */

/*         INFO = 1 Successful return.  The covariance matrix has been */
/*                  calculated and stored in the upper N by N */
/*                  submatrix of R. */

/*         INFO = 2 The Jacobian matrix is singular for the input value */
/*                  of X.  The covariance matrix cannot be calculated. */
/*                  The upper N by N submatrix of R contains the QR */
/*                  factorization of the Jacobian (probably not of */
/*                  interest to the user). */

/*     WA1 is a work array of length N. */
/*     WA2 is a work array of length N. */
/*     WA3 is a work array of length N. */
/*     WA4 is a work array of length M. */

/* ***REFERENCES  (NONE) */
/* ***ROUTINES CALLED  ENORM, FDJAC3, QRFAC, RWUPDT, XERMSG */
/* ***REVISION HISTORY  (YYMMDD) */
/*   810522  DATE WRITTEN */
/*   890505  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ) */
/*   900510  Fixed an error message.  (RWC) */
/* ***END PROLOGUE  SCOV */

/*     REVISED 820707-1100 */
/*     REVISED YYMMDD HHMM */

    /* Parameter adjustments */
    --x;
    --fvec;
    r_dim1 = *ldr;
    r_offset = 1 + r_dim1;
    r__ -= r_offset;
    --wa1;
    --wa2;
    --wa3;
    --wa4;

    /* Function Body */
/* ***FIRST EXECUTABLE STATEMENT  SCOV */
    sing = FALSE_;
    iflag = 0;
    if (*m <= 0 || *n <= 0) {
	goto L300;
    }

/*     CALCULATE SIGMA = (SUM OF THE SQUARED RESIDUALS) / (M-N) */
    iflag = 1;
    (*fcn)(&iflag, m, n, &x[1], &fvec[1], &r__[r_offset], ldr);
    if (iflag < 0) {
	goto L300;
    }
    temp = enorm_(m, &fvec[1]);
    sigma = one;
    if (*m != *n) {
	sigma = temp * temp / (*m - *n);
    }

/*     CALCULATE THE JACOBIAN */
    if (*iopt == 3) {
	goto L200;
    }

/*     STORE THE FULL JACOBIAN USING M*N STORAGE */
    if (*iopt == 1) {
	goto L100;
    }

/*     USER SUPPLIES THE JACOBIAN */
    iflag = 2;
    (*fcn)(&iflag, m, n, &x[1], &fvec[1], &r__[r_offset], ldr);
    goto L110;

/*     CODE APPROXIMATES THE JACOBIAN */
L100:
    fdjac3_((S_fp)fcn, m, n, &x[1], &fvec[1], &r__[r_offset], ldr, &iflag, &
	    zero, &wa4[1]);
L110:
    if (iflag < 0) {
	goto L300;
    }

/*     COMPUTE THE QR DECOMPOSITION */
    qrfac_(m, n, &r__[r_offset], ldr, &c_false, &idum, &c__1, &wa1[1], &wa1[1]
	    , &wa1[1]);
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L120: */
	r__[i__ + i__ * r_dim1] = wa1[i__];
    }
    goto L225;

/*     COMPUTE THE QR FACTORIZATION OF THE JACOBIAN MATRIX CALCULATED ONE */
/*     ROW AT A TIME AND STORED IN THE UPPER TRIANGLE OF R. */
/*     ( (Q TRANSPOSE)*FVEC IS ALSO CALCULATED BUT NOT USED.) */
L200:
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	wa2[j] = zero;
	i__2 = *n;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    r__[i__ + j * r_dim1] = zero;
/* L205: */
	}
/* L210: */
    }
    iflag = 3;
    i__1 = *m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	nrow = i__;
	(*fcn)(&iflag, m, n, &x[1], &fvec[1], &wa1[1], &nrow);
	if (iflag < 0) {
	    goto L300;
	}
	temp = fvec[i__];
	rwupdt_(n, &r__[r_offset], ldr, &wa1[1], &wa2[1], &temp, &wa3[1], &
		wa4[1]);
/* L220: */
    }

/*     CHECK IF R IS SINGULAR. */
L225:
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (r__[i__ + i__ * r_dim1] == zero) {
	    sing = TRUE_;
	}
/* L230: */
    }
    if (sing) {
	goto L300;
    }

/*     R IS UPPER TRIANGULAR.  CALCULATE (R TRANSPOSE) INVERSE AND STORE */
/*     IN THE UPPER TRIANGLE OF R. */
    if (*n == 1) {
	goto L275;
    }
    nm1 = *n - 1;
    i__1 = nm1;
    for (k = 1; k <= i__1; ++k) {

/*     INITIALIZE THE RIGHT-HAND SIDE (WA1(*)) AS THE K-TH COLUMN OF THE */
/*     IDENTITY MATRIX. */
	i__2 = *n;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    wa1[i__] = zero;
/* L240: */
	}
	wa1[k] = one;

	r__[k + k * r_dim1] = wa1[k] / r__[k + k * r_dim1];
	kp1 = k + 1;
	i__2 = *n;
	for (i__ = kp1; i__ <= i__2; ++i__) {

/*     SUBTRACT R(K,I-1)*R(I-1,*) FROM THE RIGHT-HAND SIDE, WA1(*). */
	    i__3 = *n;
	    for (j = i__; j <= i__3; ++j) {
		wa1[j] -= r__[k + (i__ - 1) * r_dim1] * r__[i__ - 1 + j * 
			r_dim1];
/* L250: */
	    }
	    r__[k + i__ * r_dim1] = wa1[i__] / r__[i__ + i__ * r_dim1];
/* L260: */
	}
/* L270: */
    }
L275:
    r__[*n + *n * r_dim1] = one / r__[*n + *n * r_dim1];

/*     CALCULATE R-INVERSE * (R TRANSPOSE) INVERSE AND STORE IN THE UPPER */
/*     TRIANGLE OF R. */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = *n;
	for (j = i__; j <= i__2; ++j) {
	    temp = zero;
	    i__3 = *n;
	    for (k = j; k <= i__3; ++k) {
		temp += r__[i__ + k * r_dim1] * r__[j + k * r_dim1];
/* L280: */
	    }
	    r__[i__ + j * r_dim1] = temp * sigma;
/* L290: */
	}
    }
    *info = 1;

L300:
    if (*m <= 0 || *n <= 0) {
	*info = 0;
    }
    if (iflag < 0) {
	*info = iflag;
    }
    if (sing) {
	*info = 2;
    }
    if (*info < 0) {
	xermsg_("SLATEC", "SCOV", "EXECUTION TERMINATED BECAUSE USER SET IFL"
		"AG NEGATIVE.", &c__1, &c__1, (ftnlen)6, (ftnlen)4, (ftnlen)53)
		;
    }
    if (*info == 0) {
	xermsg_("SLATEC", "SCOV", "INVALID INPUT PARAMETER.", &c__2, &c__1, (
		ftnlen)6, (ftnlen)4, (ftnlen)24);
    }
    if (*info == 2) {
	xermsg_("SLATEC", "SCOV", "SINGULAR JACOBIAN MATRIX, COVARIANCE MATR"
		"IX CANNOT BE CALCULATED.", &c__1, &c__1, (ftnlen)6, (ftnlen)4,
		 (ftnlen)65);
    }
    return 0;
} /* scov_ */

