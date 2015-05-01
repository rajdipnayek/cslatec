/* dmpar.f -- translated by f2c (version 12.02.01).
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

/* DECK DMPAR */
/* Subroutine */ int dmpar_(integer *n, doublereal *r__, integer *ldr, 
	integer *ipvt, doublereal *diag, doublereal *qtb, doublereal *delta, 
	doublereal *par, doublereal *x, doublereal *sigma, doublereal *wa1, 
	doublereal *wa2)
{
    /* Initialized data */

    static doublereal p1 = .1;
    static doublereal p001 = .001;
    static doublereal zero = 0.;

    /* System generated locals */
    integer r_dim1, r_offset, i__1, i__2;
    doublereal d__1, d__2;

    /* Local variables */
    static integer i__, j, k, l;
    static doublereal fp;
    static integer jm1, jp1;
    static doublereal sum, parc, parl;
    static integer iter;
    static doublereal temp, paru, dwarf;
    static integer nsing;
    static doublereal gnorm;
    extern doublereal d1mach_(integer *), denorm_(integer *, doublereal *);
    static doublereal dxnorm;
    extern /* Subroutine */ int dqrslv_(integer *, doublereal *, integer *, 
	    integer *, doublereal *, doublereal *, doublereal *, doublereal *,
	     doublereal *);

/* ***BEGIN PROLOGUE  DMPAR */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to DNLS1 and DNLS1E */
/* ***LIBRARY   SLATEC */
/* ***TYPE      DOUBLE PRECISION (LMPAR-S, DMPAR-D) */
/* ***AUTHOR  (UNKNOWN) */
/* ***DESCRIPTION */

/*   **** Double Precision version of LMPAR **** */

/*     Given an M by N matrix A, an N by N nonsingular DIAGONAL */
/*     matrix D, an M-vector B, and a positive number DELTA, */
/*     the problem is to determine a value for the parameter */
/*     PAR such that if X solves the system */

/*           A*X = B ,     SQRT(PAR)*D*X = 0 , */

/*     in the least squares sense, and DXNORM is the Euclidean */
/*     norm of D*X, then either PAR is zero and */

/*           (DXNORM-DELTA) .LE. 0.1*DELTA , */

/*     or PAR is positive and */

/*           ABS(DXNORM-DELTA) .LE. 0.1*DELTA . */

/*     This subroutine completes the solution of the problem */
/*     if it is provided with the necessary information from the */
/*     QR factorization, with column pivoting, of A. That is, if */
/*     A*P = Q*R, where P is a permutation matrix, Q has orthogonal */
/*     columns, and R is an upper triangular matrix with diagonal */
/*     elements of nonincreasing magnitude, then DMPAR expects */
/*     the full upper triangle of R, the permutation matrix P, */
/*     and the first N components of (Q TRANSPOSE)*B. On output */
/*     DMPAR also provides an upper triangular matrix S such that */

/*            T   T                   T */
/*           P *(A *A + PAR*D*D)*P = S *S . */

/*     S is employed within DMPAR and may be of separate interest. */

/*     Only a few iterations are generally needed for convergence */
/*     of the algorithm. If, however, the limit of 10 iterations */
/*     is reached, then the output PAR will contain the best */
/*     value obtained so far. */

/*     The subroutine statement is */

/*       SUBROUTINE DMPAR(N,R,LDR,IPVT,DIAG,QTB,DELTA,PAR,X,SIGMA, */
/*                        WA1,WA2) */

/*     where */

/*       N is a positive integer input variable set to the order of R. */

/*       R is an N by N array. On input the full upper triangle */
/*         must contain the full upper triangle of the matrix R. */
/*         On output the full upper triangle is unaltered, and the */
/*         strict lower triangle contains the strict upper triangle */
/*         (transposed) of the upper triangular matrix S. */

/*       LDR is a positive integer input variable not less than N */
/*         which specifies the leading dimension of the array R. */

/*       IPVT is an integer input array of length N which defines the */
/*         permutation matrix P such that A*P = Q*R. Column J of P */
/*         is column IPVT(J) of the identity matrix. */

/*       DIAG is an input array of length N which must contain the */
/*         diagonal elements of the matrix D. */

/*       QTB is an input array of length N which must contain the first */
/*         N elements of the vector (Q TRANSPOSE)*B. */

/*       DELTA is a positive input variable which specifies an upper */
/*         bound on the Euclidean norm of D*X. */

/*       PAR is a nonnegative variable. On input PAR contains an */
/*         initial estimate of the Levenberg-Marquardt parameter. */
/*         On output PAR contains the final estimate. */

/*       X is an output array of length N which contains the least */
/*         squares solution of the system A*X = B, SQRT(PAR)*D*X = 0, */
/*         for the output PAR. */

/*       SIGMA is an output array of length N which contains the */
/*         diagonal elements of the upper triangular matrix S. */

/*       WA1 and WA2 are work arrays of length N. */

/* ***SEE ALSO  DNLS1, DNLS1E */
/* ***ROUTINES CALLED  D1MACH, DENORM, DQRSLV */
/* ***REVISION HISTORY  (YYMMDD) */
/*   800301  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890831  Modified array declarations.  (WRB) */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900326  Removed duplicate information from DESCRIPTION section. */
/*           (WRB) */
/*   900328  Added TYPE section.  (WRB) */
/* ***END PROLOGUE  DMPAR */
    /* Parameter adjustments */
    r_dim1 = *ldr;
    r_offset = 1 + r_dim1;
    r__ -= r_offset;
    --ipvt;
    --diag;
    --qtb;
    --x;
    --sigma;
    --wa1;
    --wa2;

    /* Function Body */
/* ***FIRST EXECUTABLE STATEMENT  DMPAR */
    dwarf = d1mach_(&c__1);

/*     COMPUTE AND STORE IN X THE GAUSS-NEWTON DIRECTION. IF THE */
/*     JACOBIAN IS RANK-DEFICIENT, OBTAIN A LEAST SQUARES SOLUTION. */

    nsing = *n;
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	wa1[j] = qtb[j];
	if (r__[j + j * r_dim1] == zero && nsing == *n) {
	    nsing = j - 1;
	}
	if (nsing < *n) {
	    wa1[j] = zero;
	}
/* L10: */
    }
    if (nsing < 1) {
	goto L50;
    }
    i__1 = nsing;
    for (k = 1; k <= i__1; ++k) {
	j = nsing - k + 1;
	wa1[j] /= r__[j + j * r_dim1];
	temp = wa1[j];
	jm1 = j - 1;
	if (jm1 < 1) {
	    goto L30;
	}
	i__2 = jm1;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    wa1[i__] -= r__[i__ + j * r_dim1] * temp;
/* L20: */
	}
L30:
/* L40: */
	;
    }
L50:
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	l = ipvt[j];
	x[l] = wa1[j];
/* L60: */
    }

/*     INITIALIZE THE ITERATION COUNTER. */
/*     EVALUATE THE FUNCTION AT THE ORIGIN, AND TEST */
/*     FOR ACCEPTANCE OF THE GAUSS-NEWTON DIRECTION. */

    iter = 0;
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	wa2[j] = diag[j] * x[j];
/* L70: */
    }
    dxnorm = denorm_(n, &wa2[1]);
    fp = dxnorm - *delta;
    if (fp <= p1 * *delta) {
	goto L220;
    }

/*     IF THE JACOBIAN IS NOT RANK DEFICIENT, THE NEWTON */
/*     STEP PROVIDES A LOWER BOUND, PARL, FOR THE ZERO OF */
/*     THE FUNCTION. OTHERWISE SET THIS BOUND TO ZERO. */

    parl = zero;
    if (nsing < *n) {
	goto L120;
    }
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	l = ipvt[j];
	wa1[j] = diag[l] * (wa2[l] / dxnorm);
/* L80: */
    }
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	sum = zero;
	jm1 = j - 1;
	if (jm1 < 1) {
	    goto L100;
	}
	i__2 = jm1;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    sum += r__[i__ + j * r_dim1] * wa1[i__];
/* L90: */
	}
L100:
	wa1[j] = (wa1[j] - sum) / r__[j + j * r_dim1];
/* L110: */
    }
    temp = denorm_(n, &wa1[1]);
    parl = fp / *delta / temp / temp;
L120:

/*     CALCULATE AN UPPER BOUND, PARU, FOR THE ZERO OF THE FUNCTION. */

    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	sum = zero;
	i__2 = j;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    sum += r__[i__ + j * r_dim1] * qtb[i__];
/* L130: */
	}
	l = ipvt[j];
	wa1[j] = sum / diag[l];
/* L140: */
    }
    gnorm = denorm_(n, &wa1[1]);
    paru = gnorm / *delta;
    if (paru == zero) {
	paru = dwarf / min(*delta,p1);
    }

/*     IF THE INPUT PAR LIES OUTSIDE OF THE INTERVAL (PARL,PARU), */
/*     SET PAR TO THE CLOSER ENDPOINT. */

    *par = max(*par,parl);
    *par = min(*par,paru);
    if (*par == zero) {
	*par = gnorm / dxnorm;
    }

/*     BEGINNING OF AN ITERATION. */

L150:
    ++iter;

/*        EVALUATE THE FUNCTION AT THE CURRENT VALUE OF PAR. */

    if (*par == zero) {
/* Computing MAX */
	d__1 = dwarf, d__2 = p001 * paru;
	*par = max(d__1,d__2);
    }
    temp = sqrt(*par);
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	wa1[j] = temp * diag[j];
/* L160: */
    }
    dqrslv_(n, &r__[r_offset], ldr, &ipvt[1], &wa1[1], &qtb[1], &x[1], &sigma[
	    1], &wa2[1]);
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	wa2[j] = diag[j] * x[j];
/* L170: */
    }
    dxnorm = denorm_(n, &wa2[1]);
    temp = fp;
    fp = dxnorm - *delta;

/*        IF THE FUNCTION IS SMALL ENOUGH, ACCEPT THE CURRENT VALUE */
/*        OF PAR. ALSO TEST FOR THE EXCEPTIONAL CASES WHERE PARL */
/*        IS ZERO OR THE NUMBER OF ITERATIONS HAS REACHED 10. */

    if (abs(fp) <= p1 * *delta || parl == zero && fp <= temp && temp < zero ||
	     iter == 10) {
	goto L220;
    }

/*        COMPUTE THE NEWTON CORRECTION. */

    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	l = ipvt[j];
	wa1[j] = diag[l] * (wa2[l] / dxnorm);
/* L180: */
    }
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	wa1[j] /= sigma[j];
	temp = wa1[j];
	jp1 = j + 1;
	if (*n < jp1) {
	    goto L200;
	}
	i__2 = *n;
	for (i__ = jp1; i__ <= i__2; ++i__) {
	    wa1[i__] -= r__[i__ + j * r_dim1] * temp;
/* L190: */
	}
L200:
/* L210: */
	;
    }
    temp = denorm_(n, &wa1[1]);
    parc = fp / *delta / temp / temp;

/*        DEPENDING ON THE SIGN OF THE FUNCTION, UPDATE PARL OR PARU. */

    if (fp > zero) {
	parl = max(parl,*par);
    }
    if (fp < zero) {
	paru = min(paru,*par);
    }

/*        COMPUTE AN IMPROVED ESTIMATE FOR PAR. */

/* Computing MAX */
    d__1 = parl, d__2 = *par + parc;
    *par = max(d__1,d__2);

/*        END OF AN ITERATION. */

    goto L150;
L220:

/*     TERMINATION. */

    if (iter == 0) {
	*par = zero;
    }
    return 0;

/*     LAST CARD OF SUBROUTINE DMPAR. */

} /* dmpar_ */

