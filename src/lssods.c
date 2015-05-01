/* lssods.f -- translated by f2c (version 12.02.01).
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

static integer c__3 = 3;
static integer c__2 = 2;
static integer c__1 = 1;
static integer c__4 = 4;
static integer c__0 = 0;
static logical c_false = FALSE_;
static integer c__8 = 8;

/* DECK LSSODS */
/* Subroutine */ int lssods_(real *a, real *x, real *b, integer *m, integer *
	n, integer *nrda, integer *iflag, integer *irank, integer *iscale, 
	real *q, real *diag, integer *kpivot, integer *iter, real *resnrm, 
	real *xnorm, real *z__, real *r__, real *div, real *td, real *scales)
{
    /* System generated locals */
    integer a_dim1, a_offset, q_dim1, q_offset, i__1, i__2, i__3;
    real r__1;

    /* Local variables */
    static integer j, k, l, mj, kp, it;
    static real acc, gam;
    static integer irm, irp;
    static real uro;
    static integer nfat, mmir, nmir;
    extern doublereal sdot_(integer *, real *, integer *, real *, integer *);
    static real znrm0, gamma;
    extern /* Subroutine */ int xgetf_(integer *);
    static integer iterp;
    extern /* Subroutine */ int xsetf_(integer *);
    static real znorm;
    extern doublereal r1mach_(integer *);
    extern integer j4save_(integer *, integer *, logical *);
    static integer nfatal, maxmes;
    extern doublereal sdsdot_(integer *, real *, real *, integer *, real *, 
	    integer *);
    extern /* Subroutine */ int xermax_(integer *), xermsg_(char *, char *, 
	    char *, integer *, integer *, ftnlen, ftnlen, ftnlen), orthol_(
	    real *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, real *, integer *, real *, real *, real *), ohtror_(
	    real *, integer *, integer *, real *, integer *, real *, real *);

/* ***BEGIN PROLOGUE  LSSODS */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to BVSUP */
/* ***LIBRARY   SLATEC */
/* ***TYPE      SINGLE PRECISION (LSSODS-S) */
/* ***AUTHOR  (UNKNOWN) */
/* ***DESCRIPTION */

/*     LSSODS solves the same problem as SODS (in fact, it is called by */
/*     SODS) but is somewhat more flexible in its use. In particular, */
/*     LSSODS allows for iterative refinement of the solution, makes the */
/*     transformation and triangular reduction information more */
/*     accessible, and enables the user to avoid destruction of the */
/*     original matrix A. */

/*     Modeled after the ALGOL codes in the articles in the REFERENCES */
/*     section. */

/* ********************************************************************** */
/*   INPUT */
/* ********************************************************************** */

/*     A -- Contains the matrix of M equations in N unknowns and must */
/*          be dimensioned NRDA by N. A remains unchanged */
/*     X -- Solution array of length at least N */
/*     B -- Given constant vector of length M, B remains unchanged */
/*     M -- Number of equations, M greater or equal to 1 */
/*     N -- Number of unknowns, N not larger than M */
/*  NRDA -- Row dimension of A, NRDA greater or equal to M */
/* IFLAG -- Status indicator */
/*         = 0 for the first call (and for each new problem defined by */
/*             a new matrix A) when the matrix data is treated as exact */
/*         =-K for the first call (and for each new problem defined by */
/*             a new matrix A) when the matrix data is assumed to be */
/*             accurate to about K digits */
/*         = 1 for subsequent calls whenever the matrix A has already */
/*             been decomposed (problems with new vectors B but */
/*             same matrix a can be handled efficiently) */
/* ISCALE -- Scaling indicator */
/*         =-1 if the matrix A is to be pre-scaled by */
/*             columns when appropriate */
/*             If the scaling indicator is not equal to -1 */
/*             no scaling will be attempted */
/*             For most problems scaling will probably not be necessary */
/*   ITER -- Maximum number of iterative improvement steps to be */
/*           performed,  0 .LE. ITER .LE. 10   (SODS uses ITER=0) */
/*      Q -- Matrix used for the transformation, must be dimensioned */
/*           NRDA by N  (SODS puts A in the Q location which conserves */
/*           storage but destroys A) */
/*           When iterative improvement of the solution is requested, */
/*           ITER .GT. 0, this additional storage for Q must be */
/*           made available */
/* DIAG,KPIVOT,Z,R, -- Arrays of length N (except for R which is M) */
/*   DIV,TD,SCALES     used for internal storage */

/* ********************************************************************** */
/*   OUTPUT */
/* ********************************************************************** */

/*  IFLAG -- Status indicator */
/*            =1 if solution was obtained */
/*            =2 if improper input is detected */
/*            =3 if rank of matrix is less than N */
/*               if the minimal length least squares solution is */
/*               desired, simply reset IFLAG=1 and call the code again */

/*       The next three IFLAG values can occur only when */
/*        the iterative improvement mode is being used. */
/*            =4 if the problem is ill-conditioned and maximal */
/*               machine accuracy is not achievable */
/*            =5 if the problem is very ill-conditioned and the solution */
/*               IS likely to have no correct digits */
/*            =6 if the allowable number of iterative improvement steps */
/*               has been completed without getting convergence */
/*      X -- Least squares solution of  A X = B */
/*  IRANK -- Contains the numerically determined matrix rank */
/*           the user must not alter this value on succeeding calls */
/*           with input values of IFLAG=1 */
/*      Q -- Contains the strictly upper triangular part of the reduced */
/*           matrix and the transformation information in the lower */
/*           triangular part */
/*   DIAG -- Contains the diagonal elements of the triangular reduced */
/*           matrix */
/* KPIVOT -- Contains the pivotal information.  The column interchanges */
/*           performed on the original matrix are recorded here */
/*   ITER -- The actual number of iterative corrections used */
/* RESNRM -- The Euclidean norm of the residual vector  B - A X */
/*  XNORM -- The Euclidean norm of the solution vector */
/* DIV,TD -- Contains transformation information for rank */
/*           deficient problems */
/* SCALES -- Contains the column scaling parameters */

/* ********************************************************************** */

/* ***SEE ALSO  BVSUP */
/* ***REFERENCES  G. Golub, Numerical methods for solving linear least */
/*                 squares problems, Numerische Mathematik 7, (1965), */
/*                 pp. 206-216. */
/*               P. Businger and G. Golub, Linear least squares */
/*                 solutions by Householder transformations, Numerische */
/*                 Mathematik  7, (1965), pp. 269-276. */
/* ***ROUTINES CALLED  J4SAVE, OHTROR, ORTHOL, R1MACH, SDOT, SDSDOT, */
/*                    XERMAX, XERMSG, XGETF, XSETF */
/* ***REVISION HISTORY  (YYMMDD) */
/*   750601  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890831  Modified array declarations.  (WRB) */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ) */
/*   900402  Added TYPE section.  (WRB) */
/*   910408  Updated the REFERENCES section.  (WRB) */
/*   920501  Reformatted the REFERENCES section.  (WRB) */
/* ***END PROLOGUE  LSSODS */

/* ********************************************************************** */

/*     MACHINE PRECISION (COMPUTER UNIT ROUNDOFF VALUE) IS DEFINED */
/*     THE FUNCTION R1MACH. */

/* ***FIRST EXECUTABLE STATEMENT  LSSODS */
    /* Parameter adjustments */
    --x;
    --b;
    q_dim1 = *nrda;
    q_offset = 1 + q_dim1;
    q -= q_offset;
    a_dim1 = *nrda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --diag;
    --kpivot;
    --z__;
    --r__;
    --div;
    --td;
    --scales;

    /* Function Body */
    uro = r1mach_(&c__3);

/* ********************************************************************** */

    if (*n < 1 || *m < *n || *nrda < *m) {
	goto L1;
    }
    if (*iter < 0) {
	goto L1;
    }
    if (*iflag <= 0) {
	goto L5;
    }
    if (*iflag == 1) {
	goto L15;
    }

/*     INVALID INPUT FOR LSSODS */
L1:
    *iflag = 2;
    xermsg_("SLATEC", "LSSODS", "INVALID INPUT PARAMETERS.", &c__2, &c__1, (
	    ftnlen)6, (ftnlen)6, (ftnlen)25);
    return 0;

L5:
    xgetf_(&nfatal);
    maxmes = j4save_(&c__4, &c__0, &c_false);
    if (*iflag == 0) {
	goto L7;
    }
    nfat = -1;
    if (nfatal == 0) {
	nfat = 0;
    }
    xsetf_(&nfat);
    xermax_(&c__1);

/*     COPY MATRIX A INTO MATRIX Q */

L7:
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	i__2 = *m;
	for (k = 1; k <= i__2; ++k) {
/* L10: */
	    q[k + j * q_dim1] = a[k + j * a_dim1];
	}
    }

/*     USE ORTHOGONAL TRANSFORMATIONS TO REDUCE Q TO */
/*     UPPER TRIANGULAR FORM */

    orthol_(&q[q_offset], m, n, nrda, iflag, irank, iscale, &diag[1], &kpivot[
	    1], &scales[1], &z__[1], &td[1]);

    xsetf_(&nfatal);
    xermax_(&maxmes);
    if (*irank == *n) {
	goto L12;
    }

/*     FOR RANK DEFICIENT PROBLEMS USE ADDITIONAL ORTHOGONAL */
/*     TRANSFORMATIONS TO FURTHER REDUCE Q */

    if (*irank != 0) {
	ohtror_(&q[q_offset], n, nrda, &diag[1], irank, &div[1], &td[1]);
    }
    return 0;

/*     STORE DIVISORS FOR THE TRIANGULAR SOLUTION */

L12:
    i__2 = *n;
    for (k = 1; k <= i__2; ++k) {
/* L13: */
	div[k] = diag[k];
    }

L15:
    irm = *irank - 1;
    irp = *irank + 1;
/* Computing MIN */
    i__2 = *iter + 1;
    iterp = min(i__2,11);
    acc = uro * 10.f;

/*     ZERO OUT SOLUTION ARRAY */

    i__2 = *n;
    for (k = 1; k <= i__2; ++k) {
/* L20: */
	x[k] = 0.f;
    }

    if (*irank > 0) {
	goto L25;
    }

/*     SPECIAL CASE FOR THE NULL MATRIX */
    *iter = 0;
    *xnorm = 0.f;
    *resnrm = sqrt(sdot_(m, &b[1], &c__1, &b[1], &c__1));
    return 0;

/*     COPY CONSTANT VECTOR INTO R */

L25:
    i__2 = *m;
    for (k = 1; k <= i__2; ++k) {
/* L30: */
	r__[k] = b[k];
    }

/* ********************************************************************** */
/*     SOLUTION SECTION */
/*     ITERATIVE REFINEMENT OF THE RESIDUAL VECTOR */
/* ********************************************************************** */

    i__2 = iterp;
    for (it = 1; it <= i__2; ++it) {
	*iter = it - 1;

/*        APPLY ORTHOGONAL TRANSFORMATION TO R */

	i__1 = *irank;
	for (j = 1; j <= i__1; ++j) {
	    mj = *m - j + 1;
	    gamma = sdot_(&mj, &q[j + j * q_dim1], &c__1, &r__[j], &c__1) / (
		    diag[j] * q[j + j * q_dim1]);
	    i__3 = *m;
	    for (k = j; k <= i__3; ++k) {
/* L35: */
		r__[k] += gamma * q[k + j * q_dim1];
	    }
	}

/*        BACKWARD SUBSTITUTION FOR TRIANGULAR SYSTEM SOLUTION */

	z__[*irank] = r__[*irank] / div[*irank];
	if (irm == 0) {
	    goto L45;
	}
	i__3 = irm;
	for (l = 1; l <= i__3; ++l) {
	    k = *irank - l;
	    kp = k + 1;
/* L40: */
	    z__[k] = (r__[k] - sdot_(&l, &q[k + kp * q_dim1], nrda, &z__[kp], 
		    &c__1)) / div[k];
	}

L45:
	if (*irank == *n) {
	    goto L60;
	}

/*        FOR RANK DEFICIENT PROBLEMS OBTAIN THE */
/*        MINIMAL LENGTH SOLUTION */

	nmir = *n - *irank;
	i__3 = *n;
	for (k = irp; k <= i__3; ++k) {
/* L50: */
	    z__[k] = 0.f;
	}
	i__3 = *irank;
	for (k = 1; k <= i__3; ++k) {
	    gam = (td[k] * z__[k] + sdot_(&nmir, &q[k + irp * q_dim1], nrda, &
		    z__[irp], &c__1)) / (td[k] * div[k]);
	    z__[k] += gam * td[k];
	    i__1 = *n;
	    for (j = irp; j <= i__1; ++j) {
/* L55: */
		z__[j] += gam * q[k + j * q_dim1];
	    }
	}

/*        REORDER SOLUTION COMPONENTS ACCORDING TO PIVOTAL POINTS */
/*        AND RESCALE ANSWERS AS DICTATED */

L60:
	i__1 = *n;
	for (k = 1; k <= i__1; ++k) {
	    z__[k] *= scales[k];
	    l = kpivot[k];
/* L65: */
	    x[l] += z__[k];
	}

/*        COMPUTE CORRECTION VECTOR NORM (SOLUTION NORM) */

	znorm = sqrt(sdot_(n, &z__[1], &c__1, &z__[1], &c__1));
	if (it == 1) {
	    *xnorm = znorm;
	}
	if (iterp > 1) {
	    goto L80;
	}

/*        NO ITERATIVE CORRECTIONS TO BE PERFORMED, SO COMPUTE */
/*        THE APPROXIMATE RESIDUAL NORM DEFINED BY THE EQUATIONS */
/*        WHICH ARE NOT SATISFIED BY THE SOLUTION */
/*        THEN WE ARE DONE */

	mmir = *m - *irank;
	if (mmir == 0) {
	    goto L70;
	}
	*resnrm = sqrt(sdot_(&mmir, &r__[irp], &c__1, &r__[irp], &c__1));
	return 0;
L70:
	*resnrm = 0.f;
	return 0;

/*        COMPUTE RESIDUAL VECTOR FOR THE ITERATIVE IMPROVEMENT PROCESS */

L80:
	i__1 = *m;
	for (k = 1; k <= i__1; ++k) {
/* L85: */
	    r__1 = -b[k];
	    r__[k] = -sdsdot_(n, &r__1, &a[k + a_dim1], nrda, &x[1], &c__1);
	}
	*resnrm = sqrt(sdot_(m, &r__[1], &c__1, &r__[1], &c__1));
	if (it == 1) {
	    goto L100;
	}

/*        TEST FOR CONVERGENCE */

	if (znorm <= acc * *xnorm) {
	    return 0;
	}

/*        COMPARE SUCCESSIVE REFINEMENT VECTOR NORMS */
/*        FOR LOOP TERMINATION CRITERIA */

	if (znorm <= znrm0 * .25f) {
	    goto L100;
	}
	if (it == 2) {
	    goto L90;
	}

	*iflag = 4;
	xermsg_("SLATEC", "LSSODS", "PROBLEM MAY BE ILL-CONDITIONED.  MAXIMA"
		"L MACHINE ACCURACY IS NOT ACHIEVABLE.", &c__3, &c__1, (ftnlen)
		6, (ftnlen)6, (ftnlen)76);
	return 0;

L90:
	*iflag = 5;
	xermsg_("SLATEC", "LSSODS", "PROBLEM IS VERY ILL-CONDITIONED.  ITERA"
		"TIVE IMPROVEMENT IS INEFFECTIVE.", &c__8, &c__1, (ftnlen)6, (
		ftnlen)6, (ftnlen)71);
	return 0;

L100:
	znrm0 = znorm;
    }
/* ********************************************************************** */

/* ********************************************************************** */
    *iflag = 6;
    xermsg_("SLATEC", "LSSODS", "CONVERGENCE HAS NOT BEEN OBTAINED WITH ALLO"
	    "WABLE NUMBER OF ITERATIVE IMPROVEMENT STEPS.", &c__8, &c__1, (
	    ftnlen)6, (ftnlen)6, (ftnlen)87);

    return 0;
} /* lssods_ */

