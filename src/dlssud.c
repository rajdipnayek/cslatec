/* dlssud.f -- translated by f2c (version 12.02.01).
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

static integer c__4 = 4;
static integer c__0 = 0;
static logical c_false = FALSE_;
static integer c__1 = 1;
static integer c__2 = 2;
static doublereal c_b38 = 10.;

/* DECK DLSSUD */
/* Subroutine */ int dlssud_(doublereal *a, doublereal *x, doublereal *b, 
	integer *n, integer *m, integer *nrda, doublereal *u, integer *nrdu, 
	integer *iflag, integer *mlso, integer *irank, integer *iscale, 
	doublereal *q, doublereal *diag, integer *kpivot, doublereal *s, 
	doublereal *div, doublereal *td, integer *isflg, doublereal *scales)
{
    /* System generated locals */
    integer a_dim1, a_offset, q_dim1, q_offset, u_dim1, u_offset, i__1, i__2, 
	    i__3;
    doublereal d__1, d__2, d__3;

    /* Local variables */
    static integer i__, j, k, l, mj, kp, jr, nu;
    static doublereal ss, gam, res;
    static integer irp;
    static doublereal uro;
    static integer nfat;
    extern doublereal ddot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    static integer nmir;
    static doublereal gamma;
    extern /* Subroutine */ int xgetf_(integer *), xsetf_(integer *);
    extern doublereal d1mach_(integer *);
    extern integer j4save_(integer *, integer *, logical *);
    static integer nfatal, maxmes;
    extern /* Subroutine */ int dohtrl_(doublereal *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *), dorthr_(
	    doublereal *, integer *, integer *, integer *, integer *, integer 
	    *, integer *, doublereal *, integer *, doublereal *, doublereal *,
	     doublereal *), xermax_(integer *), xermsg_(char *, char *, char *
	    , integer *, integer *, ftnlen, ftnlen, ftnlen);

/* ***BEGIN PROLOGUE  DLSSUD */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to DBVSUP and DSUDS */
/* ***LIBRARY   SLATEC */
/* ***TYPE      DOUBLE PRECISION (LSSUDS-S, DLSSUD-D) */
/* ***AUTHOR  Watts, H. A., (SNLA) */
/* ***DESCRIPTION */

/*    DLSSUD solves the underdetermined system of equations  A Z = B, */
/*    where A is N by M and N .LE. M.  In particular, if rank A equals */
/*    IRA, a vector X and a matrix U are determined such that X is the */
/*    UNIQUE solution of smallest length, satisfying A X = B, and the */
/*    columns of U form an orthonormal basis for the null space of A, */
/*    satisfying A U = 0 .  Then all solutions Z are given by */
/*              Z = X + C(1)*U(1) + ..... + C(M-IRA)*U(M-IRA) */
/*    where U(J) represents the J-th column of U and the C(J) are */
/*    arbitrary constants. */
/*    If the system of equations are not compatible, only the least */
/*    squares solution of minimal length is computed. */

/* ********************************************************************* */
/*   INPUT */
/* ********************************************************************* */

/*     A -- Contains the matrix of N equations in M unknowns, A remains */
/*          unchanged, must be dimensioned NRDA by M. */
/*     X -- Solution array of length at least M. */
/*     B -- Given constant vector of length N, B remains unchanged. */
/*     N -- Number of equations, N greater or equal to 1. */
/*     M -- Number of unknowns, M greater or equal to N. */
/*     NRDA -- Row dimension of A, NRDA greater or equal to N. */
/*     U -- Matrix used for solution, must be dimensioned NRDU by */
/*          (M - rank of A). */
/*          (storage for U may be ignored when only the minimal length */
/*           solution X is desired) */
/*     NRDU -- Row dimension of U, NRDU greater or equal to M. */
/*             (if only the minimal length solution is wanted, */
/*              NRDU=0 is acceptable) */
/*     IFLAG -- Status indicator */
/*           =0  for the first call (and for each new problem defined by */
/*               a new matrix A) when the matrix data is treated as exact */
/*           =-K for the first call (and for each new problem defined by */
/*               a new matrix A) when the matrix data is assumed to be */
/*               accurate to about K digits. */
/*           =1  for subsequent calls whenever the matrix A has already */
/*               been decomposed (problems with new vectors B but */
/*               same matrix A can be handled efficiently). */
/*     MLSO -- =0 if only the minimal length solution is wanted. */
/*             =1 if the complete solution is wanted, includes the */
/*                linear space defined by the matrix U. */
/*     IRANK -- Variable used for the rank of A, set by the code. */
/*     ISCALE -- Scaling indicator */
/*               =-1 if the matrix A is to be pre-scaled by */
/*               columns when appropriate. */
/*               If the scaling indicator is not equal to -1 */
/*               no scaling will be attempted. */
/*            For most problems scaling will probably not be necessary. */
/*     Q -- Matrix used for the transformation, must be dimensioned */
/*            NRDA by M. */
/*     DIAG,KPIVOT,S, -- Arrays of length at least N used for internal */
/*      DIV,TD,SCALES    storage (except for SCALES which is M). */
/*     ISFLG -- Storage for an internal variable. */

/* ********************************************************************* */
/*   OUTPUT */
/* ********************************************************************* */

/*     IFLAG -- Status indicator */
/*            =1 if solution was obtained. */
/*            =2 if improper input is detected. */
/*            =3 if rank of matrix is less than N. */
/*               To continue, simply reset IFLAG=1 and call DLSSUD again. */
/*            =4 if the system of equations appears to be inconsistent. */
/*               However, the least squares solution of minimal length */
/*               was obtained. */
/*     X -- Minimal length least squares solution of A Z = B */
/*     IRANK -- Numerically determined rank of A, must not be altered */
/*              on succeeding calls with input values of IFLAG=1. */
/*     U -- Matrix whose M-IRANK columns are mutually orthogonal unit */
/*          vectors which span the null space of A. This is to be ignored */
/*          when MLSO was set to zero or IFLAG=4 on output. */
/*     Q -- Contains the strictly upper triangular part of the reduced */
/*           matrix and transformation information. */
/*     DIAG -- Contains the diagonal elements of the triangular reduced */
/*             matrix. */
/*     KPIVOT -- Contains the pivotal information.  The row interchanges */
/*               performed on the original matrix are recorded here. */
/*     S -- Contains the solution of the lower triangular system. */
/*     DIV,TD -- Contains transformation information for rank */
/*               deficient problems. */
/*     SCALES -- Contains the column scaling parameters. */

/* ********************************************************************* */

/* ***SEE ALSO  DBVSUP, DSUDS */
/* ***REFERENCES  H. A. Watts, Solving linear least squares problems */
/*                 using SODS/SUDS/CODS, Sandia Report SAND77-0683, */
/*                 Sandia Laboratories, 1977. */
/* ***ROUTINES CALLED  D1MACH, DDOT, DOHTRL, DORTHR, J4SAVE, XERMAX, */
/*                    XERMSG, XGETF, XSETF */
/* ***REVISION HISTORY  (YYMMDD) */
/*   750601  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ) */
/*   900328  Added TYPE section.  (WRB) */
/*   910408  Updated the AUTHOR and REFERENCES sections.  (WRB) */
/*   920501  Reformatted the REFERENCES section.  (WRB) */
/* ***END PROLOGUE  DLSSUD */

/*     ****************************************************************** */

/*          MACHINE PRECISION (COMPUTER UNIT ROUNDOFF VALUE) IS DEFINED */
/*          BY THE FUNCTION D1MACH. */

/*     ****************************************************************** */

/*     BEGIN BLOCK PERMITTING ...EXITS TO 310 */
/*        BEGIN BLOCK PERMITTING ...EXITS TO 80 */
/* ***FIRST EXECUTABLE STATEMENT  DLSSUD */
    /* Parameter adjustments */
    --x;
    --b;
    q_dim1 = *nrda;
    q_offset = 1 + q_dim1;
    q -= q_offset;
    a_dim1 = *nrda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    u_dim1 = *nrdu;
    u_offset = 1 + u_dim1;
    u -= u_offset;
    --diag;
    --kpivot;
    --s;
    --div;
    --td;
    --scales;

    /* Function Body */
    uro = d1mach_(&c__4);

    if (*n < 1 || *m < *n || *nrda < *n) {
	goto L70;
    }
    if (*nrdu != 0 && *nrdu < *m) {
	goto L70;
    }
    if (*iflag > 0) {
	goto L60;
    }

    xgetf_(&nfatal);
    maxmes = j4save_(&c__4, &c__0, &c_false);
    *isflg = -15;
    if (*iflag == 0) {
	goto L10;
    }
    *isflg = *iflag;
    nfat = -1;
    if (nfatal == 0) {
	nfat = 0;
    }
    xsetf_(&nfat);
    xermax_(&c__1);
L10:

/*                 COPY MATRIX A INTO MATRIX Q */

    i__1 = *m;
    for (k = 1; k <= i__1; ++k) {
	i__2 = *n;
	for (j = 1; j <= i__2; ++j) {
	    q[j + k * q_dim1] = a[j + k * a_dim1];
/* L20: */
	}
/* L30: */
    }

/*                 USE ORTHOGONAL TRANSFORMATIONS TO REDUCE Q TO LOWER */
/*                 TRIANGULAR FORM */

    dorthr_(&q[q_offset], n, m, nrda, iflag, irank, iscale, &diag[1], &kpivot[
	    1], &scales[1], &div[1], &td[1]);

    xsetf_(&nfatal);
    xermax_(&maxmes);
    if (*irank == *n) {
	goto L40;
    }

/*                    FOR RANK DEFICIENT PROBLEMS USE ADDITIONAL */
/*                    ORTHOGONAL TRANSFORMATIONS TO FURTHER REDUCE Q */

    if (*irank != 0) {
	dohtrl_(&q[q_offset], n, nrda, &diag[1], irank, &div[1], &td[1]);
    }
/*     ...............EXIT */
    goto L310;
L40:

/*                 STORE DIVISORS FOR THE TRIANGULAR SOLUTION */

    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
	div[k] = diag[k];
/* L50: */
    }
/*        .........EXIT */
    goto L80;
L60:
/*        ......EXIT */
    if (*iflag == 1) {
	goto L80;
    }
L70:

/*           INVALID INPUT FOR DLSSUD */
    *iflag = 2;
    xermsg_("SLATEC", "DLSSUD", "INVALID IMPUT PARAMETERS.", &c__2, &c__1, (
	    ftnlen)6, (ftnlen)6, (ftnlen)25);
/*     ......EXIT */
    goto L310;
L80:


    if (*irank > 0) {
	goto L130;
    }

/*           SPECIAL CASE FOR THE NULL MATRIX */
    i__1 = *m;
    for (k = 1; k <= i__1; ++k) {
	x[k] = 0.;
	if (*mlso == 0) {
	    goto L100;
	}
	u[k + k * u_dim1] = 1.;
	i__2 = *m;
	for (j = 1; j <= i__2; ++j) {
	    if (j != k) {
		u[j + k * u_dim1] = 0.;
	    }
/* L90: */
	}
L100:
/* L110: */
	;
    }
    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
	if (b[k] > 0.) {
	    *iflag = 4;
	}
/* L120: */
    }
    goto L300;
L130:
/*           BEGIN BLOCK PERMITTING ...EXITS TO 180 */

/*              COPY CONSTANT VECTOR INTO S AFTER FIRST INTERCHANGING */
/*              THE ELEMENTS ACCORDING TO THE PIVOTAL SEQUENCE */

    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
	kp = kpivot[k];
	x[k] = b[kp];
/* L140: */
    }
    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
	s[k] = x[k];
/* L150: */
    }

    irp = *irank + 1;
    nu = 1;
    if (*mlso == 0) {
	nu = 0;
    }
/*           ...EXIT */
    if (*irank == *n) {
	goto L180;
    }

/*              FOR RANK DEFICIENT PROBLEMS WE MUST APPLY THE */
/*              ORTHOGONAL TRANSFORMATION TO S */
/*              WE ALSO CHECK TO SEE IF THE SYSTEM APPEARS TO BE */
/*              INCONSISTENT */

    nmir = *n - *irank;
    ss = ddot_(n, &s[1], &c__1, &s[1], &c__1);
    i__1 = *irank;
    for (l = 1; l <= i__1; ++l) {
	k = irp - l;
	gam = (td[k] * s[k] + ddot_(&nmir, &q[irp + k * q_dim1], &c__1, &s[
		irp], &c__1)) / (td[k] * div[k]);
	s[k] += gam * td[k];
	i__2 = *n;
	for (j = irp; j <= i__2; ++j) {
	    s[j] += gam * q[j + k * q_dim1];
/* L160: */
	}
/* L170: */
    }
    res = ddot_(&nmir, &s[irp], &c__1, &s[irp], &c__1);
/*           ...EXIT */
/* Computing MAX */
    d__2 = pow_di(&c_b38, isflg), d__3 = uro * 10.;
/* Computing 2nd power */
    d__1 = max(d__2,d__3) * 10.;
    if (res <= ss * (d__1 * d__1)) {
	goto L180;
    }

/*              INCONSISTENT SYSTEM */
    *iflag = 4;
    nu = 0;
L180:

/*           APPLY FORWARD SUBSTITUTION TO SOLVE LOWER TRIANGULAR SYSTEM */

    s[1] /= div[1];
    if (*irank < 2) {
	goto L200;
    }
    i__1 = *irank;
    for (k = 2; k <= i__1; ++k) {
	i__2 = k - 1;
	s[k] = (s[k] - ddot_(&i__2, &q[k + q_dim1], nrda, &s[1], &c__1)) / 
		div[k];
/* L190: */
    }
L200:

/*           INITIALIZE X VECTOR AND THEN APPLY ORTHOGONAL TRANSFORMATION */

    i__1 = *m;
    for (k = 1; k <= i__1; ++k) {
	x[k] = 0.;
	if (k <= *irank) {
	    x[k] = s[k];
	}
/* L210: */
    }

    i__1 = *irank;
    for (jr = 1; jr <= i__1; ++jr) {
	j = irp - jr;
	mj = *m - j + 1;
	gamma = ddot_(&mj, &q[j + j * q_dim1], nrda, &x[j], &c__1) / (diag[j] 
		* q[j + j * q_dim1]);
	i__2 = *m;
	for (k = j; k <= i__2; ++k) {
	    x[k] += gamma * q[j + k * q_dim1];
/* L220: */
	}
/* L230: */
    }

/*           RESCALE ANSWERS AS DICTATED */

    i__1 = *m;
    for (k = 1; k <= i__1; ++k) {
	x[k] *= scales[k];
/* L240: */
    }

    if (nu == 0 || *m == *irank) {
	goto L290;
    }

/*              INITIALIZE U MATRIX AND THEN APPLY ORTHOGONAL */
/*              TRANSFORMATION */

    l = *m - *irank;
    i__1 = l;
    for (k = 1; k <= i__1; ++k) {
	i__2 = *m;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    u[i__ + k * u_dim1] = 0.;
	    if (i__ == *irank + k) {
		u[i__ + k * u_dim1] = 1.;
	    }
/* L250: */
	}

	i__2 = *irank;
	for (jr = 1; jr <= i__2; ++jr) {
	    j = irp - jr;
	    mj = *m - j + 1;
	    gamma = ddot_(&mj, &q[j + j * q_dim1], nrda, &u[j + k * u_dim1], &
		    c__1) / (diag[j] * q[j + j * q_dim1]);
	    i__3 = *m;
	    for (i__ = j; i__ <= i__3; ++i__) {
		u[i__ + k * u_dim1] += gamma * q[j + i__ * q_dim1];
/* L260: */
	    }
/* L270: */
	}
/* L280: */
    }
L290:
L300:
L310:

    return 0;
} /* dlssud_ */

