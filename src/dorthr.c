/* dorthr.f -- translated by f2c (version 12.02.01).
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
static integer c__2 = 2;
static integer c__1 = 1;
static doublereal c_b10 = 10.;

/* DECK DORTHR */
/* Subroutine */ int dorthr_(doublereal *a, integer *n, integer *m, integer *
	nrda, integer *iflag, integer *irank, integer *iscale, doublereal *
	diag, integer *kpivot, doublereal *scales, doublereal *rows, 
	doublereal *rs)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3;
    doublereal d__1, d__2;

    /* Local variables */
    static integer j, k, l;
    static doublereal as;
    static integer mk, kp;
    static doublereal acc, akk, sad, sig, dum, uro, rss;
    extern doublereal ddot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    static integer jrow;
    static doublereal diagk, asave, sigma, anorm;
    extern doublereal d1mach_(integer *);
    static doublereal sruro;
    extern /* Subroutine */ int dcscal_(doublereal *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, doublereal *, doublereal *,
	     doublereal *, doublereal *, integer *, integer *), xermsg_(char *
	    , char *, char *, integer *, integer *, ftnlen, ftnlen, ftnlen);

/* ***BEGIN PROLOGUE  DORTHR */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to DBVSUP and DSUDS */
/* ***LIBRARY   SLATEC */
/* ***TYPE      DOUBLE PRECISION (ORTHOR-S, DORTHR-D) */
/* ***AUTHOR  Watts, H. A., (SNLA) */
/* ***DESCRIPTION */

/*   Reduction of the matrix A to lower triangular form by a sequence of */
/*   orthogonal HOUSEHOLDER transformations post-multiplying A. */

/* ********************************************************************* */
/*   INPUT */
/* ********************************************************************* */

/*     A -- Contains the matrix to be decomposed, must be dimensioned */
/*           NRDA by N. */
/*     N -- Number of rows in the matrix, N greater or equal to 1. */
/*     M -- Number of columns in the matrix, M greater or equal to N. */
/*     IFLAG -- Indicates the uncertainty in the matrix data. */
/*             = 0 when the data is to be treated as exact. */
/*             =-K when the data is assumed to be accurate to about */
/*                 K digits. */
/*     ISCALE -- Scaling indicator. */
/*               =-1 if the matrix is to be pre-scaled by */
/*               columns when appropriate. */
/*               Otherwise no scaling will be attempted. */
/*     NRDA -- Row dimension of A, NRDA greater or equal to N. */
/*     DIAG,KPIVOT,ROWS, -- Arrays of length at least N used internally */
/*          RS,SCALES         (except for SCALES which is M). */

/* ********************************************************************* */
/*   OUTPUT */
/* ********************************************************************* */

/*     IFLAG - Status indicator */
/*            =1 for successful decomposition. */
/*            =2 if improper input is detected. */
/*            =3 if rank of the matrix is less than N. */
/*     A -- Contains the reduced matrix in the strictly lower triangular */
/*          part and transformation information. */
/*     IRANK -- Contains the numerically determined matrix rank. */
/*     DIAG -- Contains the diagonal elements of the reduced */
/*             triangular matrix. */
/*     KPIVOT -- Contains the pivotal information, the column */
/*               interchanges performed on the original matrix are */
/*               recorded here. */
/*     SCALES -- Contains the column scaling parameters. */

/* ********************************************************************* */

/* ***SEE ALSO  DBVSUP, DSUDS */
/* ***REFERENCES  G. Golub, Numerical methods for solving linear least */
/*                 squares problems, Numerische Mathematik 7, (1965), */
/*                 pp. 206-216. */
/*               P. Businger and G. Golub, Linear least squares */
/*                 solutions by Householder transformations, Numerische */
/*                 Mathematik  7, (1965), pp. 269-276. */
/* ***ROUTINES CALLED  D1MACH, DCSCAL, DDOT, XERMSG */
/* ***REVISION HISTORY  (YYMMDD) */
/*   750601  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ) */
/*   900328  Added TYPE section.  (WRB) */
/*   910408  Updated the AUTHOR and REFERENCES sections.  (WRB) */
/*   920501  Reformatted the REFERENCES section.  (WRB) */
/* ***END PROLOGUE  DORTHR */

/*     ****************************************************************** */

/*          MACHINE PRECISION (COMPUTER UNIT ROUNDOFF VALUE) IS DEFINED */
/*          BY THE FUNCTION D1MACH. */

/*     ****************************************************************** */

/* ***FIRST EXECUTABLE STATEMENT  DORTHR */
    /* Parameter adjustments */
    a_dim1 = *nrda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --diag;
    --kpivot;
    --scales;
    --rows;
    --rs;

    /* Function Body */
    uro = d1mach_(&c__4);
    if (*m >= *n && *n >= 1 && *nrda >= *n) {
	goto L10;
    }
    *iflag = 2;
    xermsg_("SLATEC", "DORTHR", "INVALID INPUT PARAMETERS.", &c__2, &c__1, (
	    ftnlen)6, (ftnlen)6, (ftnlen)25);
    goto L150;
L10:

    acc = uro * 10.;
    if (*iflag < 0) {
/* Computing MAX */
	d__1 = acc, d__2 = pow_di(&c_b10, iflag);
	acc = max(d__1,d__2);
    }
    sruro = sqrt(uro);
    *iflag = 1;
    *irank = *n;

/*        COMPUTE NORM**2 OF JTH ROW AND A MATRIX NORM */

    anorm = 0.;
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	kpivot[j] = j;
	rows[j] = ddot_(m, &a[j + a_dim1], nrda, &a[j + a_dim1], nrda);
	rs[j] = rows[j];
	anorm += rows[j];
/* L20: */
    }

/*        PERFORM COLUMN SCALING ON A WHEN SPECIFIED */

    dcscal_(&a[a_offset], nrda, n, m, &scales[1], &dum, &rows[1], &rs[1], &
	    anorm, &scales[1], iscale, &c__1);

    anorm = sqrt(anorm);


/*        CONSTRUCTION OF LOWER TRIANGULAR MATRIX AND RECORDING OF */
/*        ORTHOGONAL TRANSFORMATIONS */


    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
/*           BEGIN BLOCK PERMITTING ...EXITS TO 80 */
	mk = *m - k + 1;
/*           ...EXIT */
	if (k == *n) {
	    goto L80;
	}
	kp = k + 1;

/*              SEARCHING FOR PIVOTAL ROW */

	i__2 = *n;
	for (j = k; j <= i__2; ++j) {
/*                 BEGIN BLOCK PERMITTING ...EXITS TO 50 */
	    if (rows[j] >= sruro * rs[j]) {
		goto L30;
	    }
	    rows[j] = ddot_(&mk, &a[j + k * a_dim1], nrda, &a[j + k * a_dim1],
		     nrda);
	    rs[j] = rows[j];
L30:
	    if (j == k) {
		goto L40;
	    }
/*                 ......EXIT */
	    if (sigma >= rows[j] * .99) {
		goto L50;
	    }
L40:
	    sigma = rows[j];
	    jrow = j;
L50:
/* L60: */
	    ;
	}
/*           ...EXIT */
	if (jrow == k) {
	    goto L80;
	}

/*              PERFORM ROW INTERCHANGE */

	l = kpivot[k];
	kpivot[k] = kpivot[jrow];
	kpivot[jrow] = l;
	rows[jrow] = rows[k];
	rows[k] = sigma;
	rss = rs[k];
	rs[k] = rs[jrow];
	rs[jrow] = rss;
	i__2 = *m;
	for (l = 1; l <= i__2; ++l) {
	    asave = a[k + l * a_dim1];
	    a[k + l * a_dim1] = a[jrow + l * a_dim1];
	    a[jrow + l * a_dim1] = asave;
/* L70: */
	}
L80:

/*           CHECK RANK OF THE MATRIX */

	sig = ddot_(&mk, &a[k + k * a_dim1], nrda, &a[k + k * a_dim1], nrda);
	diagk = sqrt(sig);
	if (diagk > acc * anorm) {
	    goto L90;
	}

/*              RANK DEFICIENT PROBLEM */
	*iflag = 3;
	*irank = k - 1;
	xermsg_("SLATEC", "DORTHR", "RANK OF MATRIX IS LESS THAN THE NUMBER "
		"OF ROWS.", &c__1, &c__1, (ftnlen)6, (ftnlen)6, (ftnlen)47);
/*        ......EXIT */
	goto L140;
L90:

/*           CONSTRUCT AND APPLY TRANSFORMATION TO MATRIX A */

	akk = a[k + k * a_dim1];
	if (akk > 0.) {
	    diagk = -diagk;
	}
	diag[k] = diagk;
	a[k + k * a_dim1] = akk - diagk;
	if (k == *n) {
	    goto L120;
	}
	sad = diagk * akk - sig;
	i__2 = *n;
	for (j = kp; j <= i__2; ++j) {
	    as = ddot_(&mk, &a[k + k * a_dim1], nrda, &a[j + k * a_dim1], 
		    nrda) / sad;
	    i__3 = *m;
	    for (l = k; l <= i__3; ++l) {
		a[j + l * a_dim1] += as * a[k + l * a_dim1];
/* L100: */
	    }
/* Computing 2nd power */
	    d__1 = a[j + k * a_dim1];
	    rows[j] -= d__1 * d__1;
/* L110: */
	}
L120:
/* L130: */
	;
    }
L140:
L150:


    return 0;
} /* dorthr_ */

