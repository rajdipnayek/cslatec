/* dsisl.f -- translated by f2c (version 12.02.01).
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

/* DECK DSISL */
/* Subroutine */ int dsisl_(doublereal *a, integer *lda, integer *n, integer *
	kpvt, doublereal *b)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1;

    /* Local variables */
    static integer k;
    static doublereal ak, bk;
    static integer kp;
    static doublereal akm1, bkm1;
    extern doublereal ddot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    static doublereal temp, denom;
    extern /* Subroutine */ int daxpy_(integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *);

/* ***BEGIN PROLOGUE  DSISL */
/* ***PURPOSE  Solve a real symmetric system using the factors obtained */
/*            from SSIFA. */
/* ***LIBRARY   SLATEC (LINPACK) */
/* ***CATEGORY  D2B1A */
/* ***TYPE      DOUBLE PRECISION (SSISL-S, DSISL-D, CHISL-C, CSISL-C) */
/* ***KEYWORDS  LINEAR ALGEBRA, LINPACK, MATRIX, SOLVE, SYMMETRIC */
/* ***AUTHOR  Bunch, J., (UCSD) */
/* ***DESCRIPTION */

/*     DSISL solves the double precision symmetric system */
/*     A * X = B */
/*     using the factors computed by DSIFA. */

/*     On Entry */

/*        A       DOUBLE PRECISION(LDA,N) */
/*                the output from DSIFA. */

/*        LDA     INTEGER */
/*                the leading dimension of the array  A . */

/*        N       INTEGER */
/*                the order of the matrix  A . */

/*        KPVT    INTEGER(N) */
/*                the pivot vector from DSIFA. */

/*        B       DOUBLE PRECISION(N) */
/*                the right hand side vector. */

/*     On Return */

/*        B       the solution vector  X . */

/*     Error Condition */

/*        A division by zero may occur if  DSICO  has set RCOND .EQ. 0.0 */
/*        or  DSIFA  has set INFO .NE. 0  . */

/*     To compute  INVERSE(A) * C  where  C  is a matrix */
/*     with  P  columns */
/*           CALL DSIFA(A,LDA,N,KPVT,INFO) */
/*           IF (INFO .NE. 0) GO TO ... */
/*           DO 10 J = 1, P */
/*              CALL DSISL(A,LDA,N,KPVT,C(1,J)) */
/*        10 CONTINUE */

/* ***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W. */
/*                 Stewart, LINPACK Users' Guide, SIAM, 1979. */
/* ***ROUTINES CALLED  DAXPY, DDOT */
/* ***REVISION HISTORY  (YYMMDD) */
/*   780814  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890831  Modified array declarations.  (WRB) */
/*   891107  Modified routine equivalence list.  (WRB) */
/*   891107  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900326  Removed duplicate information from DESCRIPTION section. */
/*           (WRB) */
/*   920501  Reformatted the REFERENCES section.  (WRB) */
/* ***END PROLOGUE  DSISL */


/*     LOOP BACKWARD APPLYING THE TRANSFORMATIONS AND */
/*     D INVERSE TO B. */

/* ***FIRST EXECUTABLE STATEMENT  DSISL */
    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --kpvt;
    --b;

    /* Function Body */
    k = *n;
L10:
    if (k == 0) {
	goto L80;
    }
    if (kpvt[k] < 0) {
	goto L40;
    }

/*           1 X 1 PIVOT BLOCK. */

    if (k == 1) {
	goto L30;
    }
    kp = kpvt[k];
    if (kp == k) {
	goto L20;
    }

/*                 INTERCHANGE. */

    temp = b[k];
    b[k] = b[kp];
    b[kp] = temp;
L20:

/*              APPLY THE TRANSFORMATION. */

    i__1 = k - 1;
    daxpy_(&i__1, &b[k], &a[k * a_dim1 + 1], &c__1, &b[1], &c__1);
L30:

/*           APPLY D INVERSE. */

    b[k] /= a[k + k * a_dim1];
    --k;
    goto L70;
L40:

/*           2 X 2 PIVOT BLOCK. */

    if (k == 2) {
	goto L60;
    }
    kp = (i__1 = kpvt[k], abs(i__1));
    if (kp == k - 1) {
	goto L50;
    }

/*                 INTERCHANGE. */

    temp = b[k - 1];
    b[k - 1] = b[kp];
    b[kp] = temp;
L50:

/*              APPLY THE TRANSFORMATION. */

    i__1 = k - 2;
    daxpy_(&i__1, &b[k], &a[k * a_dim1 + 1], &c__1, &b[1], &c__1);
    i__1 = k - 2;
    daxpy_(&i__1, &b[k - 1], &a[(k - 1) * a_dim1 + 1], &c__1, &b[1], &c__1);
L60:

/*           APPLY D INVERSE. */

    ak = a[k + k * a_dim1] / a[k - 1 + k * a_dim1];
    akm1 = a[k - 1 + (k - 1) * a_dim1] / a[k - 1 + k * a_dim1];
    bk = b[k] / a[k - 1 + k * a_dim1];
    bkm1 = b[k - 1] / a[k - 1 + k * a_dim1];
    denom = ak * akm1 - 1.;
    b[k] = (akm1 * bk - bkm1) / denom;
    b[k - 1] = (ak * bkm1 - bk) / denom;
    k += -2;
L70:
    goto L10;
L80:

/*     LOOP FORWARD APPLYING THE TRANSFORMATIONS. */

    k = 1;
L90:
    if (k > *n) {
	goto L160;
    }
    if (kpvt[k] < 0) {
	goto L120;
    }

/*           1 X 1 PIVOT BLOCK. */

    if (k == 1) {
	goto L110;
    }

/*              APPLY THE TRANSFORMATION. */

    i__1 = k - 1;
    b[k] += ddot_(&i__1, &a[k * a_dim1 + 1], &c__1, &b[1], &c__1);
    kp = kpvt[k];
    if (kp == k) {
	goto L100;
    }

/*                 INTERCHANGE. */

    temp = b[k];
    b[k] = b[kp];
    b[kp] = temp;
L100:
L110:
    ++k;
    goto L150;
L120:

/*           2 X 2 PIVOT BLOCK. */

    if (k == 1) {
	goto L140;
    }

/*              APPLY THE TRANSFORMATION. */

    i__1 = k - 1;
    b[k] += ddot_(&i__1, &a[k * a_dim1 + 1], &c__1, &b[1], &c__1);
    i__1 = k - 1;
    b[k + 1] += ddot_(&i__1, &a[(k + 1) * a_dim1 + 1], &c__1, &b[1], &c__1);
    kp = (i__1 = kpvt[k], abs(i__1));
    if (kp == k) {
	goto L130;
    }

/*                 INTERCHANGE. */

    temp = b[k];
    b[k] = b[kp];
    b[kp] = temp;
L130:
L140:
    k += 2;
L150:
    goto L90;
L160:
    return 0;
} /* dsisl_ */

