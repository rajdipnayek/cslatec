/* sspsl.f -- translated by f2c (version 12.02.01).
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

/* DECK SSPSL */
/* Subroutine */ int sspsl_(real *ap, integer *n, integer *kpvt, real *b)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer k;
    static real ak, bk;
    static integer ik, kk, kp;
    static real akm1, bkm1;
    static integer ikm1, km1k, ikp1;
    static real temp;
    extern doublereal sdot_(integer *, real *, integer *, real *, integer *);
    static integer km1km1;
    static real denom;
    extern /* Subroutine */ int saxpy_(integer *, real *, real *, integer *, 
	    real *, integer *);

/* ***BEGIN PROLOGUE  SSPSL */
/* ***PURPOSE  Solve a real symmetric system using the factors obtained */
/*            from SSPFA. */
/* ***LIBRARY   SLATEC (LINPACK) */
/* ***CATEGORY  D2B1A */
/* ***TYPE      SINGLE PRECISION (SSPSL-S, DSPSL-D, CHPSL-C, CSPSL-C) */
/* ***KEYWORDS  LINEAR ALGEBRA, LINPACK, MATRIX, PACKED, SOLVE, SYMMETRIC */
/* ***AUTHOR  Bunch, J., (UCSD) */
/* ***DESCRIPTION */

/*     SSISL solves the real symmetric system */
/*     A * X = B */
/*     using the factors computed by SSPFA. */

/*     On Entry */

/*        AP      REAL(N*(N+1)/2) */
/*                the output from SSPFA. */

/*        N       INTEGER */
/*                the order of the matrix  A . */

/*        KPVT    INTEGER(N) */
/*                the pivot vector from SSPFA. */

/*        B       REAL(N) */
/*                the right hand side vector. */

/*     On Return */

/*        B       the solution vector  X . */

/*     Error Condition */

/*        A division by zero may occur if  SSPCO  has set RCOND .EQ. 0.0 */
/*        or  SSPFA  has set INFO .NE. 0  . */

/*     To compute  INVERSE(A) * C  where  C  is a matrix */
/*     with  P  columns */
/*           CALL SSPFA(AP,N,KPVT,INFO) */
/*           IF (INFO .NE. 0) GO TO ... */
/*           DO 10 J = 1, P */
/*              CALL SSPSL(AP,N,KPVT,C(1,J)) */
/*        10 CONTINUE */

/* ***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W. */
/*                 Stewart, LINPACK Users' Guide, SIAM, 1979. */
/* ***ROUTINES CALLED  SAXPY, SDOT */
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
/* ***END PROLOGUE  SSPSL */


/*     LOOP BACKWARD APPLYING THE TRANSFORMATIONS AND */
/*     D INVERSE TO B. */

/* ***FIRST EXECUTABLE STATEMENT  SSPSL */
    /* Parameter adjustments */
    --b;
    --kpvt;
    --ap;

    /* Function Body */
    k = *n;
    ik = *n * (*n - 1) / 2;
L10:
    if (k == 0) {
	goto L80;
    }
    kk = ik + k;
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
    saxpy_(&i__1, &b[k], &ap[ik + 1], &c__1, &b[1], &c__1);
L30:

/*           APPLY D INVERSE. */

    b[k] /= ap[kk];
    --k;
    ik -= k;
    goto L70;
L40:

/*           2 X 2 PIVOT BLOCK. */

    ikm1 = ik - (k - 1);
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
    saxpy_(&i__1, &b[k], &ap[ik + 1], &c__1, &b[1], &c__1);
    i__1 = k - 2;
    saxpy_(&i__1, &b[k - 1], &ap[ikm1 + 1], &c__1, &b[1], &c__1);
L60:

/*           APPLY D INVERSE. */

    km1k = ik + k - 1;
    kk = ik + k;
    ak = ap[kk] / ap[km1k];
    km1km1 = ikm1 + k - 1;
    akm1 = ap[km1km1] / ap[km1k];
    bk = b[k] / ap[km1k];
    bkm1 = b[k - 1] / ap[km1k];
    denom = ak * akm1 - 1.f;
    b[k] = (akm1 * bk - bkm1) / denom;
    b[k - 1] = (ak * bkm1 - bk) / denom;
    k += -2;
    ik = ik - (k + 1) - k;
L70:
    goto L10;
L80:

/*     LOOP FORWARD APPLYING THE TRANSFORMATIONS. */

    k = 1;
    ik = 0;
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
    b[k] += sdot_(&i__1, &ap[ik + 1], &c__1, &b[1], &c__1);
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
    ik += k;
    ++k;
    goto L150;
L120:

/*           2 X 2 PIVOT BLOCK. */

    if (k == 1) {
	goto L140;
    }

/*              APPLY THE TRANSFORMATION. */

    i__1 = k - 1;
    b[k] += sdot_(&i__1, &ap[ik + 1], &c__1, &b[1], &c__1);
    ikp1 = ik + k;
    i__1 = k - 1;
    b[k + 1] += sdot_(&i__1, &ap[ikp1 + 1], &c__1, &b[1], &c__1);
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
    ik = ik + k + k + 1;
    k += 2;
L150:
    goto L90;
L160:
    return 0;
} /* sspsl_ */

