/* csisl.f -- translated by f2c (version 12.02.01).
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

/* DECK CSISL */
/* Subroutine */ int csisl_(complex *a, integer *lda, integer *n, integer *
	kpvt, complex *b)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3;
    complex q__1, q__2, q__3;

    /* Local variables */
    static integer k;
    static complex ak, bk;
    static integer kp;
    static complex akm1, bkm1, temp, denom;
    extern /* Complex */ void cdotu_(complex *, integer *, complex *, integer 
	    *, complex *, integer *);
    extern /* Subroutine */ int caxpy_(integer *, complex *, complex *, 
	    integer *, complex *, integer *);

/* ***BEGIN PROLOGUE  CSISL */
/* ***PURPOSE  Solve a complex symmetric system using the factors obtained */
/*            from CSIFA. */
/* ***LIBRARY   SLATEC (LINPACK) */
/* ***CATEGORY  D2C1 */
/* ***TYPE      COMPLEX (SSISL-S, DSISL-D, CHISL-C, CSISL-C) */
/* ***KEYWORDS  LINEAR ALGEBRA, LINPACK, MATRIX, SOLVE, SYMMETRIC */
/* ***AUTHOR  Bunch, J., (UCSD) */
/* ***DESCRIPTION */

/*     CSISL solves the complex symmetric system */
/*     A * X = B */
/*     using the factors computed by CSIFA. */

/*     On Entry */

/*        A       COMPLEX(LDA,N) */
/*                the output from CSIFA. */

/*        LDA     INTEGER */
/*                the leading dimension of the array  A . */

/*        N       INTEGER */
/*                the order of the matrix  A . */

/*        KVPT    INTEGER(N) */
/*                the pivot vector from CSIFA. */

/*        B       COMPLEX(N) */
/*                the right hand side vector. */

/*     On Return */

/*        B       the solution vector  X . */

/*     Error Condition */

/*        A division by zero may occur if  CSICO  has set RCOND .EQ. 0.0 */
/*        or  CSIFA  has set INFO .NE. 0  . */

/*     To compute  INVERSE(A) * C  where  C  is a matrix */
/*     with  P  columns */
/*           CALL CSIFA(A,LDA,N,KVPT,INFO) */
/*           If (INFO .NE. 0) GO TO ... */
/*           DO 10 J = 1, P */
/*              CALL CSISL(A,LDA,N,KVPT,C(1,j)) */
/*        10 CONTINUE */

/* ***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W. */
/*                 Stewart, LINPACK Users' Guide, SIAM, 1979. */
/* ***ROUTINES CALLED  CAXPY, CDOTU */
/* ***REVISION HISTORY  (YYMMDD) */
/*   780814  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890831  Modified array declarations.  (WRB) */
/*   891107  Corrected category and modified routine equivalence */
/*           list.  (WRB) */
/*   891107  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900326  Removed duplicate information from DESCRIPTION section. */
/*           (WRB) */
/*   920501  Reformatted the REFERENCES section.  (WRB) */
/* ***END PROLOGUE  CSISL */


/*     LOOP BACKWARD APPLYING THE TRANSFORMATIONS AND */
/*     D INVERSE TO B. */

/* ***FIRST EXECUTABLE STATEMENT  CSISL */
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

    i__1 = k;
    temp.r = b[i__1].r, temp.i = b[i__1].i;
    i__1 = k;
    i__2 = kp;
    b[i__1].r = b[i__2].r, b[i__1].i = b[i__2].i;
    i__1 = kp;
    b[i__1].r = temp.r, b[i__1].i = temp.i;
L20:

/*              APPLY THE TRANSFORMATION. */

    i__1 = k - 1;
    caxpy_(&i__1, &b[k], &a[k * a_dim1 + 1], &c__1, &b[1], &c__1);
L30:

/*           APPLY D INVERSE. */

    i__1 = k;
    c_div(&q__1, &b[k], &a[k + k * a_dim1]);
    b[i__1].r = q__1.r, b[i__1].i = q__1.i;
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

    i__1 = k - 1;
    temp.r = b[i__1].r, temp.i = b[i__1].i;
    i__1 = k - 1;
    i__2 = kp;
    b[i__1].r = b[i__2].r, b[i__1].i = b[i__2].i;
    i__1 = kp;
    b[i__1].r = temp.r, b[i__1].i = temp.i;
L50:

/*              APPLY THE TRANSFORMATION. */

    i__1 = k - 2;
    caxpy_(&i__1, &b[k], &a[k * a_dim1 + 1], &c__1, &b[1], &c__1);
    i__1 = k - 2;
    caxpy_(&i__1, &b[k - 1], &a[(k - 1) * a_dim1 + 1], &c__1, &b[1], &c__1);
L60:

/*           APPLY D INVERSE. */

    c_div(&q__1, &a[k + k * a_dim1], &a[k - 1 + k * a_dim1]);
    ak.r = q__1.r, ak.i = q__1.i;
    c_div(&q__1, &a[k - 1 + (k - 1) * a_dim1], &a[k - 1 + k * a_dim1]);
    akm1.r = q__1.r, akm1.i = q__1.i;
    c_div(&q__1, &b[k], &a[k - 1 + k * a_dim1]);
    bk.r = q__1.r, bk.i = q__1.i;
    c_div(&q__1, &b[k - 1], &a[k - 1 + k * a_dim1]);
    bkm1.r = q__1.r, bkm1.i = q__1.i;
    q__2.r = ak.r * akm1.r - ak.i * akm1.i, q__2.i = ak.r * akm1.i + ak.i * 
	    akm1.r;
    q__1.r = q__2.r - 1.f, q__1.i = q__2.i;
    denom.r = q__1.r, denom.i = q__1.i;
    i__1 = k;
    q__3.r = akm1.r * bk.r - akm1.i * bk.i, q__3.i = akm1.r * bk.i + akm1.i * 
	    bk.r;
    q__2.r = q__3.r - bkm1.r, q__2.i = q__3.i - bkm1.i;
    c_div(&q__1, &q__2, &denom);
    b[i__1].r = q__1.r, b[i__1].i = q__1.i;
    i__1 = k - 1;
    q__3.r = ak.r * bkm1.r - ak.i * bkm1.i, q__3.i = ak.r * bkm1.i + ak.i * 
	    bkm1.r;
    q__2.r = q__3.r - bk.r, q__2.i = q__3.i - bk.i;
    c_div(&q__1, &q__2, &denom);
    b[i__1].r = q__1.r, b[i__1].i = q__1.i;
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

    i__1 = k;
    i__2 = k;
    i__3 = k - 1;
    cdotu_(&q__2, &i__3, &a[k * a_dim1 + 1], &c__1, &b[1], &c__1);
    q__1.r = b[i__2].r + q__2.r, q__1.i = b[i__2].i + q__2.i;
    b[i__1].r = q__1.r, b[i__1].i = q__1.i;
    kp = kpvt[k];
    if (kp == k) {
	goto L100;
    }

/*                 INTERCHANGE. */

    i__1 = k;
    temp.r = b[i__1].r, temp.i = b[i__1].i;
    i__1 = k;
    i__2 = kp;
    b[i__1].r = b[i__2].r, b[i__1].i = b[i__2].i;
    i__1 = kp;
    b[i__1].r = temp.r, b[i__1].i = temp.i;
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

    i__1 = k;
    i__2 = k;
    i__3 = k - 1;
    cdotu_(&q__2, &i__3, &a[k * a_dim1 + 1], &c__1, &b[1], &c__1);
    q__1.r = b[i__2].r + q__2.r, q__1.i = b[i__2].i + q__2.i;
    b[i__1].r = q__1.r, b[i__1].i = q__1.i;
    i__1 = k + 1;
    i__2 = k + 1;
    i__3 = k - 1;
    cdotu_(&q__2, &i__3, &a[(k + 1) * a_dim1 + 1], &c__1, &b[1], &c__1);
    q__1.r = b[i__2].r + q__2.r, q__1.i = b[i__2].i + q__2.i;
    b[i__1].r = q__1.r, b[i__1].i = q__1.i;
    kp = (i__1 = kpvt[k], abs(i__1));
    if (kp == k) {
	goto L130;
    }

/*                 INTERCHANGE. */

    i__1 = k;
    temp.r = b[i__1].r, temp.i = b[i__1].i;
    i__1 = k;
    i__2 = kp;
    b[i__1].r = b[i__2].r, b[i__1].i = b[i__2].i;
    i__1 = kp;
    b[i__1].r = temp.r, b[i__1].i = temp.i;
L130:
L140:
    k += 2;
L150:
    goto L90;
L160:
    return 0;
} /* csisl_ */

