/* cqrsl.f -- translated by f2c (version 12.02.01).
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

/* DECK CQRSL */
/* Subroutine */ int cqrsl_(complex *x, integer *ldx, integer *n, integer *k, 
	complex *qraux, complex *y, complex *qy, complex *qty, complex *b, 
	complex *rsd, complex *xb, integer *job, integer *info)
{
    /* System generated locals */
    integer x_dim1, x_offset, i__1, i__2, i__3;
    real r__1, r__2;
    complex q__1, q__2, q__3;

    /* Local variables */
    static integer i__, j;
    static complex t;
    static logical cb;
    static integer jj;
    static logical cr;
    static integer ju, kp1;
    static logical cxb, cqy;
    static complex temp;
    static logical cqty;
    extern /* Complex */ void cdotc_(complex *, integer *, complex *, integer 
	    *, complex *, integer *);
    extern /* Subroutine */ int ccopy_(integer *, complex *, integer *, 
	    complex *, integer *), caxpy_(integer *, complex *, complex *, 
	    integer *, complex *, integer *);

/* ***BEGIN PROLOGUE  CQRSL */
/* ***PURPOSE  Apply the output of CQRDC to compute coordinate transfor- */
/*            mations, projections, and least squares solutions. */
/* ***LIBRARY   SLATEC (LINPACK) */
/* ***CATEGORY  D9, D2C1 */
/* ***TYPE      COMPLEX (SQRSL-S, DQRSL-D, CQRSL-C) */
/* ***KEYWORDS  LINEAR ALGEBRA, LINPACK, MATRIX, ORTHOGONAL TRIANGULAR, */
/*             SOLVE */
/* ***AUTHOR  Stewart, G. W., (U. of Maryland) */
/* ***DESCRIPTION */

/*     CQRSL applies the output of CQRDC to compute coordinate */
/*     transformations, projections, and least squares solutions. */
/*     For K .LE. MIN(N,P), let XK be the matrix */

/*            XK = (X(JVPT(1)),X(JVPT(2)), ... ,X(JVPT(K))) */

/*     formed from columns JVPT(1), ... ,JVPT(K) of the original */
/*     N x P matrix X that was input to CQRDC (if no pivoting was */
/*     done, XK consists of the first K columns of X in their */
/*     original order).  CQRDC produces a factored unitary matrix Q */
/*     and an upper triangular matrix R such that */

/*              XK = Q * (R) */
/*                       (0) */

/*     This information is contained in coded form in the arrays */
/*     X and QRAUX. */

/*     On Entry */

/*        X      COMPLEX(LDX,P). */
/*               X contains the output of CQRDC. */

/*        LDX    INTEGER. */
/*               LDX is the leading dimension of the array X. */

/*        N      INTEGER. */
/*               N is the number of rows of the matrix XK.  It must */
/*               have the same value as N in CQRDC. */

/*        K      INTEGER. */
/*               K is the number of columns of the matrix XK.  K */
/*               must not be greater than (N,P), where P is the */
/*               same as in the calling sequence to CQRDC. */

/*        QRAUX  COMPLEX(P). */
/*               QRAUX contains the auxiliary output from CQRDC. */

/*        Y      COMPLEX(N) */
/*               Y contains an N-vector that is to be manipulated */
/*               by CQRSL. */

/*        JOB    INTEGER. */
/*               JOB specifies what is to be computed.  JOB has */
/*               the decimal expansion ABCDE, with the following */
/*               meaning. */

/*                    If A .NE. 0, compute QY. */
/*                    If B,C,D, or E .NE. 0, compute QTY. */
/*                    If C .NE. 0, compute B. */
/*                    If D .NE. 0, compute RSD . */
/*                    If E .NE. 0, compute  XB. */

/*               Note that a request to compute B, RSD, or XB */
/*               automatically triggers the computation of QTY, for */
/*               which an array must be provided in the calling */
/*               sequence. */

/*     On Return */

/*        QY     COMPLEX(N). */
/*               QY contains Q*Y, if its computation has been */
/*               requested. */

/*        QTY    COMPLEX(N). */
/*               QTY contains CTRANS(Q)*Y, if its computation has */
/*               been requested.  Here CTRANS(Q) is the conjugate */
/*               transpose of the matrix Q. */

/*        B      COMPLEX(K) */
/*               B contains the solution of the least squares problem */

/*                    minimize NORM2(Y - XK*B), */

/*               if its computation has been requested.  (Note that */
/*               if pivoting was requested in CQRDC, the J-th */
/*               component of B will be associated with column JVPT(J) */
/*               of the original matrix X that was input into CQRDC.) */

/*        RSD    COMPLEX(N). */
/*               RSD contains the least squares residual Y - XK*B, */
/*               if its computation has been requested.  RSD is */
/*               also the orthogonal projection of Y onto the */
/*               orthogonal complement of the column space of XK. */

/*        XB     COMPLEX(N). */
/*               XB contains the least squares approximation XK*B, */
/*               if its computation has been requested.  XB is also */
/*               the orthogonal projection of Y onto the column space */
/*               of X. */

/*        INFO   INTEGER. */
/*               INFO is zero unless the computation of B has */
/*               been requested and R is exactly singular.  In */
/*               this case, INFO is the index of the first zero */
/*               diagonal element of R and B is left unaltered. */

/*     The parameters QY, QTY, B, RSD, and XB are not referenced */
/*     if their computation is not requested and in this case */
/*     can be replaced by dummy variables in the calling program. */
/*     To save storage, the user may in some cases use the same */
/*     array for different parameters in the calling sequence.  A */
/*     frequently occurring example is when one wishes to compute */
/*     any of B, RSD, or XB and does not need Y or QTY.  In this */
/*     case one may identify Y, QTY, and one of B, RSD, or XB, while */
/*     providing separate arrays for anything else that is to be */
/*     computed.  Thus the calling sequence */

/*          CALL CQRSL(X,LDX,N,K,QRAUX,Y,DUM,Y,B,Y,DUM,110,INFO) */

/*     will result in the computation of B and RSD, with RSD */
/*     overwriting Y.  More generally, each item in the following */
/*     list contains groups of permissible identifications for */
/*     a single calling sequence. */

/*          1. (Y,QTY,B) (RSD) (XB) (QY) */

/*          2. (Y,QTY,RSD) (B) (XB) (QY) */

/*          3. (Y,QTY,XB) (B) (RSD) (QY) */

/*          4. (Y,QY) (QTY,B) (RSD) (XB) */

/*          5. (Y,QY) (QTY,RSD) (B) (XB) */

/*          6. (Y,QY) (QTY,XB) (B) (RSD) */

/*     In any group the value returned in the array allocated to */
/*     the group corresponds to the last member of the group. */

/* ***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W. */
/*                 Stewart, LINPACK Users' Guide, SIAM, 1979. */
/* ***ROUTINES CALLED  CAXPY, CCOPY, CDOTC */
/* ***REVISION HISTORY  (YYMMDD) */
/*   780814  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890831  Modified array declarations.  (WRB) */
/*   890831  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900326  Removed duplicate information from DESCRIPTION section. */
/*           (WRB) */
/*   920501  Reformatted the REFERENCES section.  (WRB) */
/* ***END PROLOGUE  CQRSL */

/* ***FIRST EXECUTABLE STATEMENT  CQRSL */

/*     SET INFO FLAG. */

    /* Parameter adjustments */
    x_dim1 = *ldx;
    x_offset = 1 + x_dim1;
    x -= x_offset;
    --qraux;
    --y;
    --qy;
    --qty;
    --b;
    --rsd;
    --xb;

    /* Function Body */
    *info = 0;

/*     DETERMINE WHAT IS TO BE COMPUTED. */

    cqy = *job / 10000 != 0;
    cqty = *job % 10000 != 0;
    cb = *job % 1000 / 100 != 0;
    cr = *job % 100 / 10 != 0;
    cxb = *job % 10 != 0;
/* Computing MIN */
    i__1 = *k, i__2 = *n - 1;
    ju = min(i__1,i__2);

/*     SPECIAL ACTION WHEN N=1. */

    if (ju != 0) {
	goto L40;
    }
    if (cqy) {
	qy[1].r = y[1].r, qy[1].i = y[1].i;
    }
    if (cqty) {
	qty[1].r = y[1].r, qty[1].i = y[1].i;
    }
    if (cxb) {
	xb[1].r = y[1].r, xb[1].i = y[1].i;
    }
    if (! cb) {
	goto L30;
    }
    i__1 = x_dim1 + 1;
    if ((r__1 = x[i__1].r, dabs(r__1)) + (r__2 = r_imag(&x[x_dim1 + 1]), dabs(
	    r__2)) != 0.f) {
	goto L10;
    }
    *info = 1;
    goto L20;
L10:
    c_div(&q__1, &y[1], &x[x_dim1 + 1]);
    b[1].r = q__1.r, b[1].i = q__1.i;
L20:
L30:
    if (cr) {
	rsd[1].r = 0.f, rsd[1].i = 0.f;
    }
    goto L250;
L40:

/*        SET UP TO COMPUTE QY OR QTY. */

    if (cqy) {
	ccopy_(n, &y[1], &c__1, &qy[1], &c__1);
    }
    if (cqty) {
	ccopy_(n, &y[1], &c__1, &qty[1], &c__1);
    }
    if (! cqy) {
	goto L70;
    }

/*           COMPUTE QY. */

    i__1 = ju;
    for (jj = 1; jj <= i__1; ++jj) {
	j = ju - jj + 1;
	i__2 = j;
	if ((r__1 = qraux[i__2].r, dabs(r__1)) + (r__2 = r_imag(&qraux[j]), 
		dabs(r__2)) == 0.f) {
	    goto L50;
	}
	i__2 = j + j * x_dim1;
	temp.r = x[i__2].r, temp.i = x[i__2].i;
	i__2 = j + j * x_dim1;
	i__3 = j;
	x[i__2].r = qraux[i__3].r, x[i__2].i = qraux[i__3].i;
	i__2 = *n - j + 1;
	cdotc_(&q__3, &i__2, &x[j + j * x_dim1], &c__1, &qy[j], &c__1);
	q__2.r = -q__3.r, q__2.i = -q__3.i;
	c_div(&q__1, &q__2, &x[j + j * x_dim1]);
	t.r = q__1.r, t.i = q__1.i;
	i__2 = *n - j + 1;
	caxpy_(&i__2, &t, &x[j + j * x_dim1], &c__1, &qy[j], &c__1);
	i__2 = j + j * x_dim1;
	x[i__2].r = temp.r, x[i__2].i = temp.i;
L50:
/* L60: */
	;
    }
L70:
    if (! cqty) {
	goto L100;
    }

/*           COMPUTE CTRANS(Q)*Y. */

    i__1 = ju;
    for (j = 1; j <= i__1; ++j) {
	i__2 = j;
	if ((r__1 = qraux[i__2].r, dabs(r__1)) + (r__2 = r_imag(&qraux[j]), 
		dabs(r__2)) == 0.f) {
	    goto L80;
	}
	i__2 = j + j * x_dim1;
	temp.r = x[i__2].r, temp.i = x[i__2].i;
	i__2 = j + j * x_dim1;
	i__3 = j;
	x[i__2].r = qraux[i__3].r, x[i__2].i = qraux[i__3].i;
	i__2 = *n - j + 1;
	cdotc_(&q__3, &i__2, &x[j + j * x_dim1], &c__1, &qty[j], &c__1);
	q__2.r = -q__3.r, q__2.i = -q__3.i;
	c_div(&q__1, &q__2, &x[j + j * x_dim1]);
	t.r = q__1.r, t.i = q__1.i;
	i__2 = *n - j + 1;
	caxpy_(&i__2, &t, &x[j + j * x_dim1], &c__1, &qty[j], &c__1);
	i__2 = j + j * x_dim1;
	x[i__2].r = temp.r, x[i__2].i = temp.i;
L80:
/* L90: */
	;
    }
L100:

/*        SET UP TO COMPUTE B, RSD, OR XB. */

    if (cb) {
	ccopy_(k, &qty[1], &c__1, &b[1], &c__1);
    }
    kp1 = *k + 1;
    if (cxb) {
	ccopy_(k, &qty[1], &c__1, &xb[1], &c__1);
    }
    if (cr && *k < *n) {
	i__1 = *n - *k;
	ccopy_(&i__1, &qty[kp1], &c__1, &rsd[kp1], &c__1);
    }
    if (! cxb || kp1 > *n) {
	goto L120;
    }
    i__1 = *n;
    for (i__ = kp1; i__ <= i__1; ++i__) {
	i__2 = i__;
	xb[i__2].r = 0.f, xb[i__2].i = 0.f;
/* L110: */
    }
L120:
    if (! cr) {
	goto L140;
    }
    i__1 = *k;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = i__;
	rsd[i__2].r = 0.f, rsd[i__2].i = 0.f;
/* L130: */
    }
L140:
    if (! cb) {
	goto L190;
    }

/*           COMPUTE B. */

    i__1 = *k;
    for (jj = 1; jj <= i__1; ++jj) {
	j = *k - jj + 1;
	i__2 = j + j * x_dim1;
	if ((r__1 = x[i__2].r, dabs(r__1)) + (r__2 = r_imag(&x[j + j * x_dim1]
		), dabs(r__2)) != 0.f) {
	    goto L150;
	}
	*info = j;
	goto L180;
L150:
	i__2 = j;
	c_div(&q__1, &b[j], &x[j + j * x_dim1]);
	b[i__2].r = q__1.r, b[i__2].i = q__1.i;
	if (j == 1) {
	    goto L160;
	}
	i__2 = j;
	q__1.r = -b[i__2].r, q__1.i = -b[i__2].i;
	t.r = q__1.r, t.i = q__1.i;
	i__2 = j - 1;
	caxpy_(&i__2, &t, &x[j * x_dim1 + 1], &c__1, &b[1], &c__1);
L160:
/* L170: */
	;
    }
L180:
L190:
    if (! cr && ! cxb) {
	goto L240;
    }

/*           COMPUTE RSD OR XB AS REQUIRED. */

    i__1 = ju;
    for (jj = 1; jj <= i__1; ++jj) {
	j = ju - jj + 1;
	i__2 = j;
	if ((r__1 = qraux[i__2].r, dabs(r__1)) + (r__2 = r_imag(&qraux[j]), 
		dabs(r__2)) == 0.f) {
	    goto L220;
	}
	i__2 = j + j * x_dim1;
	temp.r = x[i__2].r, temp.i = x[i__2].i;
	i__2 = j + j * x_dim1;
	i__3 = j;
	x[i__2].r = qraux[i__3].r, x[i__2].i = qraux[i__3].i;
	if (! cr) {
	    goto L200;
	}
	i__2 = *n - j + 1;
	cdotc_(&q__3, &i__2, &x[j + j * x_dim1], &c__1, &rsd[j], &c__1);
	q__2.r = -q__3.r, q__2.i = -q__3.i;
	c_div(&q__1, &q__2, &x[j + j * x_dim1]);
	t.r = q__1.r, t.i = q__1.i;
	i__2 = *n - j + 1;
	caxpy_(&i__2, &t, &x[j + j * x_dim1], &c__1, &rsd[j], &c__1);
L200:
	if (! cxb) {
	    goto L210;
	}
	i__2 = *n - j + 1;
	cdotc_(&q__3, &i__2, &x[j + j * x_dim1], &c__1, &xb[j], &c__1);
	q__2.r = -q__3.r, q__2.i = -q__3.i;
	c_div(&q__1, &q__2, &x[j + j * x_dim1]);
	t.r = q__1.r, t.i = q__1.i;
	i__2 = *n - j + 1;
	caxpy_(&i__2, &t, &x[j + j * x_dim1], &c__1, &xb[j], &c__1);
L210:
	i__2 = j + j * x_dim1;
	x[i__2].r = temp.r, x[i__2].i = temp.i;
L220:
/* L230: */
	;
    }
L240:
L250:
    return 0;
} /* cqrsl_ */

