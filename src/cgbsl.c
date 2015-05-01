/* cgbsl.f -- translated by f2c (version 12.02.01).
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

/* DECK CGBSL */
/* Subroutine */ int cgbsl_(complex *abd, integer *lda, integer *n, integer *
	ml, integer *mu, integer *ipvt, complex *b, integer *job)
{
    /* System generated locals */
    integer abd_dim1, abd_offset, i__1, i__2, i__3;
    complex q__1, q__2, q__3;

    /* Local variables */
    static integer k, l, m;
    static complex t;
    static integer kb, la, lb, lm, nm1;
    extern /* Complex */ void cdotc_(complex *, integer *, complex *, integer 
	    *, complex *, integer *);
    extern /* Subroutine */ int caxpy_(integer *, complex *, complex *, 
	    integer *, complex *, integer *);

/* ***BEGIN PROLOGUE  CGBSL */
/* ***PURPOSE  Solve the complex band system A*X=B or CTRANS(A)*X=B using */
/*            the factors computed by CGBCO or CGBFA. */
/* ***LIBRARY   SLATEC (LINPACK) */
/* ***CATEGORY  D2C2 */
/* ***TYPE      COMPLEX (SGBSL-S, DGBSL-D, CGBSL-C) */
/* ***KEYWORDS  BANDED, LINEAR ALGEBRA, LINPACK, MATRIX, SOLVE */
/* ***AUTHOR  Moler, C. B., (U. of New Mexico) */
/* ***DESCRIPTION */

/*     CGBSL solves the complex band system */
/*     A * X = B  or  CTRANS(A) * X = B */
/*     using the factors computed by CGBCO or CGBFA. */

/*     On Entry */

/*        ABD     COMPLEX(LDA, N) */
/*                the output from CGBCO or CGBFA. */

/*        LDA     INTEGER */
/*                the leading dimension of the array  ABD . */

/*        N       INTEGER */
/*                the order of the original matrix. */

/*        ML      INTEGER */
/*                number of diagonals below the main diagonal. */

/*        MU      INTEGER */
/*                number of diagonals above the main diagonal. */

/*        IPVT    INTEGER(N) */
/*                the pivot vector from CGBCO or CGBFA. */

/*        B       COMPLEX(N) */
/*                the right hand side vector. */

/*        JOB     INTEGER */
/*                = 0         to solve  A*X = B , */
/*                = nonzero   to solve  CTRANS(A)*X = B , where */
/*                            CTRANS(A)  is the conjugate transpose. */

/*     On Return */

/*        B       the solution vector  X . */

/*     Error Condition */

/*        A division by zero will occur if the input factor contains a */
/*        zero on the diagonal.  Technically this indicates singularity */
/*        but it is often caused by improper arguments or improper */
/*        setting of LDA .  It will not occur if the subroutines are */
/*        called correctly and if CGBCO has set RCOND .GT. 0.0 */
/*        or CGBFA has set INFO .EQ. 0 . */

/*     To compute  INVERSE(A) * C  where  C  is a matrix */
/*     with  P  columns */
/*           CALL CGBCO(ABD,LDA,N,ML,MU,IPVT,RCOND,Z) */
/*           IF (RCOND is too small) GO TO ... */
/*           DO 10 J = 1, P */
/*              CALL CGBSL(ABD,LDA,N,ML,MU,IPVT,C(1,J),0) */
/*        10 CONTINUE */

/* ***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W. */
/*                 Stewart, LINPACK Users' Guide, SIAM, 1979. */
/* ***ROUTINES CALLED  CAXPY, CDOTC */
/* ***REVISION HISTORY  (YYMMDD) */
/*   780814  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890831  Modified array declarations.  (WRB) */
/*   890831  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900326  Removed duplicate information from DESCRIPTION section. */
/*           (WRB) */
/*   920501  Reformatted the REFERENCES section.  (WRB) */
/* ***END PROLOGUE  CGBSL */

/* ***FIRST EXECUTABLE STATEMENT  CGBSL */
    /* Parameter adjustments */
    abd_dim1 = *lda;
    abd_offset = 1 + abd_dim1;
    abd -= abd_offset;
    --ipvt;
    --b;

    /* Function Body */
    m = *mu + *ml + 1;
    nm1 = *n - 1;
    if (*job != 0) {
	goto L50;
    }

/*        JOB = 0 , SOLVE  A * X = B */
/*        FIRST SOLVE L*Y = B */

    if (*ml == 0) {
	goto L30;
    }
    if (nm1 < 1) {
	goto L30;
    }
    i__1 = nm1;
    for (k = 1; k <= i__1; ++k) {
/* Computing MIN */
	i__2 = *ml, i__3 = *n - k;
	lm = min(i__2,i__3);
	l = ipvt[k];
	i__2 = l;
	t.r = b[i__2].r, t.i = b[i__2].i;
	if (l == k) {
	    goto L10;
	}
	i__2 = l;
	i__3 = k;
	b[i__2].r = b[i__3].r, b[i__2].i = b[i__3].i;
	i__2 = k;
	b[i__2].r = t.r, b[i__2].i = t.i;
L10:
	caxpy_(&lm, &t, &abd[m + 1 + k * abd_dim1], &c__1, &b[k + 1], &c__1);
/* L20: */
    }
L30:

/*        NOW SOLVE  U*X = Y */

    i__1 = *n;
    for (kb = 1; kb <= i__1; ++kb) {
	k = *n + 1 - kb;
	i__2 = k;
	c_div(&q__1, &b[k], &abd[m + k * abd_dim1]);
	b[i__2].r = q__1.r, b[i__2].i = q__1.i;
	lm = min(k,m) - 1;
	la = m - lm;
	lb = k - lm;
	i__2 = k;
	q__1.r = -b[i__2].r, q__1.i = -b[i__2].i;
	t.r = q__1.r, t.i = q__1.i;
	caxpy_(&lm, &t, &abd[la + k * abd_dim1], &c__1, &b[lb], &c__1);
/* L40: */
    }
    goto L100;
L50:

/*        JOB = NONZERO, SOLVE  CTRANS(A) * X = B */
/*        FIRST SOLVE  CTRANS(U)*Y = B */

    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
	lm = min(k,m) - 1;
	la = m - lm;
	lb = k - lm;
	cdotc_(&q__1, &lm, &abd[la + k * abd_dim1], &c__1, &b[lb], &c__1);
	t.r = q__1.r, t.i = q__1.i;
	i__2 = k;
	i__3 = k;
	q__2.r = b[i__3].r - t.r, q__2.i = b[i__3].i - t.i;
	r_cnjg(&q__3, &abd[m + k * abd_dim1]);
	c_div(&q__1, &q__2, &q__3);
	b[i__2].r = q__1.r, b[i__2].i = q__1.i;
/* L60: */
    }

/*        NOW SOLVE CTRANS(L)*X = Y */

    if (*ml == 0) {
	goto L90;
    }
    if (nm1 < 1) {
	goto L90;
    }
    i__1 = nm1;
    for (kb = 1; kb <= i__1; ++kb) {
	k = *n - kb;
/* Computing MIN */
	i__2 = *ml, i__3 = *n - k;
	lm = min(i__2,i__3);
	i__2 = k;
	i__3 = k;
	cdotc_(&q__2, &lm, &abd[m + 1 + k * abd_dim1], &c__1, &b[k + 1], &
		c__1);
	q__1.r = b[i__3].r + q__2.r, q__1.i = b[i__3].i + q__2.i;
	b[i__2].r = q__1.r, b[i__2].i = q__1.i;
	l = ipvt[k];
	if (l == k) {
	    goto L70;
	}
	i__2 = l;
	t.r = b[i__2].r, t.i = b[i__2].i;
	i__2 = l;
	i__3 = k;
	b[i__2].r = b[i__3].r, b[i__2].i = b[i__3].i;
	i__2 = k;
	b[i__2].r = t.r, b[i__2].i = t.i;
L70:
/* L80: */
	;
    }
L90:
L100:
    return 0;
} /* cgbsl_ */

