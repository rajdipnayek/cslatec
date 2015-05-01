/* cptsl.f -- translated by f2c (version 12.02.01).
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

/* DECK CPTSL */
/* Subroutine */ int cptsl_(integer *n, complex *d__, complex *e, complex *b)
{
    /* System generated locals */
    integer i__1, i__2, i__3, i__4, i__5;
    complex q__1, q__2, q__3, q__4;

    /* Local variables */
    static integer k;
    static complex t1, t2;
    static integer ke, kf, kp1, nm1, kbm1, nm1d2;

/* ***BEGIN PROLOGUE  CPTSL */
/* ***PURPOSE  Solve a positive definite tridiagonal linear system. */
/* ***LIBRARY   SLATEC (LINPACK) */
/* ***CATEGORY  D2D2A */
/* ***TYPE      COMPLEX (SPTSL-S, DPTSL-D, CPTSL-C) */
/* ***KEYWORDS  LINEAR ALGEBRA, LINPACK, MATRIX, POSITIVE DEFINITE, SOLVE, */
/*             TRIDIAGONAL */
/* ***AUTHOR  Dongarra, J., (ANL) */
/* ***DESCRIPTION */

/*     CPTSL given a positive definite tridiagonal matrix and a right */
/*     hand side will find the solution. */

/*     On Entry */

/*        N        INTEGER */
/*                 is the order of the tridiagonal matrix. */

/*        D        COMPLEX(N) */
/*                 is the diagonal of the tridiagonal matrix. */
/*                 On output D is destroyed. */

/*        E        COMPLEX(N) */
/*                 is the offdiagonal of the tridiagonal matrix. */
/*                 E(1) through E(N-1) should contain the */
/*                 offdiagonal. */

/*        B        COMPLEX(N) */
/*                 is the right hand side vector. */

/*     On Return */

/*        B        contains the solution. */

/* ***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W. */
/*                 Stewart, LINPACK Users' Guide, SIAM, 1979. */
/* ***ROUTINES CALLED  (NONE) */
/* ***REVISION HISTORY  (YYMMDD) */
/*   780814  DATE WRITTEN */
/*   890505  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900326  Removed duplicate information from DESCRIPTION section. */
/*           (WRB) */
/*   920501  Reformatted the REFERENCES section.  (WRB) */
/* ***END PROLOGUE  CPTSL */


/*     CHECK FOR 1 X 1 CASE */

/* ***FIRST EXECUTABLE STATEMENT  CPTSL */
    /* Parameter adjustments */
    --b;
    --e;
    --d__;

    /* Function Body */
    if (*n != 1) {
	goto L10;
    }
    c_div(&q__1, &b[1], &d__[1]);
    b[1].r = q__1.r, b[1].i = q__1.i;
    goto L70;
L10:
    nm1 = *n - 1;
    nm1d2 = nm1 / 2;
    if (*n == 2) {
	goto L30;
    }
    kbm1 = *n - 1;

/*           ZERO TOP HALF OF SUBDIAGONAL AND BOTTOM HALF OF */
/*           SUPERDIAGONAL */

    i__1 = nm1d2;
    for (k = 1; k <= i__1; ++k) {
	r_cnjg(&q__2, &e[k]);
	c_div(&q__1, &q__2, &d__[k]);
	t1.r = q__1.r, t1.i = q__1.i;
	i__2 = k + 1;
	i__3 = k + 1;
	i__4 = k;
	q__2.r = t1.r * e[i__4].r - t1.i * e[i__4].i, q__2.i = t1.r * e[i__4]
		.i + t1.i * e[i__4].r;
	q__1.r = d__[i__3].r - q__2.r, q__1.i = d__[i__3].i - q__2.i;
	d__[i__2].r = q__1.r, d__[i__2].i = q__1.i;
	i__2 = k + 1;
	i__3 = k + 1;
	i__4 = k;
	q__2.r = t1.r * b[i__4].r - t1.i * b[i__4].i, q__2.i = t1.r * b[i__4]
		.i + t1.i * b[i__4].r;
	q__1.r = b[i__3].r - q__2.r, q__1.i = b[i__3].i - q__2.i;
	b[i__2].r = q__1.r, b[i__2].i = q__1.i;
	c_div(&q__1, &e[kbm1], &d__[kbm1 + 1]);
	t2.r = q__1.r, t2.i = q__1.i;
	i__2 = kbm1;
	i__3 = kbm1;
	r_cnjg(&q__3, &e[kbm1]);
	q__2.r = t2.r * q__3.r - t2.i * q__3.i, q__2.i = t2.r * q__3.i + t2.i 
		* q__3.r;
	q__1.r = d__[i__3].r - q__2.r, q__1.i = d__[i__3].i - q__2.i;
	d__[i__2].r = q__1.r, d__[i__2].i = q__1.i;
	i__2 = kbm1;
	i__3 = kbm1;
	i__4 = kbm1 + 1;
	q__2.r = t2.r * b[i__4].r - t2.i * b[i__4].i, q__2.i = t2.r * b[i__4]
		.i + t2.i * b[i__4].r;
	q__1.r = b[i__3].r - q__2.r, q__1.i = b[i__3].i - q__2.i;
	b[i__2].r = q__1.r, b[i__2].i = q__1.i;
	--kbm1;
/* L20: */
    }
L30:
    kp1 = nm1d2 + 1;

/*        CLEAN UP FOR POSSIBLE 2 X 2 BLOCK AT CENTER */

    if (*n % 2 != 0) {
	goto L40;
    }
    r_cnjg(&q__2, &e[kp1]);
    c_div(&q__1, &q__2, &d__[kp1]);
    t1.r = q__1.r, t1.i = q__1.i;
    i__1 = kp1 + 1;
    i__2 = kp1 + 1;
    i__3 = kp1;
    q__2.r = t1.r * e[i__3].r - t1.i * e[i__3].i, q__2.i = t1.r * e[i__3].i + 
	    t1.i * e[i__3].r;
    q__1.r = d__[i__2].r - q__2.r, q__1.i = d__[i__2].i - q__2.i;
    d__[i__1].r = q__1.r, d__[i__1].i = q__1.i;
    i__1 = kp1 + 1;
    i__2 = kp1 + 1;
    i__3 = kp1;
    q__2.r = t1.r * b[i__3].r - t1.i * b[i__3].i, q__2.i = t1.r * b[i__3].i + 
	    t1.i * b[i__3].r;
    q__1.r = b[i__2].r - q__2.r, q__1.i = b[i__2].i - q__2.i;
    b[i__1].r = q__1.r, b[i__1].i = q__1.i;
    ++kp1;
L40:

/*        BACK SOLVE STARTING AT THE CENTER, GOING TOWARDS THE TOP */
/*        AND BOTTOM */

    i__1 = kp1;
    c_div(&q__1, &b[kp1], &d__[kp1]);
    b[i__1].r = q__1.r, b[i__1].i = q__1.i;
    if (*n == 2) {
	goto L60;
    }
    k = kp1 - 1;
    ke = kp1 + nm1d2 - 1;
    i__1 = ke;
    for (kf = kp1; kf <= i__1; ++kf) {
	i__2 = k;
	i__3 = k;
	i__4 = k;
	i__5 = k + 1;
	q__3.r = e[i__4].r * b[i__5].r - e[i__4].i * b[i__5].i, q__3.i = e[
		i__4].r * b[i__5].i + e[i__4].i * b[i__5].r;
	q__2.r = b[i__3].r - q__3.r, q__2.i = b[i__3].i - q__3.i;
	c_div(&q__1, &q__2, &d__[k]);
	b[i__2].r = q__1.r, b[i__2].i = q__1.i;
	i__2 = kf + 1;
	i__3 = kf + 1;
	r_cnjg(&q__4, &e[kf]);
	i__4 = kf;
	q__3.r = q__4.r * b[i__4].r - q__4.i * b[i__4].i, q__3.i = q__4.r * b[
		i__4].i + q__4.i * b[i__4].r;
	q__2.r = b[i__3].r - q__3.r, q__2.i = b[i__3].i - q__3.i;
	c_div(&q__1, &q__2, &d__[kf + 1]);
	b[i__2].r = q__1.r, b[i__2].i = q__1.i;
	--k;
/* L50: */
    }
L60:
    if (*n % 2 == 0) {
	q__3.r = e[1].r * b[2].r - e[1].i * b[2].i, q__3.i = e[1].r * b[2].i 
		+ e[1].i * b[2].r;
	q__2.r = b[1].r - q__3.r, q__2.i = b[1].i - q__3.i;
	c_div(&q__1, &q__2, &d__[1]);
	b[1].r = q__1.r, b[1].i = q__1.i;
    }
L70:
    return 0;
} /* cptsl_ */

