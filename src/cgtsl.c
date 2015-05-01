/* cgtsl.f -- translated by f2c (version 12.02.01).
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

/* DECK CGTSL */
/* Subroutine */ int cgtsl_(integer *n, complex *c__, complex *d__, complex *
	e, complex *b, integer *info)
{
    /* System generated locals */
    integer i__1, i__2, i__3, i__4, i__5, i__6, i__7;
    real r__1, r__2, r__3, r__4;
    complex q__1, q__2, q__3, q__4, q__5;

    /* Local variables */
    static integer k;
    static complex t;
    static integer kb, kp1, nm1, nm2;

/* ***BEGIN PROLOGUE  CGTSL */
/* ***PURPOSE  Solve a tridiagonal linear system. */
/* ***LIBRARY   SLATEC (LINPACK) */
/* ***CATEGORY  D2C2A */
/* ***TYPE      COMPLEX (SGTSL-S, DGTSL-D, CGTSL-C) */
/* ***KEYWORDS  LINEAR ALGEBRA, LINPACK, MATRIX, SOLVE, TRIDIAGONAL */
/* ***AUTHOR  Dongarra, J., (ANL) */
/* ***DESCRIPTION */

/*     CGTSL given a general tridiagonal matrix and a right hand */
/*     side will find the solution. */

/*     On Entry */

/*        N       INTEGER */
/*                is the order of the tridiagonal matrix. */

/*        C       COMPLEX(N) */
/*                is the subdiagonal of the tridiagonal matrix. */
/*                C(2) through C(N) should contain the subdiagonal. */
/*                On output C is destroyed. */

/*        D       COMPLEX(N) */
/*                is the diagonal of the tridiagonal matrix. */
/*                On output D is destroyed. */

/*        E       COMPLEX(N) */
/*                is the superdiagonal of the tridiagonal matrix. */
/*                E(1) through E(N-1) should contain the superdiagonal. */
/*                On output E is destroyed. */

/*        B       COMPLEX(N) */
/*                is the right hand side vector. */

/*     On Return */

/*        B       is the solution vector. */

/*        INFO    INTEGER */
/*                = 0 normal value. */
/*                = K if the K-th element of the diagonal becomes */
/*                    exactly zero.  The subroutine returns when */
/*                    this is detected. */

/* ***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W. */
/*                 Stewart, LINPACK Users' Guide, SIAM, 1979. */
/* ***ROUTINES CALLED  (NONE) */
/* ***REVISION HISTORY  (YYMMDD) */
/*   780814  DATE WRITTEN */
/*   890831  Modified array declarations.  (WRB) */
/*   890831  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900326  Removed duplicate information from DESCRIPTION section. */
/*           (WRB) */
/*   920501  Reformatted the REFERENCES section.  (WRB) */
/* ***END PROLOGUE  CGTSL */

/* ***FIRST EXECUTABLE STATEMENT  CGTSL */
    /* Parameter adjustments */
    --b;
    --e;
    --d__;
    --c__;

    /* Function Body */
    *info = 0;
    c__[1].r = d__[1].r, c__[1].i = d__[1].i;
    nm1 = *n - 1;
    if (nm1 < 1) {
	goto L40;
    }
    d__[1].r = e[1].r, d__[1].i = e[1].i;
    e[1].r = 0.f, e[1].i = 0.f;
    i__1 = *n;
    e[i__1].r = 0.f, e[i__1].i = 0.f;

    i__1 = nm1;
    for (k = 1; k <= i__1; ++k) {
	kp1 = k + 1;

/*              FIND THE LARGEST OF THE TWO ROWS */

	i__2 = kp1;
	i__3 = k;
	if ((r__1 = c__[i__2].r, dabs(r__1)) + (r__2 = r_imag(&c__[kp1]), 
		dabs(r__2)) < (r__3 = c__[i__3].r, dabs(r__3)) + (r__4 = 
		r_imag(&c__[k]), dabs(r__4))) {
	    goto L10;
	}

/*                 INTERCHANGE ROW */

	i__2 = kp1;
	t.r = c__[i__2].r, t.i = c__[i__2].i;
	i__2 = kp1;
	i__3 = k;
	c__[i__2].r = c__[i__3].r, c__[i__2].i = c__[i__3].i;
	i__2 = k;
	c__[i__2].r = t.r, c__[i__2].i = t.i;
	i__2 = kp1;
	t.r = d__[i__2].r, t.i = d__[i__2].i;
	i__2 = kp1;
	i__3 = k;
	d__[i__2].r = d__[i__3].r, d__[i__2].i = d__[i__3].i;
	i__2 = k;
	d__[i__2].r = t.r, d__[i__2].i = t.i;
	i__2 = kp1;
	t.r = e[i__2].r, t.i = e[i__2].i;
	i__2 = kp1;
	i__3 = k;
	e[i__2].r = e[i__3].r, e[i__2].i = e[i__3].i;
	i__2 = k;
	e[i__2].r = t.r, e[i__2].i = t.i;
	i__2 = kp1;
	t.r = b[i__2].r, t.i = b[i__2].i;
	i__2 = kp1;
	i__3 = k;
	b[i__2].r = b[i__3].r, b[i__2].i = b[i__3].i;
	i__2 = k;
	b[i__2].r = t.r, b[i__2].i = t.i;
L10:

/*              ZERO ELEMENTS */

	i__2 = k;
	if ((r__1 = c__[i__2].r, dabs(r__1)) + (r__2 = r_imag(&c__[k]), dabs(
		r__2)) != 0.f) {
	    goto L20;
	}
	*info = k;
	goto L100;
L20:
	i__2 = kp1;
	q__2.r = -c__[i__2].r, q__2.i = -c__[i__2].i;
	c_div(&q__1, &q__2, &c__[k]);
	t.r = q__1.r, t.i = q__1.i;
	i__2 = kp1;
	i__3 = kp1;
	i__4 = k;
	q__2.r = t.r * d__[i__4].r - t.i * d__[i__4].i, q__2.i = t.r * d__[
		i__4].i + t.i * d__[i__4].r;
	q__1.r = d__[i__3].r + q__2.r, q__1.i = d__[i__3].i + q__2.i;
	c__[i__2].r = q__1.r, c__[i__2].i = q__1.i;
	i__2 = kp1;
	i__3 = kp1;
	i__4 = k;
	q__2.r = t.r * e[i__4].r - t.i * e[i__4].i, q__2.i = t.r * e[i__4].i 
		+ t.i * e[i__4].r;
	q__1.r = e[i__3].r + q__2.r, q__1.i = e[i__3].i + q__2.i;
	d__[i__2].r = q__1.r, d__[i__2].i = q__1.i;
	i__2 = kp1;
	e[i__2].r = 0.f, e[i__2].i = 0.f;
	i__2 = kp1;
	i__3 = kp1;
	i__4 = k;
	q__2.r = t.r * b[i__4].r - t.i * b[i__4].i, q__2.i = t.r * b[i__4].i 
		+ t.i * b[i__4].r;
	q__1.r = b[i__3].r + q__2.r, q__1.i = b[i__3].i + q__2.i;
	b[i__2].r = q__1.r, b[i__2].i = q__1.i;
/* L30: */
    }
L40:
    i__1 = *n;
    if ((r__1 = c__[i__1].r, dabs(r__1)) + (r__2 = r_imag(&c__[*n]), dabs(
	    r__2)) != 0.f) {
	goto L50;
    }
    *info = *n;
    goto L90;
L50:

/*           BACK SOLVE */

    nm2 = *n - 2;
    i__1 = *n;
    c_div(&q__1, &b[*n], &c__[*n]);
    b[i__1].r = q__1.r, b[i__1].i = q__1.i;
    if (*n == 1) {
	goto L80;
    }
    i__1 = nm1;
    i__2 = nm1;
    i__3 = nm1;
    i__4 = *n;
    q__3.r = d__[i__3].r * b[i__4].r - d__[i__3].i * b[i__4].i, q__3.i = d__[
	    i__3].r * b[i__4].i + d__[i__3].i * b[i__4].r;
    q__2.r = b[i__2].r - q__3.r, q__2.i = b[i__2].i - q__3.i;
    c_div(&q__1, &q__2, &c__[nm1]);
    b[i__1].r = q__1.r, b[i__1].i = q__1.i;
    if (nm2 < 1) {
	goto L70;
    }
    i__1 = nm2;
    for (kb = 1; kb <= i__1; ++kb) {
	k = nm2 - kb + 1;
	i__2 = k;
	i__3 = k;
	i__4 = k;
	i__5 = k + 1;
	q__4.r = d__[i__4].r * b[i__5].r - d__[i__4].i * b[i__5].i, q__4.i = 
		d__[i__4].r * b[i__5].i + d__[i__4].i * b[i__5].r;
	q__3.r = b[i__3].r - q__4.r, q__3.i = b[i__3].i - q__4.i;
	i__6 = k;
	i__7 = k + 2;
	q__5.r = e[i__6].r * b[i__7].r - e[i__6].i * b[i__7].i, q__5.i = e[
		i__6].r * b[i__7].i + e[i__6].i * b[i__7].r;
	q__2.r = q__3.r - q__5.r, q__2.i = q__3.i - q__5.i;
	c_div(&q__1, &q__2, &c__[k]);
	b[i__2].r = q__1.r, b[i__2].i = q__1.i;
/* L60: */
    }
L70:
L80:
L90:
L100:

    return 0;
} /* cgtsl_ */

