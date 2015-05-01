/* cpzero.f -- translated by f2c (version 12.02.01).
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

static logical c_false = FALSE_;
static integer c__0 = 0;
static logical c_true = TRUE_;

/* DECK CPZERO */
/* Subroutine */ int cpzero_(integer *in, complex *a, complex *r__, complex *
	t, integer *iflg, real *s)
{
    /* System generated locals */
    integer i__1, i__2, i__3, i__4, i__5;
    real r__1, r__2, r__3, r__4, r__5, r__6;
    doublereal d__1, d__2;
    complex q__1, q__2, q__3;

    /* Local variables */
    static integer i__, j, n;
    static real u, v, x;
    static integer n1;
    static complex pn;
    static integer nr, nit, imax, nmax;
    static complex temp;
    extern /* Subroutine */ int cpevl_(integer *, integer *, complex *, 
	    complex *, complex *, complex *, logical *);

/* ***BEGIN PROLOGUE  CPZERO */
/* ***PURPOSE  Find the zeros of a polynomial with complex coefficients. */
/* ***LIBRARY   SLATEC */
/* ***CATEGORY  F1A1B */
/* ***TYPE      COMPLEX (RPZERO-S, CPZERO-C) */
/* ***KEYWORDS  POLYNOMIAL ROOTS, POLYNOMIAL ZEROS, REAL ROOTS */
/* ***AUTHOR  Kahaner, D. K., (NBS) */
/* ***DESCRIPTION */

/*      Find the zeros of the complex polynomial */
/*         P(Z)= A(1)*Z**N + A(2)*Z**(N-1) +...+ A(N+1) */

/*    Input... */
/*       IN = degree of P(Z) */
/*       A = complex vector containing coefficients of P(Z), */
/*            A(I) = coefficient of Z**(N+1-i) */
/*       R = N word complex vector containing initial estimates for zeros */
/*            if these are known. */
/*       T = 4(N+1) word array used for temporary storage */
/*       IFLG = flag to indicate if initial estimates of */
/*              zeros are input. */
/*            If IFLG .EQ. 0, no estimates are input. */
/*            If IFLG .NE. 0, the vector R contains estimates of */
/*               the zeros */
/*       ** WARNING ****** If estimates are input, they must */
/*                         be separated, that is, distinct or */
/*                         not repeated. */
/*       S = an N word array */

/*    Output... */
/*       R(I) = Ith zero, */
/*       S(I) = bound for R(I) . */
/*       IFLG = error diagnostic */
/*    Error Diagnostics... */
/*       If IFLG .EQ. 0 on return, all is well */
/*       If IFLG .EQ. 1 on return, A(1)=0.0 or N=0 on input */
/*       If IFLG .EQ. 2 on return, the program failed to converge */
/*                after 25*N iterations.  Best current estimates of the */
/*                zeros are in R(I).  Error bounds are not calculated. */

/* ***REFERENCES  (NONE) */
/* ***ROUTINES CALLED  CPEVL */
/* ***REVISION HISTORY  (YYMMDD) */
/*   810223  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890531  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/* ***END PROLOGUE  CPZERO */

/* ***FIRST EXECUTABLE STATEMENT  CPZERO */
    /* Parameter adjustments */
    --s;
    --t;
    --r__;
    --a;

    /* Function Body */
    if (*in <= 0 || c_abs(&a[1]) == 0.f) {
	goto L30;
    }

/*       CHECK FOR EASILY OBTAINED ZEROS */

    n = *in;
    n1 = n + 1;
    if (*iflg != 0) {
	goto L14;
    }
L1:
    n1 = n + 1;
    if (n > 1) {
	goto L2;
    }
    q__2.r = -a[2].r, q__2.i = -a[2].i;
    c_div(&q__1, &q__2, &a[1]);
    r__[1].r = q__1.r, r__[1].i = q__1.i;
    s[1] = 0.f;
    return 0;
L2:
    if (c_abs(&a[n1]) != 0.f) {
	goto L3;
    }
    i__1 = n;
    r__[i__1].r = 0.f, r__[i__1].i = 0.f;
    s[n] = 0.f;
    --n;
    goto L1;

/*          IF INITIAL ESTIMATES FOR ZEROS NOT GIVEN, FIND SOME */

L3:
    q__2.r = -a[2].r, q__2.i = -a[2].i;
    d__1 = (doublereal) n;
    q__3.r = d__1 * a[1].r, q__3.i = d__1 * a[1].i;
    c_div(&q__1, &q__2, &q__3);
    temp.r = q__1.r, temp.i = q__1.i;
    cpevl_(&n, &n, &a[1], &temp, &t[1], &t[1], &c_false);
    imax = n + 2;
    i__1 = n1;
    r__1 = c_abs(&t[n1]);
    t[i__1].r = r__1, t[i__1].i = 0.f;
    i__1 = n1;
    for (i__ = 2; i__ <= i__1; ++i__) {
	i__2 = n + i__;
	r__1 = -c_abs(&t[n + 2 - i__]);
	t[i__2].r = r__1, t[i__2].i = 0.f;
	i__2 = n + i__;
	i__3 = imax;
	if (t[i__2].r < t[i__3].r) {
	    imax = n + i__;
	}
/* L6: */
    }
    i__1 = imax;
    i__2 = n1;
    d__1 = (doublereal) (-t[i__1].r / t[i__2].r);
    d__2 = (doublereal) (1.f / (imax - n1));
    x = pow_dd(&d__1, &d__2);
L7:
    x *= 2.f;
    q__1.r = x, q__1.i = 0.f;
    cpevl_(&n, &c__0, &t[n1], &q__1, &pn, &pn, &c_false);
    if (pn.r < 0.f) {
	goto L7;
    }
    u = x * .5f;
    v = x;
L10:
    x = (u + v) * .5f;
    q__1.r = x, q__1.i = 0.f;
    cpevl_(&n, &c__0, &t[n1], &q__1, &pn, &pn, &c_false);
    if (pn.r > 0.f) {
	v = x;
    }
    if (pn.r <= 0.f) {
	u = x;
    }
    if (v - u > (v + 1.f) * .001f) {
	goto L10;
    }
    i__1 = n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	u = 3.14159265f / n * ((i__ << 1) - 1.5f);
/* L13: */
	i__2 = i__;
/* Computing MAX */
	r__2 = x, r__3 = c_abs(&temp) * .001f;
	r__1 = dmax(r__2,r__3);
	r__4 = cos(u);
	r__5 = sin(u);
	q__3.r = r__4, q__3.i = r__5;
	q__2.r = r__1 * q__3.r, q__2.i = r__1 * q__3.i;
	q__1.r = q__2.r + temp.r, q__1.i = q__2.i + temp.i;
	r__[i__2].r = q__1.r, r__[i__2].i = q__1.i;
    }

/*          MAIN ITERATION LOOP STARTS HERE */

L14:
    nr = 0;
    nmax = n * 25;
    i__2 = nmax;
    for (nit = 1; nit <= i__2; ++nit) {
	i__1 = n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    if (nit != 1 && c_abs(&t[i__]) == 0.f) {
		goto L18;
	    }
	    cpevl_(&n, &c__0, &a[1], &r__[i__], &pn, &temp, &c_true);
	    if ((r__1 = pn.r, dabs(r__1)) + (r__2 = r_imag(&pn), dabs(r__2)) 
		    > temp.r + r_imag(&temp)) {
		goto L16;
	    }
	    i__3 = i__;
	    t[i__3].r = 0.f, t[i__3].i = 0.f;
	    ++nr;
	    goto L18;
L16:
	    temp.r = a[1].r, temp.i = a[1].i;
	    i__3 = n;
	    for (j = 1; j <= i__3; ++j) {
/* L17: */
		if (j != i__) {
		    i__4 = i__;
		    i__5 = j;
		    q__2.r = r__[i__4].r - r__[i__5].r, q__2.i = r__[i__4].i 
			    - r__[i__5].i;
		    q__1.r = temp.r * q__2.r - temp.i * q__2.i, q__1.i = 
			    temp.r * q__2.i + temp.i * q__2.r;
		    temp.r = q__1.r, temp.i = q__1.i;
		}
	    }
	    i__4 = i__;
	    c_div(&q__1, &pn, &temp);
	    t[i__4].r = q__1.r, t[i__4].i = q__1.i;
L18:
	    ;
	}
	i__1 = n;
	for (i__ = 1; i__ <= i__1; ++i__) {
/* L15: */
	    i__4 = i__;
	    i__5 = i__;
	    i__3 = i__;
	    q__1.r = r__[i__5].r - t[i__3].r, q__1.i = r__[i__5].i - t[i__3]
		    .i;
	    r__[i__4].r = q__1.r, r__[i__4].i = q__1.i;
	}
	if (nr == n) {
	    goto L21;
	}
/* L19: */
    }
    goto L26;

/*          CALCULATE ERROR BOUNDS FOR ZEROS */

L21:
    i__2 = n;
    for (nr = 1; nr <= i__2; ++nr) {
	cpevl_(&n, &n, &a[1], &r__[nr], &t[1], &t[n + 2], &c_true);
	r__3 = (r__1 = t[1].r, dabs(r__1));
	r__4 = (r__2 = r_imag(&t[1]), dabs(r__2));
	q__2.r = r__3, q__2.i = r__4;
	i__4 = n + 2;
	q__1.r = q__2.r + t[i__4].r, q__1.i = q__2.i + t[i__4].i;
	x = c_abs(&q__1);
	s[nr] = 0.f;
	i__4 = n;
	for (i__ = 1; i__ <= i__4; ++i__) {
	    x = x * (real) (n1 - i__) / i__;
/* Computing MAX */
	    i__5 = i__ + 1;
	    i__3 = n1 + i__;
	    r__4 = (r__1 = t[i__5].r, dabs(r__1)) - t[i__3].r;
	    r__3 = dmax(r__4,0.f);
/* Computing MAX */
	    r__6 = (r__2 = r_imag(&t[i__ + 1]), dabs(r__2)) - r_imag(&t[n1 + 
		    i__]);
	    r__5 = dmax(r__6,0.f);
	    q__1.r = r__3, q__1.i = r__5;
	    temp.r = q__1.r, temp.i = q__1.i;
/* L23: */
/* Computing MAX */
	    d__1 = (doublereal) (c_abs(&temp) / x);
	    d__2 = (doublereal) (1.f / i__);
	    r__1 = s[nr], r__2 = pow_dd(&d__1, &d__2);
	    s[nr] = dmax(r__1,r__2);
	}
/* L25: */
	s[nr] = 1.f / s[nr];
    }
    return 0;
/*        ERROR EXITS */
L26:
    *iflg = 2;
    return 0;
L30:
    *iflg = 1;
    return 0;
} /* cpzero_ */

