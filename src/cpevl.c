/* cpevl.f -- translated by f2c (version 12.02.01).
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

static integer c__10 = 10;
static integer c__11 = 11;

/* DECK CPEVL */
/* Subroutine */ int cpevl_(integer *n, integer *m, complex *a, complex *z__, 
	complex *c__, complex *b, logical *kbd)
{
    /* Initialized data */

    static real d1 = 0.f;

    /* System generated locals */
    integer i__1, i__2, i__3;
    real r__1, r__2, r__3, r__4, r__5, r__6;
    complex q__1, q__2, q__3, q__4, q__5, q__6;

    /* Local variables */
    static integer i__, j;
    static real r__, s;
    static complex t, bi, ci;
    static integer np1;
    static complex bim1, cim1;
    static integer mini;
    extern integer i1mach_(integer *);

/* ***BEGIN PROLOGUE  CPEVL */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to CPZERO */
/* ***LIBRARY   SLATEC */
/* ***TYPE      SINGLE PRECISION (CPEVL-S) */
/* ***AUTHOR  (UNKNOWN) */
/* ***DESCRIPTION */

/*        Evaluate a complex polynomial and its derivatives. */
/*        Optionally compute error bounds for these values. */

/*   INPUT... */
/*        N = Degree of the polynomial */
/*        M = Number of derivatives to be calculated, */
/*            M=0 evaluates only the function */
/*            M=1 evaluates the function and first derivative, etc. */
/*             if M .GT. N+1 function and all N derivatives will be */
/*                calculated. */
/*       A = Complex vector containing the N+1 coefficients of polynomial */
/*               A(I)= coefficient of Z**(N+1-I) */
/*        Z = Complex point at which the evaluation is to take place. */
/*        C = Array of 2(M+1) words into which values are placed. */
/*        B = Array of 2(M+1) words only needed if bounds are to be */
/*              calculated.  It is not used otherwise. */
/*        KBD = A logical variable, e.g. .TRUE. or .FALSE. which is */
/*              to be set .TRUE. if bounds are to be computed. */

/*  OUTPUT... */
/*        C =  C(I+1) contains the complex value of the I-th */
/*              derivative at Z, I=0,...,M */
/*        B =  B(I) contains the bounds on the real and imaginary parts */
/*              of C(I) if they were requested. */

/* ***SEE ALSO  CPZERO */
/* ***ROUTINES CALLED  I1MACH */
/* ***REVISION HISTORY  (YYMMDD) */
/*   810223  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890831  Modified array declarations.  (WRB) */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900402  Added TYPE section.  (WRB) */
/* ***END PROLOGUE  CPEVL */

    /* Parameter adjustments */
    --b;
    --c__;
    --a;

    /* Function Body */
/* ***FIRST EXECUTABLE STATEMENT  CPEVL */
    if (d1 == 0.f) {
	r__1 = (real) i1mach_(&c__10);
	i__1 = 1 - i1mach_(&c__11);
	d1 = pow_ri(&r__1, &i__1);
    }
    np1 = *n + 1;
    i__1 = np1;
    for (j = 1; j <= i__1; ++j) {
	ci.r = 0.f, ci.i = 0.f;
	i__2 = j;
	cim1.r = a[i__2].r, cim1.i = a[i__2].i;
	bi.r = 0.f, bi.i = 0.f;
	bim1.r = 0.f, bim1.i = 0.f;
/* Computing MIN */
	i__2 = *m + 1, i__3 = *n + 2 - j;
	mini = min(i__2,i__3);
	i__2 = mini;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    if (j != 1) {
		i__3 = i__;
		ci.r = c__[i__3].r, ci.i = c__[i__3].i;
	    }
	    if (i__ != 1) {
		i__3 = i__ - 1;
		cim1.r = c__[i__3].r, cim1.i = c__[i__3].i;
	    }
	    i__3 = i__;
	    q__2.r = z__->r * ci.r - z__->i * ci.i, q__2.i = z__->r * ci.i + 
		    z__->i * ci.r;
	    q__1.r = cim1.r + q__2.r, q__1.i = cim1.i + q__2.i;
	    c__[i__3].r = q__1.r, c__[i__3].i = q__1.i;
	    if (! (*kbd)) {
		goto L1;
	    }
	    if (j != 1) {
		i__3 = i__;
		bi.r = b[i__3].r, bi.i = b[i__3].i;
	    }
	    if (i__ != 1) {
		i__3 = i__ - 1;
		bim1.r = b[i__3].r, bim1.i = b[i__3].i;
	    }
	    r__3 = d1 * 3.f + d1 * 4.f * d1;
	    r__4 = (r__1 = ci.r, dabs(r__1));
	    r__5 = (r__2 = r_imag(&ci), dabs(r__2));
	    q__3.r = r__4, q__3.i = r__5;
	    q__2.r = r__3 * q__3.r, q__2.i = r__3 * q__3.i;
	    q__1.r = bi.r + q__2.r, q__1.i = bi.i + q__2.i;
	    t.r = q__1.r, t.i = q__1.i;
	    r__3 = (r__1 = z__->r, dabs(r__1));
	    r__4 = (r__2 = r_imag(z__), dabs(r__2));
	    q__2.r = r__3, q__2.i = r__4;
	    r__5 = t.r;
	    r__6 = -r_imag(&t);
	    q__3.r = r__5, q__3.i = r__6;
	    q__1.r = q__2.r * q__3.r - q__2.i * q__3.i, q__1.i = q__2.r * 
		    q__3.i + q__2.i * q__3.r;
	    r__ = q__1.r;
	    r__3 = (r__1 = z__->r, dabs(r__1));
	    r__4 = (r__2 = r_imag(z__), dabs(r__2));
	    q__2.r = r__3, q__2.i = r__4;
	    q__1.r = q__2.r * t.r - q__2.i * t.i, q__1.i = q__2.r * t.i + 
		    q__2.i * t.r;
	    s = r_imag(&q__1);
	    i__3 = i__;
	    r__3 = d1 * 8.f + 1.f;
	    r__4 = (r__1 = cim1.r, dabs(r__1));
	    r__5 = (r__2 = r_imag(&cim1), dabs(r__2));
	    q__5.r = r__4, q__5.i = r__5;
	    q__4.r = d1 * q__5.r, q__4.i = d1 * q__5.i;
	    q__3.r = bim1.r + q__4.r, q__3.i = bim1.i + q__4.i;
	    q__6.r = r__, q__6.i = s;
	    q__2.r = q__3.r + q__6.r, q__2.i = q__3.i + q__6.i;
	    q__1.r = r__3 * q__2.r, q__1.i = r__3 * q__2.i;
	    b[i__3].r = q__1.r, b[i__3].i = q__1.i;
	    if (j == 1) {
		i__3 = i__;
		b[i__3].r = 0.f, b[i__3].i = 0.f;
	    }
L1:
	    ;
	}
    }
    return 0;
} /* cpevl_ */

