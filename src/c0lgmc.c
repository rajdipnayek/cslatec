/* c0lgmc.f -- translated by f2c (version 12.02.01).
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

static integer c__3 = 3;
static complex c_b3 = {1.f,0.f};
static integer c__2 = 2;

/* DECK C0LGMC */
/* Complex */ void c0lgmc_(complex * ret_val, complex *z__)
{
    /* Initialized data */

    static real rbig = 0.f;

    /* System generated locals */
    complex q__1, q__2, q__3, q__4, q__5, q__6, q__7;

    /* Local variables */
    static complex q;
    static real cabsz;
    extern /* Complex */ void c9ln2r_(complex *, complex *);
    extern doublereal r1mach_(integer *);

/* ***BEGIN PROLOGUE  C0LGMC */
/* ***PURPOSE  Evaluate (Z+0.5)*LOG((Z+1.)/Z) - 1.0 with relative */
/*            accuracy. */
/* ***LIBRARY   SLATEC (FNLIB) */
/* ***CATEGORY  C7A */
/* ***TYPE      COMPLEX (C0LGMC-C) */
/* ***KEYWORDS  FNLIB, GAMMA FUNCTION, SPECIAL FUNCTIONS */
/* ***AUTHOR  Fullerton, W., (LANL) */
/* ***DESCRIPTION */

/* Evaluate  (Z+0.5)*LOG((Z+1.0)/Z) - 1.0  with relative error accuracy */
/* Let Q = 1.0/Z so that */
/*     (Z+0.5)*LOG(1+1/Z) - 1 = (Z+0.5)*(LOG(1+Q) - Q + Q*Q/2) - Q*Q/4 */
/*        = (Z+0.5)*Q**3*C9LN2R(Q) - Q**2/4, */
/* where  C9LN2R  is (LOG(1+Q) - Q + 0.5*Q**2) / Q**3. */

/* ***REFERENCES  (NONE) */
/* ***ROUTINES CALLED  C9LN2R, R1MACH */
/* ***REVISION HISTORY  (YYMMDD) */
/*   780401  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890531  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/* ***END PROLOGUE  C0LGMC */
/* ***FIRST EXECUTABLE STATEMENT  C0LGMC */
    if (rbig == 0.f) {
	rbig = 1.f / r1mach_(&c__3);
    }

    cabsz = c_abs(z__);
    if (cabsz > rbig) {
	q__4.r = z__->r + .5f, q__4.i = z__->i;
	q__3.r = -q__4.r, q__3.i = -q__4.i;
	c_log(&q__5, z__);
	q__2.r = q__3.r * q__5.r - q__3.i * q__5.i, q__2.i = q__3.r * q__5.i 
		+ q__3.i * q__5.r;
	q__1.r = q__2.r - z__->r, q__1.i = q__2.i - z__->i;
	 ret_val->r = q__1.r,  ret_val->i = q__1.i;
    }
    if (cabsz > rbig) {
	return ;
    }

    c_div(&q__1, &c_b3, z__);
    q.r = q__1.r, q.i = q__1.i;
    if (cabsz <= 1.23f) {
	q__3.r = z__->r + .5f, q__3.i = z__->i;
	q__5.r = q.r + 1.f, q__5.i = q.i;
	c_log(&q__4, &q__5);
	q__2.r = q__3.r * q__4.r - q__3.i * q__4.i, q__2.i = q__3.r * q__4.i 
		+ q__3.i * q__4.r;
	q__1.r = q__2.r - 1.f, q__1.i = q__2.i;
	 ret_val->r = q__1.r,  ret_val->i = q__1.i;
    }
    if (cabsz > 1.23f) {
	q__5.r = q.r * .5f, q__5.i = q.i * .5f;
	q__4.r = q__5.r + 1.f, q__4.i = q__5.i;
	c9ln2r_(&q__6, &q);
	q__3.r = q__4.r * q__6.r - q__4.i * q__6.i, q__3.i = q__4.r * q__6.i 
		+ q__4.i * q__6.r;
	q__2.r = q__3.r - .25f, q__2.i = q__3.i;
	pow_ci(&q__7, &q, &c__2);
	q__1.r = q__2.r * q__7.r - q__2.i * q__7.i, q__1.i = q__2.r * q__7.i 
		+ q__2.i * q__7.r;
	 ret_val->r = q__1.r,  ret_val->i = q__1.i;
    }

    return ;
} /* c0lgmc_ */

