/* crotg.f -- translated by f2c (version 12.02.01).
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

/* DECK CROTG */
/* Subroutine */ int crotg_(complex *ca, complex *cb, real *c__, complex *s)
{
    /* System generated locals */
    real r__1, r__2;
    complex q__1, q__2, q__3;

    /* Local variables */
    static real norm;
    static complex alpha;
    static real scale;

/* ***BEGIN PROLOGUE  CROTG */
/* ***PURPOSE  Construct a Givens transformation. */
/* ***LIBRARY   SLATEC (BLAS) */
/* ***CATEGORY  D1B10 */
/* ***TYPE      COMPLEX (SROTG-S, DROTG-D, CROTG-C) */
/* ***KEYWORDS  BLAS, GIVENS ROTATION, GIVENS TRANSFORMATION, */
/*             LINEAR ALGEBRA, VECTOR */
/* ***AUTHOR  (UNKNOWN) */
/* ***DESCRIPTION */

/*    Complex Givens transformation */

/*    Construct the Givens transformation */

/*             (C    S) */
/*       G  =  (      ),  C**2 + ABS(S)**2 =1, */
/*             (-S   C) */

/*    which zeros the second entry of the complex 2-vector (CA,CB)**T */

/*    The quantity CA/ABS(CA)*NORM(CA,CB) overwrites CA in storage. */

/*    Input: */
/*        CA (Complex) */
/*        CB (Complex) */

/*    Output: */
/*        CA (Complex)      CA/ABS(CA)*NORM(CA,CB) */
/*        C  (Real) */
/*        S  (Complex) */

/* ***REFERENCES  (NONE) */
/* ***ROUTINES CALLED  (NONE) */
/* ***REVISION HISTORY  (YYMMDD) */
/*   790101  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890531  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/* ***END PROLOGUE  CROTG */
/* ***FIRST EXECUTABLE STATEMENT  CROTG */
    if (c_abs(ca) == 0.f) {
	*c__ = 0.f;
	s->r = 1.f, s->i = 0.f;
	ca->r = cb->r, ca->i = cb->i;
    } else {
	scale = c_abs(ca) + c_abs(cb);
	q__1.r = ca->r / scale, q__1.i = ca->i / scale;
/* Computing 2nd power */
	r__1 = c_abs(&q__1);
	q__2.r = cb->r / scale, q__2.i = cb->i / scale;
/* Computing 2nd power */
	r__2 = c_abs(&q__2);
	norm = scale * sqrt(r__1 * r__1 + r__2 * r__2);
	r__1 = c_abs(ca);
	q__1.r = ca->r / r__1, q__1.i = ca->i / r__1;
	alpha.r = q__1.r, alpha.i = q__1.i;
	*c__ = c_abs(ca) / norm;
	r_cnjg(&q__3, cb);
	q__2.r = alpha.r * q__3.r - alpha.i * q__3.i, q__2.i = alpha.r * 
		q__3.i + alpha.i * q__3.r;
	q__1.r = q__2.r / norm, q__1.i = q__2.i / norm;
	s->r = q__1.r, s->i = q__1.i;
	q__1.r = norm * alpha.r, q__1.i = norm * alpha.i;
	ca->r = q__1.r, ca->i = q__1.i;
    }
    return 0;
} /* crotg_ */

