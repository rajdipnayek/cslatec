/* cgamr.f -- translated by f2c (version 12.02.01).
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

/* DECK CGAMR */
/* Complex */ void cgamr_(complex * ret_val, complex *z__)
{
    /* System generated locals */
    complex q__1, q__2;

    /* Local variables */
    static real x;
    static integer irold;
    extern /* Subroutine */ int xgetf_(integer *), xsetf_(integer *);
    extern /* Complex */ void clngam_(complex *, complex *);
    extern /* Subroutine */ int xerclr_(void);

/* ***BEGIN PROLOGUE  CGAMR */
/* ***PURPOSE  Compute the reciprocal of the Gamma function. */
/* ***LIBRARY   SLATEC (FNLIB) */
/* ***CATEGORY  C7A */
/* ***TYPE      COMPLEX (GAMR-S, DGAMR-D, CGAMR-C) */
/* ***KEYWORDS  FNLIB, RECIPROCAL GAMMA FUNCTION, SPECIAL FUNCTIONS */
/* ***AUTHOR  Fullerton, W., (LANL) */
/* ***DESCRIPTION */

/* CGAMR(Z) calculates the reciprocal gamma function for COMPLEX */
/* argument Z.  This is a preliminary version that is not accurate. */

/* ***REFERENCES  (NONE) */
/* ***ROUTINES CALLED  CLNGAM, XERCLR, XGETF, XSETF */
/* ***REVISION HISTORY  (YYMMDD) */
/*   770701  DATE WRITTEN */
/*   861211  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/* ***END PROLOGUE  CGAMR */
/* ***FIRST EXECUTABLE STATEMENT  CGAMR */
     ret_val->r = 0.f,  ret_val->i = 0.f;
    x = z__->r;
    if (x <= 0.f && r_int(&x) == x && r_imag(z__) == 0.f) {
	return ;
    }

    xgetf_(&irold);
    xsetf_(&c__1);
    clngam_(&q__1, z__);
     ret_val->r = q__1.r,  ret_val->i = q__1.i;
    xerclr_();
    xsetf_(&irold);
    q__2.r = - ret_val->r, q__2.i = - ret_val->i;
    c_exp(&q__1, &q__2);
     ret_val->r = q__1.r,  ret_val->i = q__1.i;

    return ;
} /* cgamr_ */

