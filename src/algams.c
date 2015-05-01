/* algams.f -- translated by f2c (version 12.02.01).
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

static real c_b2 = 2.f;

/* DECK ALGAMS */
/* Subroutine */ int algams_(real *x, real *algam, real *sgngam)
{
    /* System generated locals */
    real r__1;

    /* Local variables */
    static integer int__;
    extern doublereal alngam_(real *);

/* ***BEGIN PROLOGUE  ALGAMS */
/* ***PURPOSE  Compute the logarithm of the absolute value of the Gamma */
/*            function. */
/* ***LIBRARY   SLATEC (FNLIB) */
/* ***CATEGORY  C7A */
/* ***TYPE      SINGLE PRECISION (ALGAMS-S, DLGAMS-D) */
/* ***KEYWORDS  ABSOLUTE VALUE OF THE LOGARITHM OF THE GAMMA FUNCTION, */
/*             FNLIB, SPECIAL FUNCTIONS */
/* ***AUTHOR  Fullerton, W., (LANL) */
/* ***DESCRIPTION */

/* Evaluates the logarithm of the absolute value of the gamma */
/* function. */
/*     X           - input argument */
/*     ALGAM       - result */
/*     SGNGAM      - is set to the sign of GAMMA(X) and will */
/*                   be returned at +1.0 or -1.0. */

/* ***REFERENCES  (NONE) */
/* ***ROUTINES CALLED  ALNGAM */
/* ***REVISION HISTORY  (YYMMDD) */
/*   770701  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890531  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/* ***END PROLOGUE  ALGAMS */
/* ***FIRST EXECUTABLE STATEMENT  ALGAMS */
    *algam = alngam_(x);
    *sgngam = 1.f;
    if (*x > 0.f) {
	return 0;
    }

    r__1 = -r_int(x);
    int__ = r_mod(&r__1, &c_b2) + .1f;
    if (int__ == 0) {
	*sgngam = -1.f;
    }

    return 0;
} /* algams_ */

