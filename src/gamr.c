/* gamr.f -- translated by f2c (version 12.02.01).
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

/* DECK GAMR */
doublereal gamr_(real *x)
{
    /* System generated locals */
    real ret_val;

    /* Local variables */
    extern doublereal gamma_(real *);
    static integer irold;
    static real alngx;
    extern /* Subroutine */ int xgetf_(integer *);
    static real sgngx;
    extern /* Subroutine */ int xsetf_(integer *), algams_(real *, real *, 
	    real *), xerclr_(void);

/* ***BEGIN PROLOGUE  GAMR */
/* ***PURPOSE  Compute the reciprocal of the Gamma function. */
/* ***LIBRARY   SLATEC (FNLIB) */
/* ***CATEGORY  C7A */
/* ***TYPE      SINGLE PRECISION (GAMR-S, DGAMR-D, CGAMR-C) */
/* ***KEYWORDS  FNLIB, RECIPROCAL GAMMA FUNCTION, SPECIAL FUNCTIONS */
/* ***AUTHOR  Fullerton, W., (LANL) */
/* ***DESCRIPTION */

/* GAMR is a single precision function that evaluates the reciprocal */
/* of the gamma function for single precision argument X. */

/* ***REFERENCES  (NONE) */
/* ***ROUTINES CALLED  ALGAMS, GAMMA, XERCLR, XGETF, XSETF */
/* ***REVISION HISTORY  (YYMMDD) */
/*   770701  DATE WRITTEN */
/*   861211  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900727  Added EXTERNAL statement.  (WRB) */
/* ***END PROLOGUE  GAMR */
/* ***FIRST EXECUTABLE STATEMENT  GAMR */
    ret_val = 0.f;
    if (*x <= 0.f && r_int(x) == *x) {
	return ret_val;
    }

    xgetf_(&irold);
    xsetf_(&c__1);
    if (dabs(*x) > 10.f) {
	goto L10;
    }
    ret_val = 1.f / gamma_(x);
    xerclr_();
    xsetf_(&irold);
    return ret_val;

L10:
    algams_(x, &alngx, &sgngx);
    xerclr_();
    xsetf_(&irold);
    ret_val = sgngx * exp(-alngx);
    return ret_val;

} /* gamr_ */

