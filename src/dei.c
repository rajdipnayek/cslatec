/* dei.f -- translated by f2c (version 12.02.01).
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

/* DECK DEI */
doublereal dei_(doublereal *x)
{
    /* System generated locals */
    doublereal ret_val, d__1;

    /* Local variables */
    extern doublereal de1_(doublereal *);

/* ***BEGIN PROLOGUE  DEI */
/* ***PURPOSE  Compute the exponential integral Ei(X). */
/* ***LIBRARY   SLATEC (FNLIB) */
/* ***CATEGORY  C5 */
/* ***TYPE      DOUBLE PRECISION (EI-S, DEI-D) */
/* ***KEYWORDS  EI FUNCTION, EXPONENTIAL INTEGRAL, FNLIB, */
/*             SPECIAL FUNCTIONS */
/* ***AUTHOR  Fullerton, W., (LANL) */
/* ***DESCRIPTION */

/* DEI calculates the double precision exponential integral, Ei(X), for */
/* positive double precision argument X and the Cauchy principal value */
/* for negative X.  If principal values are used everywhere, then, for */
/* all X, */

/*    Ei(X) = -E1(-X) */
/* or */
/*    E1(X) = -Ei(-X). */

/* ***REFERENCES  (NONE) */
/* ***ROUTINES CALLED  DE1 */
/* ***REVISION HISTORY  (YYMMDD) */
/*   770701  DATE WRITTEN */
/*   891115  Modified prologue description.  (WRB) */
/*   891115  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/* ***END PROLOGUE  DEI */
/* ***FIRST EXECUTABLE STATEMENT  DEI */
    d__1 = -(*x);
    ret_val = -de1_(&d__1);

    return ret_val;
} /* dei_ */

