/* dcosdg.f -- translated by f2c (version 12.02.01).
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

static doublereal c_b2 = 90.;
static doublereal c_b3 = 1.;

/* DECK DCOSDG */
doublereal dcosdg_(doublereal *x)
{
    /* Initialized data */

    static doublereal raddeg = .017453292519943295769236907684886;

    /* System generated locals */
    doublereal ret_val;

    /* Local variables */
    static integer n;

/* ***BEGIN PROLOGUE  DCOSDG */
/* ***PURPOSE  Compute the cosine of an argument in degrees. */
/* ***LIBRARY   SLATEC (FNLIB) */
/* ***CATEGORY  C4A */
/* ***TYPE      DOUBLE PRECISION (COSDG-S, DCOSDG-D) */
/* ***KEYWORDS  COSINE, DEGREES, ELEMENTARY FUNCTIONS, FNLIB, */
/*             TRIGONOMETRIC */
/* ***AUTHOR  Fullerton, W., (LANL) */
/* ***DESCRIPTION */

/* DCOSDG(X) calculates the double precision trigonometric cosine */
/* for double precision argument X in units of degrees. */

/* ***REFERENCES  (NONE) */
/* ***ROUTINES CALLED  (NONE) */
/* ***REVISION HISTORY  (YYMMDD) */
/*   770601  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890531  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/* ***END PROLOGUE  DCOSDG */
/* ***FIRST EXECUTABLE STATEMENT  DCOSDG */
    ret_val = cos(raddeg * *x);

    if (d_mod(x, &c_b2) != 0.) {
	return ret_val;
    }
    n = (integer) (abs(*x) / 90. + .5);
    n %= 2;
    if (n == 0) {
	ret_val = d_sign(&c_b3, &ret_val);
    }
    if (n == 1) {
	ret_val = 0.;
    }

    return ret_val;
} /* dcosdg_ */

