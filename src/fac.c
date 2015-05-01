/* fac.f -- translated by f2c (version 12.02.01).
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
static integer c__2 = 2;

/* DECK FAC */
doublereal fac_(integer *n)
{
    /* Initialized data */

    static real facn[26] = { 1.f,1.f,2.f,6.f,24.f,120.f,720.f,5040.f,40320.f,
	    362880.f,3628800.f,39916800.f,479001600.f,6227020800.f,
	    87178291200.f,1.307674368e12f,2.0922789888e13f,3.55687428096e14f,
	    6.402373705728e15f,1.21645100408832e17f,2.43290200817664e18f,
	    5.109094217170944e19f,1.1240007277776077e21f,
	    2.5852016738884977e22f,6.2044840173323944e23f,
	    1.5511210043330986e25f };
    static real sq2pil = .91893853320467274f;
    static integer nmax = 0;

    /* System generated locals */
    real ret_val;

    /* Local variables */
    static real x, xmin, xmax;
    extern doublereal r9lgmc_(real *);
    extern /* Subroutine */ int gamlim_(real *, real *), xermsg_(char *, char 
	    *, char *, integer *, integer *, ftnlen, ftnlen, ftnlen);

/* ***BEGIN PROLOGUE  FAC */
/* ***PURPOSE  Compute the factorial function. */
/* ***LIBRARY   SLATEC (FNLIB) */
/* ***CATEGORY  C1 */
/* ***TYPE      SINGLE PRECISION (FAC-S, DFAC-D) */
/* ***KEYWORDS  FACTORIAL, FNLIB, SPECIAL FUNCTIONS */
/* ***AUTHOR  Fullerton, W., (LANL) */
/* ***DESCRIPTION */

/* FAC(N) evaluates the factorial function of N.  FAC is single */
/* precision.  N must be an integer between 0 and 25 inclusive. */

/* ***REFERENCES  (NONE) */
/* ***ROUTINES CALLED  GAMLIM, R9LGMC, XERMSG */
/* ***REVISION HISTORY  (YYMMDD) */
/*   770601  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890531  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ) */
/* ***END PROLOGUE  FAC */
/* ***FIRST EXECUTABLE STATEMENT  FAC */
    if (nmax != 0) {
	goto L10;
    }
    gamlim_(&xmin, &xmax);
    nmax = xmax - 1.f;

L10:
    if (*n < 0) {
	xermsg_("SLATEC", "FAC", "FACTORIAL OF NEGATIVE INTEGER UNDEFINED", &
		c__1, &c__2, (ftnlen)6, (ftnlen)3, (ftnlen)39);
    }

    if (*n <= 25) {
	ret_val = facn[*n];
    }
    if (*n <= 25) {
	return ret_val;
    }

    if (*n > nmax) {
	xermsg_("SLATEC", "FAC", "N SO BIG FACTORIAL(N) OVERFLOWS", &c__2, &
		c__2, (ftnlen)6, (ftnlen)3, (ftnlen)31);
    }

    x = (real) (*n + 1);
    ret_val = exp((x - .5f) * log(x) - x + sq2pil + r9lgmc_(&x));

    return ret_val;
} /* fac_ */

