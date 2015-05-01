/* dfac.f -- translated by f2c (version 12.02.01).
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

/* DECK DFAC */
doublereal dfac_(integer *n)
{
    /* Initialized data */

    static doublereal facn[31] = { 1.,1.,2.,6.,24.,120.,720.,5040.,40320.,
	    362880.,3628800.,39916800.,479001600.,6227020800.,87178291200.,
	    1.307674368e12,2.0922789888e13,3.55687428096e14,6.402373705728e15,
	    1.21645100408832e17,2.43290200817664e18,5.109094217170944e19,
	    1.12400072777760768e21,2.585201673888497664e22,
	    6.2044840173323943936e23,1.5511210043330985984e25,
	    4.03291461126605635584e26,1.0888869450418352160768e28,
	    3.04888344611713860501504e29,8.841761993739701954543616e30,
	    2.6525285981219105863630848e32 };
    static doublereal sq2pil = .91893853320467274178032973640562;
    static integer nmax = 0;

    /* System generated locals */
    doublereal ret_val;

    /* Local variables */
    static doublereal x, xmin, xmax;
    extern doublereal d9lgmc_(doublereal *);
    extern /* Subroutine */ int dgamlm_(doublereal *, doublereal *), xermsg_(
	    char *, char *, char *, integer *, integer *, ftnlen, ftnlen, 
	    ftnlen);

/* ***BEGIN PROLOGUE  DFAC */
/* ***PURPOSE  Compute the factorial function. */
/* ***LIBRARY   SLATEC (FNLIB) */
/* ***CATEGORY  C1 */
/* ***TYPE      DOUBLE PRECISION (FAC-S, DFAC-D) */
/* ***KEYWORDS  FACTORIAL, FNLIB, SPECIAL FUNCTIONS */
/* ***AUTHOR  Fullerton, W., (LANL) */
/* ***DESCRIPTION */

/* DFAC(N) calculates the double precision factorial for integer */
/* argument N. */

/* ***REFERENCES  (NONE) */
/* ***ROUTINES CALLED  D9LGMC, DGAMLM, XERMSG */
/* ***REVISION HISTORY  (YYMMDD) */
/*   770601  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890531  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ) */
/* ***END PROLOGUE  DFAC */
/* ***FIRST EXECUTABLE STATEMENT  DFAC */
    if (nmax != 0) {
	goto L10;
    }
    dgamlm_(&xmin, &xmax);
    nmax = (integer) (xmax - 1.);

L10:
    if (*n < 0) {
	xermsg_("SLATEC", "DFAC", "FACTORIAL OF NEGATIVE INTEGER UNDEFINED", &
		c__1, &c__2, (ftnlen)6, (ftnlen)4, (ftnlen)39);
    }

    if (*n <= 30) {
	ret_val = facn[*n];
    }
    if (*n <= 30) {
	return ret_val;
    }

    if (*n > nmax) {
	xermsg_("SLATEC", "DFAC", "N SO BIG FACTORIAL(N) OVERFLOWS", &c__2, &
		c__2, (ftnlen)6, (ftnlen)4, (ftnlen)31);
    }

    x = (doublereal) (*n + 1);
    ret_val = exp((x - .5) * log(x) - x + sq2pil + d9lgmc_(&x));

    return ret_val;
} /* dfac_ */

