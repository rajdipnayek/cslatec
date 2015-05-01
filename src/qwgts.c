/* qwgts.f -- translated by f2c (version 12.02.01).
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

/* DECK QWGTS */
doublereal qwgts_(real *x, real *a, real *b, real *alfa, real *beta, integer *
	integr)
{
    /* System generated locals */
    real ret_val;
    doublereal d__1, d__2, d__3, d__4;

    /* Local variables */
    static real xma, bmx;

/* ***BEGIN PROLOGUE  QWGTS */
/* ***SUBSIDIARY */
/* ***PURPOSE  This function subprogram is used together with the */
/*            routine QAWS and defines the WEIGHT function. */
/* ***LIBRARY   SLATEC */
/* ***TYPE      SINGLE PRECISION (QWGTS-S, DQWGTS-D) */
/* ***KEYWORDS  ALGEBRAICO-LOGARITHMIC, END POINT SINGULARITIES, */
/*             WEIGHT FUNCTION */
/* ***AUTHOR  Piessens, Robert */
/*             Applied Mathematics and Programming Division */
/*             K. U. Leuven */
/*           de Doncker, Elise */
/*             Applied Mathematics and Programming Division */
/*             K. U. Leuven */
/* ***SEE ALSO  QK15W */
/* ***ROUTINES CALLED  (NONE) */
/* ***REVISION HISTORY  (YYMMDD) */
/*   810101  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890531  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900328  Added TYPE section.  (WRB) */
/* ***END PROLOGUE  QWGTS */

/* ***FIRST EXECUTABLE STATEMENT  QWGTS */
    xma = *x - *a;
    bmx = *b - *x;
    d__1 = (doublereal) xma;
    d__2 = (doublereal) (*alfa);
    d__3 = (doublereal) bmx;
    d__4 = (doublereal) (*beta);
    ret_val = pow_dd(&d__1, &d__2) * pow_dd(&d__3, &d__4);
    switch (*integr) {
	case 1:  goto L40;
	case 2:  goto L10;
	case 3:  goto L20;
	case 4:  goto L30;
    }
L10:
    ret_val *= log(xma);
    goto L40;
L20:
    ret_val *= log(bmx);
    goto L40;
L30:
    ret_val = ret_val * log(xma) * log(bmx);
L40:
    return ret_val;
} /* qwgts_ */

