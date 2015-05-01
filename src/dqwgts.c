/* dqwgts.f -- translated by f2c (version 12.02.01).
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

/* DECK DQWGTS */
doublereal dqwgts_(doublereal *x, doublereal *a, doublereal *b, doublereal *
	alfa, doublereal *beta, integer *integr)
{
    /* System generated locals */
    doublereal ret_val;

    /* Local variables */
    static doublereal xma, bmx;

/* ***BEGIN PROLOGUE  DQWGTS */
/* ***SUBSIDIARY */
/* ***PURPOSE  This function subprogram is used together with the */
/*            routine DQAWS and defines the WEIGHT function. */
/* ***LIBRARY   SLATEC */
/* ***TYPE      DOUBLE PRECISION (QWGTS-S, DQWGTS-D) */
/* ***KEYWORDS  ALGEBRAICO-LOGARITHMIC, END POINT SINGULARITIES, */
/*             WEIGHT FUNCTION */
/* ***AUTHOR  Piessens, Robert */
/*             Applied Mathematics and Programming Division */
/*             K. U. Leuven */
/*           de Doncker, Elise */
/*             Applied Mathematics and Programming Division */
/*             K. U. Leuven */
/* ***SEE ALSO  DQK15W */
/* ***ROUTINES CALLED  (NONE) */
/* ***REVISION HISTORY  (YYMMDD) */
/*   810101  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890531  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900328  Added TYPE section.  (WRB) */
/* ***END PROLOGUE  DQWGTS */

/* ***FIRST EXECUTABLE STATEMENT  DQWGTS */
    xma = *x - *a;
    bmx = *b - *x;
    ret_val = pow_dd(&xma, alfa) * pow_dd(&bmx, beta);
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
} /* dqwgts_ */

