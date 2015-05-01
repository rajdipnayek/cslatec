/* dqwgtc.f -- translated by f2c (version 12.02.01).
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

/* DECK DQWGTC */
doublereal dqwgtc_(doublereal *x, doublereal *c__, doublereal *p2, doublereal 
	*p3, doublereal *p4, integer *kp)
{
    /* System generated locals */
    doublereal ret_val;

/* ***BEGIN PROLOGUE  DQWGTC */
/* ***SUBSIDIARY */
/* ***PURPOSE  This function subprogram is used together with the */
/*            routine DQAWC and defines the WEIGHT function. */
/* ***LIBRARY   SLATEC */
/* ***TYPE      DOUBLE PRECISION (QWGTC-S, DQWGTC-D) */
/* ***KEYWORDS  CAUCHY PRINCIPAL VALUE, WEIGHT FUNCTION */
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
/*   830518  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900328  Added TYPE section.  (WRB) */
/* ***END PROLOGUE  DQWGTC */

/* ***FIRST EXECUTABLE STATEMENT  DQWGTC */
    ret_val = 1. / (*x - *c__);
    return ret_val;
} /* dqwgtc_ */

