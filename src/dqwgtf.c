/* dqwgtf.f -- translated by f2c (version 12.02.01).
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

/* DECK DQWGTF */
doublereal dqwgtf_(doublereal *x, doublereal *omega, doublereal *p2, 
	doublereal *p3, doublereal *p4, integer *integr)
{
    /* System generated locals */
    doublereal ret_val;

    /* Local variables */
    static doublereal omx;

/* ***BEGIN PROLOGUE  DQWGTF */
/* ***SUBSIDIARY */
/* ***PURPOSE  This function subprogram is used together with the */
/*            routine DQAWF and defines the WEIGHT function. */
/* ***LIBRARY   SLATEC */
/* ***TYPE      DOUBLE PRECISION (QWGTF-S, DQWGTF-D) */
/* ***KEYWORDS  COS OR SIN IN WEIGHT FUNCTION */
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
/* ***END PROLOGUE  DQWGTF */

/* ***FIRST EXECUTABLE STATEMENT  DQWGTF */
    omx = *omega * *x;
    switch (*integr) {
	case 1:  goto L10;
	case 2:  goto L20;
    }
L10:
    ret_val = cos(omx);
    goto L30;
L20:
    ret_val = sin(omx);
L30:
    return ret_val;
} /* dqwgtf_ */

