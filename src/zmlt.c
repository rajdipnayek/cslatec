/* zmlt.f -- translated by f2c (version 12.02.01).
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

/* DECK ZMLT */
/* Subroutine */ int zmlt_(doublereal *ar, doublereal *ai, doublereal *br, 
	doublereal *bi, doublereal *cr, doublereal *ci)
{
    static doublereal ca, cb;

/* ***BEGIN PROLOGUE  ZMLT */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to ZBESH, ZBESI, ZBESJ, ZBESK, ZBESY, ZAIRY and */
/*            ZBIRY */
/* ***LIBRARY   SLATEC */
/* ***TYPE      ALL (ZMLT-A) */
/* ***AUTHOR  Amos, D. E., (SNL) */
/* ***DESCRIPTION */

/*     DOUBLE PRECISION COMPLEX MULTIPLY, C=A*B. */

/* ***SEE ALSO  ZAIRY, ZBESH, ZBESI, ZBESJ, ZBESK, ZBESY, ZBIRY */
/* ***ROUTINES CALLED  (NONE) */
/* ***REVISION HISTORY  (YYMMDD) */
/*   830501  DATE WRITTEN */
/*   910415  Prologue converted to Version 4.0 format.  (BAB) */
/* ***END PROLOGUE  ZMLT */
/* ***FIRST EXECUTABLE STATEMENT  ZMLT */
    ca = *ar * *br - *ai * *bi;
    cb = *ar * *bi + *ai * *br;
    *cr = ca;
    *ci = cb;
    return 0;
} /* zmlt_ */

