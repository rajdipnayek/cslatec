/* zabs.f -- translated by f2c (version 12.02.01).
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

/* DECK ZABS */
doublereal zabs_(doublereal *zr, doublereal *zi)
{
    /* System generated locals */
    doublereal ret_val;

    /* Local variables */
    static doublereal q, s, u, v;

/* ***BEGIN PROLOGUE  ZABS */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to ZBESH, ZBESI, ZBESJ, ZBESK, ZBESY, ZAIRY and */
/*            ZBIRY */
/* ***LIBRARY   SLATEC */
/* ***TYPE      ALL (ZABS-A) */
/* ***AUTHOR  Amos, D. E., (SNL) */
/* ***DESCRIPTION */

/*     ZABS COMPUTES THE ABSOLUTE VALUE OR MAGNITUDE OF A DOUBLE */
/*     PRECISION COMPLEX VARIABLE CMPLX(ZR,ZI) */

/* ***SEE ALSO  ZAIRY, ZBESH, ZBESI, ZBESJ, ZBESK, ZBESY, ZBIRY */
/* ***ROUTINES CALLED  (NONE) */
/* ***REVISION HISTORY  (YYMMDD) */
/*   830501  DATE WRITTEN */
/*   910415  Prologue converted to Version 4.0 format.  (BAB) */
/* ***END PROLOGUE  ZABS */
/* ***FIRST EXECUTABLE STATEMENT  ZABS */
    u = abs(*zr);
    v = abs(*zi);
    s = u + v;
/* ----------------------------------------------------------------------- */
/*     S*1.0D0 MAKES AN UNNORMALIZED UNDERFLOW ON CDC MACHINES INTO A */
/*     TRUE FLOATING ZERO */
/* ----------------------------------------------------------------------- */
    s *= 1.;
    if (s == 0.) {
	goto L20;
    }
    if (u > v) {
	goto L10;
    }
    q = u / v;
    ret_val = v * sqrt(q * q + 1.);
    return ret_val;
L10:
    q = v / u;
    ret_val = u * sqrt(q * q + 1.);
    return ret_val;
L20:
    ret_val = 0.;
    return ret_val;
} /* zabs_ */

