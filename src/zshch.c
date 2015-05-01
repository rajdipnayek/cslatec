/* zshch.f -- translated by f2c (version 12.02.01).
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

/* DECK ZSHCH */
/* Subroutine */ int zshch_(doublereal *zr, doublereal *zi, doublereal *cshr, 
	doublereal *cshi, doublereal *cchr, doublereal *cchi)
{
    /* Local variables */
    static doublereal ch, cn, sh, sn;

/* ***BEGIN PROLOGUE  ZSHCH */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to ZBESH and ZBESK */
/* ***LIBRARY   SLATEC */
/* ***TYPE      ALL (CSHCH-A, ZSHCH-A) */
/* ***AUTHOR  Amos, D. E., (SNL) */
/* ***DESCRIPTION */

/*     ZSHCH COMPUTES THE COMPLEX HYPERBOLIC FUNCTIONS CSH=SINH(X+I*Y) */
/*     AND CCH=COSH(X+I*Y), WHERE I**2=-1. */

/* ***SEE ALSO  ZBESH, ZBESK */
/* ***ROUTINES CALLED  (NONE) */
/* ***REVISION HISTORY  (YYMMDD) */
/*   830501  DATE WRITTEN */
/*   910415  Prologue converted to Version 4.0 format.  (BAB) */
/* ***END PROLOGUE  ZSHCH */

/* ***FIRST EXECUTABLE STATEMENT  ZSHCH */
    sh = sinh(*zr);
    ch = cosh(*zr);
    sn = sin(*zi);
    cn = cos(*zi);
    *cshr = sh * cn;
    *cshi = ch * sn;
    *cchr = ch * cn;
    *cchi = sh * sn;
    return 0;
} /* zshch_ */

