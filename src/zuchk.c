/* zuchk.f -- translated by f2c (version 12.02.01).
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

/* DECK ZUCHK */
/* Subroutine */ int zuchk_(doublereal *yr, doublereal *yi, integer *nz, 
	doublereal *ascle, doublereal *tol)
{
    static doublereal wi, ss, st, wr;

/* ***BEGIN PROLOGUE  ZUCHK */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to SERI, ZUOIK, ZUNK1, ZUNK2, ZUNI1, ZUNI2 and */
/*            ZKSCL */
/* ***LIBRARY   SLATEC */
/* ***TYPE      ALL (CUCHK-A, ZUCHK-A) */
/* ***AUTHOR  Amos, D. E., (SNL) */
/* ***DESCRIPTION */

/*      Y ENTERS AS A SCALED QUANTITY WHOSE MAGNITUDE IS GREATER THAN */
/*      EXP(-ALIM)=ASCLE=1.0E+3*D1MACH(1)/TOL. THE TEST IS MADE TO SEE */
/*      IF THE MAGNITUDE OF THE REAL OR IMAGINARY PART WOULD UNDERFLOW */
/*      WHEN Y IS SCALED (BY TOL) TO ITS PROPER VALUE. Y IS ACCEPTED */
/*      IF THE UNDERFLOW IS AT LEAST ONE PRECISION BELOW THE MAGNITUDE */
/*      OF THE LARGEST COMPONENT; OTHERWISE THE PHASE ANGLE DOES NOT HAVE */
/*      ABSOLUTE ACCURACY AND AN UNDERFLOW IS ASSUMED. */

/* ***SEE ALSO  SERI, ZKSCL, ZUNI1, ZUNI2, ZUNK1, ZUNK2, ZUOIK */
/* ***ROUTINES CALLED  (NONE) */
/* ***REVISION HISTORY  (YYMMDD) */
/*   ??????  DATE WRITTEN */
/*   910415  Prologue converted to Version 4.0 format.  (BAB) */
/* ***END PROLOGUE  ZUCHK */

/*     COMPLEX Y */
/* ***FIRST EXECUTABLE STATEMENT  ZUCHK */
    *nz = 0;
    wr = abs(*yr);
    wi = abs(*yi);
    st = min(wr,wi);
    if (st > *ascle) {
	return 0;
    }
    ss = max(wr,wi);
    st /= *tol;
    if (ss < st) {
	*nz = 1;
    }
    return 0;
} /* zuchk_ */

