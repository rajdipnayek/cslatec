/* cuchk.f -- translated by f2c (version 12.02.01).
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

/* DECK CUCHK */
/* Subroutine */ int cuchk_(complex *y, integer *nz, real *ascle, real *tol)
{
    /* Local variables */
    static real yi, ss, st, yr;

/* ***BEGIN PROLOGUE  CUCHK */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to SERI, CUOIK, CUNK1, CUNK2, CUNI1, CUNI2 and */
/*            CKSCL */
/* ***LIBRARY   SLATEC */
/* ***TYPE      ALL (CUCHK-A, ZUCHK-A) */
/* ***AUTHOR  Amos, D. E., (SNL) */
/* ***DESCRIPTION */

/*      Y ENTERS AS A SCALED QUANTITY WHOSE MAGNITUDE IS GREATER THAN */
/*      EXP(-ALIM)=ASCLE=1.0E+3*R1MACH(1)/TOL. THE TEST IS MADE TO SEE */
/*      IF THE MAGNITUDE OF THE REAL OR IMAGINARY PART WOULD UNDER FLOW */
/*      WHEN Y IS SCALED (BY TOL) TO ITS PROPER VALUE. Y IS ACCEPTED */
/*      IF THE UNDERFLOW IS AT LEAST ONE PRECISION BELOW THE MAGNITUDE */
/*      OF THE LARGEST COMPONENT; OTHERWISE THE PHASE ANGLE DOES NOT HAVE */
/*      ABSOLUTE ACCURACY AND AN UNDERFLOW IS ASSUMED. */

/* ***SEE ALSO  CKSCL, CUNI1, CUNI2, CUNK1, CUNK2, CUOIK, SERI */
/* ***ROUTINES CALLED  (NONE) */
/* ***REVISION HISTORY  (YYMMDD) */
/*   ??????  DATE WRITTEN */
/*   910415  Prologue converted to Version 4.0 format.  (BAB) */
/* ***END PROLOGUE  CUCHK */

/* ***FIRST EXECUTABLE STATEMENT  CUCHK */
    *nz = 0;
    yr = y->r;
    yi = r_imag(y);
    yr = dabs(yr);
    yi = dabs(yi);
    st = dmin(yr,yi);
    if (st > *ascle) {
	return 0;
    }
    ss = dmax(yr,yi);
    st /= *tol;
    if (ss < st) {
	*nz = 1;
    }
    return 0;
} /* cuchk_ */

