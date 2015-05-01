/* cbunk.f -- translated by f2c (version 12.02.01).
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

/* DECK CBUNK */
/* Subroutine */ int cbunk_(complex *z__, real *fnu, integer *kode, integer *
	mr, integer *n, complex *y, integer *nz, real *tol, real *elim, real *
	alim)
{
    /* Local variables */
    static real ax, ay, xx, yy;
    extern /* Subroutine */ int cunk1_(complex *, real *, integer *, integer *
	    , integer *, complex *, integer *, real *, real *, real *), 
	    cunk2_(complex *, real *, integer *, integer *, integer *, 
	    complex *, integer *, real *, real *, real *);

/* ***BEGIN PROLOGUE  CBUNK */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to CBESH and CBESK */
/* ***LIBRARY   SLATEC */
/* ***TYPE      ALL (CBUNK-A, ZBUNK-A) */
/* ***AUTHOR  Amos, D. E., (SNL) */
/* ***DESCRIPTION */

/*     CBUNK COMPUTES THE K BESSEL FUNCTION FOR FNU.GT.FNUL. */
/*     ACCORDING TO THE UNIFORM ASYMPTOTIC EXPANSION FOR K(FNU,Z) */
/*     IN CUNK1 AND THE EXPANSION FOR H(2,FNU,Z) IN CUNK2 */

/* ***SEE ALSO  CBESH, CBESK */
/* ***ROUTINES CALLED  CUNK1, CUNK2 */
/* ***REVISION HISTORY  (YYMMDD) */
/*   830501  DATE WRITTEN */
/*   910415  Prologue converted to Version 4.0 format.  (BAB) */
/* ***END PROLOGUE  CBUNK */
/* ***FIRST EXECUTABLE STATEMENT  CBUNK */
    /* Parameter adjustments */
    --y;

    /* Function Body */
    *nz = 0;
    xx = z__->r;
    yy = r_imag(z__);
    ax = dabs(xx) * 1.7321f;
    ay = dabs(yy);
    if (ay > ax) {
	goto L10;
    }
/* ----------------------------------------------------------------------- */
/*     ASYMPTOTIC EXPANSION FOR K(FNU,Z) FOR LARGE FNU APPLIED IN */
/*     -PI/3.LE.ARG(Z).LE.PI/3 */
/* ----------------------------------------------------------------------- */
    cunk1_(z__, fnu, kode, mr, n, &y[1], nz, tol, elim, alim);
    goto L20;
L10:
/* ----------------------------------------------------------------------- */
/*     ASYMPTOTIC EXPANSION FOR H(2,FNU,Z*EXP(M*HPI)) FOR LARGE FNU */
/*     APPLIED IN PI/3.LT.ABS(ARG(Z)).LE.PI/2 WHERE M=+I OR -I */
/*     AND HPI=PI/2 */
/* ----------------------------------------------------------------------- */
    cunk2_(z__, fnu, kode, mr, n, &y[1], nz, tol, elim, alim);
L20:
    return 0;
} /* cbunk_ */

