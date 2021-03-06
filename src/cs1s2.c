/* cs1s2.f -- translated by f2c (version 12.02.01).
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

/* DECK CS1S2 */
/* Subroutine */ int cs1s2_(complex *zr, complex *s1, complex *s2, integer *
	nz, real *ascle, real *alim, integer *iuf)
{
    /* Initialized data */

    static complex czero = {0.f,0.f};

    /* System generated locals */
    complex q__1, q__2, q__3;

    /* Local variables */
    static complex c1;
    static real aa, xx, as1, as2;
    static complex s1d;
    static real aln;

/* ***BEGIN PROLOGUE  CS1S2 */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to CAIRY and CBESK */
/* ***LIBRARY   SLATEC */
/* ***TYPE      ALL (CS1S2-A, ZS1S2-A) */
/* ***AUTHOR  Amos, D. E., (SNL) */
/* ***DESCRIPTION */

/*     CS1S2 TESTS FOR A POSSIBLE UNDERFLOW RESULTING FROM THE */
/*     ADDITION OF THE I AND K FUNCTIONS IN THE ANALYTIC CON- */
/*     TINUATION FORMULA WHERE S1=K FUNCTION AND S2=I FUNCTION. */
/*     ON KODE=1 THE I AND K FUNCTIONS ARE DIFFERENT ORDERS OF */
/*     MAGNITUDE, BUT FOR KODE=2 THEY CAN BE OF THE SAME ORDER */
/*     OF MAGNITUDE AND THE MAXIMUM MUST BE AT LEAST ONE */
/*     PRECISION ABOVE THE UNDERFLOW LIMIT. */

/* ***SEE ALSO  CAIRY, CBESK */
/* ***ROUTINES CALLED  (NONE) */
/* ***REVISION HISTORY  (YYMMDD) */
/*   830501  DATE WRITTEN */
/*   910415  Prologue converted to Version 4.0 format.  (BAB) */
/* ***END PROLOGUE  CS1S2 */
/* ***FIRST EXECUTABLE STATEMENT  CS1S2 */
    *nz = 0;
    as1 = c_abs(s1);
    as2 = c_abs(s2);
    aa = s1->r;
    aln = r_imag(s1);
    if (aa == 0.f && aln == 0.f) {
	goto L10;
    }
    if (as1 == 0.f) {
	goto L10;
    }
    xx = zr->r;
    aln = -xx - xx + log(as1);
    s1d.r = s1->r, s1d.i = s1->i;
    s1->r = czero.r, s1->i = czero.i;
    as1 = 0.f;
    if (aln < -(*alim)) {
	goto L10;
    }
    c_log(&q__3, &s1d);
    q__2.r = q__3.r - zr->r, q__2.i = q__3.i - zr->i;
    q__1.r = q__2.r - zr->r, q__1.i = q__2.i - zr->i;
    c1.r = q__1.r, c1.i = q__1.i;
    c_exp(&q__1, &c1);
    s1->r = q__1.r, s1->i = q__1.i;
    as1 = c_abs(s1);
    ++(*iuf);
L10:
    aa = dmax(as1,as2);
    if (aa > *ascle) {
	return 0;
    }
    s1->r = czero.r, s1->i = czero.i;
    s2->r = czero.r, s2->i = czero.i;
    *nz = 1;
    *iuf = 0;
    return 0;
} /* cs1s2_ */

