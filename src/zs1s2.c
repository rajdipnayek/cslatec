/* zs1s2.f -- translated by f2c (version 12.02.01).
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

/* DECK ZS1S2 */
/* Subroutine */ int zs1s2_(doublereal *zrr, doublereal *zri, doublereal *s1r,
	 doublereal *s1i, doublereal *s2r, doublereal *s2i, integer *nz, 
	doublereal *ascle, doublereal *alim, integer *iuf)
{
    /* Initialized data */

    static doublereal zeror = 0.;
    static doublereal zeroi = 0.;

    /* Local variables */
    static doublereal aa, c1i, as1, as2, c1r, aln, s1di, s1dr;
    static integer idum;
    extern doublereal zabs_(doublereal *, doublereal *);
    extern /* Subroutine */ int zlog_(doublereal *, doublereal *, doublereal *
	    , doublereal *, integer *), zexp_(doublereal *, doublereal *, 
	    doublereal *, doublereal *);

/* ***BEGIN PROLOGUE  ZS1S2 */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to ZAIRY and ZBESK */
/* ***LIBRARY   SLATEC */
/* ***TYPE      ALL (CS1S2-A, ZS1S2-A) */
/* ***AUTHOR  Amos, D. E., (SNL) */
/* ***DESCRIPTION */

/*     ZS1S2 TESTS FOR A POSSIBLE UNDERFLOW RESULTING FROM THE */
/*     ADDITION OF THE I AND K FUNCTIONS IN THE ANALYTIC CON- */
/*     TINUATION FORMULA WHERE S1=K FUNCTION AND S2=I FUNCTION. */
/*     ON KODE=1 THE I AND K FUNCTIONS ARE DIFFERENT ORDERS OF */
/*     MAGNITUDE, BUT FOR KODE=2 THEY CAN BE OF THE SAME ORDER */
/*     OF MAGNITUDE AND THE MAXIMUM MUST BE AT LEAST ONE */
/*     PRECISION ABOVE THE UNDERFLOW LIMIT. */

/* ***SEE ALSO  ZAIRY, ZBESK */
/* ***ROUTINES CALLED  ZABS, ZEXP, ZLOG */
/* ***REVISION HISTORY  (YYMMDD) */
/*   830501  DATE WRITTEN */
/*   910415  Prologue converted to Version 4.0 format.  (BAB) */
/*   930122  Added ZEXP and ZLOG to EXTERNAL statement.  (RWC) */
/* ***END PROLOGUE  ZS1S2 */
/*     COMPLEX CZERO,C1,S1,S1D,S2,ZR */
/* ***FIRST EXECUTABLE STATEMENT  ZS1S2 */
    *nz = 0;
    as1 = zabs_(s1r, s1i);
    as2 = zabs_(s2r, s2i);
    if (*s1r == 0. && *s1i == 0.) {
	goto L10;
    }
    if (as1 == 0.) {
	goto L10;
    }
    aln = -(*zrr) - *zrr + log(as1);
    s1dr = *s1r;
    s1di = *s1i;
    *s1r = zeror;
    *s1i = zeroi;
    as1 = zeror;
    if (aln < -(*alim)) {
	goto L10;
    }
    zlog_(&s1dr, &s1di, &c1r, &c1i, &idum);
    c1r = c1r - *zrr - *zrr;
    c1i = c1i - *zri - *zri;
    zexp_(&c1r, &c1i, s1r, s1i);
    as1 = zabs_(s1r, s1i);
    ++(*iuf);
L10:
    aa = max(as1,as2);
    if (aa > *ascle) {
	return 0;
    }
    *s1r = zeror;
    *s1i = zeroi;
    *s2r = zeror;
    *s2i = zeroi;
    *nz = 1;
    *iuf = 0;
    return 0;
} /* zs1s2_ */

