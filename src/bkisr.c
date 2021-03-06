/* bkisr.f -- translated by f2c (version 12.02.01).
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

/* Table of constant values */

static integer c__4 = 4;
static integer c__3 = 3;
static integer c__2 = 2;
static integer c__1 = 1;

/* DECK BKISR */
/* Subroutine */ int bkisr_(real *x, integer *n, real *sum, integer *ierr)
{
    /* Initialized data */

    static real c__[2] = { 1.57079632679489662f,1.f };

    /* System generated locals */
    integer i__1;
    real r__1;

    /* Local variables */
    static integer i__, k, k1;
    static real ak, bk, fk, fn;
    static integer kk, np;
    static real hx, pr;
    static integer kkn;
    static real pol, tkp, tol, xln, hxs, trm, atol;
    extern doublereal psixn_(integer *), r1mach_(integer *);

/* ***BEGIN PROLOGUE  BKISR */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to BSKIN */
/* ***LIBRARY   SLATEC */
/* ***TYPE      SINGLE PRECISION (BKISR-S, DBKISR-D) */
/* ***AUTHOR  Amos, D. E., (SNLA) */
/* ***DESCRIPTION */

/*     BKISR computes repeated integrals of the K0 Bessel function */
/*     by the series for N=0,1, and 2. */

/* ***SEE ALSO  BSKIN */
/* ***ROUTINES CALLED  PSIXN, R1MACH */
/* ***REVISION HISTORY  (YYMMDD) */
/*   820601  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900328  Added TYPE section.  (WRB) */
/*   910722  Updated AUTHOR section.  (ALS) */
/* ***END PROLOGUE  BKISR */

/* ***FIRST EXECUTABLE STATEMENT  BKISR */
    *ierr = 0;
/* Computing MAX */
    r__1 = r1mach_(&c__4);
    tol = dmax(r__1,1e-18f);
    if (*x < tol) {
	goto L50;
    }
    pr = 1.f;
    pol = 0.f;
    if (*n == 0) {
	goto L20;
    }
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	pol = -pol * *x + c__[i__ - 1];
	pr = pr * *x / i__;
/* L10: */
    }
L20:
    hx = *x * .5f;
    hxs = hx * hx;
    xln = log(hx);
    np = *n + 1;
    tkp = 3.f;
    fk = 2.f;
    fn = (real) (*n);
    bk = 4.f;
    ak = 2.f / ((fn + 1.f) * (fn + 2.f));
    i__1 = *n + 3;
    *sum = ak * (psixn_(&i__1) - psixn_(&c__3) + psixn_(&c__2) - xln);
    atol = *sum * tol * .75f;
    for (k = 2; k <= 20; ++k) {
	ak = ak * (hxs / bk) * ((tkp + 1.f) / (tkp + fn + 1.f)) * (tkp / (tkp 
		+ fn));
	k1 = k + 1;
	kk = k1 + k;
	kkn = kk + *n;
	trm = (psixn_(&k1) + psixn_(&kkn) - psixn_(&kk) - xln) * ak;
	*sum += trm;
	if (dabs(trm) <= atol) {
	    goto L40;
	}
	tkp += 2.f;
	bk += tkp;
	fk += 1.f;
/* L30: */
    }
    goto L80;
L40:
    *sum = (*sum * hxs + psixn_(&np) - xln) * pr;
    if (*n == 1) {
	*sum = -(*sum);
    }
    *sum = pol + *sum;
    return 0;
/* ----------------------------------------------------------------------- */
/*     SMALL X CASE, X.LT.WORD TOLERANCE */
/* ----------------------------------------------------------------------- */
L50:
    if (*n > 0) {
	goto L60;
    }
    hx = *x * .5f;
    *sum = psixn_(&c__1) - log(hx);
    return 0;
L60:
    *sum = c__[*n - 1];
    return 0;
L80:
    *ierr = 2;
    return 0;
} /* bkisr_ */

