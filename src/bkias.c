/* bkias.f -- translated by f2c (version 12.02.01).
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

/* DECK BKIAS */
/* Subroutine */ int bkias_(real *x, integer *n, integer *ktrms, real *t, 
	real *ans, integer *ind, integer *ms, real *gmrn, real *h__, integer *
	ierr)
{
    /* Initialized data */

    static real b[120] = { 1.f,1.f,-2.f,1.f,-8.f,6.f,1.f,-22.f,58.f,-24.f,1.f,
	    -52.f,328.f,-444.f,120.f,1.f,-114.f,1452.f,-4400.f,3708.f,-720.f,
	    1.f,-240.f,5610.f,-32120.f,58140.f,-33984.f,5040.f,1.f,-494.f,
	    19950.f,-195800.f,644020.f,-785304.f,341136.f,-40320.f,1.f,
	    -1004.f,67260.f,-1062500.f,5765500.f,-12440064.f,11026296.f,
	    -3733920.f,362880.f,1.f,-2026.f,218848.f,-5326160.f,4.4765e7f,
	    -155357384.f,238904904.f,-162186912.f,44339040.f,-3628800.f,1.f,
	    -4072.f,695038.f,-25243904.f,314369720.f,-1648384304.f,
	    4002695088.f,-4642163952.f,2507481216.f,-568356480.f,39916800.f,
	    1.f,-8166.f,2170626.f,-114876376.f,2051482776.f,-15548960784.f,
	    56041398784.f,-101180433024.f,92199790224.f,-40788301824.f,
	    7827719040.f,-479001600.f,1.f,-16356.f,6699696.f,-507259276.f,
	    12669817776.f,-134323420224.f,687720046384.f,-1818188642304.f,
	    2549865473424.f,-1883079661824.f,697929436800.f,-115336085760.f,
	    6227020800.f,1.f,-32738.f,20507988.f,-2189829808.f,75016052228.f,
	    -1084676512416.f,7634832149392.f,-28299910066112.f,
	    57494373464592.f,-64728375139872.f,39689578055808.f,
	    -12550904017920.f,1810992556800.f,-87178291200.f,1.f,-65504.f,
	    62407890.f,-9292526920.f,429826006340.f,-8308444327968.f,
	    78391384831312.f,-394365587815520.f,1111747472569680.f,
	    -1797171220690560.f,1666424486271456.f,-865023253219584.f,
	    236908271543040.f,-30196376985600.f,1.307674368e12f };
    static real bnd[15] = { 1.f,1.f,1.f,1.f,3.1f,5.18f,11.7f,29.8f,90.4f,
	    297.f,1070.f,4290.f,18100.f,84700.f,4.08e5f };
    static real hrtpi = .886226925452758014f;

    /* System generated locals */
    integer i__1, i__2, i__3;
    real r__1;

    /* Local variables */
    static integer i__, j, k;
    static real s[31], v[52], w[52], z__, g1, fj, fk;
    static integer ii, kk;
    static real er;
    static integer jn, km, mm;
    static real gs, hn;
    static integer mp;
    static real ss, xp[16], rz, fm1, rg1;
    static integer jmi;
    static real fln, rat, err, tol, rxp, rzx, den1, den2, den3, sumi, sumj;
    extern /* Subroutine */ int bdiff_(integer *, real *);
    extern doublereal gamrn_(real *);
    extern /* Subroutine */ int hkseq_(real *, integer *, real *, integer *);
    extern doublereal r1mach_(integer *);

/* ***BEGIN PROLOGUE  BKIAS */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to BSKIN */
/* ***LIBRARY   SLATEC */
/* ***TYPE      SINGLE PRECISION (BKIAS-S, DBKIAS-D) */
/* ***AUTHOR  Amos, D. E., (SNLA) */
/* ***DESCRIPTION */

/*     BKIAS computes repeated integrals of the K0 Bessel function */
/*     by the asymptotic expansion */

/* ***SEE ALSO  BSKIN */
/* ***ROUTINES CALLED  BDIFF, GAMRN, HKSEQ, R1MACH */
/* ***REVISION HISTORY  (YYMMDD) */
/*   820601  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900328  Added TYPE section.  (WRB) */
/*   910722  Updated AUTHOR section.  (ALS) */
/* ***END PROLOGUE  BKIAS */
/* ----------------------------------------------------------------------- */
/*             COEFFICIENTS OF POLYNOMIAL P(J-1,X), J=1,15 */
/* ----------------------------------------------------------------------- */
    /* Parameter adjustments */
    --h__;
    --t;

    /* Function Body */
/* ----------------------------------------------------------------------- */
/*             BOUNDS B(M,K) , K=M-3 */
/* ----------------------------------------------------------------------- */

/* ***FIRST EXECUTABLE STATEMENT  BKIAS */
    *ierr = 0;
/* Computing MAX */
    r__1 = r1mach_(&c__4);
    tol = dmax(r__1,1e-18f);
    fln = (real) (*n);
    rz = 1.f / (*x + fln);
    rzx = *x * rz;
    z__ = (*x + fln) * .5f;
    if (*ind > 1) {
	goto L10;
    }
    *gmrn = gamrn_(&z__);
L10:
    gs = hrtpi * *gmrn;
    g1 = gs + gs;
    rg1 = 1.f / g1;
    *gmrn = (rz + rz) / *gmrn;
    if (*ind > 1) {
	goto L70;
    }
/* ----------------------------------------------------------------------- */
/*     EVALUATE ERROR FOR M=MS */
/* ----------------------------------------------------------------------- */
    hn = fln * .5f;
    den2 = (real) (*ktrms + *ktrms + *n);
    den3 = den2 - 2.f;
    den1 = *x + den2;
    err = rg1 * (*x + *x) / (den1 - 1.f);
    if (*n == 0) {
	goto L20;
    }
    rat = 1.f / (fln * fln);
L20:
    if (*ktrms == 0) {
	goto L30;
    }
    fj = (real) (*ktrms);
    rat = .25f / (hrtpi * den3 * sqrt(fj));
L30:
    err *= rat;
    fj = -3.f;
    for (j = 1; j <= 15; ++j) {
	if (j <= 5) {
	    err /= den1;
	}
	fm1 = dmax(1.f,fj);
	fj += 1.f;
	er = bnd[j - 1] * err;
	if (*ktrms == 0) {
	    goto L40;
	}
	er /= fm1;
	if (er < tol) {
	    goto L60;
	}
	if (j >= 5) {
	    err /= den3;
	}
	goto L50;
L40:
	er *= hn / fm1 + 1.f;
	if (er < tol) {
	    goto L60;
	}
	if (j >= 5) {
	    err /= fln;
	}
L50:
	;
    }
    goto L200;
L60:
    *ms = j;
L70:
    mm = *ms + *ms;
    mp = mm + 1;
/* ----------------------------------------------------------------------- */
/*     H(K)=(-Z)**(K)*(PSI(K-1,Z)-PSI(K-1,Z+0.5))/GAMMA(K) , K=1,2,...,MM */
/* ----------------------------------------------------------------------- */
    if (*ind > 1) {
	goto L80;
    }
    hkseq_(&z__, &mm, &h__[1], ierr);
    goto L100;
L80:
    rat = z__ / (z__ - .5f);
    rxp = rat;
    i__1 = mm;
    for (i__ = 1; i__ <= i__1; ++i__) {
	h__[i__] = rxp * (1.f - h__[i__]);
	rxp *= rat;
/* L90: */
    }
L100:
/* ----------------------------------------------------------------------- */
/*     SCALED S SEQUENCE */
/* ----------------------------------------------------------------------- */
    s[0] = 1.f;
    fk = 1.f;
    i__1 = mp;
    for (k = 2; k <= i__1; ++k) {
	ss = 0.f;
	km = k - 1;
	i__ = km;
	i__2 = km;
	for (j = 1; j <= i__2; ++j) {
	    ss += s[j - 1] * h__[i__];
	    --i__;
/* L110: */
	}
	s[k - 1] = ss / fk;
	fk += 1.f;
/* L120: */
    }
/* ----------------------------------------------------------------------- */
/*     SCALED S-TILDA SEQUENCE */
/* ----------------------------------------------------------------------- */
    if (*ktrms == 0) {
	goto L160;
    }
    fk = 0.f;
    ss = 0.f;
    rg1 /= z__;
    i__1 = *ktrms;
    for (k = 1; k <= i__1; ++k) {
	v[k - 1] = z__ / (z__ + fk);
	w[k - 1] = t[k] * v[k - 1];
	ss += w[k - 1];
	fk += 1.f;
/* L130: */
    }
    s[0] -= ss * rg1;
    i__1 = mp;
    for (i__ = 2; i__ <= i__1; ++i__) {
	ss = 0.f;
	i__2 = *ktrms;
	for (k = 1; k <= i__2; ++k) {
	    w[k - 1] *= v[k - 1];
	    ss += w[k - 1];
/* L140: */
	}
	s[i__ - 1] -= ss * rg1;
/* L150: */
    }
L160:
/* ----------------------------------------------------------------------- */
/*     SUM ON J */
/* ----------------------------------------------------------------------- */
    sumj = 0.f;
    jn = 1;
    rxp = 1.f;
    xp[0] = 1.f;
    i__1 = *ms;
    for (j = 1; j <= i__1; ++j) {
	jn = jn + j - 1;
	xp[j] = xp[j - 1] * rzx;
	rxp *= rz;
/* ----------------------------------------------------------------------- */
/*     SUM ON I */
/* ----------------------------------------------------------------------- */
	sumi = 0.f;
	ii = jn;
	i__2 = j;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    jmi = j - i__ + 1;
	    kk = j + i__ + 1;
	    i__3 = jmi;
	    for (k = 1; k <= i__3; ++k) {
		v[k - 1] = s[kk - 1] * xp[k - 1];
		++kk;
/* L170: */
	    }
	    bdiff_(&jmi, v);
	    sumi += b[ii - 1] * v[jmi - 1] * xp[i__];
	    ++ii;
/* L180: */
	}
	sumj += sumi * rxp;
/* L190: */
    }
    *ans = gs * (s[0] - sumj);
    return 0;
L200:
    *ierr = 2;
    return 0;
} /* bkias_ */

