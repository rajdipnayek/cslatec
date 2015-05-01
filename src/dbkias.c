/* dbkias.f -- translated by f2c (version 12.02.01).
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

/* DECK DBKIAS */
/* Subroutine */ int dbkias_(doublereal *x, integer *n, integer *ktrms, 
	doublereal *t, doublereal *ans, integer *ind, integer *ms, doublereal 
	*gmrn, doublereal *h__, integer *ierr)
{
    /* Initialized data */

    static doublereal b[120] = { 1.,1.,-2.,1.,-8.,6.,1.,-22.,58.,-24.,1.,-52.,
	    328.,-444.,120.,1.,-114.,1452.,-4400.,3708.,-720.,1.,-240.,5610.,
	    -32120.,58140.,-33984.,5040.,1.,-494.,19950.,-195800.,644020.,
	    -785304.,341136.,-40320.,1.,-1004.,67260.,-1062500.,5765500.,
	    -12440064.,11026296.,-3733920.,362880.,1.,-2026.,218848.,
	    -5326160.,4.4765e7,-155357384.,238904904.,-162186912.,44339040.,
	    -3628800.,1.,-4072.,695038.,-25243904.,314369720.,-1648384304.,
	    4002695088.,-4642163952.,2507481216.,-568356480.,39916800.,1.,
	    -8166.,2170626.,-114876376.,2051482776.,-15548960784.,
	    56041398784.,-101180433024.,92199790224.,-40788301824.,
	    7827719040.,-479001600.,1.,-16356.,6699696.,-507259276.,
	    12669817776.,-134323420224.,687720046384.,-1818188642304.,
	    2549865473424.,-1883079661824.,697929436800.,-115336085760.,
	    6227020800.,1.,-32738.,20507988.,-2189829808.,75016052228.,
	    -1084676512416.,7634832149392.,-28299910066112.,57494373464592.,
	    -64728375139872.,39689578055808.,-12550904017920.,1810992556800.,
	    -87178291200.,1.,-65504.,62407890.,-9292526920.,429826006340.,
	    -8308444327968.,78391384831312.,-394365587815520.,
	    1111747472569680.,-1797171220690560.,1666424486271456.,
	    -865023253219584.,236908271543040.,-30196376985600.,
	    1.307674368e12 };
    static doublereal bnd[15] = { 1.,1.,1.,1.,3.1,5.18,11.7,29.8,90.4,297.,
	    1070.,4290.,18100.,84700.,4.08e5 };
    static doublereal hrtpi = .886226925452758014;

    /* System generated locals */
    integer i__1, i__2, i__3;
    doublereal d__1;

    /* Local variables */
    static integer i__, j, k;
    static doublereal s[31], v[52], w[52], z__, g1, fj, fk;
    static integer ii, kk;
    static doublereal er;
    static integer jn, km, mm;
    static doublereal gs, hn;
    static integer mp;
    static doublereal ss, xp[16], rz, fm1, rg1;
    static integer jmi;
    static doublereal fln, rat, err, tol, rxp, rzx, den1, den2, den3, sumi, 
	    sumj;
    extern doublereal d1mach_(integer *);
    extern /* Subroutine */ int dbdiff_(integer *, doublereal *);
    extern doublereal dgamrn_(doublereal *);
    extern /* Subroutine */ int dhkseq_(doublereal *, integer *, doublereal *,
	     integer *);

/* ***BEGIN PROLOGUE  DBKIAS */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to DBSKIN */
/* ***LIBRARY   SLATEC */
/* ***TYPE      DOUBLE PRECISION (BKIAS-S, DBKIAS-D) */
/* ***AUTHOR  Amos, D. E., (SNLA) */
/* ***DESCRIPTION */

/*     DBKIAS computes repeated integrals of the K0 Bessel function */
/*     by the asymptotic expansion */

/* ***SEE ALSO  DBSKIN */
/* ***ROUTINES CALLED  D1MACH, DBDIFF, DGAMRN, DHKSEQ */
/* ***REVISION HISTORY  (YYMMDD) */
/*   820601  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890911  Removed unnecessary intrinsics.  (WRB) */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900328  Added TYPE section.  (WRB) */
/*   910722  Updated AUTHOR section.  (ALS) */
/* ***END PROLOGUE  DBKIAS */
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

/* ***FIRST EXECUTABLE STATEMENT  DBKIAS */
    *ierr = 0;
/* Computing MAX */
    d__1 = d1mach_(&c__4);
    tol = max(d__1,1e-18);
    fln = (doublereal) (*n);
    rz = 1. / (*x + fln);
    rzx = *x * rz;
    z__ = (*x + fln) * .5;
    if (*ind > 1) {
	goto L10;
    }
    *gmrn = dgamrn_(&z__);
L10:
    gs = hrtpi * *gmrn;
    g1 = gs + gs;
    rg1 = 1. / g1;
    *gmrn = (rz + rz) / *gmrn;
    if (*ind > 1) {
	goto L70;
    }
/* ----------------------------------------------------------------------- */
/*     EVALUATE ERROR FOR M=MS */
/* ----------------------------------------------------------------------- */
    hn = fln * .5;
    den2 = (doublereal) (*ktrms + *ktrms + *n);
    den3 = den2 - 2.;
    den1 = *x + den2;
    err = rg1 * (*x + *x) / (den1 - 1.);
    if (*n == 0) {
	goto L20;
    }
    rat = 1. / (fln * fln);
L20:
    if (*ktrms == 0) {
	goto L30;
    }
    fj = (doublereal) (*ktrms);
    rat = .25 / (hrtpi * den3 * sqrt(fj));
L30:
    err *= rat;
    fj = -3.;
    for (j = 1; j <= 15; ++j) {
	if (j <= 5) {
	    err /= den1;
	}
	fm1 = max(1.,fj);
	fj += 1.;
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
	er *= hn / fm1 + 1.;
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
    dhkseq_(&z__, &mm, &h__[1], ierr);
    goto L100;
L80:
    rat = z__ / (z__ - .5);
    rxp = rat;
    i__1 = mm;
    for (i__ = 1; i__ <= i__1; ++i__) {
	h__[i__] = rxp * (1. - h__[i__]);
	rxp *= rat;
/* L90: */
    }
L100:
/* ----------------------------------------------------------------------- */
/*     SCALED S SEQUENCE */
/* ----------------------------------------------------------------------- */
    s[0] = 1.;
    fk = 1.;
    i__1 = mp;
    for (k = 2; k <= i__1; ++k) {
	ss = 0.;
	km = k - 1;
	i__ = km;
	i__2 = km;
	for (j = 1; j <= i__2; ++j) {
	    ss += s[j - 1] * h__[i__];
	    --i__;
/* L110: */
	}
	s[k - 1] = ss / fk;
	fk += 1.;
/* L120: */
    }
/* ----------------------------------------------------------------------- */
/*     SCALED S-TILDA SEQUENCE */
/* ----------------------------------------------------------------------- */
    if (*ktrms == 0) {
	goto L160;
    }
    fk = 0.;
    ss = 0.;
    rg1 /= z__;
    i__1 = *ktrms;
    for (k = 1; k <= i__1; ++k) {
	v[k - 1] = z__ / (z__ + fk);
	w[k - 1] = t[k] * v[k - 1];
	ss += w[k - 1];
	fk += 1.;
/* L130: */
    }
    s[0] -= ss * rg1;
    i__1 = mp;
    for (i__ = 2; i__ <= i__1; ++i__) {
	ss = 0.;
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
    sumj = 0.;
    jn = 1;
    rxp = 1.;
    xp[0] = 1.;
    i__1 = *ms;
    for (j = 1; j <= i__1; ++j) {
	jn = jn + j - 1;
	xp[j] = xp[j - 1] * rzx;
	rxp *= rz;
/* ----------------------------------------------------------------------- */
/*     SUM ON I */
/* ----------------------------------------------------------------------- */
	sumi = 0.;
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
	    dbdiff_(&jmi, v);
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
} /* dbkias_ */

