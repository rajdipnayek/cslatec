/* besj.f -- translated by f2c (version 12.02.01).
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

static integer c__3 = 3;
static integer c__11 = 11;
static integer c__12 = 12;
static integer c__5 = 5;
static integer c__1 = 1;
static integer c__2 = 2;

/* DECK BESJ */
/* Subroutine */ int besj_(real *x, real *alpha, integer *n, real *y, integer 
	*nz)
{
    /* Initialized data */

    static real rtwo = 1.34839972492648f;
    static real pdf = .785398163397448f;
    static real rttp = .797884560802865f;
    static real pidt = 1.5707963267949f;
    static real pp[4] = { 8.72909153935547f,.26569393226503f,
	    .124578576865586f,7.70133747430388e-4f };
    static integer inlim = 150;
    static real fnulim[2] = { 100.f,60.f };

    /* System generated locals */
    integer i__1;
    real r__1;

    /* Local variables */
    static integer i__, k;
    static real s, t;
    static integer i1, i2;
    static real s1, s2, t1, t2, ak, ap, fn, sa;
    static integer kk, in, km;
    static real sb, ta, tb;
    static integer is, nn, kt, ns;
    static real tm, wk[7], tx, xo2, dfn, akm, arg, fnf, fni, gln, ans, dtm, 
	    tfn, fnu, tau, tol, etx, rtx, trx, fnp1, xo2l, sxo2, coef, earg, 
	    relb;
    static integer ialp;
    static real rden;
    static integer iflw;
    static real slim, temp[3], rtol, elim1, fidal;
    static integer idalp;
    static real flgjy;
    extern /* Subroutine */ int jairy_();
    static real rzden, tolln;
    extern /* Subroutine */ int asyjy_(U_fp, real *, real *, real *, integer *
	    , real *, real *, integer *);
    extern integer i1mach_(integer *);
    extern doublereal r1mach_(integer *);
    static real dalpha;
    extern doublereal alngam_(real *);
    extern /* Subroutine */ int xermsg_(char *, char *, char *, integer *, 
	    integer *, ftnlen, ftnlen, ftnlen);

/* ***BEGIN PROLOGUE  BESJ */
/* ***PURPOSE  Compute an N member sequence of J Bessel functions */
/*            J/SUB(ALPHA+K-1)/(X), K=1,...,N for non-negative ALPHA */
/*            and X. */
/* ***LIBRARY   SLATEC */
/* ***CATEGORY  C10A3 */
/* ***TYPE      SINGLE PRECISION (BESJ-S, DBESJ-D) */
/* ***KEYWORDS  J BESSEL FUNCTION, SPECIAL FUNCTIONS */
/* ***AUTHOR  Amos, D. E., (SNLA) */
/*           Daniel, S. L., (SNLA) */
/*           Weston, M. K., (SNLA) */
/* ***DESCRIPTION */

/*     Abstract */
/*         BESJ computes an N member sequence of J Bessel functions */
/*         J/sub(ALPHA+K-1)/(X), K=1,...,N for non-negative ALPHA and X. */
/*         A combination of the power series, the asymptotic expansion */
/*         for X to infinity and the uniform asymptotic expansion for */
/*         NU to infinity are applied over subdivisions of the (NU,X) */
/*         plane.  For values of (NU,X) not covered by one of these */
/*         formulae, the order is incremented or decremented by integer */
/*         values into a region where one of the formulae apply. Backward */
/*         recursion is applied to reduce orders by integer values except */
/*         where the entire sequence lies in the oscillatory region.  In */
/*         this case forward recursion is stable and values from the */
/*         asymptotic expansion for X to infinity start the recursion */
/*         when it is efficient to do so.  Leading terms of the series */
/*         and uniform expansion are tested for underflow.  If a sequence */
/*         is requested and the last member would underflow, the result */
/*         is set to zero and the next lower order tried, etc., until a */
/*         member comes on scale or all members are set to zero. */
/*         Overflow cannot occur. */

/*     Description of Arguments */

/*         Input */
/*           X      - X .GE. 0.0E0 */
/*           ALPHA  - order of first member of the sequence, */
/*                    ALPHA .GE. 0.0E0 */
/*           N      - number of members in the sequence, N .GE. 1 */

/*         Output */
/*           Y      - a vector whose first  N components contain */
/*                    values for J/sub(ALPHA+K-1)/(X), K=1,...,N */
/*           NZ     - number of components of Y set to zero due to */
/*                    underflow, */
/*                    NZ=0   , normal return, computation completed */
/*                    NZ .NE. 0, last NZ components of Y set to zero, */
/*                             Y(K)=0.0E0, K=N-NZ+1,...,N. */

/*     Error Conditions */
/*         Improper input arguments - a fatal error */
/*         Underflow  - a non-fatal error (NZ .NE. 0) */

/* ***REFERENCES  D. E. Amos, S. L. Daniel and M. K. Weston, CDC 6600 */
/*                 subroutines IBESS and JBESS for Bessel functions */
/*                 I(NU,X) and J(NU,X), X .GE. 0, NU .GE. 0, ACM */
/*                 Transactions on Mathematical Software 3, (1977), */
/*                 pp. 76-92. */
/*               F. W. J. Olver, Tables of Bessel Functions of Moderate */
/*                 or Large Orders, NPL Mathematical Tables 6, Her */
/*                 Majesty's Stationery Office, London, 1962. */
/* ***ROUTINES CALLED  ALNGAM, ASYJY, I1MACH, JAIRY, R1MACH, XERMSG */
/* ***REVISION HISTORY  (YYMMDD) */
/*   750101  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890531  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ) */
/*   900326  Removed duplicate information from DESCRIPTION section. */
/*           (WRB) */
/*   920501  Reformatted the REFERENCES section.  (WRB) */
/* ***END PROLOGUE  BESJ */
    /* Parameter adjustments */
    --y;

    /* Function Body */
/* ***FIRST EXECUTABLE STATEMENT  BESJ */
    *nz = 0;
    kt = 1;
    ns = 0;
/*     I1MACH(14) REPLACES I1MACH(11) IN A DOUBLE PRECISION CODE */
/*     I1MACH(15) REPLACES I1MACH(12) IN A DOUBLE PRECISION CODE */
    ta = r1mach_(&c__3);
    tol = dmax(ta,1e-15f);
    i1 = i1mach_(&c__11) + 1;
    i2 = i1mach_(&c__12);
    tb = r1mach_(&c__5);
    elim1 = (i2 * tb + 3.f) * -2.303f;
    rtol = 1.f / tol;
    slim = r1mach_(&c__1) * 1e3f * rtol;
/*     TOLLN = -LN(TOL) */
    tolln = tb * 2.303f * i1;
    tolln = dmin(tolln,34.5388f);
    if ((i__1 = *n - 1) < 0) {
	goto L720;
    } else if (i__1 == 0) {
	goto L10;
    } else {
	goto L20;
    }
L10:
    kt = 2;
L20:
    nn = *n;
    if (*x < 0.f) {
	goto L730;
    } else if (*x == 0) {
	goto L30;
    } else {
	goto L80;
    }
L30:
    if (*alpha < 0.f) {
	goto L710;
    } else if (*alpha == 0) {
	goto L40;
    } else {
	goto L50;
    }
L40:
    y[1] = 1.f;
    if (*n == 1) {
	return 0;
    }
    i1 = 2;
    goto L60;
L50:
    i1 = 1;
L60:
    i__1 = *n;
    for (i__ = i1; i__ <= i__1; ++i__) {
	y[i__] = 0.f;
/* L70: */
    }
    return 0;
L80:
    if (*alpha < 0.f) {
	goto L710;
    }

    ialp = (integer) (*alpha);
    fni = (real) (ialp + *n - 1);
    fnf = *alpha - ialp;
    dfn = fni + fnf;
    fnu = dfn;
    xo2 = *x * .5f;
    sxo2 = xo2 * xo2;

/*     DECISION TREE FOR REGION WHERE SERIES, ASYMPTOTIC EXPANSION FOR X */
/*     TO INFINITY AND ASYMPTOTIC EXPANSION FOR NU TO INFINITY ARE */
/*     APPLIED. */

    if (sxo2 <= fnu + 1.f) {
	goto L90;
    }
    ta = dmax(20.f,fnu);
    if (*x > ta) {
	goto L120;
    }
    if (*x > 12.f) {
	goto L110;
    }
    xo2l = log(xo2);
    ns = (integer) (sxo2 - fnu) + 1;
    goto L100;
L90:
    fn = fnu;
    fnp1 = fn + 1.f;
    xo2l = log(xo2);
    is = kt;
    if (*x <= .5f) {
	goto L330;
    }
    ns = 0;
L100:
    fni += ns;
    dfn = fni + fnf;
    fn = dfn;
    fnp1 = fn + 1.f;
    is = kt;
    if (*n - 1 + ns > 0) {
	is = 3;
    }
    goto L330;
L110:
/* Computing MAX */
    r__1 = 36.f - fnu;
    ans = dmax(r__1,0.f);
    ns = (integer) ans;
    fni += ns;
    dfn = fni + fnf;
    fn = dfn;
    is = kt;
    if (*n - 1 + ns > 0) {
	is = 3;
    }
    goto L130;
L120:
    rtx = sqrt(*x);
    tau = rtwo * rtx;
    ta = tau + fnulim[kt - 1];
    if (fnu <= ta) {
	goto L480;
    }
    fn = fnu;
    is = kt;

/*     UNIFORM ASYMPTOTIC EXPANSION FOR NU TO INFINITY */

L130:
    i1 = (i__1 = 3 - is, abs(i__1));
    i1 = max(i1,1);
    flgjy = 1.f;
    asyjy_((U_fp)jairy_, x, &fn, &flgjy, &i1, &temp[is - 1], wk, &iflw);
    if (iflw != 0) {
	goto L380;
    }
    switch (is) {
	case 1:  goto L320;
	case 2:  goto L450;
	case 3:  goto L620;
    }
L310:
    temp[0] = temp[2];
    kt = 1;
L320:
    is = 2;
    fni += -1.f;
    dfn = fni + fnf;
    fn = dfn;
    if (i1 == 2) {
	goto L450;
    }
    goto L130;

/*     SERIES FOR (X/2)**2.LE.NU+1 */

L330:
    gln = alngam_(&fnp1);
    arg = fn * xo2l - gln;
    if (arg < -elim1) {
	goto L400;
    }
    earg = exp(arg);
L340:
    s = 1.f;
    if (*x < tol) {
	goto L360;
    }
    ak = 3.f;
    t2 = 1.f;
    t = 1.f;
    s1 = fn;
    for (k = 1; k <= 17; ++k) {
	s2 = t2 + s1;
	t = -t * sxo2 / s2;
	s += t;
	if (dabs(t) < tol) {
	    goto L360;
	}
	t2 += ak;
	ak += 2.f;
	s1 += fn;
/* L350: */
    }
L360:
    temp[is - 1] = s * earg;
    switch (is) {
	case 1:  goto L370;
	case 2:  goto L450;
	case 3:  goto L610;
    }
L370:
    earg = earg * fn / xo2;
    fni += -1.f;
    dfn = fni + fnf;
    fn = dfn;
    is = 2;
    goto L340;

/*     SET UNDERFLOW VALUE AND UPDATE PARAMETERS */
/*     UNDERFLOW CAN ONLY OCCUR FOR NS=0 SINCE THE ORDER MUST BE */
/*     LARGER THAN 36. THEREFORE, NS NEED NOT BE CONSIDERED. */

L380:
    y[nn] = 0.f;
    --nn;
    fni += -1.f;
    dfn = fni + fnf;
    fn = dfn;
    if ((i__1 = nn - 1) < 0) {
	goto L440;
    } else if (i__1 == 0) {
	goto L390;
    } else {
	goto L130;
    }
L390:
    kt = 2;
    is = 2;
    goto L130;
L400:
    y[nn] = 0.f;
    --nn;
    fnp1 = fn;
    fni += -1.f;
    dfn = fni + fnf;
    fn = dfn;
    if ((i__1 = nn - 1) < 0) {
	goto L440;
    } else if (i__1 == 0) {
	goto L410;
    } else {
	goto L420;
    }
L410:
    kt = 2;
    is = 2;
L420:
    if (sxo2 <= fnp1) {
	goto L430;
    }
    goto L130;
L430:
    arg = arg - xo2l + log(fnp1);
    if (arg < -elim1) {
	goto L400;
    }
    goto L330;
L440:
    *nz = *n - nn;
    return 0;

/*     BACKWARD RECURSION SECTION */

L450:
    if (ns != 0) {
	goto L451;
    }
    *nz = *n - nn;
    if (kt == 2) {
	goto L470;
    }
/*     BACKWARD RECUR FROM INDEX ALPHA+NN-1 TO ALPHA */
    y[nn] = temp[0];
    y[nn - 1] = temp[1];
    if (nn == 2) {
	return 0;
    }
L451:
    trx = 2.f / *x;
    dtm = fni;
    tm = (dtm + fnf) * trx;
    ak = 1.f;
    ta = temp[0];
    tb = temp[1];
    if (dabs(ta) > slim) {
	goto L455;
    }
    ta *= rtol;
    tb *= rtol;
    ak = tol;
L455:
    kk = 2;
    in = ns - 1;
    if (in == 0) {
	goto L690;
    }
    if (ns != 0) {
	goto L670;
    }
    k = nn - 2;
    i__1 = nn;
    for (i__ = 3; i__ <= i__1; ++i__) {
	s = tb;
	tb = tm * tb - ta;
	ta = s;
	y[k] = tb * ak;
	--k;
	dtm += -1.f;
	tm = (dtm + fnf) * trx;
/* L460: */
    }
    return 0;
L470:
    y[1] = temp[1];
    return 0;

/*     ASYMPTOTIC EXPANSION FOR X TO INFINITY WITH FORWARD RECURSION IN */
/*     OSCILLATORY REGION X.GT.MAX(20, NU), PROVIDED THE LAST MEMBER */
/*     OF THE SEQUENCE IS ALSO IN THE REGION. */

L480:
    in = (integer) (*alpha - tau + 2.f);
    if (in <= 0) {
	goto L490;
    }
    idalp = ialp - in - 1;
    kt = 1;
    goto L500;
L490:
    idalp = ialp;
    in = 0;
L500:
    is = kt;
    fidal = (real) idalp;
    dalpha = fidal + fnf;
    arg = *x - pidt * dalpha - pdf;
    sa = sin(arg);
    sb = cos(arg);
    coef = rttp / rtx;
    etx = *x * 8.f;
L510:
    dtm = fidal + fidal;
    dtm *= dtm;
    tm = 0.f;
    if (fidal == 0.f && dabs(fnf) < tol) {
	goto L520;
    }
    tm = fnf * 4.f * (fidal + fidal + fnf);
L520:
    trx = dtm - 1.f;
    t2 = (trx + tm) / etx;
    s2 = t2;
    relb = tol * dabs(t2);
    t1 = etx;
    s1 = 1.f;
    fn = 1.f;
    ak = 8.f;
    for (k = 1; k <= 13; ++k) {
	t1 += etx;
	fn += ak;
	trx = dtm - fn;
	ap = trx + tm;
	t2 = -t2 * ap / t1;
	s1 += t2;
	t1 += etx;
	ak += 8.f;
	fn += ak;
	trx = dtm - fn;
	ap = trx + tm;
	t2 = t2 * ap / t1;
	s2 += t2;
	if (dabs(t2) <= relb) {
	    goto L540;
	}
	ak += 8.f;
/* L530: */
    }
L540:
    temp[is - 1] = coef * (s1 * sb - s2 * sa);
    if (is == 2) {
	goto L560;
    }
    fidal += 1.f;
    dalpha = fidal + fnf;
    is = 2;
    tb = sa;
    sa = -sb;
    sb = tb;
    goto L510;

/*     FORWARD RECURSION SECTION */

L560:
    if (kt == 2) {
	goto L470;
    }
    s1 = temp[0];
    s2 = temp[1];
    tx = 2.f / *x;
    tm = dalpha * tx;
    if (in == 0) {
	goto L580;
    }

/*     FORWARD RECUR TO INDEX ALPHA */

    i__1 = in;
    for (i__ = 1; i__ <= i__1; ++i__) {
	s = s2;
	s2 = tm * s2 - s1;
	tm += tx;
	s1 = s;
/* L570: */
    }
    if (nn == 1) {
	goto L600;
    }
    s = s2;
    s2 = tm * s2 - s1;
    tm += tx;
    s1 = s;
L580:

/*     FORWARD RECUR FROM INDEX ALPHA TO ALPHA+N-1 */

    y[1] = s1;
    y[2] = s2;
    if (nn == 2) {
	return 0;
    }
    i__1 = nn;
    for (i__ = 3; i__ <= i__1; ++i__) {
	y[i__] = tm * y[i__ - 1] - y[i__ - 2];
	tm += tx;
/* L590: */
    }
    return 0;
L600:
    y[1] = s2;
    return 0;

/*     BACKWARD RECURSION WITH NORMALIZATION BY */
/*     ASYMPTOTIC EXPANSION FOR NU TO INFINITY OR POWER SERIES. */

L610:
/*     COMPUTATION OF LAST ORDER FOR SERIES NORMALIZATION */
/* Computing MAX */
    r__1 = 3.f - fn;
    akm = dmax(r__1,0.f);
    km = (integer) akm;
    tfn = fn + km;
    ta = (gln + tfn - .9189385332f - .0833333333f / tfn) / (tfn + .5f);
    ta = xo2l - ta;
    tb = -(1.f - 1.5f / tfn) / tfn;
    akm = tolln / (-ta + sqrt(ta * ta - tolln * tb)) + 1.5f;
    in = km + (integer) akm;
    goto L660;
L620:
/*     COMPUTATION OF LAST ORDER FOR ASYMPTOTIC EXPANSION NORMALIZATION */
    gln = wk[2] + wk[1];
    if (wk[5] > 30.f) {
	goto L640;
    }
    rden = (pp[3] * wk[5] + pp[2]) * wk[5] + 1.f;
    rzden = pp[0] + pp[1] * wk[5];
    ta = rzden / rden;
    if (wk[0] < .1f) {
	goto L630;
    }
    tb = gln / wk[4];
    goto L650;
L630:
    tb = ((wk[0] * .0887944358f + .167989473f) * wk[0] + 1.259921049f) / wk[6]
	    ;
    goto L650;
L640:
    ta = tolln * .5f / wk[3];
    ta = ((ta * .049382716f - .1111111111f) * ta + .6666666667f) * ta * wk[5];
    if (wk[0] < .1f) {
	goto L630;
    }
    tb = gln / wk[4];
L650:
    in = (integer) (ta / tb + 1.5f);
    if (in > inlim) {
	goto L310;
    }
L660:
    dtm = fni + in;
    trx = 2.f / *x;
    tm = (dtm + fnf) * trx;
    ta = 0.f;
    tb = tol;
    kk = 1;
    ak = 1.f;
L670:

/*     BACKWARD RECUR UNINDEXED AND SCALE WHEN MAGNITUDES ARE CLOSE TO */
/*     UNDERFLOW LIMITS (LESS THAN SLIM=R1MACH(1)*1.0E+3/TOL) */

    i__1 = in;
    for (i__ = 1; i__ <= i__1; ++i__) {
	s = tb;
	tb = tm * tb - ta;
	ta = s;
	dtm += -1.f;
	tm = (dtm + fnf) * trx;
/* L680: */
    }
/*     NORMALIZATION */
    if (kk != 1) {
	goto L690;
    }
    s = temp[2];
    sa = ta / tb;
    ta = s;
    tb = s;
    if (dabs(s) > slim) {
	goto L685;
    }
    ta *= rtol;
    tb *= rtol;
    ak = tol;
L685:
    ta *= sa;
    kk = 2;
    in = ns;
    if (ns != 0) {
	goto L670;
    }
L690:
    y[nn] = tb * ak;
    *nz = *n - nn;
    if (nn == 1) {
	return 0;
    }
    k = nn - 1;
    s = tb;
    tb = tm * tb - ta;
    ta = s;
    y[k] = tb * ak;
    if (nn == 2) {
	return 0;
    }
    dtm += -1.f;
    tm = (dtm + fnf) * trx;
    k = nn - 2;

/*     BACKWARD RECUR INDEXED */

    i__1 = nn;
    for (i__ = 3; i__ <= i__1; ++i__) {
	s = tb;
	tb = tm * tb - ta;
	ta = s;
	y[k] = tb * ak;
	dtm += -1.f;
	tm = (dtm + fnf) * trx;
	--k;
/* L700: */
    }
    return 0;



L710:
    xermsg_("SLATEC", "BESJ", "ORDER, ALPHA, LESS THAN ZERO.", &c__2, &c__1, (
	    ftnlen)6, (ftnlen)4, (ftnlen)29);
    return 0;
L720:
    xermsg_("SLATEC", "BESJ", "N LESS THAN ONE.", &c__2, &c__1, (ftnlen)6, (
	    ftnlen)4, (ftnlen)16);
    return 0;
L730:
    xermsg_("SLATEC", "BESJ", "X LESS THAN ZERO.", &c__2, &c__1, (ftnlen)6, (
	    ftnlen)4, (ftnlen)17);
    return 0;
} /* besj_ */

