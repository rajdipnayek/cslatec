/* asyjy.f -- translated by f2c (version 12.02.01).
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
static integer c__5 = 5;
static integer c__12 = 12;
static integer c__11 = 11;
static integer c__1 = 1;

/* DECK ASYJY */
/* Subroutine */ int asyjy_(S_fp funjy, real *x, real *fnu, real *flgjy, 
	integer *in, real *y, real *wk, integer *iflw)
{
    /* Initialized data */

    static real con2 = .333333333333333f;
    static real con548 = .104166666666667f;
    static real ar[8] = { .0835503472222222f,.128226574556327f,
	    .29184902646414f,.881627267443758f,3.32140828186277f,
	    14.9957629868626f,78.9230130115865f,474.451538868264f };
    static real br[10] = { -.145833333333333f,-.0987413194444444f,
	    -.143312053915895f,-.317227202678414f,-.94242914795712f,
	    -3.51120304082635f,-15.727263620368f,-82.2814390971859f,
	    -492.355370523671f,-3316.21856854797f };
    static real c__[65] = { -.208333333333333f,.125f,.334201388888889f,
	    -.401041666666667f,.0703125f,-1.02581259645062f,1.84646267361111f,
	    -.8912109375f,.0732421875f,4.66958442342625f,-11.207002616223f,
	    8.78912353515625f,-2.3640869140625f,.112152099609375f,
	    -28.2120725582002f,84.6362176746007f,-91.81824154324f,
	    42.5349987453885f,-7.36879435947963f,.227108001708984f,
	    212.570130039217f,-765.252468141182f,1059.990452528f,
	    -699.579627376133f,218.190511744212f,-26.4914304869516f,
	    .572501420974731f,-1919.45766231841f,8061.72218173731f,
	    -13586.5500064341f,11655.3933368645f,-5305.6469786134f,
	    1200.90291321635f,-108.090919788395f,1.72772750258446f,
	    20204.2913309661f,-96980.5983886375f,192547.001232532f,
	    -203400.177280416f,122200.464983017f,-41192.6549688976f,
	    7109.51430248936f,-493.915304773088f,6.07404200127348f,
	    -242919.187900551f,1311763.61466298f,-2998015.91853811f,
	    3763271.2976564f,-2813563.22658653f,1268365.27332162f,
	    -331645.172484564f,45218.7689813627f,-2499.83048181121f,
	    24.3805296995561f,3284469.85307204f,-19706819.1184322f,
	    50952602.4926646f,-74105148.2115327f,66344512.274729f,
	    -37567176.6607634f,13288767.1664218f,-2785618.12808645f,
	    308186.404612662f,-13886.089753717f,110.017140269247f };
    static real gama[26] = { .629960524947437f,.251984209978975f,
	    .154790300415656f,.110713062416159f,.0857309395527395f,
	    .0697161316958684f,.0586085671893714f,.0504698873536311f,
	    .0442600580689155f,.039372066154351f,.0354283195924455f,
	    .0321818857502098f,.0294646240791158f,.0271581677112934f,
	    .0251768272973862f,.0234570755306079f,.0219508390134907f,
	    .0206210828235646f,.0194388240897881f,.0183810633800683f,
	    .0174293213231963f,.0165685837786612f,.0157865285987918f,
	    .0150729501494096f,.0144193250839955f,.0138184805735342f };
    static real tols = -6.90775527898214f;
    static real con1 = .666666666666667f;
    static struct {
	real e_1[104];
	} equiv_1 = { -.00444444444444444f, -9.22077922077922e-4f, 
		-8.84892884892885e-5f, 1.6592768783245e-4f, 
		2.46691372741793e-4f, 2.65995589346255e-4f, 
		2.61824297061501e-4f, 2.48730437344656e-4f, 
		2.32721040083232e-4f, 2.16362485712365e-4f, 
		2.00738858762752e-4f, 1.86267636637545e-4f, 
		1.73060775917876e-4f, 1.61091705929016e-4f, 
		1.50274774160908e-4f, 1.4050349739127e-4f, 
		1.31668816545923e-4f, 1.23667445598253e-4f, 
		1.16405271474738e-4f, 1.09798298372713e-4f, 
		1.03772410422993e-4f, 9.82626078369363e-5f, 
		9.32120517249503e-5f, 8.85710852478712e-5f, 
		8.429631057157e-5f, 8.03497548407791e-5f, 
		6.93735541354589e-4f, 2.32241745182922e-4f, 
		-1.41986273556691e-5f, -1.16444931672049e-4f, 
		-1.50803558053049e-4f, -1.55121924918096e-4f, 
		-1.46809756646466e-4f, -1.33815503867491e-4f, 
		-1.19744975684254e-4f, -1.06184319207974e-4f, 
		-9.37699549891194e-5f, -8.26923045588193e-5f, 
		-7.29374348155221e-5f, -6.44042357721016e-5f, 
		-5.69611566009369e-5f, -5.04731044303562e-5f, 
		-4.48134868008883e-5f, -3.98688727717599e-5f, 
		-3.55400532972042e-5f, -3.17414256609022e-5f, 
		-2.83996793904175e-5f, -2.54522720634871e-5f, 
		-2.28459297164725e-5f, -2.05352753106481e-5f, 
		-1.84816217627666e-5f, -1.66519330021394e-5f, 
		-3.54211971457744e-4f, -1.56161263945159e-4f, 
		3.04465503594936e-5f, 1.30198655773243e-4f, 
		1.67471106699712e-4f, 1.70222587683593e-4f, 
		1.56501427608595e-4f, 1.36339170977445e-4f, 
		1.14886692029825e-4f, 9.45869093034688e-5f, 
		7.64498419250898e-5f, 6.07570334965197e-5f, 
		4.74394299290509e-5f, 3.62757512005344e-5f, 
		2.69939714979225e-5f, 1.93210938247939e-5f, 
		1.30056674793963e-5f, 7.82620866744497e-6f, 
		3.59257485819352e-6f, 1.44040049814252e-7f, 
		-2.65396769697939e-6f, -4.91346867098486e-6f, 
		-6.72739296091248e-6f, -8.17269379678658e-6f, 
		-9.31304715093561e-6f, -1.02011418798016e-5f, 
		3.78194199201773e-4f, 2.02471952761816e-4f, 
		-6.37938506318862e-5f, -2.38598230603006e-4f, 
		-3.10916256027362e-4f, -3.13680115247576e-4f, 
		-2.78950273791323e-4f, -2.28564082619141e-4f, 
		-1.75245280340847e-4f, -1.2554406306069e-4f, 
		-8.22982872820208e-5f, -4.62860730588116e-5f, 
		-1.72334302366962e-5f, 5.60690482304602e-6f, 
		2.31395443148287e-5f, 3.62642745856794e-5f, 
		4.58006124490189e-5f, 5.24595294959114e-5f, 
		5.68396208545815e-5f, 5.94349820393104e-5f, 
		6.06478527578422e-5f, 6.08023907788436e-5f, 
		6.0157789453946e-5f, 5.89199657344698e-5f, 
		5.72515823777593e-5f, 5.52804375585853e-5f };

    static struct {
	real e_1[130];
	} equiv_4 = { .0179988721413553f, .00559964911064388f, 
		.00288501402231133f, .00180096606761054f, .00124753110589199f,
		 9.22878876572938e-4f, 7.14430421727287e-4f, 
		5.71787281789705e-4f, 4.69431007606482e-4f, 
		3.93232835462917e-4f, 3.34818889318298e-4f, 
		2.88952148495752e-4f, 2.52211615549573e-4f, 
		2.22280580798883e-4f, 1.97541838033063e-4f, 
		1.76836855019718e-4f, 1.59316899661821e-4f, 
		1.44347930197334e-4f, 1.31448068119965e-4f, 
		1.20245444949303e-4f, 1.10449144504599e-4f, 
		1.01828770740567e-4f, 9.41998224204238e-5f, 
		8.74130545753834e-5f, 8.13466262162801e-5f, 
		7.59002269646219e-5f, -.00149282953213429f, 
		-8.78204709546389e-4f, -5.02916549572035e-4f, 
		-2.94822138512746e-4f, -1.75463996970783e-4f, 
		-1.04008550460816e-4f, -5.96141953046458e-5f, 
		-3.12038929076098e-5f, -1.2608973598023e-5f, 
		-2.4289260857573e-7f, 8.05996165414274e-6f, 
		1.36507009262147e-5f, 1.73964125472926e-5f, 
		1.98672978842134e-5f, 2.14463263790823e-5f, 
		2.23954659232457e-5f, 2.28967783814713e-5f, 
		2.30785389811178e-5f, 2.30321976080909e-5f, 
		2.28236073720349e-5f, 2.25005881105292e-5f, 
		2.20981015361991e-5f, 2.16418427448104e-5f, 
		2.11507649256221e-5f, 2.06388749782171e-5f, 
		2.01165241997082e-5f, 5.52213076721293e-4f, 
		4.47932581552385e-4f, 2.79520653992021e-4f, 
		1.52468156198447e-4f, 6.93271105657044e-5f, 
		1.76258683069991e-5f, -1.35744996343269e-5f, 
		-3.17972413350427e-5f, -4.18861861696693e-5f, 
		-4.69004889379141e-5f, -4.87665447413787e-5f, 
		-4.87010031186735e-5f, -4.74755620890087e-5f, 
		-4.55813058138628e-5f, -4.33309644511266e-5f, 
		-4.0923019315775e-5f, -3.84822638603221e-5f, 
		-3.60857167535411e-5f, -3.37793306123367e-5f, 
		-3.1588856077211e-5f, -2.95269561750807e-5f, 
		-2.75978914828336e-5f, -2.58006174666884e-5f, 
		-2.4130835676128e-5f, -2.25823509518346e-5f, 
		-2.11479656768913e-5f, -4.7461779655996e-4f, 
		-4.77864567147321e-4f, -3.20390228067038e-4f, 
		-1.61105016119962e-4f, -4.25778101285435e-5f, 
		3.44571294294968e-5f, 7.97092684075675e-5f, 
		1.03138236708272e-4f, 1.12466775262204e-4f, 
		1.13103642108481e-4f, 1.08651634848774e-4f, 
		1.01437951597662e-4f, 9.29298396593364e-5f, 
		8.4029313301609e-5f, 7.52727991349134e-5f, 
		6.69632521975731e-5f, 5.92564547323195e-5f, 
		5.22169308826976e-5f, 4.58539485165361e-5f, 
		4.01445513891487e-5f, 3.50481730031328e-5f, 
		3.05157995034347e-5f, 2.64956119950516e-5f, 
		2.29363633690998e-5f, 1.97893056664022e-5f, 
		1.70091984636413e-5f, 7.36465810572578e-4f, 
		8.72790805146194e-4f, 6.22614862573135e-4f, 
		2.85998154194304e-4f, 3.84737672879366e-6f, 
		-1.87906003636972e-4f, -2.97603646594555e-4f, 
		-3.45998126832656e-4f, -3.53382470916038e-4f, 
		-3.35715635775049e-4f, -3.0432112478904e-4f, 
		-2.66722723047613e-4f, -2.2765421412282e-4f, 
		-1.89922611854562e-4f, -1.55058918599094e-4f, 
		-1.23778240761874e-4f, -9.62926147717644e-5f, 
		-7.25178327714425e-5f, -5.22070028895634e-5f, 
		-3.50347750511901e-5f, -2.06489761035552e-5f, 
		-8.70106096849767e-6f, 1.136986866751e-6f, 
		9.16426474122779e-6f, 1.56477785428873e-5f, 
		2.08223629482467e-5f };


    /* System generated locals */
    integer i__1, i__2, i__3;
    real r__1;
    doublereal d__1, d__2;

    /* Local variables */
    static integer i__, j, k, l;
    static real z__, s1, t2;
    static integer kb;
    static real fi, ap, cr[10], dr[10];
    static integer jn;
    static real fn, sa, az;
    static integer jr;
    static real sb;
    static integer ks, ju, lr;
    static real ta, tb, z32, xx, fn2;
    static integer kp1;
    static real dfi, akm, phi, tfn, tau, rcz, tol, rtz, abw2, rfn2;
    static integer ksp1, lrp1;
#define alfa ((real *)&equiv_1)
#define beta ((real *)&equiv_4)
    static real relb, elim, rden;
    static integer kmax[5];
    static real crz32, asum, bsum, suma, sumb, upol[10];
#define alfa1 ((real *)&equiv_1)
#define alfa2 ((real *)&equiv_1 + 52)
#define beta1 ((real *)&equiv_4)
#define beta2 ((real *)&equiv_4 + 52)
#define beta3 ((real *)&equiv_4 + 104)
    static integer iseta, isetb, klast;
    static real rzden;
    extern integer i1mach_(integer *);
    extern doublereal r1mach_(integer *);
    static integer kstemp;

/* ***BEGIN PROLOGUE  ASYJY */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to BESJ and BESY */
/* ***LIBRARY   SLATEC */
/* ***TYPE      SINGLE PRECISION (ASYJY-S, DASYJY-D) */
/* ***AUTHOR  Amos, D. E., (SNLA) */
/* ***DESCRIPTION */

/*                 ASYJY computes Bessel functions J and Y */
/*               for arguments X.GT.0.0 and orders FNU.GE.35.0 */
/*               on FLGJY = 1 and FLGJY = -1 respectively */

/*                                  INPUT */

/*      FUNJY - external function JAIRY or YAIRY */
/*          X - argument, X.GT.0.0E0 */
/*        FNU - order of the first Bessel function */
/*      FLGJY - selection flag */
/*              FLGJY =  1.0E0 gives the J function */
/*              FLGJY = -1.0E0 gives the Y function */
/*         IN - number of functions desired, IN = 1 or 2 */

/*                                  OUTPUT */

/*         Y  - a vector whose first in components contain the sequence */
/*       IFLW - a flag indicating underflow or overflow */
/*                    return variables for BESJ only */
/*      WK(1) = 1 - (X/FNU)**2 = W**2 */
/*      WK(2) = SQRT(ABS(WK(1))) */
/*      WK(3) = ABS(WK(2) - ATAN(WK(2)))  or */
/*              ABS(LN((1 + WK(2))/(X/FNU)) - WK(2)) */
/*            = ABS((2/3)*ZETA**(3/2)) */
/*      WK(4) = FNU*WK(3) */
/*      WK(5) = (1.5*WK(3)*FNU)**(1/3) = SQRT(ZETA)*FNU**(1/3) */
/*      WK(6) = SIGN(1.,W**2)*WK(5)**2 = SIGN(1.,W**2)*ZETA*FNU**(2/3) */
/*      WK(7) = FNU**(1/3) */

/*     Abstract */
/*         ASYJY implements the uniform asymptotic expansion of */
/*         the J and Y Bessel functions for FNU.GE.35 and real */
/*         X.GT.0.0E0. The forms are identical except for a change */
/*         in sign of some of the terms. This change in sign is */
/*         accomplished by means of the flag FLGJY = 1 or -1. On */
/*         FLGJY = 1 the AIRY functions AI(X) and DAI(X) are */
/*         supplied by the external function JAIRY, and on */
/*         FLGJY = -1 the AIRY functions BI(X) and DBI(X) are */
/*         supplied by the external function YAIRY. */

/* ***SEE ALSO  BESJ, BESY */
/* ***ROUTINES CALLED  I1MACH, R1MACH */
/* ***REVISION HISTORY  (YYMMDD) */
/*   750101  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   891009  Removed unreferenced variable.  (WRB) */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900328  Added TYPE section.  (WRB) */
/*   910408  Updated the AUTHOR section.  (WRB) */
/* ***END PROLOGUE  ASYJY */
    /* Parameter adjustments */
    --wk;
    --y;

    /* Function Body */
/* ***FIRST EXECUTABLE STATEMENT  ASYJY */
    ta = r1mach_(&c__3);
    tol = dmax(ta,1e-15f);
    tb = r1mach_(&c__5);
    ju = i1mach_(&c__12);
    if (*flgjy == 1.f) {
	goto L6;
    }
    jr = i1mach_(&c__11);
    elim = tb * -2.303f * (ju + jr);
    goto L7;
L6:
    elim = (tb * ju + 3.f) * -2.303f;
L7:
    fn = *fnu;
    *iflw = 0;
    i__1 = *in;
    for (jn = 1; jn <= i__1; ++jn) {
	xx = *x / fn;
	wk[1] = 1.f - xx * xx;
	abw2 = dabs(wk[1]);
	wk[2] = sqrt(abw2);
	d__1 = (doublereal) fn;
	d__2 = (doublereal) con2;
	wk[7] = pow_dd(&d__1, &d__2);
	if (abw2 > .2775f) {
	    goto L80;
	}

/*     ASYMPTOTIC EXPANSION */
/*     CASES NEAR X=FN, ABS(1.-(X/FN)**2).LE.0.2775 */
/*     COEFFICIENTS OF ASYMPTOTIC EXPANSION BY SERIES */

/*     ZETA AND TRUNCATION FOR A(ZETA) AND B(ZETA) SERIES */

/*     KMAX IS TRUNCATION INDEX FOR A(ZETA) AND B(ZETA) SERIES=MAX(2,SA) */

	sa = 0.f;
	if (abw2 == 0.f) {
	    goto L10;
	}
	sa = tols / log(abw2);
L10:
	sb = sa;
	for (i__ = 1; i__ <= 5; ++i__) {
	    akm = dmax(sa,2.f);
	    kmax[i__ - 1] = (integer) akm;
	    sa += sb;
/* L20: */
	}
	kb = kmax[4];
	klast = kb - 1;
	sa = gama[kb - 1];
	i__2 = klast;
	for (k = 1; k <= i__2; ++k) {
	    --kb;
	    sa = sa * wk[1] + gama[kb - 1];
/* L30: */
	}
	z__ = wk[1] * sa;
	az = dabs(z__);
	rtz = sqrt(az);
	wk[3] = con1 * az * rtz;
	wk[4] = wk[3] * fn;
	wk[5] = rtz * wk[7];
	wk[6] = -wk[5] * wk[5];
	if (z__ <= 0.f) {
	    goto L35;
	}
	if (wk[4] > elim) {
	    goto L75;
	}
	wk[6] = -wk[6];
L35:
	phi = sqrt(sqrt(sa + sa + sa + sa));

/*     B(ZETA) FOR S=0 */

	kb = kmax[4];
	klast = kb - 1;
	sb = beta[kb - 1];
	i__2 = klast;
	for (k = 1; k <= i__2; ++k) {
	    --kb;
	    sb = sb * wk[1] + beta[kb - 1];
/* L40: */
	}
	ksp1 = 1;
	fn2 = fn * fn;
	rfn2 = 1.f / fn2;
	rden = 1.f;
	asum = 1.f;
	relb = tol * dabs(sb);
	bsum = sb;
	for (ks = 1; ks <= 4; ++ks) {
	    ++ksp1;
	    rden *= rfn2;

/*     A(ZETA) AND B(ZETA) FOR S=1,2,3,4 */

	    kstemp = 5 - ks;
	    kb = kmax[kstemp - 1];
	    klast = kb - 1;
	    sa = alfa[kb + ks * 26 - 27];
	    sb = beta[kb + ksp1 * 26 - 27];
	    i__2 = klast;
	    for (k = 1; k <= i__2; ++k) {
		--kb;
		sa = sa * wk[1] + alfa[kb + ks * 26 - 27];
		sb = sb * wk[1] + beta[kb + ksp1 * 26 - 27];
/* L50: */
	    }
	    ta = sa * rden;
	    tb = sb * rden;
	    asum += ta;
	    bsum += tb;
	    if (dabs(ta) <= tol && dabs(tb) <= relb) {
		goto L70;
	    }
/* L60: */
	}
L70:
	bsum /= fn * wk[7];
	goto L160;

L75:
	*iflw = 1;
	return 0;

L80:
	upol[0] = 1.f;
	tau = 1.f / wk[2];
	t2 = 1.f / wk[1];
	if (wk[1] >= 0.f) {
	    goto L90;
	}

/*     CASES FOR (X/FN).GT.SQRT(1.2775) */

	wk[3] = (r__1 = wk[2] - atan(wk[2]), dabs(r__1));
	wk[4] = wk[3] * fn;
	rcz = -con1 / wk[4];
	z32 = wk[3] * 1.5f;
	d__1 = (doublereal) z32;
	d__2 = (doublereal) con2;
	rtz = pow_dd(&d__1, &d__2);
	wk[5] = rtz * wk[7];
	wk[6] = -wk[5] * wk[5];
	goto L100;
L90:

/*     CASES FOR (X/FN).LT.SQRT(0.7225) */

	wk[3] = (r__1 = log((wk[2] + 1.f) / xx) - wk[2], dabs(r__1));
	wk[4] = wk[3] * fn;
	rcz = con1 / wk[4];
	if (wk[4] > elim) {
	    goto L75;
	}
	z32 = wk[3] * 1.5f;
	d__1 = (doublereal) z32;
	d__2 = (doublereal) con2;
	rtz = pow_dd(&d__1, &d__2);
	d__1 = (doublereal) fn;
	d__2 = (doublereal) con2;
	wk[7] = pow_dd(&d__1, &d__2);
	wk[5] = rtz * wk[7];
	wk[6] = wk[5] * wk[5];
L100:
	phi = sqrt((rtz + rtz) * tau);
	tb = 1.f;
	asum = 1.f;
	tfn = tau / fn;
	rden = 1.f / fn;
	rfn2 = rden * rden;
	rden = 1.f;
	upol[1] = (c__[0] * t2 + c__[1]) * tfn;
	crz32 = con548 * rcz;
	bsum = upol[1] + crz32;
	relb = tol * dabs(bsum);
	ap = tfn;
	ks = 0;
	kp1 = 2;
	rzden = rcz;
	l = 2;
	iseta = 0;
	isetb = 0;
	for (lr = 2; lr <= 8; lr += 2) {

/*     COMPUTE TWO U POLYNOMIALS FOR NEXT A(ZETA) AND B(ZETA) */

	    lrp1 = lr + 1;
	    i__2 = lrp1;
	    for (k = lr; k <= i__2; ++k) {
		++ks;
		++kp1;
		++l;
		s1 = c__[l - 1];
		i__3 = kp1;
		for (j = 2; j <= i__3; ++j) {
		    ++l;
		    s1 = s1 * t2 + c__[l - 1];
/* L110: */
		}
		ap *= tfn;
		upol[kp1 - 1] = ap * s1;
		cr[ks - 1] = br[ks - 1] * rzden;
		rzden *= rcz;
		dr[ks - 1] = ar[ks - 1] * rzden;
/* L120: */
	    }
	    suma = upol[lrp1 - 1];
	    sumb = upol[lr + 1] + upol[lrp1 - 1] * crz32;
	    ju = lrp1;
	    i__2 = lr;
	    for (jr = 1; jr <= i__2; ++jr) {
		--ju;
		suma += cr[jr - 1] * upol[ju - 1];
		sumb += dr[jr - 1] * upol[ju - 1];
/* L130: */
	    }
	    rden *= rfn2;
	    tb = -tb;
	    if (wk[1] > 0.f) {
		tb = dabs(tb);
	    }
	    if (rden < tol) {
		goto L131;
	    }
	    asum += suma * tb;
	    bsum += sumb * tb;
	    goto L140;
L131:
	    if (iseta == 1) {
		goto L132;
	    }
	    if (dabs(suma) < tol) {
		iseta = 1;
	    }
	    asum += suma * tb;
L132:
	    if (isetb == 1) {
		goto L133;
	    }
	    if (dabs(sumb) < relb) {
		isetb = 1;
	    }
	    bsum += sumb * tb;
L133:
	    if (iseta == 1 && isetb == 1) {
		goto L150;
	    }
L140:
	    ;
	}
L150:
	tb = wk[5];
	if (wk[1] > 0.f) {
	    tb = -tb;
	}
	bsum /= tb;

L160:
	(*funjy)(&wk[6], &wk[5], &wk[4], &fi, &dfi);
	ta = 1.f / tol;
	tb = r1mach_(&c__1) * ta * 1e3f;
	if (dabs(fi) > tb) {
	    goto L165;
	}
	fi *= ta;
	dfi *= ta;
	phi *= tol;
L165:
	y[jn] = *flgjy * phi * (fi * asum + dfi * bsum) / wk[7];
	fn -= *flgjy;
/* L170: */
    }
    return 0;
} /* asyjy_ */

#undef beta3
#undef beta2
#undef beta1
#undef alfa2
#undef alfa1
#undef beta
#undef alfa


