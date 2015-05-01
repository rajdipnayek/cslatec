/* gamln.f -- translated by f2c (version 12.02.01).
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
static integer c__11 = 11;
static integer c__5 = 5;
static integer c__2 = 2;

/* DECK GAMLN */
doublereal gamln_(real *z__, integer *ierr)
{
    /* Initialized data */

    static real gln[100] = { 0.f,0.f,.693147180559945309f,1.791759469228055f,
	    3.17805383034794562f,4.78749174278204599f,6.579251212010101f,
	    8.5251613610654143f,10.6046029027452502f,12.8018274800814696f,
	    15.1044125730755153f,17.5023078458738858f,19.9872144956618861f,
	    22.5521638531234229f,25.1912211827386815f,27.8992713838408916f,
	    30.6718601060806728f,33.5050734501368889f,36.3954452080330536f,
	    39.339884187199494f,42.335616460753485f,45.380138898476908f,
	    48.4711813518352239f,51.6066755677643736f,54.7847293981123192f,
	    58.0036052229805199f,61.261701761002002f,64.5575386270063311f,
	    67.889743137181535f,71.257038967168009f,74.6582363488301644f,
	    78.0922235533153106f,81.5579594561150372f,85.0544670175815174f,
	    88.5808275421976788f,92.1361756036870925f,95.7196945421432025f,
	    99.3306124547874269f,102.968198614513813f,106.631760260643459f,
	    110.320639714757395f,114.034211781461703f,117.771881399745072f,
	    121.533081515438634f,125.317271149356895f,129.123933639127215f,
	    132.95257503561631f,136.802722637326368f,140.673923648234259f,
	    144.565743946344886f,148.477766951773032f,152.409592584497358f,
	    156.360836303078785f,160.331128216630907f,164.320112263195181f,
	    168.327445448427652f,172.352797139162802f,176.395848406997352f,
	    180.456291417543771f,184.533828861449491f,188.628173423671591f,
	    192.739047287844902f,196.866181672889994f,201.009316399281527f,
	    205.168199482641199f,209.342586752536836f,213.532241494563261f,
	    217.736934113954227f,221.956441819130334f,226.190548323727593f,
	    230.439043565776952f,234.701723442818268f,238.978389561834323f,
	    243.268849002982714f,247.572914096186884f,251.890402209723194f,
	    256.221135550009525f,260.564940971863209f,264.921649798552801f,
	    269.291097651019823f,273.673124285693704f,278.067573440366143f,
	    282.474292687630396f,286.893133295426994f,291.323950094270308f,
	    295.766601350760624f,300.220948647014132f,304.686856765668715f,
	    309.164193580146922f,313.652829949879062f,318.152639620209327f,
	    322.663499126726177f,327.185287703775217f,331.717887196928473f,
	    336.261181979198477f,340.815058870799018f,345.379407062266854f,
	    349.954118040770237f,354.539085519440809f,359.134205369575399f };
    static real cf[22] = { .0833333333333333333f,-.00277777777777777778f,
	    7.93650793650793651e-4f,-5.95238095238095238e-4f,
	    8.41750841750841751e-4f,-.00191752691752691753f,
	    .00641025641025641026f,-.0295506535947712418f,
	    .179644372368830573f,-1.39243221690590112f,13.402864044168392f,
	    -156.848284626002017f,2193.10333333333333f,-36108.7712537249894f,
	    691472.268851313067f,-15238221.5394074162f,382900751.391414141f,
	    -10882266035.7843911f,347320283765.002252f,-12369602142269.2745f,
	    488788064793079.335f,-21320333960919373.9f };
    static real con = 1.83787706640934548f;

    /* System generated locals */
    integer i__1;
    real ret_val;

    /* Local variables */
    static integer i__, k;
    static real s, t1, fz;
    static integer mz, nz;
    static real zm, zp;
    static integer i1m;
    static real fln, tlg, rln, trm, tst, zsq, zinc, zmin, zdmy, wdtol;
    extern integer i1mach_(integer *);
    extern doublereal r1mach_(integer *);

/* ***BEGIN PROLOGUE  GAMLN */
/* ***SUBSIDIARY */
/* ***PURPOSE  Compute the logarithm of the Gamma function */
/* ***LIBRARY   SLATEC */
/* ***CATEGORY  C7A */
/* ***TYPE      SINGLE PRECISION (GAMLN-S, DGAMLN-D) */
/* ***KEYWORDS  LOGARITHM OF GAMMA FUNCTION */
/* ***AUTHOR  Amos, D. E., (SNL) */
/* ***DESCRIPTION */

/*         GAMLN COMPUTES THE NATURAL LOG OF THE GAMMA FUNCTION FOR */
/*         Z.GT.0.  THE ASYMPTOTIC EXPANSION IS USED TO GENERATE VALUES */
/*         GREATER THAN ZMIN WHICH ARE ADJUSTED BY THE RECURSION */
/*         G(Z+1)=Z*G(Z) FOR Z.LE.ZMIN.  THE FUNCTION WAS MADE AS */
/*         PORTABLE AS POSSIBLE BY COMPUTING ZMIN FROM THE NUMBER OF BASE */
/*         10 DIGITS IN A WORD, RLN=MAX(-ALOG10(R1MACH(4)),0.5E-18) */
/*         LIMITED TO 18 DIGITS OF (RELATIVE) ACCURACY. */

/*         SINCE INTEGER ARGUMENTS ARE COMMON, A TABLE LOOK UP ON 100 */
/*         VALUES IS USED FOR SPEED OF EXECUTION. */

/*     DESCRIPTION OF ARGUMENTS */

/*         INPUT */
/*           Z      - REAL ARGUMENT, Z.GT.0.0E0 */

/*         OUTPUT */
/*           GAMLN  - NATURAL LOG OF THE GAMMA FUNCTION AT Z */
/*           IERR   - ERROR FLAG */
/*                    IERR=0, NORMAL RETURN, COMPUTATION COMPLETED */
/*                    IERR=1, Z.LE.0.0E0,    NO COMPUTATION */

/* ***REFERENCES  COMPUTATION OF BESSEL FUNCTIONS OF COMPLEX ARGUMENT */
/*                 BY D. E. AMOS, SAND83-0083, MAY, 1983. */
/* ***ROUTINES CALLED  I1MACH, R1MACH */
/* ***REVISION HISTORY  (YYMMDD) */
/*   830501  DATE WRITTEN */
/*   830501  REVISION DATE from Version 3.2 */
/*   910415  Prologue converted to Version 4.0 format.  (BAB) */
/*   920128  Category corrected.  (WRB) */
/*   921215  GAMLN defined for Z negative.  (WRB) */
/* ***END PROLOGUE  GAMLN */

/*           LNGAMMA(N), N=1,100 */
/*             COEFFICIENTS OF ASYMPTOTIC EXPANSION */

/*             LN(2*PI) */

/* ***FIRST EXECUTABLE STATEMENT  GAMLN */
    *ierr = 0;
    if (*z__ <= 0.f) {
	goto L70;
    }
    if (*z__ > 101.f) {
	goto L10;
    }
    nz = *z__;
    fz = *z__ - nz;
    if (fz > 0.f) {
	goto L10;
    }
    if (nz > 100) {
	goto L10;
    }
    ret_val = gln[nz - 1];
    return ret_val;
L10:
    wdtol = r1mach_(&c__4);
    wdtol = dmax(wdtol,5e-19f);
    i1m = i1mach_(&c__11);
    rln = r1mach_(&c__5) * i1m;
    fln = dmin(rln,20.f);
    fln = dmax(fln,3.f);
    fln += -3.f;
    zm = fln * .3875f + 1.8f;
    mz = zm + 1;
    zmin = (real) mz;
    zdmy = *z__;
    zinc = 0.f;
    if (*z__ >= zmin) {
	goto L20;
    }
    zinc = zmin - nz;
    zdmy = *z__ + zinc;
L20:
    zp = 1.f / zdmy;
    t1 = cf[0] * zp;
    s = t1;
    if (zp < wdtol) {
	goto L40;
    }
    zsq = zp * zp;
    tst = t1 * wdtol;
    for (k = 2; k <= 22; ++k) {
	zp *= zsq;
	trm = cf[k - 1] * zp;
	if (dabs(trm) < tst) {
	    goto L40;
	}
	s += trm;
/* L30: */
    }
L40:
    if (zinc != 0.f) {
	goto L50;
    }
    tlg = log(*z__);
    ret_val = *z__ * (tlg - 1.f) + (con - tlg) * .5f + s;
    return ret_val;
L50:
    zp = 1.f;
    nz = zinc;
    i__1 = nz;
    for (i__ = 1; i__ <= i__1; ++i__) {
	zp *= *z__ + (i__ - 1);
/* L60: */
    }
    tlg = log(zdmy);
    ret_val = zdmy * (tlg - 1.f) - log(zp) + (con - tlg) * .5f + s;
    return ret_val;


L70:
    ret_val = r1mach_(&c__2);
    *ierr = 1;
    return ret_val;
} /* gamln_ */

