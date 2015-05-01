/* cunik.f -- translated by f2c (version 12.02.01).
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

static integer c__1 = 1;

/* DECK CUNIK */
/* Subroutine */ int cunik_(complex *zr, real *fnu, integer *ikflg, integer *
	ipmtr, real *tol, integer *init, complex *phi, complex *zeta1, 
	complex *zeta2, complex *sum, complex *cwrk)
{
    /* Initialized data */

    static complex czero = {0.f,0.f};
    static complex cone = {1.f,0.f};
    static complex con[2] = { {.398942280401432678f,0.f},{
	    1.25331413731550025f,0.f} };
    static real c__[120] = { 1.f,-.208333333333333333f,.125f,
	    .334201388888888889f,-.401041666666666667f,.0703125f,
	    -1.02581259645061728f,1.84646267361111111f,-.8912109375f,
	    .0732421875f,4.66958442342624743f,-11.2070026162229938f,
	    8.78912353515625f,-2.3640869140625f,.112152099609375f,
	    -28.2120725582002449f,84.6362176746007346f,-91.8182415432400174f,
	    42.5349987453884549f,-7.3687943594796317f,.227108001708984375f,
	    212.570130039217123f,-765.252468141181642f,1059.99045252799988f,
	    -699.579627376132541f,218.19051174421159f,-26.4914304869515555f,
	    .572501420974731445f,-1919.457662318407f,8061.72218173730938f,
	    -13586.5500064341374f,11655.3933368645332f,-5305.64697861340311f,
	    1200.90291321635246f,-108.090919788394656f,1.7277275025844574f,
	    20204.2913309661486f,-96980.5983886375135f,192547.001232531532f,
	    -203400.177280415534f,122200.46498301746f,-41192.6549688975513f,
	    7109.51430248936372f,-493.915304773088012f,6.07404200127348304f,
	    -242919.187900551333f,1311763.6146629772f,-2998015.91853810675f,
	    3763271.297656404f,-2813563.22658653411f,1268365.27332162478f,
	    -331645.172484563578f,45218.7689813627263f,-2499.83048181120962f,
	    24.3805296995560639f,3284469.85307203782f,-19706819.1184322269f,
	    50952602.4926646422f,-74105148.2115326577f,66344512.2747290267f,
	    -37567176.6607633513f,13288767.1664218183f,-2785618.12808645469f,
	    308186.404612662398f,-13886.0897537170405f,110.017140269246738f,
	    -49329253.664509962f,325573074.185765749f,-939462359.681578403f,
	    1553596899.57058006f,-1621080552.10833708f,1106842816.82301447f,
	    -495889784.275030309f,142062907.797533095f,-24474062.7257387285f,
	    2243768.17792244943f,-84005.4336030240853f,551.335896122020586f,
	    814789096.118312115f,-5866481492.05184723f,18688207509.2958249f,
	    -34632043388.1587779f,41280185579.753974f,-33026599749.8007231f,
	    17954213731.1556001f,-6563293792.61928433f,1559279864.87925751f,
	    -225105661.889415278f,17395107.5539781645f,-549842.327572288687f,
	    3038.09051092238427f,-14679261247.6956167f,114498237732.02581f,
	    -399096175224.466498f,819218669548.577329f,-1098375156081.22331f,
	    1008158106865.38209f,-645364869245.376503f,287900649906.150589f,
	    -87867072178.0232657f,17634730606.8349694f,-2167164983.22379509f,
	    143157876.718888981f,-3871833.44257261262f,18257.7554742931747f,
	    286464035717.679043f,-2406297900028.50396f,9109341185239.89896f,
	    -20516899410934.4374f,30565125519935.3206f,-31667088584785.1584f,
	    23348364044581.8409f,-12320491305598.2872f,4612725780849.13197f,
	    -1196552880196.1816f,205914503232.410016f,-21822927757.5292237f,
	    1247009293.51271032f,-29188388.1222208134f,118838.426256783253f };

    /* System generated locals */
    integer i__1, i__2;
    real r__1;
    complex q__1, q__2, q__3;

    /* Local variables */
    static integer i__, j, k, l;
    static complex s, t, t2;
    static real ac;
    static complex sr, zn, cfn;
    static real rfn;
    static complex crfn;
    static real test, tsti, tstr;
    extern doublereal r1mach_(integer *);

/* ***BEGIN PROLOGUE  CUNIK */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to CBESI and CBESK */
/* ***LIBRARY   SLATEC */
/* ***TYPE      ALL (CUNIK-A, ZUNIK-A) */
/* ***AUTHOR  Amos, D. E., (SNL) */
/* ***DESCRIPTION */

/*        CUNIK COMPUTES PARAMETERS FOR THE UNIFORM ASYMPTOTIC */
/*        EXPANSIONS OF THE I AND K FUNCTIONS ON IKFLG= 1 OR 2 */
/*        RESPECTIVELY BY */

/*        W(FNU,ZR) = PHI*EXP(ZETA)*SUM */

/*        WHERE       ZETA=-ZETA1 + ZETA2       OR */
/*                          ZETA1 - ZETA2 */

/*        THE FIRST CALL MUST HAVE INIT=0. SUBSEQUENT CALLS WITH THE */
/*        SAME ZR AND FNU WILL RETURN THE I OR K FUNCTION ON IKFLG= */
/*        1 OR 2 WITH NO CHANGE IN INIT. CWRK IS A COMPLEX WORK */
/*        ARRAY. IPMTR=0 COMPUTES ALL PARAMETERS. IPMTR=1 COMPUTES PHI, */
/*        ZETA1,ZETA2. */

/* ***SEE ALSO  CBESI, CBESK */
/* ***ROUTINES CALLED  R1MACH */
/* ***REVISION HISTORY  (YYMMDD) */
/*   830501  DATE WRITTEN */
/*   910415  Prologue converted to Version 4.0 format.  (BAB) */
/* ***END PROLOGUE  CUNIK */
    /* Parameter adjustments */
    --cwrk;

    /* Function Body */
/* ***FIRST EXECUTABLE STATEMENT  CUNIK */
    if (*init != 0) {
	goto L40;
    }
/* ----------------------------------------------------------------------- */
/*     INITIALIZE ALL VARIABLES */
/* ----------------------------------------------------------------------- */
    rfn = 1.f / *fnu;
    q__1.r = rfn, q__1.i = 0.f;
    crfn.r = q__1.r, crfn.i = q__1.i;
/*     T = ZR*CRFN */
/* ----------------------------------------------------------------------- */
/*     OVERFLOW TEST (ZR/FNU TOO SMALL) */
/* ----------------------------------------------------------------------- */
    tstr = zr->r;
    tsti = r_imag(zr);
    test = r1mach_(&c__1) * 1e3f;
    ac = *fnu * test;
    if (dabs(tstr) > ac || dabs(tsti) > ac) {
	goto L15;
    }
    ac = (r__1 = log(test), dabs(r__1)) * 2.f + *fnu;
    q__1.r = ac, q__1.i = 0.f;
    zeta1->r = q__1.r, zeta1->i = q__1.i;
    q__1.r = *fnu, q__1.i = 0.f;
    zeta2->r = q__1.r, zeta2->i = q__1.i;
    phi->r = cone.r, phi->i = cone.i;
    return 0;
L15:
    q__1.r = zr->r * crfn.r - zr->i * crfn.i, q__1.i = zr->r * crfn.i + zr->i 
	    * crfn.r;
    t.r = q__1.r, t.i = q__1.i;
    q__2.r = t.r * t.r - t.i * t.i, q__2.i = t.r * t.i + t.i * t.r;
    q__1.r = cone.r + q__2.r, q__1.i = cone.i + q__2.i;
    s.r = q__1.r, s.i = q__1.i;
    c_sqrt(&q__1, &s);
    sr.r = q__1.r, sr.i = q__1.i;
    q__1.r = *fnu, q__1.i = 0.f;
    cfn.r = q__1.r, cfn.i = q__1.i;
    q__2.r = cone.r + sr.r, q__2.i = cone.i + sr.i;
    c_div(&q__1, &q__2, &t);
    zn.r = q__1.r, zn.i = q__1.i;
    c_log(&q__2, &zn);
    q__1.r = cfn.r * q__2.r - cfn.i * q__2.i, q__1.i = cfn.r * q__2.i + cfn.i 
	    * q__2.r;
    zeta1->r = q__1.r, zeta1->i = q__1.i;
    q__1.r = cfn.r * sr.r - cfn.i * sr.i, q__1.i = cfn.r * sr.i + cfn.i * 
	    sr.r;
    zeta2->r = q__1.r, zeta2->i = q__1.i;
    c_div(&q__1, &cone, &sr);
    t.r = q__1.r, t.i = q__1.i;
    q__1.r = t.r * crfn.r - t.i * crfn.i, q__1.i = t.r * crfn.i + t.i * 
	    crfn.r;
    sr.r = q__1.r, sr.i = q__1.i;
    c_sqrt(&q__1, &sr);
    cwrk[16].r = q__1.r, cwrk[16].i = q__1.i;
    i__1 = *ikflg - 1;
    q__1.r = cwrk[16].r * con[i__1].r - cwrk[16].i * con[i__1].i, q__1.i = 
	    cwrk[16].r * con[i__1].i + cwrk[16].i * con[i__1].r;
    phi->r = q__1.r, phi->i = q__1.i;
    if (*ipmtr != 0) {
	return 0;
    }
    c_div(&q__1, &cone, &s);
    t2.r = q__1.r, t2.i = q__1.i;
    cwrk[1].r = cone.r, cwrk[1].i = cone.i;
    crfn.r = cone.r, crfn.i = cone.i;
    ac = 1.f;
    l = 1;
    for (k = 2; k <= 15; ++k) {
	s.r = czero.r, s.i = czero.i;
	i__1 = k;
	for (j = 1; j <= i__1; ++j) {
	    ++l;
	    q__2.r = s.r * t2.r - s.i * t2.i, q__2.i = s.r * t2.i + s.i * 
		    t2.r;
	    i__2 = l - 1;
	    q__3.r = c__[i__2], q__3.i = 0.f;
	    q__1.r = q__2.r + q__3.r, q__1.i = q__2.i + q__3.i;
	    s.r = q__1.r, s.i = q__1.i;
/* L10: */
	}
	q__1.r = crfn.r * sr.r - crfn.i * sr.i, q__1.i = crfn.r * sr.i + 
		crfn.i * sr.r;
	crfn.r = q__1.r, crfn.i = q__1.i;
	i__1 = k;
	q__1.r = crfn.r * s.r - crfn.i * s.i, q__1.i = crfn.r * s.i + crfn.i *
		 s.r;
	cwrk[i__1].r = q__1.r, cwrk[i__1].i = q__1.i;
	ac *= rfn;
	i__1 = k;
	tstr = cwrk[i__1].r;
	tsti = r_imag(&cwrk[k]);
	test = dabs(tstr) + dabs(tsti);
	if (ac < *tol && test < *tol) {
	    goto L30;
	}
/* L20: */
    }
    k = 15;
L30:
    *init = k;
L40:
    if (*ikflg == 2) {
	goto L60;
    }
/* ----------------------------------------------------------------------- */
/*     COMPUTE SUM FOR THE I FUNCTION */
/* ----------------------------------------------------------------------- */
    s.r = czero.r, s.i = czero.i;
    i__1 = *init;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = i__;
	q__1.r = s.r + cwrk[i__2].r, q__1.i = s.i + cwrk[i__2].i;
	s.r = q__1.r, s.i = q__1.i;
/* L50: */
    }
    sum->r = s.r, sum->i = s.i;
    q__1.r = cwrk[16].r * con[0].r - cwrk[16].i * con[0].i, q__1.i = cwrk[16]
	    .r * con[0].i + cwrk[16].i * con[0].r;
    phi->r = q__1.r, phi->i = q__1.i;
    return 0;
L60:
/* ----------------------------------------------------------------------- */
/*     COMPUTE SUM FOR THE K FUNCTION */
/* ----------------------------------------------------------------------- */
    s.r = czero.r, s.i = czero.i;
    t.r = cone.r, t.i = cone.i;
    i__1 = *init;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = i__;
	q__2.r = t.r * cwrk[i__2].r - t.i * cwrk[i__2].i, q__2.i = t.r * cwrk[
		i__2].i + t.i * cwrk[i__2].r;
	q__1.r = s.r + q__2.r, q__1.i = s.i + q__2.i;
	s.r = q__1.r, s.i = q__1.i;
	q__1.r = -t.r, q__1.i = -t.i;
	t.r = q__1.r, t.i = q__1.i;
/* L70: */
    }
    sum->r = s.r, sum->i = s.i;
    q__1.r = cwrk[16].r * con[1].r - cwrk[16].i * con[1].i, q__1.i = cwrk[16]
	    .r * con[1].i + cwrk[16].i * con[1].r;
    phi->r = q__1.r, phi->i = q__1.i;
    return 0;
} /* cunik_ */

