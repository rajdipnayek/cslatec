/* daie.f -- translated by f2c (version 12.02.01).
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
static integer c__13 = 13;
static integer c__57 = 57;
static integer c__37 = 37;
static doublereal c_b7 = .3333;
static integer c__2 = 2;
static doublereal c_b9 = .6666;

/* DECK DAIE */
doublereal daie_(doublereal *x)
{
    /* Initialized data */

    static doublereal aifcs[13] = { -.037971358496669997496197089469414,
	    .059191888537263638574319728013777,
	    9.862928057727997536560389104406e-4,
	    6.8488438190765667554854830182412e-6,
	    2.5942025962194713019489279081403e-8,
	    6.1766127740813750329445749697236e-11,
	    1.0092454172466117901429556224601e-13,
	    1.2014792511179938141288033225333e-16,
	    1.0882945588716991878525295466666e-19,
	    7.75137721966848870392384e-23,
	    4.4548112037175638391466666666666e-26,
	    2.1092845231692343466666666666666e-29,
	    8.3701735910741333333333333333333e-33 };
    static doublereal aigcs[13] = { .018152365581161273011556209957864,
	    .021572563166010755534030638819968,
	    2.5678356987483249659052428090133e-4,
	    1.4265214119792403898829496921721e-6,
	    4.5721149200180426070434097558191e-9,
	    9.5251708435647098607392278840592e-12,
	    1.392563460577139905115042068619e-14,
	    1.5070999142762379592306991138666e-17,
	    1.2559148312567778822703205333333e-20,
	    8.3063073770821340343829333333333e-24,
	    4.4657538493718567445333333333333e-27,
	    1.9900855034518869333333333333333e-30,
	    7.4702885256533333333333333333333e-34 };
    static doublereal aip1cs[57] = { -.02146951858910538455460863467778,
	    -.007535382535043301166219720865565,
	    5.971527949026380852035388881994e-4,
	    -7.283251254207610648502368291548e-5,
	    1.11029713073929966651738182114e-5,
	    -1.950386152284405710346930314033e-6,
	    3.786973885159515193885319670057e-7,
	    -7.929675297350978279039072879154e-8,
	    1.762247638674256075568420122202e-8,
	    -4.110767539667195045029896593893e-9,
	    9.984770057857892247183414107544e-10,
	    -2.510093251387122211349867730034e-10,
	    6.500501929860695409272038601725e-11,
	    -1.727818405393616515478877107366e-11,
	    4.699378842824512578362292872307e-12,
	    -1.304675656297743914491241246272e-12,
	    3.689698478462678810473948382282e-13,
	    -1.061087206646806173650359679035e-13,
	    3.09841438487818743866021007011e-14,
	    -9.174908079824139307833423547851e-15,
	    2.752049140347210895693579062271e-15,
	    -8.35375011592204655809139330188e-16,
	    2.563931129357934947568636168612e-16,
	    -7.950633762598854983273747289822e-17,
	    2.489283634603069977437281175644e-17,
	    -7.864326933928735569664626221296e-18,
	    2.505687311439975672324470645019e-18,
	    -8.047420364163909524537958682241e-19,
	    2.604097118952053964443401104392e-19,
	    -8.486954164056412259482488834184e-20,
	    2.784706882142337843359429186027e-20,
	    -9.195858953498612913687224151354e-21,
	    3.055304318374238742247668225583e-21,
	    -1.021035455479477875902177048439e-21,
	    3.431118190743757844000555680836e-22,
	    -1.159129341797749513376922463109e-22,
	    3.935772844200255610836268229154e-23,
	    -1.342880980296717611956718989038e-23,
	    4.603287883520002741659190305314e-24,
	    -1.585043927004064227810772499387e-24,
	    5.481275667729675908925523755008e-25,
	    -1.903349371855047259064017948945e-25,
	    6.635682302374008716777612115968e-26,
	    -2.322311650026314307975200986453e-26,
	    8.157640113429179313142743695359e-27,
	    -2.875824240632900490057489929557e-27,
	    1.017329450942901435079714319018e-27,
	    -3.610879108742216446575703490559e-28,
	    1.285788540363993421256640342698e-28,
	    -4.592901037378547425160693022719e-29,
	    1.645597033820713725812102485333e-29,
	    -5.91342129984350184208792027136e-30,
	    2.131057006604993303479369509546e-30,
	    -7.701158157787598216982761745066e-31,
	    2.79053330796893041758178377728e-31,
	    -1.013807715111284006452241367039e-31,
	    3.692580158719624093658286216533e-32 };
    static doublereal aip2cs[37] = { -.00174314496929375513390355844011,
	    -.0016789385432554167163219061348,
	    3.59653403352166035885983858114e-5,
	    -1.380818602739228354573993831e-6,
	    7.41122807731505298848699095233e-8,
	    -5.00238203900133013130422866325e-9,
	    4.00693917417184240675446866355e-10,
	    -3.67331242795905044199318496207e-11,
	    3.76034439592373852439592002918e-12,
	    -4.22321332718747538026564938968e-13,
	    5.1350945403365707091961875412e-14,
	    -6.69095850390477595651681356676e-15,
	    9.26667545641290648239550724382e-16,
	    -1.35514382416070576333397356591e-16,
	    2.08115496312830995299006549335e-17,
	    -3.34116499159176856871277570256e-18,
	    5.58578584585924316868032946585e-19,
	    -9.69219040152365247518658209109e-20,
	    1.74045700128893206465696557738e-20,
	    -3.22640979731130400247846333098e-21,
	    6.16074471106625258533259618986e-22,
	    -1.20936347982490059076420676266e-22,
	    2.43632763310138108261570095786e-23,
	    -5.02914221497457468943403144533e-24,
	    1.06224175543635689495470626133e-24,
	    -2.29284284895989241509856324266e-25,
	    5.05181733929503744986884778666e-26,
	    -1.1349812371441240497979392e-26,2.59765565985606980698374144e-27,
	    -6.05124621542939506172231679999e-28,
	    1.43359777966772800720295253333e-28,
	    -3.45147757060899986280721066666e-29,
	    8.43875190213646740427025066666e-30,
	    -2.09396142298188169434453333333e-30,
	    5.27008873478945503182848e-31,
	    -1.34457433014553385789030399999e-31,
	    3.47570964526601147340117333333e-32 };
    static logical first = TRUE_;

    /* System generated locals */
    doublereal ret_val, d__1;

    /* Local variables */
    static doublereal z__, xm;
    static real eta;
    static integer naif, naig;
    static doublereal xbig;
    static integer naip1, naip2;
    static doublereal x3sml, theta, x32sml;
    extern doublereal d1mach_(integer *);
    static doublereal sqrtx;
    extern /* Subroutine */ int d9aimp_(doublereal *, doublereal *, 
	    doublereal *);
    extern doublereal dcsevl_(doublereal *, doublereal *, integer *);
    extern integer initds_(doublereal *, integer *, real *);

/* ***BEGIN PROLOGUE  DAIE */
/* ***PURPOSE  Calculate the Airy function for a negative argument and an */
/*            exponentially scaled Airy function for a non-negative */
/*            argument. */
/* ***LIBRARY   SLATEC (FNLIB) */
/* ***CATEGORY  C10D */
/* ***TYPE      DOUBLE PRECISION (AIE-S, DAIE-D) */
/* ***KEYWORDS  EXPONENTIALLY SCALED AIRY FUNCTION, FNLIB, */
/*             SPECIAL FUNCTIONS */
/* ***AUTHOR  Fullerton, W., (LANL) */
/* ***DESCRIPTION */

/* DAIE(X) calculates the Airy function or the exponentially scaled */
/* Airy function depending on the value of the argument.  The function */
/* and argument are both double precision. */

/* Evaluate AI(X) for X .LE. 0.0 and AI(X)*EXP(ZETA) where */
/* ZETA = 2/3 * X**(3/2)  for X .GE. 0.0 */

/* Series for AIF        on the interval -1.00000E+00 to  1.00000E+00 */
/*                                        with weighted error   8.37E-33 */
/*                                         log weighted error  32.08 */
/*                               significant figures required  30.87 */
/*                                    decimal places required  32.63 */

/* Series for AIG        on the interval -1.00000E+00 to  1.00000E+00 */
/*                                        with weighted error   7.47E-34 */
/*                                         log weighted error  33.13 */
/*                               significant figures required  31.50 */
/*                                    decimal places required  33.68 */

/* Series for AIP1       on the interval  1.25000E-01 to  1.00000E+00 */
/*                                        with weighted error   3.69E-32 */
/*                                         log weighted error  31.43 */
/*                               significant figures required  29.55 */
/*                                    decimal places required  32.31 */

/* Series for AIP2       on the interval  0.          to  1.25000E-01 */
/*                                        with weighted error   3.48E-32 */
/*                                         log weighted error  31.46 */
/*                               significant figures required  28.74 */
/*                                    decimal places required  32.24 */

/* ***REFERENCES  (NONE) */
/* ***ROUTINES CALLED  D1MACH, D9AIMP, DCSEVL, INITDS */
/* ***REVISION HISTORY  (YYMMDD) */
/*   770701  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890531  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   920618  Removed space from variable names.  (RWC, WRB) */
/* ***END PROLOGUE  DAIE */
/* ***FIRST EXECUTABLE STATEMENT  DAIE */
    if (first) {
	eta = (real) d1mach_(&c__3) * .1f;
	naif = initds_(aifcs, &c__13, &eta);
	naig = initds_(aigcs, &c__13, &eta);
	naip1 = initds_(aip1cs, &c__57, &eta);
	naip2 = initds_(aip2cs, &c__37, &eta);

	d__1 = (doublereal) eta;
	x3sml = pow_dd(&d__1, &c_b7);
/* Computing 2nd power */
	d__1 = x3sml;
	x32sml = d__1 * d__1 * 1.3104;
	d__1 = d1mach_(&c__2);
	xbig = pow_dd(&d__1, &c_b9);
    }
    first = FALSE_;

    if (*x >= -1.) {
	goto L20;
    }
    d9aimp_(x, &xm, &theta);
    ret_val = xm * cos(theta);
    return ret_val;

L20:
    if (*x > 1.) {
	goto L30;
    }
    z__ = 0.;
    if (abs(*x) > x3sml) {
/* Computing 3rd power */
	d__1 = *x;
	z__ = d__1 * (d__1 * d__1);
    }
    ret_val = dcsevl_(&z__, aifcs, &naif) - *x * (dcsevl_(&z__, aigcs, &naig) 
	    + .25) + .375;
    if (*x > x32sml) {
	ret_val *= exp(*x * 2. * sqrt(*x) / 3.);
    }
    return ret_val;

L30:
    if (*x > 4.) {
	goto L40;
    }
    sqrtx = sqrt(*x);
    z__ = (16. / (*x * sqrtx) - 9.) / 7.;
    ret_val = (dcsevl_(&z__, aip1cs, &naip1) + .28125) / sqrt(sqrtx);
    return ret_val;

L40:
    sqrtx = sqrt(*x);
    z__ = -1.;
    if (*x < xbig) {
	z__ = 16. / (*x * sqrtx) - 1.;
    }
    ret_val = (dcsevl_(&z__, aip2cs, &naip2) + .28125) / sqrt(sqrtx);
    return ret_val;

} /* daie_ */

