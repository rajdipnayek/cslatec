/* derfc.f -- translated by f2c (version 12.02.01).
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
static integer c__21 = 21;
static integer c__59 = 59;
static integer c__49 = 49;
static integer c__1 = 1;

/* DECK DERFC */
doublereal derfc_(doublereal *x)
{
    /* Initialized data */

    static doublereal erfcs[21] = { -.049046121234691808039984544033376,
	    -.14226120510371364237824741899631,
	    .010035582187599795575754676712933,
	    -5.7687646997674847650827025509167e-4,
	    2.7419931252196061034422160791471e-5,
	    -1.1043175507344507604135381295905e-6,
	    3.8488755420345036949961311498174e-8,
	    -1.1808582533875466969631751801581e-9,
	    3.2334215826050909646402930953354e-11,
	    -7.9910159470045487581607374708595e-13,
	    1.7990725113961455611967245486634e-14,
	    -3.7186354878186926382316828209493e-16,
	    7.1035990037142529711689908394666e-18,
	    -1.2612455119155225832495424853333e-19,
	    2.0916406941769294369170500266666e-21,
	    -3.253973102931407298236416e-23,
	    4.7668672097976748332373333333333e-25,
	    -6.5980120782851343155199999999999e-27,
	    8.6550114699637626197333333333333e-29,
	    -1.0788925177498064213333333333333e-30,
	    1.2811883993017002666666666666666e-32 };
    static doublereal erc2cs[49] = { -.06960134660230950112739150826197,
	    -.04110133936262089348982212084666,
	    .003914495866689626881561143705244,
	    -4.906395650548979161280935450774e-4,
	    7.157479001377036380760894141825e-5,
	    -1.153071634131232833808232847912e-5,
	    1.994670590201997635052314867709e-6,
	    -3.642666471599222873936118430711e-7,
	    6.944372610005012589931277214633e-8,
	    -1.37122090210436601953460514121e-8,
	    2.788389661007137131963860348087e-9,
	    -5.814164724331161551864791050316e-10,
	    1.23892049175275318118016881795e-10,
	    -2.690639145306743432390424937889e-11,
	    5.94261435084791098244470968384e-12,
	    -1.33238673575811957928775442057e-12,
	    3.028046806177132017173697243304e-13,
	    -6.966648814941032588795867588954e-14,
	    1.620854541053922969812893227628e-14,
	    -3.809934465250491999876913057729e-15,
	    9.040487815978831149368971012975e-16,
	    -2.164006195089607347809812047003e-16,
	    5.222102233995854984607980244172e-17,
	    -1.26972960236455533637241552778e-17,
	    3.109145504276197583836227412951e-18,
	    -7.663762920320385524009566714811e-19,
	    1.90081925136274520253692973329e-19,
	    -4.742207279069039545225655999965e-20,
	    1.189649200076528382880683078451e-20,
	    -3.000035590325780256845271313066e-21,
	    7.602993453043246173019385277098e-22,
	    -1.93590944760687288156981104913e-22,
	    4.951399124773337881000042386773e-23,
	    -1.271807481336371879608621989888e-23,
	    3.280049600469513043315841652053e-24,
	    -8.492320176822896568924792422399e-25,
	    2.206917892807560223519879987199e-25,
	    -5.755617245696528498312819507199e-26,
	    1.506191533639234250354144051199e-26,
	    -3.954502959018796953104285695999e-27,
	    1.041529704151500979984645051733e-27,
	    -2.751487795278765079450178901333e-28,
	    7.29005820549755740899770368e-29,
	    -1.936939645915947804077501098666e-29,
	    5.160357112051487298370054826666e-30,
	    -1.3784193221930940993896448e-30,
	    3.691326793107069042251093333333e-31,
	    -9.909389590624365420653226666666e-32,
	    2.666491705195388413323946666666e-32 };
    static doublereal erfccs[59] = { .0715179310202924774503697709496,
	    -.0265324343376067157558893386681,
	    .00171115397792085588332699194606,
	    -1.63751663458517884163746404749e-4,
	    1.98712935005520364995974806758e-5,
	    -2.84371241276655508750175183152e-6,
	    4.60616130896313036969379968464e-7,
	    -8.22775302587920842057766536366e-8,
	    1.59214187277090112989358340826e-8,
	    -3.29507136225284321486631665072e-9,
	    7.2234397604005554658126115389e-10,
	    -1.66485581339872959344695966886e-10,
	    4.01039258823766482077671768814e-11,
	    -1.00481621442573113272170176283e-11,
	    2.60827591330033380859341009439e-12,
	    -6.99111056040402486557697812476e-13,
	    1.92949233326170708624205749803e-13,
	    -5.47013118875433106490125085271e-14,
	    1.58966330976269744839084032762e-14,
	    -4.7268939801975548392036958429e-15,
	    1.4358733767849847867287399784e-15,
	    -4.44951056181735839417250062829e-16,
	    1.40481088476823343737305537466e-16,
	    -4.51381838776421089625963281623e-17,
	    1.47452154104513307787018713262e-17,
	    -4.89262140694577615436841552532e-18,
	    1.64761214141064673895301522827e-18,
	    -5.62681717632940809299928521323e-19,
	    1.94744338223207851429197867821e-19,
	    -6.82630564294842072956664144723e-20,
	    2.42198888729864924018301125438e-20,
	    -8.69341413350307042563800861857e-21,
	    3.15518034622808557122363401262e-21,
	    -1.15737232404960874261239486742e-21,
	    4.28894716160565394623737097442e-22,
	    -1.60503074205761685005737770964e-22,
	    6.06329875745380264495069923027e-23,
	    -2.31140425169795849098840801367e-23,
	    8.88877854066188552554702955697e-24,
	    -3.44726057665137652230718495566e-24,
	    1.34786546020696506827582774181e-24,
	    -5.31179407112502173645873201807e-25,
	    2.10934105861978316828954734537e-25,
	    -8.43836558792378911598133256738e-26,
	    3.39998252494520890627359576337e-26,
	    -1.3794523880732420900223837711e-26,
	    5.63449031183325261513392634811e-27,
	    -2.316490434477065448234277527e-27,
	    9.58446284460181015263158381226e-28,
	    -3.99072288033010972624224850193e-28,
	    1.67212922594447736017228709669e-28,
	    -7.04599152276601385638803782587e-29,
	    2.97976840286420635412357989444e-29,
	    -1.26252246646061929722422632994e-29,
	    5.39543870454248793985299653154e-30,
	    -2.38099288253145918675346190062e-30,
	    1.0990528301027615735972668375e-30,
	    -4.86771374164496572732518677435e-31,
	    1.52587726411035756763200828211e-31 };
    static doublereal sqrtpi = 1.77245385090551602729816748334115;
    static logical first = TRUE_;

    /* System generated locals */
    doublereal ret_val, d__1;

    /* Local variables */
    static doublereal y;
    static real eta;
    static doublereal xmax, xsml;
    static integer nterf;
    static doublereal sqeps;
    extern doublereal d1mach_(integer *);
    static doublereal txmax;
    static integer nterc2;
    extern doublereal dcsevl_(doublereal *, doublereal *, integer *);
    static integer nterfc;
    extern integer initds_(doublereal *, integer *, real *);
    extern /* Subroutine */ int xermsg_(char *, char *, char *, integer *, 
	    integer *, ftnlen, ftnlen, ftnlen);

/* ***BEGIN PROLOGUE  DERFC */
/* ***PURPOSE  Compute the complementary error function. */
/* ***LIBRARY   SLATEC (FNLIB) */
/* ***CATEGORY  C8A, L5A1E */
/* ***TYPE      DOUBLE PRECISION (ERFC-S, DERFC-D) */
/* ***KEYWORDS  COMPLEMENTARY ERROR FUNCTION, ERFC, FNLIB, */
/*             SPECIAL FUNCTIONS */
/* ***AUTHOR  Fullerton, W., (LANL) */
/* ***DESCRIPTION */

/* DERFC(X) calculates the double precision complementary error function */
/* for double precision argument X. */

/* Series for ERF        on the interval  0.          to  1.00000E+00 */
/*                                        with weighted Error   1.28E-32 */
/*                                         log weighted Error  31.89 */
/*                               significant figures required  31.05 */
/*                                    decimal places required  32.55 */

/* Series for ERC2       on the interval  2.50000E-01 to  1.00000E+00 */
/*                                        with weighted Error   2.67E-32 */
/*                                         log weighted Error  31.57 */
/*                               significant figures required  30.31 */
/*                                    decimal places required  32.42 */

/* Series for ERFC       on the interval  0.          to  2.50000E-01 */
/*                                        with weighted error   1.53E-31 */
/*                                         log weighted error  30.82 */
/*                               significant figures required  29.47 */
/*                                    decimal places required  31.70 */

/* ***REFERENCES  (NONE) */
/* ***ROUTINES CALLED  D1MACH, DCSEVL, INITDS, XERMSG */
/* ***REVISION HISTORY  (YYMMDD) */
/*   770701  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890531  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ) */
/*   920618  Removed space from variable names.  (RWC, WRB) */
/* ***END PROLOGUE  DERFC */
/* ***FIRST EXECUTABLE STATEMENT  DERFC */
    if (first) {
	eta = (real) d1mach_(&c__3) * .1f;
	nterf = initds_(erfcs, &c__21, &eta);
	nterfc = initds_(erfccs, &c__59, &eta);
	nterc2 = initds_(erc2cs, &c__49, &eta);

	xsml = -sqrt(-log(sqrtpi * d1mach_(&c__3)));
	txmax = sqrt(-log(sqrtpi * d1mach_(&c__1)));
	xmax = txmax - log(txmax) * .5 / txmax - .01;
	sqeps = sqrt(d1mach_(&c__3) * 2.);
    }
    first = FALSE_;

    if (*x > xsml) {
	goto L20;
    }

/* ERFC(X) = 1.0 - ERF(X)  FOR  X .LT. XSML */

    ret_val = 2.;
    return ret_val;

L20:
    if (*x > xmax) {
	goto L40;
    }
    y = abs(*x);
    if (y > 1.) {
	goto L30;
    }

/* ERFC(X) = 1.0 - ERF(X)  FOR ABS(X) .LE. 1.0 */

    if (y < sqeps) {
	ret_val = 1. - *x * 2. / sqrtpi;
    }
    if (y >= sqeps) {
	d__1 = *x * 2. * *x - 1.;
	ret_val = 1. - *x * (dcsevl_(&d__1, erfcs, &nterf) + 1.);
    }
    return ret_val;

/* ERFC(X) = 1.0 - ERF(X)  FOR  1.0 .LT. ABS(X) .LE. XMAX */

L30:
    y *= y;
    if (y <= 4.) {
	d__1 = (8. / y - 5.) / 3.;
	ret_val = exp(-y) / abs(*x) * (dcsevl_(&d__1, erc2cs, &nterc2) + .5);
    }
    if (y > 4.) {
	d__1 = 8. / y - 1.;
	ret_val = exp(-y) / abs(*x) * (dcsevl_(&d__1, erfccs, &nterfc) + .5);
    }
    if (*x < 0.) {
	ret_val = 2. - ret_val;
    }
    return ret_val;

L40:
    xermsg_("SLATEC", "DERFC", "X SO BIG ERFC UNDERFLOWS", &c__1, &c__1, (
	    ftnlen)6, (ftnlen)5, (ftnlen)24);
    ret_val = 0.;
    return ret_val;

} /* derfc_ */

