/* d9ln2r.f -- translated by f2c (version 12.02.01).
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
static integer c__50 = 50;
static integer c__37 = 37;
static integer c__4 = 4;
static integer c__1 = 1;
static integer c__2 = 2;

/* DECK D9LN2R */
doublereal d9ln2r_(doublereal *x)
{
    /* Initialized data */

    static doublereal ln21cs[50] = { .18111962513478809875894953043071,
	    -.15627123192872462669625155541078,
	    .028676305361557275209540627102051,
	    -.0055586996559481398781157725126781,
	    .0011178976652299837657335666279727,
	    -2.3080508982327947182299279585705e-4,
	    4.859885334110017587468155806875e-5,
	    -1.0390127388903210765514242633338e-5,
	    2.2484563707390128494621804946408e-6,
	    -4.9140592739266484875327802597091e-7,
	    1.0828256507077483336620152971597e-7,
	    -2.4025872763420701435976675416719e-8,
	    5.3624600472708133762984443250163e-9,
	    -1.2029951362138772264671646424377e-9,
	    2.7107889277591860785622551632266e-10,
	    -6.132356261831901006879672843069e-11,
	    1.3920858369159469857436908543978e-11,
	    -3.1699300330223494015283057260883e-12,
	    7.2383754044307505335214326197011e-13,
	    -1.6570017184764411391498805506268e-13,
	    3.8018428663117424257364422631876e-14,
	    -8.7411189296972700259724429899137e-15,
	    2.0135619845055748302118751028154e-15,
	    -4.6464456409033907031102008154477e-16,
	    1.0739282147018339453453338554925e-16,
	    -2.485853461993779475553402183396e-17,
	    5.7620197950800189813888142628181e-18,
	    -1.337306376980439470140219995805e-18,
	    3.1074653227331824966533807166805e-19,
	    -7.2288104083040539906901957917627e-20,
	    1.6833783788037385103313258186888e-20,
	    -3.9239463312069958052519372739925e-21,
	    9.1551468387536789746385528640853e-22,
	    -2.1378895321320159520982095801002e-22,
	    4.9964507479047864699828564568746e-23,
	    -1.1686240636080170135360806147413e-23,
	    2.7353123470391863775628686786559e-24,
	    -6.4068025084792111965050345881599e-25,
	    1.5016293204334124162949071940266e-25,
	    -3.5217372410398479759497145002666e-26,
	    8.2643901014814767012482733397333e-27,
	    -1.9404930275943401918036617898666e-27,
	    4.5587880018841283562451588437333e-28,
	    -1.0715492087545202154378625023999e-28,
	    2.5199408007927592978096674133333e-29,
	    -5.92890884001209693417504768e-30,
	    1.3955864061057513058237153279999e-30,
	    -3.2864578813478583431436697599999e-31,
	    7.7424967950478166247254698666666e-32,
	    -1.8247735667260887638125226666666e-32 };
    static doublereal ln22cs[37] = { -.2224253253502046082986015223552,
	    -.06104710010807862398680104755764,
	    .007427235009750394590519629755729,
	    -9.335018261636970565612779606397e-4,
	    1.200499076872601283350731287359e-4,
	    -1.570472295282004112823352608243e-5,
	    2.081874781051271096050783592759e-6,
	    -2.789195577646713654057213051375e-7,
	    3.769355823760132058422895135447e-8,
	    -5.130902896527711258240589938003e-9,
	    7.027141178150694738206218215392e-10,
	    -9.674859550134342389243972005137e-11,
	    1.338104645924887306588496449748e-11,
	    -1.858102603534063981628453846591e-12,
	    2.58929442252791974930860012307e-13,
	    -3.619568316141588674466025382172e-14,
	    5.074037398016623088006858917396e-15,
	    -7.131012977031127302700938748927e-16,
	    1.004490328554567481853386784126e-16,
	    -1.417906532184025791904405075285e-17,
	    2.005297034743326117891086396074e-18,
	    -2.840996662339803305365396717567e-19,
	    4.031469883969079899599878662826e-20,
	    -5.729325241832207320455498956799e-21,
	    8.153488253890010675848928733866e-22,
	    -1.161825588549721787606027468799e-22,
	    1.657516611662538343659339775999e-23,
	    -2.36733670471080519011401728e-24,
	    3.384670367975521386076569599999e-25,
	    -4.843940829215718204296396799999e-26,
	    6.938759162514273718676138666666e-27,
	    -9.948142607031436571923797333333e-28,
	    1.427440611211698610634752e-28,
	    -2.049794721898234911566506666666e-29,
	    2.945648756401362222885546666666e-30,
	    -4.235973185184957027669333333333e-31,
	    6.095532614003832040106666666666e-32 };
    static logical first = TRUE_;

    /* System generated locals */
    real r__1;
    doublereal ret_val, d__1;

    /* Local variables */
    static real eps;
    static doublereal xbig, xmin, xmax, txbig;
    static integer ntln21, ntln22;
    static real sqeps;
    extern doublereal d1mach_(integer *);
    static doublereal txmax;
    extern doublereal dcsevl_(doublereal *, doublereal *, integer *);
    extern integer initds_(doublereal *, integer *, real *);
    extern /* Subroutine */ int xermsg_(char *, char *, char *, integer *, 
	    integer *, ftnlen, ftnlen, ftnlen);

/* ***BEGIN PROLOGUE  D9LN2R */
/* ***SUBSIDIARY */
/* ***PURPOSE  Evaluate LOG(1+X) from second order relative accuracy so */
/*            that LOG(1+X) = X - X**2/2 + X**3*D9LN2R(X) */
/* ***LIBRARY   SLATEC (FNLIB) */
/* ***CATEGORY  C4B */
/* ***TYPE      DOUBLE PRECISION (R9LN2R-S, D9LN2R-D, C9LN2R-C) */
/* ***KEYWORDS  ELEMENTARY FUNCTIONS, FNLIB, LOGARITHM, SECOND ORDER */
/* ***AUTHOR  Fullerton, W., (LANL) */
/* ***DESCRIPTION */

/* Evaluate  LOG(1+X)  from 2-nd order with relative error accuracy so */
/* that    LOG(1+X) = X - X**2/2 + X**3*D9LN2R(X) */

/* Series for LN21       on the interval -6.25000E-01 to  0. */
/*                                        with weighted error   1.82E-32 */
/*                                         log weighted error  31.74 */
/*                               significant figures required  31.00 */
/*                                    decimal places required  32.59 */

/* Series for LN22       on the interval  0.          to  8.12500E-01 */
/*                                        with weighted error   6.10E-32 */
/*                                         log weighted error  31.21 */
/*                               significant figures required  30.32 */
/*                                    decimal places required  32.00 */

/* ***REFERENCES  (NONE) */
/* ***ROUTINES CALLED  D1MACH, DCSEVL, INITDS, XERMSG */
/* ***REVISION HISTORY  (YYMMDD) */
/*   780401  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890911  Removed unnecessary intrinsics.  (WRB) */
/*   890911  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ) */
/*   900720  Routine changed from user-callable to subsidiary.  (WRB) */
/* ***END PROLOGUE  D9LN2R */
/* ***FIRST EXECUTABLE STATEMENT  D9LN2R */
    if (first) {
	eps = d1mach_(&c__3);
	r__1 = eps * .1f;
	ntln21 = initds_(ln21cs, &c__50, &r__1);
	r__1 = eps * .1f;
	ntln22 = initds_(ln22cs, &c__37, &r__1);

	xmin = sqrt(d1mach_(&c__4)) - 1.;
	sqeps = sqrt(eps);
	txmax = 8.f / sqeps;
/* Computing 2nd power */
	d__1 = txmax;
	xmax = txmax - (eps * (d__1 * d__1) - log(txmax) * 2.) / (eps * 2. * 
		txmax);
	txbig = 6.f / sqrt(sqeps);
/* Computing 2nd power */
	d__1 = txbig;
	xbig = txbig - (sqeps * (d__1 * d__1) - log(txbig) * 2.) / (sqeps * 
		2. * txbig);
    }
    first = FALSE_;

    if (*x < -.625 || *x > .8125) {
	goto L20;
    }

    if (*x < 0.) {
	d__1 = *x * 16. / 5. + 1.;
	ret_val = dcsevl_(&d__1, ln21cs, &ntln21) + .375;
    }
    if (*x >= 0.) {
	d__1 = *x * 32. / 13. - 1.;
	ret_val = dcsevl_(&d__1, ln22cs, &ntln22) + .375;
    }
    return ret_val;

L20:
    if (*x < xmin) {
	xermsg_("SLATEC", "D9LN2R", "ANSWER LT HALF PRECISION BECAUSE X IS T"
		"OO NEAR -1", &c__1, &c__1, (ftnlen)6, (ftnlen)6, (ftnlen)49);
    }
    if (*x > xmax) {
	xermsg_("SLATEC", "D9LN2R", "NO PRECISION IN ANSWER BECAUSE X IS TOO"
		" BIG", &c__3, &c__2, (ftnlen)6, (ftnlen)6, (ftnlen)43);
    }
    if (*x > xbig) {
	xermsg_("SLATEC", "D9LN2R", "ANSWER LT HALF PRECISION BECAUSE X IS T"
		"OO BIG", &c__2, &c__1, (ftnlen)6, (ftnlen)6, (ftnlen)45);
    }

/* Computing 3rd power */
    d__1 = *x;
    ret_val = (log(*x + 1.) - *x * (1. - *x * .5)) / (d__1 * (d__1 * d__1));
    return ret_val;

} /* d9ln2r_ */

