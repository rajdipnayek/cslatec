/* dqng.f -- translated by f2c (version 12.02.01).
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
static integer c__1 = 1;
static doublereal c_b17 = 1.5;
static integer c__0 = 0;

/* DECK DQNG */
/* Subroutine */ int dqng_(D_fp f, doublereal *a, doublereal *b, doublereal *
	epsabs, doublereal *epsrel, doublereal *result, doublereal *abserr, 
	integer *neval, integer *ier)
{
    /* Initialized data */

    static doublereal x1[5] = { .973906528517171720077964012084452,
	    .865063366688984510732096688423493,
	    .679409568299024406234327365114874,
	    .433395394129247190799265943165784,
	    .14887433898163121088482600112972 };
    static doublereal w87a[21] = { .00814837738414917290000287844819,
	    .018761438201562822243935059003794,
	    .027347451050052286161582829741283,
	    .033677707311637930046581056957588,
	    .036935099820427907614589586742499,
	    .002884872430211530501334156248695,
	    .013685946022712701888950035273128,
	    .023280413502888311123409291030404,
	    .030872497611713358675466394126442,
	    .035693633639418770719351355457044,
	    9.15283345202241360843392549948e-4,
	    .005399280219300471367738743391053,
	    .010947679601118931134327826856808,
	    .01629873169678733526266570322328,
	    .02108156888920383511243306018819,
	    .02537096976925382724346799983171,
	    .02918969775647575250144615408492,
	    .032373202467202789685788194889595,
	    .034783098950365142750781997949596,
	    .036412220731351787562801163687577,
	    .037253875503047708539592001191226 };
    static doublereal w87b[23] = { 2.74145563762072350016527092881e-4,
	    .001807124155057942948341311753254,
	    .00409686928275916486445807068348,
	    .006758290051847378699816577897424,
	    .009549957672201646536053581325377,
	    .01232944765224485369462663996378,
	    .015010447346388952376697286041943,
	    .0175489679862431910996653529259,
	    .019938037786440888202278192730714,
	    .022194935961012286796332102959499,
	    .024339147126000805470360647041454,
	    .026374505414839207241503786552615,
	    .02828691078877120065996800298796,
	    .030052581128092695322521110347341,
	    .031646751371439929404586051078883,
	    .033050413419978503290785944862689,
	    .034255099704226061787082821046821,
	    .035262412660156681033782717998428,
	    .036076989622888701185500318003895,
	    .036698604498456094498018047441094,
	    .037120549269832576114119958413599,
	    .037334228751935040321235449094698,
	    .037361073762679023410321241766599 };
    static doublereal w10[5] = { .066671344308688137593568809893332,
	    .149451349150580593145776339657697,
	    .219086362515982043995534934228163,
	    .269266719309996355091226921569469,
	    .295524224714752870173892994651338 };
    static doublereal x2[5] = { .995657163025808080735527280689003,
	    .930157491355708226001207180059508,
	    .780817726586416897063717578345042,
	    .562757134668604683339000099272694,
	    .294392862701460198131126603103866 };
    static doublereal w21a[5] = { .03255816230796472747881897245939,
	    .07503967481091995276704314091619,
	    .109387158802297641899210590325805,
	    .134709217311473325928054001771707,
	    .147739104901338491374841515972068 };
    static doublereal w21b[6] = { .011694638867371874278064396062192,
	    .05475589657435199603138130024458,
	    .093125454583697605535065465083366,
	    .123491976262065851077958109831074,
	    .142775938577060080797094273138717,
	    .149445554002916905664936468389821 };
    static doublereal x3[11] = { .999333360901932081394099323919911,
	    .987433402908088869795961478381209,
	    .954807934814266299257919200290473,
	    .900148695748328293625099494069092,
	    .82519831498311415084706673258852,
	    .732148388989304982612354848755461,
	    .622847970537725238641159120344323,
	    .499479574071056499952214885499755,
	    .364901661346580768043989548502644,
	    .222254919776601296498260928066212,
	    .074650617461383322043914435796506 };
    static doublereal w43a[10] = { .016296734289666564924281974617663,
	    .037522876120869501461613795898115,
	    .054694902058255442147212685465005,
	    .067355414609478086075553166302174,
	    .073870199632393953432140695251367,
	    .005768556059769796184184327908655,
	    .027371890593248842081276069289151,
	    .046560826910428830743339154433824,
	    .061744995201442564496240336030883,
	    .071387267268693397768559114425516 };
    static doublereal w43b[12] = { .001844477640212414100389106552965,
	    .010798689585891651740465406741293,
	    .021895363867795428102523123075149,
	    .032597463975345689443882222526137,
	    .042163137935191811847627924327955,
	    .050741939600184577780189020092084,
	    .058379395542619248375475369330206,
	    .064746404951445885544689259517511,
	    .069566197912356484528633315038405,
	    .072824441471833208150939535192842,
	    .074507751014175118273571813842889,
	    .074722147517403005594425168280423 };
    static doublereal x4[22] = { .999902977262729234490529830591582,
	    .99798989598667874542749632236596,
	    .992175497860687222808523352251425,
	    .981358163572712773571916941623894,
	    .965057623858384619128284110607926,
	    .943167613133670596816416634507426,
	    .91580641468550720959182643072005,
	    .883221657771316501372117548744163,
	    .845710748462415666605902011504855,
	    .803557658035230982788739474980964,
	    .75700573068549555832894279343202,
	    .70627320978732181982409427474084,
	    .651589466501177922534422205016736,
	    .593223374057961088875273770349144,
	    .531493605970831932285268948562671,
	    .46676362304202284487196678165927,
	    .399424847859218804732101665817923,
	    .329874877106188288265053371824597,
	    .258503559202161551802280975429025,
	    .185695396568346652015917141167606,
	    .111842213179907468172398359241362,
	    .037352123394619870814998165437704 };

    /* System generated locals */
    doublereal d__1, d__2, d__3, d__4;

    /* Local variables */
    static integer k, l;
    static doublereal fv1[5], fv2[5], fv3[5], fv4[5];
    static integer ipx;
    static doublereal absc, fval, res10, res21, res43, res87, fval1, fval2, 
	    hlgth, centr, reskh, uflow;
    extern doublereal d1mach_(integer *);
    static doublereal epmach, dhlgth, resabs, resasc, fcentr, savfun[21];
    extern /* Subroutine */ int xermsg_(char *, char *, char *, integer *, 
	    integer *, ftnlen, ftnlen, ftnlen);

/* ***BEGIN PROLOGUE  DQNG */
/* ***PURPOSE  The routine calculates an approximation result to a */
/*            given definite integral I = integral of F over (A,B), */
/*            hopefully satisfying following claim for accuracy */
/*            ABS(I-RESULT).LE.MAX(EPSABS,EPSREL*ABS(I)). */
/* ***LIBRARY   SLATEC (QUADPACK) */
/* ***CATEGORY  H2A1A1 */
/* ***TYPE      DOUBLE PRECISION (QNG-S, DQNG-D) */
/* ***KEYWORDS  AUTOMATIC INTEGRATOR, GAUSS-KRONROD(PATTERSON) RULES, */
/*             NONADAPTIVE, QUADPACK, QUADRATURE, SMOOTH INTEGRAND */
/* ***AUTHOR  Piessens, Robert */
/*             Applied Mathematics and Programming Division */
/*             K. U. Leuven */
/*           de Doncker, Elise */
/*             Applied Mathematics and Programming Division */
/*             K. U. Leuven */
/* ***DESCRIPTION */

/* NON-ADAPTIVE INTEGRATION */
/* STANDARD FORTRAN SUBROUTINE */
/* DOUBLE PRECISION VERSION */

/*           F      - Double precision */
/*                    Function subprogram defining the integrand function */
/*                    F(X). The actual name for F needs to be declared */
/*                    E X T E R N A L in the driver program. */

/*           A      - Double precision */
/*                    Lower limit of integration */

/*           B      - Double precision */
/*                    Upper limit of integration */

/*           EPSABS - Double precision */
/*                    Absolute accuracy requested */
/*           EPSREL - Double precision */
/*                    Relative accuracy requested */
/*                    If  EPSABS.LE.0 */
/*                    And EPSREL.LT.MAX(50*REL.MACH.ACC.,0.5D-28), */
/*                    The routine will end with IER = 6. */

/*         ON RETURN */
/*           RESULT - Double precision */
/*                    Approximation to the integral I */
/*                    Result is obtained by applying the 21-POINT */
/*                    GAUSS-KRONROD RULE (RES21) obtained by optimal */
/*                    addition of abscissae to the 10-POINT GAUSS RULE */
/*                    (RES10), or by applying the 43-POINT RULE (RES43) */
/*                    obtained by optimal addition of abscissae to the */
/*                    21-POINT GAUSS-KRONROD RULE, or by applying the */
/*                    87-POINT RULE (RES87) obtained by optimal addition */
/*                    of abscissae to the 43-POINT RULE. */

/*           ABSERR - Double precision */
/*                    Estimate of the modulus of the absolute error, */
/*                    which should EQUAL or EXCEED ABS(I-RESULT) */

/*           NEVAL  - Integer */
/*                    Number of integrand evaluations */

/*           IER    - IER = 0 normal and reliable termination of the */
/*                            routine. It is assumed that the requested */
/*                            accuracy has been achieved. */
/*                    IER.GT.0 Abnormal termination of the routine. It is */
/*                            assumed that the requested accuracy has */
/*                            not been achieved. */
/*           ERROR MESSAGES */
/*                    IER = 1 The maximum number of steps has been */
/*                            executed. The integral is probably too */
/*                            difficult to be calculated by DQNG. */
/*                        = 6 The input is invalid, because */
/*                            EPSABS.LE.0 AND */
/*                            EPSREL.LT.MAX(50*REL.MACH.ACC.,0.5D-28). */
/*                            RESULT, ABSERR and NEVAL are set to zero. */

/* ***REFERENCES  (NONE) */
/* ***ROUTINES CALLED  D1MACH, XERMSG */
/* ***REVISION HISTORY  (YYMMDD) */
/*   800101  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890531  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ) */
/* ***END PROLOGUE  DQNG */



/*           THE FOLLOWING DATA STATEMENTS CONTAIN THE */
/*           ABSCISSAE AND WEIGHTS OF THE INTEGRATION RULES USED. */

/*           X1      ABSCISSAE COMMON TO THE 10-, 21-, 43- AND 87- */
/*                   POINT RULE */
/*           X2      ABSCISSAE COMMON TO THE 21-, 43- AND 87-POINT RULE */
/*           X3      ABSCISSAE COMMON TO THE 43- AND 87-POINT RULE */
/*           X4      ABSCISSAE OF THE 87-POINT RULE */
/*           W10     WEIGHTS OF THE 10-POINT FORMULA */
/*           W21A    WEIGHTS OF THE 21-POINT FORMULA FOR ABSCISSAE X1 */
/*           W21B    WEIGHTS OF THE 21-POINT FORMULA FOR ABSCISSAE X2 */
/*           W43A    WEIGHTS OF THE 43-POINT FORMULA FOR ABSCISSAE X1, X3 */
/*           W43B    WEIGHTS OF THE 43-POINT FORMULA FOR ABSCISSAE X3 */
/*           W87A    WEIGHTS OF THE 87-POINT FORMULA FOR ABSCISSAE X1, */
/*                   X2, X3 */
/*           W87B    WEIGHTS OF THE 87-POINT FORMULA FOR ABSCISSAE X4 */


/* GAUSS-KRONROD-PATTERSON QUADRATURE COEFFICIENTS FOR USE IN */
/* QUADPACK ROUTINE QNG.  THESE COEFFICIENTS WERE CALCULATED WITH */
/* 101 DECIMAL DIGIT ARITHMETIC BY L. W. FULLERTON, BELL LABS, NOV 1981. */





/*           LIST OF MAJOR VARIABLES */
/*           ----------------------- */

/*           CENTR  - MID POINT OF THE INTEGRATION INTERVAL */
/*           HLGTH  - HALF-LENGTH OF THE INTEGRATION INTERVAL */
/*           FCENTR - FUNCTION VALUE AT MID POINT */
/*           ABSC   - ABSCISSA */
/*           FVAL   - FUNCTION VALUE */
/*           SAVFUN - ARRAY OF FUNCTION VALUES WHICH HAVE ALREADY BEEN */
/*                    COMPUTED */
/*           RES10  - 10-POINT GAUSS RESULT */
/*           RES21  - 21-POINT KRONROD RESULT */
/*           RES43  - 43-POINT RESULT */
/*           RES87  - 87-POINT RESULT */
/*           RESABS - APPROXIMATION TO THE INTEGRAL OF ABS(F) */
/*           RESASC - APPROXIMATION TO THE INTEGRAL OF ABS(F-I/(B-A)) */

/*           MACHINE DEPENDENT CONSTANTS */
/*           --------------------------- */

/*           EPMACH IS THE LARGEST RELATIVE SPACING. */
/*           UFLOW IS THE SMALLEST POSITIVE MAGNITUDE. */

/* ***FIRST EXECUTABLE STATEMENT  DQNG */
    epmach = d1mach_(&c__4);
    uflow = d1mach_(&c__1);

/*           TEST ON VALIDITY OF PARAMETERS */
/*           ------------------------------ */

    *result = 0.;
    *abserr = 0.;
    *neval = 0;
    *ier = 6;
/* Computing MAX */
    d__1 = epmach * 50.;
    if (*epsabs <= 0. && *epsrel < max(d__1,5e-29)) {
	goto L80;
    }
    hlgth = (*b - *a) * .5;
    dhlgth = abs(hlgth);
    centr = (*b + *a) * .5;
    fcentr = (*f)(&centr);
    *neval = 21;
    *ier = 1;

/*          COMPUTE THE INTEGRAL USING THE 10- AND 21-POINT FORMULA. */

    for (l = 1; l <= 3; ++l) {
	switch (l) {
	    case 1:  goto L5;
	    case 2:  goto L25;
	    case 3:  goto L45;
	}
L5:
	res10 = 0.;
	res21 = w21b[5] * fcentr;
	resabs = w21b[5] * abs(fcentr);
	for (k = 1; k <= 5; ++k) {
	    absc = hlgth * x1[k - 1];
	    d__1 = centr + absc;
	    fval1 = (*f)(&d__1);
	    d__1 = centr - absc;
	    fval2 = (*f)(&d__1);
	    fval = fval1 + fval2;
	    res10 += w10[k - 1] * fval;
	    res21 += w21a[k - 1] * fval;
	    resabs += w21a[k - 1] * (abs(fval1) + abs(fval2));
	    savfun[k - 1] = fval;
	    fv1[k - 1] = fval1;
	    fv2[k - 1] = fval2;
/* L10: */
	}
	ipx = 5;
	for (k = 1; k <= 5; ++k) {
	    ++ipx;
	    absc = hlgth * x2[k - 1];
	    d__1 = centr + absc;
	    fval1 = (*f)(&d__1);
	    d__1 = centr - absc;
	    fval2 = (*f)(&d__1);
	    fval = fval1 + fval2;
	    res21 += w21b[k - 1] * fval;
	    resabs += w21b[k - 1] * (abs(fval1) + abs(fval2));
	    savfun[ipx - 1] = fval;
	    fv3[k - 1] = fval1;
	    fv4[k - 1] = fval2;
/* L15: */
	}

/*          TEST FOR CONVERGENCE. */

	*result = res21 * hlgth;
	resabs *= dhlgth;
	reskh = res21 * .5;
	resasc = w21b[5] * (d__1 = fcentr - reskh, abs(d__1));
	for (k = 1; k <= 5; ++k) {
	    resasc = resasc + w21a[k - 1] * ((d__1 = fv1[k - 1] - reskh, abs(
		    d__1)) + (d__2 = fv2[k - 1] - reskh, abs(d__2))) + w21b[k 
		    - 1] * ((d__3 = fv3[k - 1] - reskh, abs(d__3)) + (d__4 = 
		    fv4[k - 1] - reskh, abs(d__4)));
/* L20: */
	}
	*abserr = (d__1 = (res21 - res10) * hlgth, abs(d__1));
	resasc *= dhlgth;
	goto L65;

/*          COMPUTE THE INTEGRAL USING THE 43-POINT FORMULA. */

L25:
	res43 = w43b[11] * fcentr;
	*neval = 43;
	for (k = 1; k <= 10; ++k) {
	    res43 += savfun[k - 1] * w43a[k - 1];
/* L30: */
	}
	for (k = 1; k <= 11; ++k) {
	    ++ipx;
	    absc = hlgth * x3[k - 1];
	    d__1 = absc + centr;
	    d__2 = centr - absc;
	    fval = (*f)(&d__1) + (*f)(&d__2);
	    res43 += fval * w43b[k - 1];
	    savfun[ipx - 1] = fval;
/* L40: */
	}

/*          TEST FOR CONVERGENCE. */

	*result = res43 * hlgth;
	*abserr = (d__1 = (res43 - res21) * hlgth, abs(d__1));
	goto L65;

/*          COMPUTE THE INTEGRAL USING THE 87-POINT FORMULA. */

L45:
	res87 = w87b[22] * fcentr;
	*neval = 87;
	for (k = 1; k <= 21; ++k) {
	    res87 += savfun[k - 1] * w87a[k - 1];
/* L50: */
	}
	for (k = 1; k <= 22; ++k) {
	    absc = hlgth * x4[k - 1];
	    d__1 = absc + centr;
	    d__2 = centr - absc;
	    res87 += w87b[k - 1] * ((*f)(&d__1) + (*f)(&d__2));
/* L60: */
	}
	*result = res87 * hlgth;
	*abserr = (d__1 = (res87 - res43) * hlgth, abs(d__1));
L65:
	if (resasc != 0. && *abserr != 0.) {
/* Computing MIN */
	    d__3 = *abserr * 200. / resasc;
	    d__1 = 1., d__2 = pow_dd(&d__3, &c_b17);
	    *abserr = resasc * min(d__1,d__2);
	}
	if (resabs > uflow / (epmach * 50.)) {
/* Computing MAX */
	    d__1 = epmach * 50. * resabs;
	    *abserr = max(d__1,*abserr);
	}
/* Computing MAX */
	d__1 = *epsabs, d__2 = *epsrel * abs(*result);
	if (*abserr <= max(d__1,d__2)) {
	    *ier = 0;
	}
/* ***JUMP OUT OF DO-LOOP */
	if (*ier == 0) {
	    goto L999;
	}
/* L70: */
    }
L80:
    xermsg_("SLATEC", "DQNG", "ABNORMAL RETURN", ier, &c__0, (ftnlen)6, (
	    ftnlen)4, (ftnlen)15);
L999:
    return 0;
} /* dqng_ */

