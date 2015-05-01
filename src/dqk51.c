/* dqk51.f -- translated by f2c (version 12.02.01).
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
static doublereal c_b7 = 1.5;

/* DECK DQK51 */
/* Subroutine */ int dqk51_(D_fp f, doublereal *a, doublereal *b, doublereal *
	result, doublereal *abserr, doublereal *resabs, doublereal *resasc)
{
    /* Initialized data */

    static doublereal wg[13] = { .011393798501026287947902964113235,
	    .026354986615032137261901815295299,
	    .040939156701306312655623487711646,
	    .054904695975835191925936891540473,
	    .068038333812356917207187185656708,
	    .080140700335001018013234959669111,
	    .091028261982963649811497220702892,
	    .100535949067050644202206890392686,
	    .108519624474263653116093957050117,
	    .114858259145711648339325545869556,
	    .119455763535784772228178126512901,
	    .122242442990310041688959518945852,
	    .12317605372671545120390287307905 };
    static doublereal xgk[26] = { .999262104992609834193457486540341,
	    .995556969790498097908784946893902,
	    .988035794534077247637331014577406,
	    .976663921459517511498315386479594,
	    .961614986425842512418130033660167,
	    .942974571228974339414011169658471,
	    .920747115281701561746346084546331,
	    .894991997878275368851042006782805,
	    .86584706529327559544899696958834,
	    .83344262876083400142102110869357,
	    .797873797998500059410410904994307,
	    .759259263037357630577282865204361,
	    .717766406813084388186654079773298,
	    .673566368473468364485120633247622,
	    .626810099010317412788122681624518,
	    .577662930241222967723689841612654,
	    .52632528433471918259962377815801,
	    .473002731445714960522182115009192,
	    .417885382193037748851814394594572,
	    .361172305809387837735821730127641,
	    .303089538931107830167478909980339,
	    .243866883720988432045190362797452,
	    .183718939421048892015969888759528,
	    .122864692610710396387359818808037,
	    .061544483005685078886546392366797,0. };
    static doublereal wgk[26] = { .001987383892330315926507851882843,
	    .005561932135356713758040236901066,
	    .009473973386174151607207710523655,
	    .013236229195571674813656405846976,
	    .016847817709128298231516667536336,
	    .020435371145882835456568292235939,
	    .024009945606953216220092489164881,
	    .027475317587851737802948455517811,
	    .030792300167387488891109020215229,
	    .034002130274329337836748795229551,
	    .03711627148341554356033062536762,
	    .040083825504032382074839284467076,
	    .042872845020170049476895792439495,
	    .04550291304992178890987058475266,
	    .047982537138836713906392255756915,
	    .05027767908071567196332525943344,
	    .052362885806407475864366712137873,
	    .054251129888545490144543370459876,
	    .055950811220412317308240686382747,
	    .057437116361567832853582693939506,
	    .058689680022394207961974175856788,
	    .059720340324174059979099291932562,
	    .060539455376045862945360267517565,
	    .061128509717053048305859030416293,
	    .061471189871425316661544131965264,
	    .061580818067832935078759824240055 };

    /* System generated locals */
    doublereal d__1, d__2, d__3;

    /* Local variables */
    static integer j;
    static doublereal fc, fv1[25], fv2[25];
    static integer jtw;
    static doublereal absc, resg, resk, fsum, fval1, fval2;
    static integer jtwm1;
    static doublereal hlgth, centr, reskh, uflow;
    extern doublereal d1mach_(integer *);
    static doublereal epmach, dhlgth;

/* ***BEGIN PROLOGUE  DQK51 */
/* ***PURPOSE  To compute I = Integral of F over (A,B) with error */
/*                           estimate */
/*                       J = Integral of ABS(F) over (A,B) */
/* ***LIBRARY   SLATEC (QUADPACK) */
/* ***CATEGORY  H2A1A2 */
/* ***TYPE      DOUBLE PRECISION (QK51-S, DQK51-D) */
/* ***KEYWORDS  51-POINT GAUSS-KRONROD RULES, QUADPACK, QUADRATURE */
/* ***AUTHOR  Piessens, Robert */
/*             Applied Mathematics and Programming Division */
/*             K. U. Leuven */
/*           de Doncker, Elise */
/*             Applied Mathematics and Programming Division */
/*             K. U. Leuven */
/* ***DESCRIPTION */

/*           Integration rules */
/*           Standard fortran subroutine */
/*           Double precision version */

/*           PARAMETERS */
/*            ON ENTRY */
/*              F      - Double precision */
/*                       Function subroutine defining the integrand */
/*                       function F(X). The actual name for F needs to be */
/*                       declared E X T E R N A L in the calling program. */

/*              A      - Double precision */
/*                       Lower limit of integration */

/*              B      - Double precision */
/*                       Upper limit of integration */

/*            ON RETURN */
/*              RESULT - Double precision */
/*                       Approximation to the integral I */
/*                       RESULT is computed by applying the 51-point */
/*                       Kronrod rule (RESK) obtained by optimal addition */
/*                       of abscissae to the 25-point Gauss rule (RESG). */

/*              ABSERR - Double precision */
/*                       Estimate of the modulus of the absolute error, */
/*                       which should not exceed ABS(I-RESULT) */

/*              RESABS - Double precision */
/*                       Approximation to the integral J */

/*              RESASC - Double precision */
/*                       Approximation to the integral of ABS(F-I/(B-A)) */
/*                       over (A,B) */

/* ***REFERENCES  (NONE) */
/* ***ROUTINES CALLED  D1MACH */
/* ***REVISION HISTORY  (YYMMDD) */
/*   800101  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890531  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   910819  Added WGK(26) to code.  (WRB) */
/* ***END PROLOGUE  DQK51 */



/*           THE ABSCISSAE AND WEIGHTS ARE GIVEN FOR THE INTERVAL (-1,1). */
/*           BECAUSE OF SYMMETRY ONLY THE POSITIVE ABSCISSAE AND THEIR */
/*           CORRESPONDING WEIGHTS ARE GIVEN. */

/*           XGK    - ABSCISSAE OF THE 51-POINT KRONROD RULE */
/*                    XGK(2), XGK(4), ...  ABSCISSAE OF THE 25-POINT */
/*                    GAUSS RULE */
/*                    XGK(1), XGK(3), ...  ABSCISSAE WHICH ARE OPTIMALLY */
/*                    ADDED TO THE 25-POINT GAUSS RULE */

/*           WGK    - WEIGHTS OF THE 51-POINT KRONROD RULE */

/*           WG     - WEIGHTS OF THE 25-POINT GAUSS RULE */


/* GAUSS QUADRATURE WEIGHTS AND KRONROD QUADRATURE ABSCISSAE AND WEIGHTS */
/* AS EVALUATED WITH 80 DECIMAL DIGIT ARITHMETIC BY L. W. FULLERTON, */
/* BELL LABS, NOV. 1981. */





/*           LIST OF MAJOR VARIABLES */
/*           ----------------------- */

/*           CENTR  - MID POINT OF THE INTERVAL */
/*           HLGTH  - HALF-LENGTH OF THE INTERVAL */
/*           ABSC   - ABSCISSA */
/*           FVAL*  - FUNCTION VALUE */
/*           RESG   - RESULT OF THE 25-POINT GAUSS FORMULA */
/*           RESK   - RESULT OF THE 51-POINT KRONROD FORMULA */
/*           RESKH  - APPROXIMATION TO THE MEAN VALUE OF F OVER (A,B), */
/*                    I.E. TO I/(B-A) */

/*           MACHINE DEPENDENT CONSTANTS */
/*           --------------------------- */

/*           EPMACH IS THE LARGEST RELATIVE SPACING. */
/*           UFLOW IS THE SMALLEST POSITIVE MAGNITUDE. */

/* ***FIRST EXECUTABLE STATEMENT  DQK51 */
    epmach = d1mach_(&c__4);
    uflow = d1mach_(&c__1);

    centr = (*a + *b) * .5;
    hlgth = (*b - *a) * .5;
    dhlgth = abs(hlgth);

/*           COMPUTE THE 51-POINT KRONROD APPROXIMATION TO */
/*           THE INTEGRAL, AND ESTIMATE THE ABSOLUTE ERROR. */

    fc = (*f)(&centr);
    resg = wg[12] * fc;
    resk = wgk[25] * fc;
    *resabs = abs(resk);
    for (j = 1; j <= 12; ++j) {
	jtw = j << 1;
	absc = hlgth * xgk[jtw - 1];
	d__1 = centr - absc;
	fval1 = (*f)(&d__1);
	d__1 = centr + absc;
	fval2 = (*f)(&d__1);
	fv1[jtw - 1] = fval1;
	fv2[jtw - 1] = fval2;
	fsum = fval1 + fval2;
	resg += wg[j - 1] * fsum;
	resk += wgk[jtw - 1] * fsum;
	*resabs += wgk[jtw - 1] * (abs(fval1) + abs(fval2));
/* L10: */
    }
    for (j = 1; j <= 13; ++j) {
	jtwm1 = (j << 1) - 1;
	absc = hlgth * xgk[jtwm1 - 1];
	d__1 = centr - absc;
	fval1 = (*f)(&d__1);
	d__1 = centr + absc;
	fval2 = (*f)(&d__1);
	fv1[jtwm1 - 1] = fval1;
	fv2[jtwm1 - 1] = fval2;
	fsum = fval1 + fval2;
	resk += wgk[jtwm1 - 1] * fsum;
	*resabs += wgk[jtwm1 - 1] * (abs(fval1) + abs(fval2));
/* L15: */
    }
    reskh = resk * .5;
    *resasc = wgk[25] * (d__1 = fc - reskh, abs(d__1));
    for (j = 1; j <= 25; ++j) {
	*resasc += wgk[j - 1] * ((d__1 = fv1[j - 1] - reskh, abs(d__1)) + (
		d__2 = fv2[j - 1] - reskh, abs(d__2)));
/* L20: */
    }
    *result = resk * hlgth;
    *resabs *= dhlgth;
    *resasc *= dhlgth;
    *abserr = (d__1 = (resk - resg) * hlgth, abs(d__1));
    if (*resasc != 0. && *abserr != 0.) {
/* Computing MIN */
	d__3 = *abserr * 200. / *resasc;
	d__1 = 1., d__2 = pow_dd(&d__3, &c_b7);
	*abserr = *resasc * min(d__1,d__2);
    }
    if (*resabs > uflow / (epmach * 50.)) {
/* Computing MAX */
	d__1 = epmach * 50. * *resabs;
	*abserr = max(d__1,*abserr);
    }
    return 0;
} /* dqk51_ */

