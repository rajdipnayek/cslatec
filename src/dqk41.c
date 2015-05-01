/* dqk41.f -- translated by f2c (version 12.02.01).
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

/* DECK DQK41 */
/* Subroutine */ int dqk41_(D_fp f, doublereal *a, doublereal *b, doublereal *
	result, doublereal *abserr, doublereal *resabs, doublereal *resasc)
{
    /* Initialized data */

    static doublereal wg[10] = { .017614007139152118311861962351853,
	    .040601429800386941331039952274932,
	    .062672048334109063569506535187042,
	    .083276741576704748724758143222046,
	    .10193011981724043503675013548035,
	    .118194531961518417312377377711382,
	    .131688638449176626898494499748163,
	    .142096109318382051329298325067165,
	    .149172986472603746787828737001969,
	    .152753387130725850698084331955098 };
    static doublereal xgk[21] = { .998859031588277663838315576545863,
	    .99312859918509492478612238847132,
	    .981507877450250259193342994720217,
	    .963971927277913791267666131197277,
	    .940822633831754753519982722212443,
	    .912234428251325905867752441203298,
	    .878276811252281976077442995113078,
	    .839116971822218823394529061701521,
	    .795041428837551198350638833272788,
	    .746331906460150792614305070355642,
	    .693237656334751384805490711845932,
	    .636053680726515025452836696226286,
	    .575140446819710315342946036586425,
	    .510867001950827098004364050955251,
	    .44359317523872510319999221349264,
	    .373706088715419560672548177024927,
	    .301627868114913004320555356858592,
	    .227785851141645078080496195368575,
	    .152605465240922675505220241022678,
	    .076526521133497333754640409398838,0. };
    static doublereal wgk[21] = { .003073583718520531501218293246031,
	    .008600269855642942198661787950102,
	    .014626169256971252983787960308868,
	    .020388373461266523598010231432755,
	    .025882133604951158834505067096153,
	    .031287306777032798958543119323801,
	    .036600169758200798030557240707211,
	    .041668873327973686263788305936895,
	    .046434821867497674720231880926108,
	    .050944573923728691932707670050345,
	    .055195105348285994744832372419777,
	    .059111400880639572374967220648594,
	    .062653237554781168025870122174255,
	    .065834597133618422111563556969398,
	    .068648672928521619345623411885368,
	    .07105442355344406830579036172321,
	    .073030690332786667495189417658913,
	    .074582875400499188986581418362488,
	    .075704497684556674659542775376617,
	    .076377867672080736705502835038061,
	    .076600711917999656445049901530102 };

    /* System generated locals */
    doublereal d__1, d__2, d__3;

    /* Local variables */
    static integer j;
    static doublereal fc, fv1[20], fv2[20];
    static integer jtw;
    static doublereal absc, resg, resk, fsum, fval1, fval2;
    static integer jtwm1;
    static doublereal hlgth, centr, reskh, uflow;
    extern doublereal d1mach_(integer *);
    static doublereal epmach, dhlgth;

/* ***BEGIN PROLOGUE  DQK41 */
/* ***PURPOSE  To compute I = Integral of F over (A,B), with error */
/*                           estimate */
/*                       J = Integral of ABS(F) over (A,B) */
/* ***LIBRARY   SLATEC (QUADPACK) */
/* ***CATEGORY  H2A1A2 */
/* ***TYPE      DOUBLE PRECISION (QK41-S, DQK41-D) */
/* ***KEYWORDS  41-POINT GAUSS-KRONROD RULES, QUADPACK, QUADRATURE */
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
/*                       Function subprogram defining the integrand */
/*                       FUNCTION F(X). The actual name for F needs to be */
/*                       declared E X T E R N A L in the calling program. */

/*              A      - Double precision */
/*                       Lower limit of integration */

/*              B      - Double precision */
/*                       Upper limit of integration */

/*            ON RETURN */
/*              RESULT - Double precision */
/*                       Approximation to the integral I */
/*                       RESULT is computed by applying the 41-POINT */
/*                       GAUSS-KRONROD RULE (RESK) obtained by optimal */
/*                       addition of abscissae to the 20-POINT GAUSS */
/*                       RULE (RESG). */

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
/* ***END PROLOGUE  DQK41 */



/*           THE ABSCISSAE AND WEIGHTS ARE GIVEN FOR THE INTERVAL (-1,1). */
/*           BECAUSE OF SYMMETRY ONLY THE POSITIVE ABSCISSAE AND THEIR */
/*           CORRESPONDING WEIGHTS ARE GIVEN. */

/*           XGK    - ABSCISSAE OF THE 41-POINT GAUSS-KRONROD RULE */
/*                    XGK(2), XGK(4), ...  ABSCISSAE OF THE 20-POINT */
/*                    GAUSS RULE */
/*                    XGK(1), XGK(3), ...  ABSCISSAE WHICH ARE OPTIMALLY */
/*                    ADDED TO THE 20-POINT GAUSS RULE */

/*           WGK    - WEIGHTS OF THE 41-POINT GAUSS-KRONROD RULE */

/*           WG     - WEIGHTS OF THE 20-POINT GAUSS RULE */


/* GAUSS QUADRATURE WEIGHTS AND KRONROD QUADRATURE ABSCISSAE AND WEIGHTS */
/* AS EVALUATED WITH 80 DECIMAL DIGIT ARITHMETIC BY L. W. FULLERTON, */
/* BELL LABS, NOV. 1981. */





/*           LIST OF MAJOR VARIABLES */
/*           ----------------------- */

/*           CENTR  - MID POINT OF THE INTERVAL */
/*           HLGTH  - HALF-LENGTH OF THE INTERVAL */
/*           ABSC   - ABSCISSA */
/*           FVAL*  - FUNCTION VALUE */
/*           RESG   - RESULT OF THE 20-POINT GAUSS FORMULA */
/*           RESK   - RESULT OF THE 41-POINT KRONROD FORMULA */
/*           RESKH  - APPROXIMATION TO MEAN VALUE OF F OVER (A,B), I.E. */
/*                    TO I/(B-A) */

/*           MACHINE DEPENDENT CONSTANTS */
/*           --------------------------- */

/*           EPMACH IS THE LARGEST RELATIVE SPACING. */
/*           UFLOW IS THE SMALLEST POSITIVE MAGNITUDE. */

/* ***FIRST EXECUTABLE STATEMENT  DQK41 */
    epmach = d1mach_(&c__4);
    uflow = d1mach_(&c__1);

    centr = (*a + *b) * .5;
    hlgth = (*b - *a) * .5;
    dhlgth = abs(hlgth);

/*           COMPUTE THE 41-POINT GAUSS-KRONROD APPROXIMATION TO */
/*           THE INTEGRAL, AND ESTIMATE THE ABSOLUTE ERROR. */

    resg = 0.;
    fc = (*f)(&centr);
    resk = wgk[20] * fc;
    *resabs = abs(resk);
    for (j = 1; j <= 10; ++j) {
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
    for (j = 1; j <= 10; ++j) {
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
    *resasc = wgk[20] * (d__1 = fc - reskh, abs(d__1));
    for (j = 1; j <= 20; ++j) {
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
} /* dqk41_ */

