/* dqk31.f -- translated by f2c (version 12.02.01).
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

/* DECK DQK31 */
/* Subroutine */ int dqk31_(D_fp f, doublereal *a, doublereal *b, doublereal *
	result, doublereal *abserr, doublereal *resabs, doublereal *resasc)
{
    /* Initialized data */

    static doublereal wg[8] = { .030753241996117268354628393577204,
	    .070366047488108124709267416450667,
	    .107159220467171935011869546685869,
	    .139570677926154314447804794511028,
	    .166269205816993933553200860481209,
	    .186161000015562211026800561866423,
	    .198431485327111576456118326443839,
	    .202578241925561272880620199967519 };
    static doublereal xgk[16] = { .998002298693397060285172840152271,
	    .987992518020485428489565718586613,
	    .967739075679139134257347978784337,
	    .937273392400705904307758947710209,
	    .897264532344081900882509656454496,
	    .848206583410427216200648320774217,
	    .790418501442465932967649294817947,
	    .724417731360170047416186054613938,
	    .650996741297416970533735895313275,
	    .570972172608538847537226737253911,
	    .485081863640239680693655740232351,
	    .394151347077563369897207370981045,
	    .299180007153168812166780024266389,
	    .201194093997434522300628303394596,
	    .101142066918717499027074231447392,0. };
    static doublereal wgk[16] = { .005377479872923348987792051430128,
	    .015007947329316122538374763075807,
	    .025460847326715320186874001019653,
	    .03534636079137584622203794847836,
	    .04458975132476487660822729937328,
	    .05348152469092808726534314723943,
	    .062009567800670640285139230960803,
	    .069854121318728258709520077099147,
	    .076849680757720378894432777482659,
	    .083080502823133021038289247286104,
	    .088564443056211770647275443693774,
	    .093126598170825321225486872747346,
	    .096642726983623678505179907627589,
	    .099173598721791959332393173484603,
	    .10076984552387559504494666261757,
	    .101330007014791549017374792767493 };

    /* System generated locals */
    doublereal d__1, d__2, d__3;

    /* Local variables */
    static integer j;
    static doublereal fc, fv1[15], fv2[15];
    static integer jtw;
    static doublereal absc, resg, resk, fsum, fval1, fval2;
    static integer jtwm1;
    static doublereal hlgth, centr, reskh, uflow;
    extern doublereal d1mach_(integer *);
    static doublereal epmach, dhlgth;

/* ***BEGIN PROLOGUE  DQK31 */
/* ***PURPOSE  To compute I = Integral of F over (A,B) with error */
/*                           estimate */
/*                       J = Integral of ABS(F) over (A,B) */
/* ***LIBRARY   SLATEC (QUADPACK) */
/* ***CATEGORY  H2A1A2 */
/* ***TYPE      DOUBLE PRECISION (QK31-S, DQK31-D) */
/* ***KEYWORDS  31-POINT GAUSS-KRONROD RULES, QUADPACK, QUADRATURE */
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
/*                       Declared E X T E R N A L in the calling program. */

/*              A      - Double precision */
/*                       Lower limit of integration */

/*              B      - Double precision */
/*                       Upper limit of integration */

/*            ON RETURN */
/*              RESULT - Double precision */
/*                       Approximation to the integral I */
/*                       RESULT is computed by applying the 31-POINT */
/*                       GAUSS-KRONROD RULE (RESK), obtained by optimal */
/*                       addition of abscissae to the 15-POINT GAUSS */
/*                       RULE (RESG). */

/*              ABSERR - Double precision */
/*                       Estimate of the modulus of the modulus, */
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
/* ***END PROLOGUE  DQK31 */


/*           THE ABSCISSAE AND WEIGHTS ARE GIVEN FOR THE INTERVAL (-1,1). */
/*           BECAUSE OF SYMMETRY ONLY THE POSITIVE ABSCISSAE AND THEIR */
/*           CORRESPONDING WEIGHTS ARE GIVEN. */

/*           XGK    - ABSCISSAE OF THE 31-POINT KRONROD RULE */
/*                    XGK(2), XGK(4), ...  ABSCISSAE OF THE 15-POINT */
/*                    GAUSS RULE */
/*                    XGK(1), XGK(3), ...  ABSCISSAE WHICH ARE OPTIMALLY */
/*                    ADDED TO THE 15-POINT GAUSS RULE */

/*           WGK    - WEIGHTS OF THE 31-POINT KRONROD RULE */

/*           WG     - WEIGHTS OF THE 15-POINT GAUSS RULE */


/* GAUSS QUADRATURE WEIGHTS AND KRONROD QUADRATURE ABSCISSAE AND WEIGHTS */
/* AS EVALUATED WITH 80 DECIMAL DIGIT ARITHMETIC BY L. W. FULLERTON, */
/* BELL LABS, NOV. 1981. */





/*           LIST OF MAJOR VARIABLES */
/*           ----------------------- */
/*           CENTR  - MID POINT OF THE INTERVAL */
/*           HLGTH  - HALF-LENGTH OF THE INTERVAL */
/*           ABSC   - ABSCISSA */
/*           FVAL*  - FUNCTION VALUE */
/*           RESG   - RESULT OF THE 15-POINT GAUSS FORMULA */
/*           RESK   - RESULT OF THE 31-POINT KRONROD FORMULA */
/*           RESKH  - APPROXIMATION TO THE MEAN VALUE OF F OVER (A,B), */
/*                    I.E. TO I/(B-A) */

/*           MACHINE DEPENDENT CONSTANTS */
/*           --------------------------- */
/*           EPMACH IS THE LARGEST RELATIVE SPACING. */
/*           UFLOW IS THE SMALLEST POSITIVE MAGNITUDE. */
/* ***FIRST EXECUTABLE STATEMENT  DQK31 */
    epmach = d1mach_(&c__4);
    uflow = d1mach_(&c__1);

    centr = (*a + *b) * .5;
    hlgth = (*b - *a) * .5;
    dhlgth = abs(hlgth);

/*           COMPUTE THE 31-POINT KRONROD APPROXIMATION TO */
/*           THE INTEGRAL, AND ESTIMATE THE ABSOLUTE ERROR. */

    fc = (*f)(&centr);
    resg = wg[7] * fc;
    resk = wgk[15] * fc;
    *resabs = abs(resk);
    for (j = 1; j <= 7; ++j) {
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
    for (j = 1; j <= 8; ++j) {
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
    *resasc = wgk[15] * (d__1 = fc - reskh, abs(d__1));
    for (j = 1; j <= 15; ++j) {
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
} /* dqk31_ */

