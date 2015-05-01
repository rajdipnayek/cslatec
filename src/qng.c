/* qng.f -- translated by f2c (version 12.02.01).
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

/* DECK QNG */
/* Subroutine */ int qng_(E_fp f, real *a, real *b, real *epsabs, real *
	epsrel, real *result, real *abserr, integer *neval, integer *ier)
{
    /* Initialized data */

    static real x1[5] = { .9739065285171717f,.8650633666889845f,
	    .6794095682990244f,.4333953941292472f,.1488743389816312f };
    static real w87a[21] = { .008148377384149173f,.01876143820156282f,
	    .02734745105005229f,.03367770731163793f,.03693509982042791f,
	    .002884872430211531f,.0136859460227127f,.02328041350288831f,
	    .03087249761171336f,.03569363363941877f,9.152833452022414e-4f,
	    .005399280219300471f,.01094767960111893f,.01629873169678734f,
	    .02108156888920384f,.02537096976925383f,.02918969775647575f,
	    .03237320246720279f,.03478309895036514f,.03641222073135179f,
	    .03725387550304771f };
    static real w87b[23] = { 2.741455637620724e-4f,.001807124155057943f,
	    .004096869282759165f,.006758290051847379f,.009549957672201647f,
	    .01232944765224485f,.01501044734638895f,.01754896798624319f,
	    .01993803778644089f,.02219493596101229f,.02433914712600081f,
	    .02637450541483921f,.0282869107887712f,.0300525811280927f,
	    .03164675137143993f,.0330504134199785f,.03425509970422606f,
	    .03526241266015668f,.0360769896228887f,.03669860449845609f,
	    .03712054926983258f,.03733422875193504f,.03736107376267902f };
    static real x2[5] = { .9956571630258081f,.9301574913557082f,
	    .7808177265864169f,.5627571346686047f,.2943928627014602f };
    static real x3[11] = { .9993333609019321f,.9874334029080889f,
	    .9548079348142663f,.9001486957483283f,.8251983149831142f,
	    .732148388989305f,.6228479705377252f,.4994795740710565f,
	    .3649016613465808f,.2222549197766013f,.07465061746138332f };
    static real x4[22] = { .9999029772627292f,.9979898959866787f,
	    .9921754978606872f,.9813581635727128f,.9650576238583846f,
	    .9431676131336706f,.9158064146855072f,.8832216577713165f,
	    .8457107484624157f,.803557658035231f,.7570057306854956f,
	    .7062732097873218f,.6515894665011779f,.5932233740579611f,
	    .5314936059708319f,.4667636230420228f,.3994248478592188f,
	    .3298748771061883f,.2585035592021616f,.1856953965683467f,
	    .1118422131799075f,.03735212339461987f };
    static real w10[5] = { .06667134430868814f,.1494513491505806f,
	    .219086362515982f,.2692667193099964f,.2955242247147529f };
    static real w21a[5] = { .03255816230796473f,.07503967481091995f,
	    .1093871588022976f,.1347092173114733f,.1477391049013385f };
    static real w21b[6] = { .01169463886737187f,.054755896574352f,
	    .09312545458369761f,.1234919762620659f,.1427759385770601f,
	    .1494455540029169f };
    static real w43a[10] = { .01629673428966656f,.0375228761208695f,
	    .05469490205825544f,.06735541460947809f,.07387019963239395f,
	    .005768556059769796f,.02737189059324884f,.04656082691042883f,
	    .06174499520144256f,.0713872672686934f };
    static real w43b[12] = { .001844477640212414f,.01079868958589165f,
	    .02189536386779543f,.03259746397534569f,.04216313793519181f,
	    .05074193960018458f,.05837939554261925f,.06474640495144589f,
	    .06956619791235648f,.07282444147183321f,.07450775101417512f,
	    .07472214751740301f };

    /* System generated locals */
    real r__1, r__2, r__3, r__4;
    doublereal d__1;

    /* Local variables */
    static integer k, l;
    static real fv1[5], fv2[5], fv3[5], fv4[5];
    static integer ipx;
    static real absc, fval, res10, res21, res43, res87, fval1, fval2, hlgth, 
	    centr, reskh, uflow;
    extern doublereal r1mach_(integer *);
    static real epmach, dhlgth, resabs, resasc, fcentr, savfun[21];
    extern /* Subroutine */ int xermsg_(char *, char *, char *, integer *, 
	    integer *, ftnlen, ftnlen, ftnlen);

/* ***BEGIN PROLOGUE  QNG */
/* ***PURPOSE  The routine calculates an approximation result to a */
/*            given definite integral I = integral of F over (A,B), */
/*            hopefully satisfying following claim for accuracy */
/*            ABS(I-RESULT).LE.MAX(EPSABS,EPSREL*ABS(I)). */
/* ***LIBRARY   SLATEC (QUADPACK) */
/* ***CATEGORY  H2A1A1 */
/* ***TYPE      SINGLE PRECISION (QNG-S, DQNG-D) */
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
/* REAL VERSION */

/*           F      - Real version */
/*                    Function subprogram defining the integrand function */
/*                    F(X). The actual name for F needs to be declared */
/*                    E X T E R N A L in the driver program. */

/*           A      - Real version */
/*                    Lower limit of integration */

/*           B      - Real version */
/*                    Upper limit of integration */

/*           EPSABS - Real */
/*                    Absolute accuracy requested */
/*           EPSREL - Real */
/*                    Relative accuracy requested */
/*                    If  EPSABS.LE.0 */
/*                    And EPSREL.LT.MAX(50*REL.MACH.ACC.,0.5D-28), */
/*                    The routine will end with IER = 6. */

/*         ON RETURN */
/*           RESULT - Real */
/*                    Approximation to the integral I */
/*                    Result is obtained by applying the 21-POINT */
/*                    GAUSS-KRONROD RULE (RES21) obtained by optimal */
/*                    addition of abscissae to the 10-POINT GAUSS RULE */
/*                    (RES10), or by applying the 43-POINT RULE (RES43) */
/*                    obtained by optimal addition of abscissae to the */
/*                    21-POINT GAUSS-KRONROD RULE, or by applying the */
/*                    87-POINT RULE (RES87) obtained by optimal addition */
/*                    of abscissae to the 43-POINT RULE. */

/*           ABSERR - Real */
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
/* ***ROUTINES CALLED  R1MACH, XERMSG */
/* ***REVISION HISTORY  (YYMMDD) */
/*   800101  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890531  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ) */
/* ***END PROLOGUE  QNG */



/*           THE FOLLOWING DATA STATEMENTS CONTAIN THE */
/*           ABSCISSAE AND WEIGHTS OF THE INTEGRATION RULES USED. */

/*           X1      ABSCISSAE COMMON TO THE 10-, 21-, 43- */
/*                   AND 87-POINT RULE */
/*           X2      ABSCISSAE COMMON TO THE 21-, 43- AND */
/*                   87-POINT RULE */
/*           X3      ABSCISSAE COMMON TO THE 43- AND 87-POINT */
/*                   RULE */
/*           X4      ABSCISSAE OF THE 87-POINT RULE */
/*           W10     WEIGHTS OF THE 10-POINT FORMULA */
/*           W21A    WEIGHTS OF THE 21-POINT FORMULA FOR */
/*                   ABSCISSAE X1 */
/*           W21B    WEIGHTS OF THE 21-POINT FORMULA FOR */
/*                   ABSCISSAE X2 */
/*           W43A    WEIGHTS OF THE 43-POINT FORMULA FOR */
/*                   ABSCISSAE X1, X3 */
/*           W43B    WEIGHTS OF THE 43-POINT FORMULA FOR */
/*                   ABSCISSAE X3 */
/*           W87A    WEIGHTS OF THE 87-POINT FORMULA FOR */
/*                   ABSCISSAE X1, X2, X3 */
/*           W87B    WEIGHTS OF THE 87-POINT FORMULA FOR */
/*                   ABSCISSAE X4 */


/*           LIST OF MAJOR VARIABLES */
/*           ----------------------- */

/*           CENTR  - MID POINT OF THE INTEGRATION INTERVAL */
/*           HLGTH  - HALF-LENGTH OF THE INTEGRATION INTERVAL */
/*           FCENTR - FUNCTION VALUE AT MID POINT */
/*           ABSC   - ABSCISSA */
/*           FVAL   - FUNCTION VALUE */
/*           SAVFUN - ARRAY OF FUNCTION VALUES WHICH */
/*                    HAVE ALREADY BEEN COMPUTED */
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

/* ***FIRST EXECUTABLE STATEMENT  QNG */
    epmach = r1mach_(&c__4);
    uflow = r1mach_(&c__1);

/*           TEST ON VALIDITY OF PARAMETERS */
/*           ------------------------------ */

    *result = 0.f;
    *abserr = 0.f;
    *neval = 0;
    *ier = 6;
/* Computing MAX */
    r__1 = 5e-15f, r__2 = epmach * 50.f;
    if (*epsabs <= 0.f && *epsrel < dmax(r__1,r__2)) {
	goto L80;
    }
    hlgth = (*b - *a) * .5f;
    dhlgth = dabs(hlgth);
    centr = (*b + *a) * .5f;
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
	res10 = 0.f;
	res21 = w21b[5] * fcentr;
	resabs = w21b[5] * dabs(fcentr);
	for (k = 1; k <= 5; ++k) {
	    absc = hlgth * x1[k - 1];
	    r__1 = centr + absc;
	    fval1 = (*f)(&r__1);
	    r__1 = centr - absc;
	    fval2 = (*f)(&r__1);
	    fval = fval1 + fval2;
	    res10 += w10[k - 1] * fval;
	    res21 += w21a[k - 1] * fval;
	    resabs += w21a[k - 1] * (dabs(fval1) + dabs(fval2));
	    savfun[k - 1] = fval;
	    fv1[k - 1] = fval1;
	    fv2[k - 1] = fval2;
/* L10: */
	}
	ipx = 5;
	for (k = 1; k <= 5; ++k) {
	    ++ipx;
	    absc = hlgth * x2[k - 1];
	    r__1 = centr + absc;
	    fval1 = (*f)(&r__1);
	    r__1 = centr - absc;
	    fval2 = (*f)(&r__1);
	    fval = fval1 + fval2;
	    res21 += w21b[k - 1] * fval;
	    resabs += w21b[k - 1] * (dabs(fval1) + dabs(fval2));
	    savfun[ipx - 1] = fval;
	    fv3[k - 1] = fval1;
	    fv4[k - 1] = fval2;
/* L15: */
	}

/*          TEST FOR CONVERGENCE. */

	*result = res21 * hlgth;
	resabs *= dhlgth;
	reskh = res21 * .5f;
	resasc = w21b[5] * (r__1 = fcentr - reskh, dabs(r__1));
	for (k = 1; k <= 5; ++k) {
	    resasc = resasc + w21a[k - 1] * ((r__1 = fv1[k - 1] - reskh, dabs(
		    r__1)) + (r__2 = fv2[k - 1] - reskh, dabs(r__2))) + w21b[
		    k - 1] * ((r__3 = fv3[k - 1] - reskh, dabs(r__3)) + (r__4 
		    = fv4[k - 1] - reskh, dabs(r__4)));
/* L20: */
	}
	*abserr = (r__1 = (res21 - res10) * hlgth, dabs(r__1));
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
	    r__1 = absc + centr;
	    r__2 = centr - absc;
	    fval = (*f)(&r__1) + (*f)(&r__2);
	    res43 += fval * w43b[k - 1];
	    savfun[ipx - 1] = fval;
/* L40: */
	}

/*          TEST FOR CONVERGENCE. */

	*result = res43 * hlgth;
	*abserr = (r__1 = (res43 - res21) * hlgth, dabs(r__1));
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
	    r__1 = absc + centr;
	    r__2 = centr - absc;
	    res87 += w87b[k - 1] * ((*f)(&r__1) + (*f)(&r__2));
/* L60: */
	}
	*result = res87 * hlgth;
	*abserr = (r__1 = (res87 - res43) * hlgth, dabs(r__1));
L65:
	if (resasc != 0.f && *abserr != 0.f) {
/* Computing MIN */
	    d__1 = (doublereal) (*abserr * 200.f / resasc);
	    r__1 = 1.f, r__2 = pow_dd(&d__1, &c_b17);
	    *abserr = resasc * dmin(r__1,r__2);
	}
	if (resabs > uflow / (epmach * 50.f)) {
/* Computing MAX */
	    r__1 = epmach * 50.f * resabs;
	    *abserr = dmax(r__1,*abserr);
	}
/* Computing MAX */
	r__1 = *epsabs, r__2 = *epsrel * dabs(*result);
	if (*abserr <= dmax(r__1,r__2)) {
	    *ier = 0;
	}
/* ***JUMP OUT OF DO-LOOP */
	if (*ier == 0) {
	    goto L999;
	}
/* L70: */
    }
L80:
    xermsg_("SLATEC", "QNG", "ABNORMAL RETURN", ier, &c__0, (ftnlen)6, (
	    ftnlen)3, (ftnlen)15);
L999:
    return 0;
} /* qng_ */

