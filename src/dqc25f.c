/* dqc25f.f -- translated by f2c (version 12.02.01).
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

static integer c__2 = 2;

/* DECK DQC25F */
/* Subroutine */ int dqc25f_(D_fp f, doublereal *a, doublereal *b, doublereal 
	*omega, integer *integr, integer *nrmom, integer *maxp1, integer *
	ksave, doublereal *result, doublereal *abserr, integer *neval, 
	doublereal *resabs, doublereal *resasc, integer *momcom, doublereal *
	chebmo)
{
    /* Initialized data */

    static doublereal x[11] = { .9914448613738104,.9659258262890683,
	    .9238795325112868,.8660254037844386,.7933533402912352,
	    .7071067811865475,.6087614290087206,.5,.3826834323650898,
	    .2588190451025208,.1305261922200516 };

    /* System generated locals */
    integer chebmo_dim1, chebmo_offset, i__1;
    doublereal d__1, d__2;

    /* Local variables */
    static doublereal d__[25];
    static integer i__, j, k, m;
    static doublereal v[28], d1[25], d2[25], p2, p3, p4, ac, an, as, an2, ass,
	     par2, conc, asap, par22, fval[25], estc, cons;
    static integer iers;
    static doublereal ests;
    static integer isym, noeq1;
    static doublereal cheb12[13], cheb24[25], resc12, resc24, hlgth, centr;
    extern /* Subroutine */ int dqk15w_(D_fp, D_fp, doublereal *, doublereal *
	    , doublereal *, doublereal *, integer *, doublereal *, doublereal 
	    *, doublereal *, doublereal *, doublereal *, doublereal *), 
	    dgtsl_(integer *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, integer *);
    static doublereal ress12, ress24, oflow;
    static integer noequ;
    extern doublereal d1mach_(integer *);
    extern /* Subroutine */ int dqcheb_(doublereal *, doublereal *, 
	    doublereal *, doublereal *);
    static doublereal cospar;
    extern doublereal dqwgtf_();
    static doublereal parint, sinpar;

/* ***BEGIN PROLOGUE  DQC25F */
/* ***PURPOSE  To compute the integral I=Integral of F(X) over (A,B) */
/*            Where W(X) = COS(OMEGA*X) or W(X)=SIN(OMEGA*X) and to */
/*            compute J = Integral of ABS(F) over (A,B). For small value */
/*            of OMEGA or small intervals (A,B) the 15-point GAUSS-KRONRO */
/*            Rule is used. Otherwise a generalized CLENSHAW-CURTIS */
/*            method is used. */
/* ***LIBRARY   SLATEC (QUADPACK) */
/* ***CATEGORY  H2A2A2 */
/* ***TYPE      DOUBLE PRECISION (QC25F-S, DQC25F-D) */
/* ***KEYWORDS  CLENSHAW-CURTIS METHOD, GAUSS-KRONROD RULES, */
/*             INTEGRATION RULES FOR FUNCTIONS WITH COS OR SIN FACTOR, */
/*             QUADPACK, QUADRATURE */
/* ***AUTHOR  Piessens, Robert */
/*             Applied Mathematics and Programming Division */
/*             K. U. Leuven */
/*           de Doncker, Elise */
/*             Applied Mathematics and Programming Division */
/*             K. U. Leuven */
/* ***DESCRIPTION */

/*        Integration rules for functions with COS or SIN factor */
/*        Standard fortran subroutine */
/*        Double precision version */

/*        PARAMETERS */
/*         ON ENTRY */
/*           F      - Double precision */
/*                    Function subprogram defining the integrand */
/*                    function F(X). The actual name for F needs to */
/*                    be declared E X T E R N A L in the calling program. */

/*           A      - Double precision */
/*                    Lower limit of integration */

/*           B      - Double precision */
/*                    Upper limit of integration */

/*           OMEGA  - Double precision */
/*                    Parameter in the WEIGHT function */

/*           INTEGR - Integer */
/*                    Indicates which WEIGHT function is to be used */
/*                       INTEGR = 1   W(X) = COS(OMEGA*X) */
/*                       INTEGR = 2   W(X) = SIN(OMEGA*X) */

/*           NRMOM  - Integer */
/*                    The length of interval (A,B) is equal to the length */
/*                    of the original integration interval divided by */
/*                    2**NRMOM (we suppose that the routine is used in an */
/*                    adaptive integration process, otherwise set */
/*                    NRMOM = 0). NRMOM must be zero at the first call. */

/*           MAXP1  - Integer */
/*                    Gives an upper bound on the number of Chebyshev */
/*                    moments which can be stored, i.e. for the */
/*                    intervals of lengths ABS(BB-AA)*2**(-L), */
/*                    L = 0,1,2, ..., MAXP1-2. */

/*           KSAVE  - Integer */
/*                    Key which is one when the moments for the */
/*                    current interval have been computed */

/*         ON RETURN */
/*           RESULT - Double precision */
/*                    Approximation to the integral I */

/*           ABSERR - Double precision */
/*                    Estimate of the modulus of the absolute */
/*                    error, which should equal or exceed ABS(I-RESULT) */

/*           NEVAL  - Integer */
/*                    Number of integrand evaluations */

/*           RESABS - Double precision */
/*                    Approximation to the integral J */

/*           RESASC - Double precision */
/*                    Approximation to the integral of ABS(F-I/(B-A)) */

/*         ON ENTRY AND RETURN */
/*           MOMCOM - Integer */
/*                    For each interval length we need to compute the */
/*                    Chebyshev moments. MOMCOM counts the number of */
/*                    intervals for which these moments have already been */
/*                    computed. If NRMOM.LT.MOMCOM or KSAVE = 1, the */
/*                    Chebyshev moments for the interval (A,B) have */
/*                    already been computed and stored, otherwise we */
/*                    compute them and we increase MOMCOM. */

/*           CHEBMO - Double precision */
/*                    Array of dimension at least (MAXP1,25) containing */
/*                    the modified Chebyshev moments for the first MOMCOM */
/*                    MOMCOM interval lengths */

/* ...................................................................... */

/* ***REFERENCES  (NONE) */
/* ***ROUTINES CALLED  D1MACH, DGTSL, DQCHEB, DQK15W, DQWGTF */
/* ***REVISION HISTORY  (YYMMDD) */
/*   810101  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890531  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/* ***END PROLOGUE  DQC25F */




/*           THE VECTOR X CONTAINS THE VALUES COS(K*PI/24) */
/*           K = 1, ...,11, TO BE USED FOR THE CHEBYSHEV EXPANSION OF F */

    /* Parameter adjustments */
    chebmo_dim1 = *maxp1;
    chebmo_offset = 1 + chebmo_dim1;
    chebmo -= chebmo_offset;

    /* Function Body */

/*           LIST OF MAJOR VARIABLES */
/*           ----------------------- */

/*           CENTR  - MID POINT OF THE INTEGRATION INTERVAL */
/*           HLGTH  - HALF-LENGTH OF THE INTEGRATION INTERVAL */
/*           FVAL   - VALUE OF THE FUNCTION F AT THE POINTS */
/*                    (B-A)*0.5*COS(K*PI/12) + (B+A)*0.5, K = 0, ..., 24 */
/*           CHEB12 - COEFFICIENTS OF THE CHEBYSHEV SERIES EXPANSION */
/*                    OF DEGREE 12, FOR THE FUNCTION F, IN THE */
/*                    INTERVAL (A,B) */
/*           CHEB24 - COEFFICIENTS OF THE CHEBYSHEV SERIES EXPANSION */
/*                    OF DEGREE 24, FOR THE FUNCTION F, IN THE */
/*                    INTERVAL (A,B) */
/*           RESC12 - APPROXIMATION TO THE INTEGRAL OF */
/*                    COS(0.5*(B-A)*OMEGA*X)*F(0.5*(B-A)*X+0.5*(B+A)) */
/*                    OVER (-1,+1), USING THE CHEBYSHEV SERIES */
/*                    EXPANSION OF DEGREE 12 */
/*           RESC24 - APPROXIMATION TO THE SAME INTEGRAL, USING THE */
/*                    CHEBYSHEV SERIES EXPANSION OF DEGREE 24 */
/*           RESS12 - THE ANALOGUE OF RESC12 FOR THE SINE */
/*           RESS24 - THE ANALOGUE OF RESC24 FOR THE SINE */


/*           MACHINE DEPENDENT CONSTANT */
/*           -------------------------- */

/*           OFLOW IS THE LARGEST POSITIVE MAGNITUDE. */

/* ***FIRST EXECUTABLE STATEMENT  DQC25F */
    oflow = d1mach_(&c__2);

    centr = (*b + *a) * .5;
    hlgth = (*b - *a) * .5;
    parint = *omega * hlgth;

/*           COMPUTE THE INTEGRAL USING THE 15-POINT GAUSS-KRONROD */
/*           FORMULA IF THE VALUE OF THE PARAMETER IN THE INTEGRAND */
/*           IS SMALL. */

    if (abs(parint) > 2.) {
	goto L10;
    }
    dqk15w_((D_fp)f, (D_fp)dqwgtf_, omega, &p2, &p3, &p4, integr, a, b, 
	    result, abserr, resabs, resasc);
    *neval = 15;
    goto L170;

/*           COMPUTE THE INTEGRAL USING THE GENERALIZED CLENSHAW- */
/*           CURTIS METHOD. */

L10:
    conc = hlgth * cos(centr * *omega);
    cons = hlgth * sin(centr * *omega);
    *resasc = oflow;
    *neval = 25;

/*           CHECK WHETHER THE CHEBYSHEV MOMENTS FOR THIS INTERVAL */
/*           HAVE ALREADY BEEN COMPUTED. */

    if (*nrmom < *momcom || *ksave == 1) {
	goto L120;
    }

/*           COMPUTE A NEW SET OF CHEBYSHEV MOMENTS. */

    m = *momcom + 1;
    par2 = parint * parint;
    par22 = par2 + 2.;
    sinpar = sin(parint);
    cospar = cos(parint);

/*           COMPUTE THE CHEBYSHEV MOMENTS WITH RESPECT TO COSINE. */

    v[0] = sinpar * 2. / parint;
    v[1] = (cospar * 8. + (par2 + par2 - 8.) * sinpar / parint) / par2;
    v[2] = ((par2 - 12.) * 32. * cospar + ((par2 - 80.) * par2 + 192.) * 2. * 
	    sinpar / parint) / (par2 * par2);
    ac = cospar * 8.;
    as = parint * 24. * sinpar;
    if (abs(parint) > 24.) {
	goto L30;
    }

/*           COMPUTE THE CHEBYSHEV MOMENTS AS THE SOLUTIONS OF A */
/*           BOUNDARY VALUE PROBLEM WITH 1 INITIAL VALUE (V(3)) AND 1 */
/*           END VALUE (COMPUTED USING AN ASYMPTOTIC FORMULA). */

    noequ = 25;
    noeq1 = noequ - 1;
    an = 6.;
    i__1 = noeq1;
    for (k = 1; k <= i__1; ++k) {
	an2 = an * an;
	d__[k - 1] = (an2 - 4.) * -2. * (par22 - an2 - an2);
	d2[k - 1] = (an - 1.) * (an - 2.) * par2;
	d1[k] = (an + 3.) * (an + 4.) * par2;
	v[k + 2] = as - (an2 - 4.) * ac;
	an += 2.;
/* L20: */
    }
    an2 = an * an;
    d__[noequ - 1] = (an2 - 4.) * -2. * (par22 - an2 - an2);
    v[noequ + 2] = as - (an2 - 4.) * ac;
    v[3] -= par2 * 56. * v[2];
    ass = parint * sinpar;
    asap = (((((par2 * 210. - 1.) * cospar - (par2 * 105. - 63.) * ass) / an2 
	    - (1. - par2 * 15.) * cospar + ass * 15.) / an2 - cospar + ass * 
	    3.) / an2 - cospar) / an2;
    v[noequ + 2] -= asap * 2. * par2 * (an - 1.) * (an - 2.);

/*           SOLVE THE TRIDIAGONAL SYSTEM BY MEANS OF GAUSSIAN */
/*           ELIMINATION WITH PARTIAL PIVOTING. */

/* ***       CALL TO DGTSL MUST BE REPLACED BY CALL TO */
/* ***       DOUBLE PRECISION VERSION OF LINPACK ROUTINE SGTSL */

    dgtsl_(&noequ, d1, d__, d2, &v[3], &iers);
    goto L50;

/*           COMPUTE THE CHEBYSHEV MOMENTS BY MEANS OF FORWARD */
/*           RECURSION. */

L30:
    an = 4.;
    for (i__ = 4; i__ <= 13; ++i__) {
	an2 = an * an;
	v[i__ - 1] = ((an2 - 4.) * ((par22 - an2 - an2) * 2. * v[i__ - 2] - 
		ac) + as - par2 * (an + 1.) * (an + 2.) * v[i__ - 3]) / (par2 
		* (an - 1.) * (an - 2.));
	an += 2.;
/* L40: */
    }
L50:
    for (j = 1; j <= 13; ++j) {
	chebmo[m + ((j << 1) - 1) * chebmo_dim1] = v[j - 1];
/* L60: */
    }

/*           COMPUTE THE CHEBYSHEV MOMENTS WITH RESPECT TO SINE. */

    v[0] = (sinpar - parint * cospar) * 2. / par2;
    v[1] = (18. - 48. / par2) * sinpar / par2 + (48. / par2 - 2.) * cospar / 
	    parint;
    ac = parint * -24. * cospar;
    as = sinpar * -8.;
    if (abs(parint) > 24.) {
	goto L80;
    }

/*           COMPUTE THE CHEBYSHEV MOMENTS AS THE SOLUTIONS OF A BOUNDARY */
/*           VALUE PROBLEM WITH 1 INITIAL VALUE (V(2)) AND 1 END VALUE */
/*           (COMPUTED USING AN ASYMPTOTIC FORMULA). */

    an = 5.;
    i__1 = noeq1;
    for (k = 1; k <= i__1; ++k) {
	an2 = an * an;
	d__[k - 1] = (an2 - 4.) * -2. * (par22 - an2 - an2);
	d2[k - 1] = (an - 1.) * (an - 2.) * par2;
	d1[k] = (an + 3.) * (an + 4.) * par2;
	v[k + 1] = ac + (an2 - 4.) * as;
	an += 2.;
/* L70: */
    }
    an2 = an * an;
    d__[noequ - 1] = (an2 - 4.) * -2. * (par22 - an2 - an2);
    v[noequ + 1] = ac + (an2 - 4.) * as;
    v[2] -= par2 * 42. * v[1];
    ass = parint * cospar;
    asap = (((((par2 * 105. - 63.) * ass + (par2 * 210. - 1.) * sinpar) / an2 
	    + (par2 * 15. - 1.) * sinpar - ass * 15.) / an2 - ass * 3. - 
	    sinpar) / an2 - sinpar) / an2;
    v[noequ + 1] -= asap * 2. * par2 * (an - 1.) * (an - 2.);

/*           SOLVE THE TRIDIAGONAL SYSTEM BY MEANS OF GAUSSIAN */
/*           ELIMINATION WITH PARTIAL PIVOTING. */

/* ***       CALL TO DGTSL MUST BE REPLACED BY CALL TO */
/* ***       DOUBLE PRECISION VERSION OF LINPACK ROUTINE SGTSL */

    dgtsl_(&noequ, d1, d__, d2, &v[2], &iers);
    goto L100;

/*           COMPUTE THE CHEBYSHEV MOMENTS BY MEANS OF FORWARD RECURSION. */

L80:
    an = 3.;
    for (i__ = 3; i__ <= 12; ++i__) {
	an2 = an * an;
	v[i__ - 1] = ((an2 - 4.) * ((par22 - an2 - an2) * 2. * v[i__ - 2] + 
		as) + ac - par2 * (an + 1.) * (an + 2.) * v[i__ - 3]) / (par2 
		* (an - 1.) * (an - 2.));
	an += 2.;
/* L90: */
    }
L100:
    for (j = 1; j <= 12; ++j) {
	chebmo[m + (j << 1) * chebmo_dim1] = v[j - 1];
/* L110: */
    }
L120:
    if (*nrmom < *momcom) {
	m = *nrmom + 1;
    }
    if (*momcom < *maxp1 - 1 && *nrmom >= *momcom) {
	++(*momcom);
    }

/*           COMPUTE THE COEFFICIENTS OF THE CHEBYSHEV EXPANSIONS */
/*           OF DEGREES 12 AND 24 OF THE FUNCTION F. */

    d__1 = centr + hlgth;
    fval[0] = (*f)(&d__1) * .5;
    fval[12] = (*f)(&centr);
    d__1 = centr - hlgth;
    fval[24] = (*f)(&d__1) * .5;
    for (i__ = 2; i__ <= 12; ++i__) {
	isym = 26 - i__;
	d__1 = hlgth * x[i__ - 2] + centr;
	fval[i__ - 1] = (*f)(&d__1);
	d__1 = centr - hlgth * x[i__ - 2];
	fval[isym - 1] = (*f)(&d__1);
/* L130: */
    }
    dqcheb_(x, fval, cheb12, cheb24);

/*           COMPUTE THE INTEGRAL AND ERROR ESTIMATES. */

    resc12 = cheb12[12] * chebmo[m + chebmo_dim1 * 13];
    ress12 = 0.;
    k = 11;
    for (j = 1; j <= 6; ++j) {
	resc12 += cheb12[k - 1] * chebmo[m + k * chebmo_dim1];
	ress12 += cheb12[k] * chebmo[m + (k + 1) * chebmo_dim1];
	k += -2;
/* L140: */
    }
    resc24 = cheb24[24] * chebmo[m + chebmo_dim1 * 25];
    ress24 = 0.;
    *resabs = abs(cheb24[24]);
    k = 23;
    for (j = 1; j <= 12; ++j) {
	resc24 += cheb24[k - 1] * chebmo[m + k * chebmo_dim1];
	ress24 += cheb24[k] * chebmo[m + (k + 1) * chebmo_dim1];
	*resabs = (d__1 = cheb24[k - 1], abs(d__1)) + (d__2 = cheb24[k], abs(
		d__2));
	k += -2;
/* L150: */
    }
    estc = (d__1 = resc24 - resc12, abs(d__1));
    ests = (d__1 = ress24 - ress12, abs(d__1));
    *resabs *= abs(hlgth);
    if (*integr == 2) {
	goto L160;
    }
    *result = conc * resc24 - cons * ress24;
    *abserr = (d__1 = conc * estc, abs(d__1)) + (d__2 = cons * ests, abs(d__2)
	    );
    goto L170;
L160:
    *result = conc * ress24 + cons * resc24;
    *abserr = (d__1 = conc * ests, abs(d__1)) + (d__2 = cons * estc, abs(d__2)
	    );
L170:
    return 0;
} /* dqc25f_ */

