/* qc25f.f -- translated by f2c (version 12.02.01).
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

/* DECK QC25F */
/* Subroutine */ int qc25f_(E_fp f, real *a, real *b, real *omega, integer *
	integr, integer *nrmom, integer *maxp1, integer *ksave, real *result, 
	real *abserr, integer *neval, real *resabs, real *resasc, integer *
	momcom, real *chebmo)
{
    /* Initialized data */

    static real x[11] = { .9914448613738104f,.9659258262890683f,
	    .9238795325112868f,.8660254037844386f,.7933533402912352f,
	    .7071067811865475f,.6087614290087206f,.5f,.3826834323650898f,
	    .2588190451025208f,.1305261922200516f };

    /* System generated locals */
    integer chebmo_dim1, chebmo_offset, i__1;
    real r__1, r__2;

    /* Local variables */
    static real d__[25];
    static integer i__, j, k, m;
    static real v[28], d1[25], d2[25], p2, p3, p4, ac, an, as, an2, ass, par2,
	     conc, asap, par22, fval[25], estc, cons;
    static integer iers;
    extern /* Subroutine */ int qk15w_(E_fp, E_fp, real *, real *, real *, 
	    real *, integer *, real *, real *, real *, real *, real *, real *)
	    ;
    static real ests;
    static integer isym, noeq1;
    static real cheb12[13], cheb24[25];
    extern /* Subroutine */ int qcheb_(real *, real *, real *, real *);
    static real resc12, resc24, hlgth, centr, ress12, ress24, oflow;
    static integer noequ;
    extern doublereal qwgtf_();
    extern /* Subroutine */ int sgtsl_(integer *, real *, real *, real *, 
	    real *, integer *);
    extern doublereal r1mach_(integer *);
    static real cospar, sinpar, parint;

/* ***BEGIN PROLOGUE  QC25F */
/* ***PURPOSE  To compute the integral I=Integral of F(X) over (A,B) */
/*            Where W(X) = COS(OMEGA*X) Or (WX)=SIN(OMEGA*X) */
/*            and to compute J=Integral of ABS(F) over (A,B). For small */
/*            value of OMEGA or small intervals (A,B) 15-point GAUSS- */
/*            KRONROD Rule used. Otherwise generalized CLENSHAW-CURTIS us */
/* ***LIBRARY   SLATEC (QUADPACK) */
/* ***CATEGORY  H2A2A2 */
/* ***TYPE      SINGLE PRECISION (QC25F-S, DQC25F-D) */
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
/*        Real version */

/*        PARAMETERS */
/*         ON ENTRY */
/*           F      - Real */
/*                    Function subprogram defining the integrand */
/*                    function F(X). The actual name for F needs to */
/*                    be declared E X T E R N A L in the calling program. */

/*           A      - Real */
/*                    Lower limit of integration */

/*           B      - Real */
/*                    Upper limit of integration */

/*           OMEGA  - Real */
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
/*           RESULT - Real */
/*                    Approximation to the integral I */

/*           ABSERR - Real */
/*                    Estimate of the modulus of the absolute */
/*                    error, which should equal or exceed ABS(I-RESULT) */

/*           NEVAL  - Integer */
/*                    Number of integrand evaluations */

/*           RESABS - Real */
/*                    Approximation to the integral J */

/*           RESASC - Real */
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

/*           CHEBMO - Real */
/*                    Array of dimension at least (MAXP1,25) containing */
/*                    the modified Chebyshev moments for the first MOMCOM */
/*                    MOMCOM interval lengths */

/* ***REFERENCES  (NONE) */
/* ***ROUTINES CALLED  QCHEB, QK15W, QWGTF, R1MACH, SGTSL */
/* ***REVISION HISTORY  (YYMMDD) */
/*   810101  DATE WRITTEN */
/*   861211  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/* ***END PROLOGUE  QC25F */




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
/*                    (B-A)*0.5*COS(K*PI/12) + (B+A)*0.5, */
/*                    K = 0, ..., 24 */
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

/* ***FIRST EXECUTABLE STATEMENT  QC25F */
    oflow = r1mach_(&c__2);

    centr = (*b + *a) * .5f;
    hlgth = (*b - *a) * .5f;
    parint = *omega * hlgth;

/*           COMPUTE THE INTEGRAL USING THE 15-POINT GAUSS-KRONROD */
/*           FORMULA IF THE VALUE OF THE PARAMETER IN THE INTEGRAND */
/*           IS SMALL. */

    if (dabs(parint) > 2.f) {
	goto L10;
    }
    qk15w_((E_fp)f, (E_fp)qwgtf_, omega, &p2, &p3, &p4, integr, a, b, result, 
	    abserr, resabs, resasc);
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
    par22 = par2 + 2.f;
    sinpar = sin(parint);
    cospar = cos(parint);

/*           COMPUTE THE CHEBYSHEV MOMENTS WITH RESPECT TO COSINE. */

    v[0] = sinpar * 2.f / parint;
    v[1] = (cospar * 8.f + (par2 + par2 - 8.f) * sinpar / parint) / par2;
    v[2] = ((par2 - 12.f) * 32.f * cospar + ((par2 - 80.f) * par2 + 192.f) * 
	    2.f * sinpar / parint) / (par2 * par2);
    ac = cospar * 8.f;
    as = parint * 24.f * sinpar;
    if (dabs(parint) > 24.f) {
	goto L30;
    }

/*           COMPUTE THE CHEBYSHEV MOMENTS AS THE */
/*           SOLUTIONS OF A BOUNDARY VALUE PROBLEM WITH 1 */
/*           INITIAL VALUE (V(3)) AND 1 END VALUE (COMPUTED */
/*           USING AN ASYMPTOTIC FORMULA). */

    noequ = 25;
    noeq1 = noequ - 1;
    an = 6.f;
    i__1 = noeq1;
    for (k = 1; k <= i__1; ++k) {
	an2 = an * an;
	d__[k - 1] = (an2 - 4.f) * -2.f * (par22 - an2 - an2);
	d2[k - 1] = (an - 1.f) * (an - 2.f) * par2;
	d1[k] = (an + 3.f) * (an + 4.f) * par2;
	v[k + 2] = as - (an2 - 4.f) * ac;
	an += 2.f;
/* L20: */
    }
    an2 = an * an;
    d__[noequ - 1] = (an2 - 4.f) * -2.f * (par22 - an2 - an2);
    v[noequ + 2] = as - (an2 - 4.f) * ac;
    v[3] -= par2 * 56.f * v[2];
    ass = parint * sinpar;
    asap = (((((par2 * 210.f - 1.f) * cospar - (par2 * 105.f - 63.f) * ass) / 
	    an2 - (1.f - par2 * 15.f) * cospar + ass * 15.f) / an2 - cospar + 
	    ass * 3.f) / an2 - cospar) / an2;
    v[noequ + 2] -= asap * 2.f * par2 * (an - 1.f) * (an - 2.f);

/*           SOLVE THE TRIDIAGONAL SYSTEM BY MEANS OF GAUSSIAN */
/*           ELIMINATION WITH PARTIAL PIVOTING. */

    sgtsl_(&noequ, d1, d__, d2, &v[3], &iers);
    goto L50;

/*           COMPUTE THE CHEBYSHEV MOMENTS BY MEANS OF FORWARD */
/*           RECURSION. */

L30:
    an = 4.f;
    for (i__ = 4; i__ <= 13; ++i__) {
	an2 = an * an;
	v[i__ - 1] = ((an2 - 4.f) * ((par22 - an2 - an2) * 2.f * v[i__ - 2] - 
		ac) + as - par2 * (an + 1.f) * (an + 2.f) * v[i__ - 3]) / (
		par2 * (an - 1.f) * (an - 2.f));
	an += 2.f;
/* L40: */
    }
L50:
    for (j = 1; j <= 13; ++j) {
	chebmo[m + ((j << 1) - 1) * chebmo_dim1] = v[j - 1];
/* L60: */
    }

/*           COMPUTE THE CHEBYSHEV MOMENTS WITH RESPECT TO SINE. */

    v[0] = (sinpar - parint * cospar) * 2.f / par2;
    v[1] = (18.f - 48.f / par2) * sinpar / par2 + (48.f / par2 - 2.f) * 
	    cospar / parint;
    ac = parint * -24.f * cospar;
    as = sinpar * -8.f;
    if (dabs(parint) > 24.f) {
	goto L80;
    }

/*           COMPUTE THE CHEBYSHEV MOMENTS AS THE */
/*           SOLUTIONS OF A BOUNDARY VALUE PROBLEM WITH 1 */
/*           INITIAL VALUE (V(2)) AND 1 END VALUE (COMPUTED */
/*           USING AN ASYMPTOTIC FORMULA). */

    an = 5.f;
    i__1 = noeq1;
    for (k = 1; k <= i__1; ++k) {
	an2 = an * an;
	d__[k - 1] = (an2 - 4.f) * -2.f * (par22 - an2 - an2);
	d2[k - 1] = (an - 1.f) * (an - 2.f) * par2;
	d1[k] = (an + 3.f) * (an + 4.f) * par2;
	v[k + 1] = ac + (an2 - 4.f) * as;
	an += 2.f;
/* L70: */
    }
    an2 = an * an;
    d__[noequ - 1] = (an2 - 4.f) * -2.f * (par22 - an2 - an2);
    v[noequ + 1] = ac + (an2 - 4.f) * as;
    v[2] -= par2 * 42.f * v[1];
    ass = parint * cospar;
    asap = (((((par2 * 105.f - 63.f) * ass + (par2 * 210.f - 1.f) * sinpar) / 
	    an2 + (par2 * 15.f - 1.f) * sinpar - ass * 15.f) / an2 - ass * 
	    3.f - sinpar) / an2 - sinpar) / an2;
    v[noequ + 1] -= asap * 2.f * par2 * (an - 1.f) * (an - 2.f);

/*           SOLVE THE TRIDIAGONAL SYSTEM BY MEANS OF GAUSSIAN */
/*           ELIMINATION WITH PARTIAL PIVOTING. */

    sgtsl_(&noequ, d1, d__, d2, &v[2], &iers);
    goto L100;

/*           COMPUTE THE CHEBYSHEV MOMENTS BY MEANS OF */
/*           FORWARD RECURSION. */

L80:
    an = 3.f;
    for (i__ = 3; i__ <= 12; ++i__) {
	an2 = an * an;
	v[i__ - 1] = ((an2 - 4.f) * ((par22 - an2 - an2) * 2.f * v[i__ - 2] + 
		as) + ac - par2 * (an + 1.f) * (an + 2.f) * v[i__ - 3]) / (
		par2 * (an - 1.f) * (an - 2.f));
	an += 2.f;
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

    r__1 = centr + hlgth;
    fval[0] = (*f)(&r__1) * .5f;
    fval[12] = (*f)(&centr);
    r__1 = centr - hlgth;
    fval[24] = (*f)(&r__1) * .5f;
    for (i__ = 2; i__ <= 12; ++i__) {
	isym = 26 - i__;
	r__1 = hlgth * x[i__ - 2] + centr;
	fval[i__ - 1] = (*f)(&r__1);
	r__1 = centr - hlgth * x[i__ - 2];
	fval[isym - 1] = (*f)(&r__1);
/* L130: */
    }
    qcheb_(x, fval, cheb12, cheb24);

/*           COMPUTE THE INTEGRAL AND ERROR ESTIMATES. */

    resc12 = cheb12[12] * chebmo[m + chebmo_dim1 * 13];
    ress12 = 0.f;
    k = 11;
    for (j = 1; j <= 6; ++j) {
	resc12 += cheb12[k - 1] * chebmo[m + k * chebmo_dim1];
	ress12 += cheb12[k] * chebmo[m + (k + 1) * chebmo_dim1];
	k += -2;
/* L140: */
    }
    resc24 = cheb24[24] * chebmo[m + chebmo_dim1 * 25];
    ress24 = 0.f;
    *resabs = dabs(cheb24[24]);
    k = 23;
    for (j = 1; j <= 12; ++j) {
	resc24 += cheb24[k - 1] * chebmo[m + k * chebmo_dim1];
	ress24 += cheb24[k] * chebmo[m + (k + 1) * chebmo_dim1];
	*resabs = (r__1 = cheb24[k - 1], dabs(r__1)) + (r__2 = cheb24[k], 
		dabs(r__2));
	k += -2;
/* L150: */
    }
    estc = (r__1 = resc24 - resc12, dabs(r__1));
    ests = (r__1 = ress24 - ress12, dabs(r__1));
    *resabs *= dabs(hlgth);
    if (*integr == 2) {
	goto L160;
    }
    *result = conc * resc24 - cons * ress24;
    *abserr = (r__1 = conc * estc, dabs(r__1)) + (r__2 = cons * ests, dabs(
	    r__2));
    goto L170;
L160:
    *result = conc * ress24 + cons * resc24;
    *abserr = (r__1 = conc * ests, dabs(r__1)) + (r__2 = cons * estc, dabs(
	    r__2));
L170:
    return 0;
} /* qc25f_ */

