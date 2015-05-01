/* dbi.f -- translated by f2c (version 12.02.01).
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
static integer c__15 = 15;
static doublereal c_b7 = .3333;
static integer c__2 = 2;
static doublereal c_b9 = .6666;
static integer c__1 = 1;

/* DECK DBI */
doublereal dbi_(doublereal *x)
{
    /* Initialized data */

    static doublereal bifcs[13] = { -.016730216471986649483537423928176,
	    .10252335834249445611426362777757,
	    .0017083092507381516539429650242013,
	    1.186254546774468117921645921004e-5,
	    4.4932907017792133694531887927242e-8,
	    1.0698207143387889067567767663628e-10,
	    1.7480643399771824706010517628573e-13,
	    2.0810231071761711025881891834399e-16,
	    1.8849814695665416509927971733333e-19,
	    1.3425779173097804625882666666666e-22,
	    7.7159593429658887893333333333333e-26,
	    3.6533879617478566399999999999999e-29,
	    1.4497565927953066666666666666666e-32 };
    static doublereal bigcs[13] = { .022466223248574522283468220139024,
	    .037364775453019545441727561666752,
	    4.4476218957212285696215294326639e-4,
	    2.4708075636329384245494591948882e-6,
	    7.9191353395149635134862426285596e-9,
	    1.6498079851827779880887872402706e-11,
	    2.4119906664835455909247501122841e-14,
	    2.6103736236091436985184781269333e-17,
	    2.1753082977160323853123792e-20,
	    1.4386946400390433219483733333333e-23,
	    7.7349125612083468629333333333333e-27,
	    3.4469292033849002666666666666666e-30,1.2938919273216e-33 };
    static doublereal bif2cs[15] = { .0998457269381604104468284257993,
	    .47862497786300553772211467318231,
	    .025155211960433011771324415436675,
	    5.8206938852326456396515697872216e-4,
	    7.4997659644377865943861457378217e-6,
	    6.1346028703493836681403010356474e-8,
	    3.4627538851480632900434268733359e-10,
	    1.4288910080270254287770846748931e-12,
	    4.49627042983346418950564721792e-15,
	    1.1142323065833011708428300106666e-17,
	    2.2304791066175002081517866666666e-20,
	    3.6815778736393142842922666666666e-23,
	    5.0960868449338261333333333333333e-26,
	    6.0003386926288554666666666666666e-29,
	    6.0827497446570666666666666666666e-32 };
    static doublereal big2cs[15] = { .033305662145514340465176188111647,
	    .161309215123197067613287532084943,
	    .00631900730961342869121615634921173,
	    1.18790456816251736389780192304567e-4,
	    1.30453458862002656147116485012843e-6,
	    9.37412599553521729546809615508936e-9,
	    4.74580188674725153788510169834595e-11,
	    1.78310726509481399800065667560946e-13,
	    5.1675919278495818037427635664e-16,
	    1.19004508386827125129496251733333e-18,
	    2.22982880666403517277063466666666e-21,
	    3.46551923027689419722666666666666e-24,
	    4.53926336320504514133333333333333e-27,
	    5.07884996513522346666666666666666e-30,
	    4.91020674696533333333333333333333e-33 };
    static logical first = TRUE_;

    /* System generated locals */
    doublereal ret_val, d__1;

    /* Local variables */
    static doublereal z__, xm;
    static real eta;
    extern doublereal dbie_(doublereal *);
    static integer nbif, nbig;
    static doublereal xmax;
    static integer nbif2, nbig2;
    static doublereal x3sml, theta;
    extern doublereal d1mach_(integer *);
    extern /* Subroutine */ int d9aimp_(doublereal *, doublereal *, 
	    doublereal *);
    extern doublereal dcsevl_(doublereal *, doublereal *, integer *);
    extern integer initds_(doublereal *, integer *, real *);
    extern /* Subroutine */ int xermsg_(char *, char *, char *, integer *, 
	    integer *, ftnlen, ftnlen, ftnlen);

/* ***BEGIN PROLOGUE  DBI */
/* ***PURPOSE  Evaluate the Bairy function (the Airy function of the */
/*            second kind). */
/* ***LIBRARY   SLATEC (FNLIB) */
/* ***CATEGORY  C10D */
/* ***TYPE      DOUBLE PRECISION (BI-S, DBI-D) */
/* ***KEYWORDS  BAIRY FUNCTION, FNLIB, SPECIAL FUNCTIONS */
/* ***AUTHOR  Fullerton, W., (LANL) */
/* ***DESCRIPTION */

/* DBI(X) calculates the double precision Airy function of the */
/* second kind for double precision argument X. */

/* Series for BIF        on the interval -1.00000E+00 to  1.00000E+00 */
/*                                        with weighted error   1.45E-32 */
/*                                         log weighted error  31.84 */
/*                               significant figures required  30.85 */
/*                                    decimal places required  32.40 */

/* Series for BIG        on the interval -1.00000E+00 to  1.00000E+00 */
/*                                        with weighted error   1.29E-33 */
/*                                         log weighted error  32.89 */
/*                               significant figures required  31.48 */
/*                                    decimal places required  33.45 */

/* Series for BIF2       on the interval  1.00000E+00 to  8.00000E+00 */
/*                                        with weighted error   6.08E-32 */
/*                                         log weighted error  31.22 */
/*                        approx significant figures required  30.8 */
/*                                    decimal places required  31.80 */

/* Series for BIG2       on the interval  1.00000E+00 to  8.00000E+00 */
/*                                        with weighted error   4.91E-33 */
/*                                         log weighted error  32.31 */
/*                        approx significant figures required  31.6 */
/*                                    decimal places required  32.90 */

/* ***REFERENCES  (NONE) */
/* ***ROUTINES CALLED  D1MACH, D9AIMP, DBIE, DCSEVL, INITDS, XERMSG */
/* ***REVISION HISTORY  (YYMMDD) */
/*   770701  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890531  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ) */
/* ***END PROLOGUE  DBI */
/* ***FIRST EXECUTABLE STATEMENT  DBI */
    if (first) {
	eta = (real) d1mach_(&c__3) * .1f;
	nbif = initds_(bifcs, &c__13, &eta);
	nbig = initds_(bigcs, &c__13, &eta);
	nbif2 = initds_(bif2cs, &c__15, &eta);
	nbig2 = initds_(big2cs, &c__15, &eta);

	d__1 = (doublereal) eta;
	x3sml = pow_dd(&d__1, &c_b7);
	d__1 = log(d1mach_(&c__2)) * 1.5f;
	xmax = pow_dd(&d__1, &c_b9);
    }
    first = FALSE_;

    if (*x >= -1.) {
	goto L20;
    }
    d9aimp_(x, &xm, &theta);
    ret_val = xm * sin(theta);
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
    ret_val = dcsevl_(&z__, bifcs, &nbif) + .625f + *x * (dcsevl_(&z__, bigcs,
	     &nbig) + .4375);
    return ret_val;

L30:
    if (*x > 2.) {
	goto L40;
    }
/* Computing 3rd power */
    d__1 = *x;
    z__ = (d__1 * (d__1 * d__1) * 2. - 9.) / 7.;
    ret_val = dcsevl_(&z__, bif2cs, &nbif2) + 1.125 + *x * (dcsevl_(&z__, 
	    big2cs, &nbig2) + .625);
    return ret_val;

L40:
    if (*x > xmax) {
	xermsg_("SLATEC", "DBI", "X SO BIG THAT BI OVERFLOWS", &c__1, &c__2, (
		ftnlen)6, (ftnlen)3, (ftnlen)26);
    }

    ret_val = dbie_(x) * exp(*x * 2. * sqrt(*x) / 3.);
    return ret_val;

} /* dbi_ */

