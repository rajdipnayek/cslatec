/* d9atn1.f -- translated by f2c (version 12.02.01).
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
static integer c__40 = 40;
static integer c__2 = 2;
static integer c__1 = 1;

/* DECK D9ATN1 */
doublereal d9atn1_(doublereal *x)
{
    /* Initialized data */

    static doublereal atn1cs[40] = { -.0328399753535520235690793992299,
	    .05833432343172412449951669914907,
	    -.007400369696719646463809011551413,
	    .001009784199337288083590357511639,
	    -1.4397871635652056214713036977e-4,
	    2.114512648992107572072112243439e-5,
	    -3.172321074254667167402564996757e-6,
	    4.8366203654607108253778593848e-7,
	    -7.467746546814112670437614322776e-8,
	    1.164800896824429830620998641342e-8,
	    -1.832088370847201392699956242452e-9,
	    2.901908277966063313175351230455e-10,
	    -4.623885312106326738351805721512e-11,
	    7.405528668775736917992197048286e-12,
	    -1.191354457845136682370820373417e-12,
	    1.924090144391772599867855692518e-13,
	    -3.118271051076194272254476155327e-14,
	    5.069240036567731789694520593032e-15,
	    -8.263694719802866053818284405964e-16,
	    1.350486709817079420526506123029e-16,
	    -2.212023650481746045840137823191e-17,
	    3.630654747381356783829047647709e-18,
	    -5.970345328847154052451215859165e-19,
	    9.834816050077133119448329005738e-20,
	    -1.62265507585506233614438760448e-20,
	    2.681186176945436796301320301226e-21,
	    -4.436309706785255479636243688106e-22,
	    7.3496918976524969450724655104e-23,
	    -1.219077508350052588289401378133e-23,
	    2.024298836805215403184540876799e-24,
	    -3.364871555797354579925576362666e-25,
	    5.598673968346988749492933973333e-26,
	    -9.323939267272320229628532053333e-27,
	    1.554133116995970222934807893333e-27,
	    -2.592569534179745922757427199999e-28,
	    4.328193466245734685037909333333e-29,
	    -7.231013125595437471192405333333e-30,
	    1.208902859830494772942165333333e-30,
	    -2.022404543449897579315199999999e-31,
	    3.385428713046493843073706666666e-32 };
    static logical first = TRUE_;

    /* System generated locals */
    real r__1;
    doublereal ret_val, d__1;

    /* Local variables */
    static doublereal y, eps, xbig, xmax, xsml;
    extern doublereal d1mach_(integer *);
    static integer ntatn1;
    extern doublereal dcsevl_(doublereal *, doublereal *, integer *);
    extern integer initds_(doublereal *, integer *, real *);
    extern /* Subroutine */ int xermsg_(char *, char *, char *, integer *, 
	    integer *, ftnlen, ftnlen, ftnlen);

/* ***BEGIN PROLOGUE  D9ATN1 */
/* ***SUBSIDIARY */
/* ***PURPOSE  Evaluate DATAN(X) from first order relative accuracy so */
/*            that DATAN(X) = X + X**3*D9ATN1(X). */
/* ***LIBRARY   SLATEC (FNLIB) */
/* ***CATEGORY  C4A */
/* ***TYPE      DOUBLE PRECISION (R9ATN1-S, D9ATN1-D) */
/* ***KEYWORDS  ARC TANGENT, ELEMENTARY FUNCTIONS, FIRST ORDER, FNLIB, */
/*             TRIGONOMETRIC */
/* ***AUTHOR  Fullerton, W., (LANL) */
/* ***DESCRIPTION */

/* Evaluate  DATAN(X)  from first order, that is, evaluate */
/* (DATAN(X)-X)/X**3  with relative error accuracy so that */
/*        DATAN(X) = X + X**3*D9ATN1(X). */

/* Series for ATN1       on the interval  0.          to  1.00000E+00 */
/*                                        with weighted error   3.39E-32 */
/*                                         log weighted error  31.47 */
/*                               significant figures required  30.26 */
/*                                    decimal places required  32.27 */

/* ***REFERENCES  (NONE) */
/* ***ROUTINES CALLED  D1MACH, DCSEVL, INITDS, XERMSG */
/* ***REVISION HISTORY  (YYMMDD) */
/*   780401  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   891115  Corrected third argument in reference to INITDS.  (WRB) */
/*   891115  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ) */
/*   900720  Routine changed from user-callable to subsidiary.  (WRB) */
/* ***END PROLOGUE  D9ATN1 */
/* ***FIRST EXECUTABLE STATEMENT  D9ATN1 */
    if (first) {
	eps = d1mach_(&c__3);
	r__1 = (real) eps * .1f;
	ntatn1 = initds_(atn1cs, &c__40, &r__1);

	xsml = sqrt(eps * .1);
	xbig = 1.571 / sqrt(eps);
	xmax = 1.571 / eps;
    }
    first = FALSE_;

    y = abs(*x);
    if (y > 1.) {
	goto L20;
    }

    if (y <= xsml) {
	ret_val = -.33333333333333331;
    }
    if (y <= xsml) {
	return ret_val;
    }

    d__1 = y * 2. * y - 1.;
    ret_val = dcsevl_(&d__1, atn1cs, &ntatn1) - .25;
    return ret_val;

L20:
    if (y > xmax) {
	xermsg_("SLATEC", "D9ATN1", "NO PRECISION IN ANSWER BECAUSE X IS TOO"
		" BIG", &c__2, &c__2, (ftnlen)6, (ftnlen)6, (ftnlen)43);
    }
    if (y > xbig) {
	xermsg_("SLATEC", "D9ATN1", "ANSWER LT HALF PRECISION BECAUSE X IS T"
		"OO BIG", &c__1, &c__1, (ftnlen)6, (ftnlen)6, (ftnlen)45);
    }

/* Computing 3rd power */
    d__1 = *x;
    ret_val = (atan(*x) - *x) / (d__1 * (d__1 * d__1));
    return ret_val;

} /* d9atn1_ */

