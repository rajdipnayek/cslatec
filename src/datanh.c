/* datanh.f -- translated by f2c (version 12.02.01).
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
static integer c__27 = 27;
static integer c__4 = 4;
static integer c__2 = 2;
static integer c__1 = 1;

/* DECK DATANH */
doublereal datanh_(doublereal *x)
{
    /* Initialized data */

    static doublereal atnhcs[27] = { .09439510239319549230842892218633,
	    .04919843705578615947200034576668,
	    .002102593522455432763479327331752,
	    1.073554449776116584640731045276e-4,
	    5.978267249293031478642787517872e-6,
	    3.5050620308891348459668348862e-7,
	    2.126374343765340350896219314431e-8,
	    1.321694535715527192129801723055e-9,
	    8.365875501178070364623604052959e-11,
	    5.370503749311002163881434587772e-12,
	    3.48665947015710792297124578429e-13,
	    2.284549509603433015524024119722e-14,
	    1.508407105944793044874229067558e-15,
	    1.002418816804109126136995722837e-16,
	    6.698674738165069539715526882986e-18,
	    4.497954546494931083083327624533e-19,
	    3.032954474279453541682367146666e-20,
	    2.052702064190936826463861418666e-21,
	    1.393848977053837713193014613333e-22,
	    9.492580637224576971958954666666e-24,
	    6.481915448242307604982442666666e-25,
	    4.43673020572361527263232e-26,3.043465618543161638912e-27,
	    2.091881298792393474047999999999e-28,
	    1.440445411234050561365333333333e-29,
	    9.935374683141640465066666666666e-31,
	    6.863462444358260053333333333333e-32 };
    static logical first = TRUE_;

    /* System generated locals */
    real r__1;
    doublereal ret_val, d__1;

    /* Local variables */
    static doublereal y, dxrel, sqeps;
    extern doublereal d1mach_(integer *), dcsevl_(doublereal *, doublereal *, 
	    integer *);
    extern integer initds_(doublereal *, integer *, real *);
    extern /* Subroutine */ int xermsg_(char *, char *, char *, integer *, 
	    integer *, ftnlen, ftnlen, ftnlen);
    static integer nterms;

/* ***BEGIN PROLOGUE  DATANH */
/* ***PURPOSE  Compute the arc hyperbolic tangent. */
/* ***LIBRARY   SLATEC (FNLIB) */
/* ***CATEGORY  C4C */
/* ***TYPE      DOUBLE PRECISION (ATANH-S, DATANH-D, CATANH-C) */
/* ***KEYWORDS  ARC HYPERBOLIC TANGENT, ATANH, ELEMENTARY FUNCTIONS, */
/*             FNLIB, INVERSE HYPERBOLIC TANGENT */
/* ***AUTHOR  Fullerton, W., (LANL) */
/* ***DESCRIPTION */

/* DATANH(X) calculates the double precision arc hyperbolic */
/* tangent for double precision argument X. */

/* Series for ATNH       on the interval  0.          to  2.50000E-01 */
/*                                        with weighted error   6.86E-32 */
/*                                         log weighted error  31.16 */
/*                               significant figures required  30.00 */
/*                                    decimal places required  31.88 */

/* ***REFERENCES  (NONE) */
/* ***ROUTINES CALLED  D1MACH, DCSEVL, INITDS, XERMSG */
/* ***REVISION HISTORY  (YYMMDD) */
/*   770601  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890531  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ) */
/* ***END PROLOGUE  DATANH */
/* ***FIRST EXECUTABLE STATEMENT  DATANH */
    if (first) {
	r__1 = (real) d1mach_(&c__3) * .1f;
	nterms = initds_(atnhcs, &c__27, &r__1);
	dxrel = sqrt(d1mach_(&c__4));
	sqeps = sqrt(d1mach_(&c__3) * 3.);
    }
    first = FALSE_;

    y = abs(*x);
    if (y >= 1.) {
	xermsg_("SLATEC", "DATANH", "ABS(X) GE 1", &c__2, &c__2, (ftnlen)6, (
		ftnlen)6, (ftnlen)11);
    }

    if (1. - y < dxrel) {
	xermsg_("SLATEC", "DATANH", "ANSWER LT HALF PRECISION BECAUSE ABS(X)"
		" TOO NEAR 1", &c__1, &c__1, (ftnlen)6, (ftnlen)6, (ftnlen)50);
    }

    ret_val = *x;
    if (y > sqeps && y <= .5) {
	d__1 = *x * 8. * *x - 1.;
	ret_val = *x * (dcsevl_(&d__1, atnhcs, &nterms) + 1.);
    }
    if (y > .5) {
	ret_val = log((*x + 1.) / (1. - *x)) * .5;
    }

    return ret_val;
} /* datanh_ */

