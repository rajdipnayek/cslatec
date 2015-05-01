/* dasinh.f -- translated by f2c (version 12.02.01).
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
static integer c__39 = 39;

/* DECK DASINH */
doublereal dasinh_(doublereal *x)
{
    /* Initialized data */

    static doublereal asnhcs[39] = { -.12820039911738186343372127359268,
	    -.058811761189951767565211757138362,
	    .0047274654322124815640725249756029,
	    -4.9383631626536172101360174790273e-4,
	    5.8506207058557412287494835259321e-5,
	    -7.4669983289313681354755069217188e-6,
	    1.0011693583558199265966192015812e-6,
	    -1.3903543858708333608616472258886e-7,
	    1.9823169483172793547317360237148e-8,
	    -2.8847468417848843612747272800317e-9,
	    4.2672965467159937953457514995907e-10,
	    -6.3976084654366357868752632309681e-11,
	    9.6991686089064704147878293131179e-12,
	    -1.4844276972043770830246658365696e-12,
	    2.2903737939027447988040184378983e-13,
	    -3.558839513273264515997894265131e-14,
	    5.5639694080056789953374539088554e-15,
	    -8.7462509599624678045666593520162e-16,
	    1.3815248844526692155868802298129e-16,
	    -2.1916688282900363984955142264149e-17,
	    3.490465852482756563831392370688e-18,
	    -5.5785788400895742439630157032106e-19,
	    8.9445146617134012551050882798933e-20,
	    -1.4383426346571317305551845239466e-20,
	    2.3191811872169963036326144682666e-21,
	    -3.7487007953314343674570604543999e-22,
	    6.073210982206427940454924288e-23,
	    -9.859940276463358317737017344e-24,
	    1.6039217452788496315232638293333e-24,
	    -2.6138847350287686596716134399999e-25,
	    4.2670849606857390833358165333333e-26,
	    -6.9770217039185243299730773333333e-27,
	    1.1425088336806858659812693333333e-27,
	    -1.8735292078860968933021013333333e-28,
	    3.076358441446492279406592e-29,
	    -5.0577364031639824787046399999999e-30,
	    8.3250754712689142224213333333333e-31,
	    -1.3718457282501044163925333333333e-31,
	    2.2629868426552784104106666666666e-32 };
    static doublereal aln2 = .69314718055994530941723212145818;
    static logical first = TRUE_;

    /* System generated locals */
    real r__1;
    doublereal ret_val, d__1;

    /* Local variables */
    static doublereal y, xmax, sqeps;
    extern doublereal d1mach_(integer *), dcsevl_(doublereal *, doublereal *, 
	    integer *);
    extern integer initds_(doublereal *, integer *, real *);
    static integer nterms;

/* ***BEGIN PROLOGUE  DASINH */
/* ***PURPOSE  Compute the arc hyperbolic sine. */
/* ***LIBRARY   SLATEC (FNLIB) */
/* ***CATEGORY  C4C */
/* ***TYPE      DOUBLE PRECISION (ASINH-S, DASINH-D, CASINH-C) */
/* ***KEYWORDS  ARC HYPERBOLIC SINE, ASINH, ELEMENTARY FUNCTIONS, FNLIB, */
/*             INVERSE HYPERBOLIC SINE */
/* ***AUTHOR  Fullerton, W., (LANL) */
/* ***DESCRIPTION */

/* DASINH(X) calculates the double precision arc hyperbolic */
/* sine for double precision argument X. */

/* ***REFERENCES  (NONE) */
/* ***ROUTINES CALLED  D1MACH, DCSEVL, INITDS */
/* ***REVISION HISTORY  (YYMMDD) */
/*   770601  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890531  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/* ***END PROLOGUE  DASINH */
/* ***FIRST EXECUTABLE STATEMENT  DASINH */
    if (first) {
	r__1 = (real) d1mach_(&c__3) * .1f;
	nterms = initds_(asnhcs, &c__39, &r__1);
	sqeps = sqrt(d1mach_(&c__3));
	xmax = 1. / sqeps;
    }
    first = FALSE_;

    y = abs(*x);
    if (y > 1.) {
	goto L20;
    }

    ret_val = *x;
    if (y > sqeps) {
	d__1 = *x * 2. * *x - 1.;
	ret_val = *x * (dcsevl_(&d__1, asnhcs, &nterms) + 1.);
    }
    return ret_val;
L20:
    if (y < xmax) {
	ret_val = log(y + sqrt(y * y + 1.));
    }
    if (y >= xmax) {
	ret_val = aln2 + log(y);
    }
    ret_val = d_sign(&ret_val, x);
    return ret_val;

} /* dasinh_ */

