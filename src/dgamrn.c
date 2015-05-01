/* dgamrn.f -- translated by f2c (version 12.02.01).
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
static integer c__14 = 14;
static integer c__5 = 5;

/* DECK DGAMRN */
doublereal dgamrn_(doublereal *x)
{
    /* Initialized data */

    static doublereal gr[12] = { 1.,-.015625,.0025634765625,
	    -.0012798309326171875,.00134351104497909546,
	    -.00243289663922041655,.00675423753364157164,
	    -.0266369606131178216,.141527455519564332,-.974384543032201613,
	    8.43686251229783675,-89.7258321640552515 };

    /* System generated locals */
    integer i__1;
    doublereal ret_val, d__1;

    /* Local variables */
    static integer i__, k;
    static doublereal s;
    static integer mx, nx;
    static doublereal xm, xp, fln, rln, tol, trm, xsq;
    static integer i1m11;
    static doublereal xinc, xmin, xdmy;
    extern doublereal d1mach_(integer *);
    extern integer i1mach_(integer *);

/* ***BEGIN PROLOGUE  DGAMRN */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to DBSKIN */
/* ***LIBRARY   SLATEC */
/* ***TYPE      DOUBLE PRECISION (GAMRN-S, DGAMRN-D) */
/* ***AUTHOR  Amos, D. E., (SNLA) */
/* ***DESCRIPTION */

/*     Abstract   * A Double Precision Routine * */
/*         DGAMRN computes the GAMMA function ratio GAMMA(X)/GAMMA(X+0.5) */
/*         for real X.gt.0. If X.ge.XMIN, an asymptotic expansion is */
/*         evaluated. If X.lt.XMIN, an integer is added to X to form a */
/*         new value of X.ge.XMIN and the asymptotic expansion is eval- */
/*         uated for this new value of X. Successive application of the */
/*         recurrence relation */

/*                      W(X)=W(X+1)*(1+0.5/X) */

/*         reduces the argument to its original value. XMIN and comp- */
/*         utational tolerances are computed as a function of the number */
/*         of digits carried in a word by calls to I1MACH and D1MACH. */
/*         However, the computational accuracy is limited to the max- */
/*         imum of unit roundoff (=D1MACH(4)) and 1.0D-18 since critical */
/*         constants are given to only 18 digits. */

/*         Input      X is Double Precision */
/*           X      - Argument, X.gt.0.0D0 */

/*         Output      DGAMRN is DOUBLE PRECISION */
/*           DGAMRN  - Ratio  GAMMA(X)/GAMMA(X+0.5) */

/* ***SEE ALSO  DBSKIN */
/* ***REFERENCES  Y. L. Luke, The Special Functions and Their */
/*                 Approximations, Vol. 1, Math In Sci. And */
/*                 Eng. Series 53, Academic Press, New York, 1969, */
/*                 pp. 34-35. */
/* ***ROUTINES CALLED  D1MACH, I1MACH */
/* ***REVISION HISTORY  (YYMMDD) */
/*   820601  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890911  Removed unnecessary intrinsics.  (WRB) */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900328  Added TYPE section.  (WRB) */
/*   910722  Updated AUTHOR section.  (ALS) */
/*   920520  Added REFERENCES section.  (WRB) */
/* ***END PROLOGUE  DGAMRN */


/* ***FIRST EXECUTABLE STATEMENT  DGAMRN */
    nx = (integer) (*x);
/* Computing MAX */
    d__1 = d1mach_(&c__4);
    tol = max(d__1,1e-18);
    i1m11 = i1mach_(&c__14);
    rln = d1mach_(&c__5) * i1m11;
    fln = min(rln,20.);
    fln = max(fln,3.);
    fln += -3.;
    xm = fln * (fln * .01723 + .2366) + 2.;
    mx = (integer) xm + 1;
    xmin = (doublereal) mx;
    xdmy = *x - .25;
    xinc = 0.;
    if (*x >= xmin) {
	goto L10;
    }
    xinc = xmin - nx;
    xdmy += xinc;
L10:
    s = 1.;
    if (xdmy * tol > 1.) {
	goto L30;
    }
    xsq = 1. / (xdmy * xdmy);
    xp = xsq;
    for (k = 2; k <= 12; ++k) {
	trm = gr[k - 1] * xp;
	if (abs(trm) < tol) {
	    goto L30;
	}
	s += trm;
	xp *= xsq;
/* L20: */
    }
L30:
    s /= sqrt(xdmy);
    if (xinc != 0.) {
	goto L40;
    }
    ret_val = s;
    return ret_val;
L40:
    nx = (integer) xinc;
    xp = 0.;
    i__1 = nx;
    for (i__ = 1; i__ <= i__1; ++i__) {
	s *= .5 / (*x + xp) + 1.;
	xp += 1.;
/* L50: */
    }
    ret_val = s;
    return ret_val;
} /* dgamrn_ */

