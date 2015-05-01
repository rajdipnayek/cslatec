/* gamrn.f -- translated by f2c (version 12.02.01).
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
static integer c__11 = 11;
static integer c__5 = 5;

/* DECK GAMRN */
doublereal gamrn_(real *x)
{
    /* Initialized data */

    static real gr[12] = { 1.f,-.015625f,.0025634765625f,
	    -.0012798309326171875f,.00134351104497909546f,
	    -.00243289663922041655f,.00675423753364157164f,
	    -.0266369606131178216f,.141527455519564332f,-.974384543032201613f,
	    8.43686251229783675f,-89.7258321640552515f };

    /* System generated locals */
    integer i__1;
    real ret_val, r__1;

    /* Local variables */
    static integer i__, k;
    static real s;
    static integer mx, nx;
    static real xm, xp, fln, rln, tol, trm, xsq;
    static integer i1m11;
    static real xinc, xmin, xdmy;
    extern integer i1mach_(integer *);
    extern doublereal r1mach_(integer *);

/* ***BEGIN PROLOGUE  GAMRN */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to BSKIN */
/* ***LIBRARY   SLATEC */
/* ***TYPE      SINGLE PRECISION (GAMRN-S, DGAMRN-D) */
/* ***AUTHOR  Amos, D. E., (SNLA) */
/* ***DESCRIPTION */

/*     Abstract */
/*         GAMRN computes the GAMMA function ratio GAMMA(X)/GAMMA(X+0.5) */
/*         for real X.gt.0. If X.ge.XMIN, an asymptotic expansion is */
/*         evaluated. If X.lt.XMIN, an integer is added to X to form a */
/*         new value of X.ge.XMIN and the asymptotic expansion is eval- */
/*         uated for this new value of X. Successive application of the */
/*         recurrence relation */

/*                      W(X)=W(X+1)*(1+0.5/X) */

/*         reduces the argument to its original value. XMIN and comp- */
/*         utational tolerances are computed as a function of the number */
/*         of digits carried in a word by calls to I1MACH and R1MACH. */
/*         However, the computational accuracy is limited to the max- */
/*         imum of unit roundoff (=R1MACH(4)) and 1.0E-18 since critical */
/*         constants are given to only 18 digits. */

/*         Input */
/*           X      - Argument, X.gt.0.0 */

/*         OUTPUT */
/*           GAMRN  - Ratio  GAMMA(X)/GAMMA(X+0.5) */

/* ***SEE ALSO  BSKIN */
/* ***REFERENCES  Y. L. Luke, The Special Functions and Their */
/*                 Approximations, Vol. 1, Math In Sci. And */
/*                 Eng. Series 53, Academic Press, New York, 1969, */
/*                 pp. 34-35. */
/* ***ROUTINES CALLED  I1MACH, R1MACH */
/* ***REVISION HISTORY  (YYMMDD) */
/*   820601  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900328  Added TYPE section.  (WRB) */
/*   910722  Updated AUTHOR section.  (ALS) */
/*   920520  Added REFERENCES section.  (WRB) */
/* ***END PROLOGUE  GAMRN */


/* ***FIRST EXECUTABLE STATEMENT  GAMRN */
    nx = (integer) (*x);
/* Computing MAX */
    r__1 = r1mach_(&c__4);
    tol = dmax(r__1,1e-18f);
    i1m11 = i1mach_(&c__11);
    rln = r1mach_(&c__5) * i1m11;
    fln = dmin(rln,20.f);
    fln = dmax(fln,3.f);
    fln += -3.f;
    xm = fln * (fln * .01723f + .2366f) + 2.f;
    mx = (integer) xm + 1;
    xmin = (real) mx;
    xdmy = *x - .25f;
    xinc = 0.f;
    if (*x >= xmin) {
	goto L10;
    }
    xinc = xmin - nx;
    xdmy += xinc;
L10:
    s = 1.f;
    if (xdmy * tol > 1.f) {
	goto L30;
    }
    xsq = 1.f / (xdmy * xdmy);
    xp = xsq;
    for (k = 2; k <= 12; ++k) {
	trm = gr[k - 1] * xp;
	if (dabs(trm) < tol) {
	    goto L30;
	}
	s += trm;
	xp *= xsq;
/* L20: */
    }
L30:
    s /= sqrt(xdmy);
    if (xinc != 0.f) {
	goto L40;
    }
    ret_val = s;
    return ret_val;
L40:
    nx = (integer) xinc;
    xp = 0.f;
    i__1 = nx;
    for (i__ = 1; i__ <= i__1; ++i__) {
	s *= .5f / (*x + xp) + 1.f;
	xp += 1.f;
/* L50: */
    }
    ret_val = s;
    return ret_val;
} /* gamrn_ */

