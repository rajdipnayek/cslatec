/* betai.f -- translated by f2c (version 12.02.01).
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
static integer c__1 = 1;
static integer c__2 = 2;

/* DECK BETAI */
doublereal betai_(real *x, real *pin, real *qin)
{
    /* Initialized data */

    static logical first = TRUE_;

    /* System generated locals */
    integer i__1;
    real ret_val, r__1;

    /* Local variables */
    static real c__;
    static integer i__, n;
    static real p, q, y, p1;
    static integer ib;
    static real xb, ps, eps, sml, term;
    extern doublereal r1mach_(integer *), albeta_(real *, real *);
    static real alneps, alnsml, finsum;
    extern /* Subroutine */ int xermsg_(char *, char *, char *, integer *, 
	    integer *, ftnlen, ftnlen, ftnlen);

/* ***BEGIN PROLOGUE  BETAI */
/* ***PURPOSE  Calculate the incomplete Beta function. */
/* ***LIBRARY   SLATEC (FNLIB) */
/* ***CATEGORY  C7F */
/* ***TYPE      SINGLE PRECISION (BETAI-S, DBETAI-D) */
/* ***KEYWORDS  FNLIB, INCOMPLETE BETA FUNCTION, SPECIAL FUNCTIONS */
/* ***AUTHOR  Fullerton, W., (LANL) */
/* ***DESCRIPTION */

/*   BETAI calculates the REAL incomplete beta function. */

/*   The incomplete beta function ratio is the probability that a */
/*   random variable from a beta distribution having parameters PIN and */
/*   QIN will be less than or equal to X. */

/*     -- Input Arguments -- All arguments are REAL. */
/*   X      upper limit of integration.  X must be in (0,1) inclusive. */
/*   PIN    first beta distribution parameter.  PIN must be .GT. 0.0. */
/*   QIN    second beta distribution parameter.  QIN must be .GT. 0.0. */

/* ***REFERENCES  Nancy E. Bosten and E. L. Battiste, Remark on Algorithm */
/*                 179, Communications of the ACM 17, 3 (March 1974), */
/*                 pp. 156. */
/* ***ROUTINES CALLED  ALBETA, R1MACH, XERMSG */
/* ***REVISION HISTORY  (YYMMDD) */
/*   770401  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890531  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ) */
/*   900326  Removed duplicate information from DESCRIPTION section. */
/*           (WRB) */
/*   920528  DESCRIPTION and REFERENCES sections revised.  (WRB) */
/* ***END PROLOGUE  BETAI */
/* ***FIRST EXECUTABLE STATEMENT  BETAI */
    if (first) {
	eps = r1mach_(&c__3);
	alneps = log(eps);
	sml = r1mach_(&c__1);
	alnsml = log(sml);
    }
    first = FALSE_;

    if (*x < 0.f || *x > 1.f) {
	xermsg_("SLATEC", "BETAI", "X IS NOT IN THE RANGE (0,1)", &c__1, &
		c__2, (ftnlen)6, (ftnlen)5, (ftnlen)27);
    }
    if (*pin <= 0.f || *qin <= 0.f) {
	xermsg_("SLATEC", "BETAI", "P AND/OR Q IS LE ZERO", &c__2, &c__2, (
		ftnlen)6, (ftnlen)5, (ftnlen)21);
    }

    y = *x;
    p = *pin;
    q = *qin;
    if (q <= p && *x < .8f) {
	goto L20;
    }
    if (*x < .2f) {
	goto L20;
    }
    y = 1.f - y;
    p = *qin;
    q = *pin;

L20:
    if ((p + q) * y / (p + 1.f) < eps) {
	goto L80;
    }

/* EVALUATE THE INFINITE SUM FIRST. */
/* TERM WILL EQUAL Y**P/BETA(PS,P) * (1.-PS)I * Y**I / FAC(I) */

    ps = q - r_int(&q);
    if (ps == 0.f) {
	ps = 1.f;
    }
    xb = p * log(y) - albeta_(&ps, &p) - log(p);
    ret_val = 0.f;
    if (xb < alnsml) {
	goto L40;
    }

    ret_val = exp(xb);
    term = ret_val * p;
    if (ps == 1.f) {
	goto L40;
    }

/* Computing MAX */
    r__1 = alneps / log(y);
    n = dmax(r__1,4.f);
    i__1 = n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	term = term * (i__ - ps) * y / i__;
	ret_val += term / (p + i__);
/* L30: */
    }

/* NOW EVALUATE THE FINITE SUM, MAYBE. */

L40:
    if (q <= 1.f) {
	goto L70;
    }

    xb = p * log(y) + q * log(1.f - y) - albeta_(&p, &q) - log(q);
/* Computing MAX */
    r__1 = xb / alnsml;
    ib = dmax(r__1,0.f);
    term = exp(xb - ib * alnsml);
    c__ = 1.f / (1.f - y);
    p1 = q * c__ / (p + q - 1.f);

    finsum = 0.f;
    n = q;
    if (q == (real) n) {
	--n;
    }
    i__1 = n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (p1 <= 1.f && term / eps <= finsum) {
	    goto L60;
	}
	term = (q - i__ + 1) * c__ * term / (p + q - i__);

	if (term > 1.f) {
	    --ib;
	}
	if (term > 1.f) {
	    term *= sml;
	}

	if (ib == 0) {
	    finsum += term;
	}
/* L50: */
    }

L60:
    ret_val += finsum;
L70:
    if (y != *x || p != *pin) {
	ret_val = 1.f - ret_val;
    }
/* Computing MAX */
    r__1 = dmin(ret_val,1.f);
    ret_val = dmax(r__1,0.f);
    return ret_val;

L80:
    ret_val = 0.f;
    xb = p * log((dmax(y,sml))) - log(p) - albeta_(&p, &q);
    if (xb > alnsml && y != 0.f) {
	ret_val = exp(xb);
    }
    if (y != *x || p != *pin) {
	ret_val = 1.f - ret_val;
    }
    return ret_val;

} /* betai_ */

