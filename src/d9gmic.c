/* d9gmic.f -- translated by f2c (version 12.02.01).
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
static integer c__4 = 4;

/* DECK D9GMIC */
doublereal d9gmic_(doublereal *a, doublereal *x, doublereal *alx)
{
    /* Initialized data */

    static doublereal euler = .5772156649015328606065120900824;
    static logical first = TRUE_;

    /* System generated locals */
    integer i__1;
    doublereal ret_val, d__1;

    /* Local variables */
    static integer k, m;
    static doublereal s, t, fk, fm, te;
    static integer mm1;
    static doublereal bot, eps, fkp1, alng, sgng;
    extern doublereal d1mach_(integer *), dlngam_(doublereal *);
    extern /* Subroutine */ int xermsg_(char *, char *, char *, integer *, 
	    integer *, ftnlen, ftnlen, ftnlen);

/* ***BEGIN PROLOGUE  D9GMIC */
/* ***SUBSIDIARY */
/* ***PURPOSE  Compute the complementary incomplete Gamma function for A */
/*            near a negative integer and X small. */
/* ***LIBRARY   SLATEC (FNLIB) */
/* ***CATEGORY  C7E */
/* ***TYPE      DOUBLE PRECISION (R9GMIC-S, D9GMIC-D) */
/* ***KEYWORDS  COMPLEMENTARY INCOMPLETE GAMMA FUNCTION, FNLIB, SMALL X, */
/*             SPECIAL FUNCTIONS */
/* ***AUTHOR  Fullerton, W., (LANL) */
/* ***DESCRIPTION */

/* Compute the complementary incomplete gamma function for A near */
/* a negative integer and for small X. */

/* ***REFERENCES  (NONE) */
/* ***ROUTINES CALLED  D1MACH, DLNGAM, XERMSG */
/* ***REVISION HISTORY  (YYMMDD) */
/*   770701  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890911  Removed unnecessary intrinsics.  (WRB) */
/*   890911  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ) */
/*   900720  Routine changed from user-callable to subsidiary.  (WRB) */
/* ***END PROLOGUE  D9GMIC */
/* ***FIRST EXECUTABLE STATEMENT  D9GMIC */
    if (first) {
	eps = d1mach_(&c__3) * .5;
	bot = log(d1mach_(&c__1));
    }
    first = FALSE_;

    if (*a > 0.) {
	xermsg_("SLATEC", "D9GMIC", "A MUST BE NEAR A NEGATIVE INTEGER", &
		c__2, &c__2, (ftnlen)6, (ftnlen)6, (ftnlen)33);
    }
    if (*x <= 0.) {
	xermsg_("SLATEC", "D9GMIC", "X MUST BE GT ZERO", &c__3, &c__2, (
		ftnlen)6, (ftnlen)6, (ftnlen)17);
    }

    m = (integer) (-(*a - .5));
    fm = (doublereal) m;

    te = 1.;
    t = 1.;
    s = t;
    for (k = 1; k <= 200; ++k) {
	fkp1 = (doublereal) (k + 1);
	te = -(*x) * te / (fm + fkp1);
	t = te / fkp1;
	s += t;
	if (abs(t) < eps * s) {
	    goto L30;
	}
/* L20: */
    }
    xermsg_("SLATEC", "D9GMIC", "NO CONVERGENCE IN 200 TERMS OF CONTINUED FR"
	    "ACTION", &c__4, &c__2, (ftnlen)6, (ftnlen)6, (ftnlen)49);

L30:
    ret_val = -(*alx) - euler + *x * s / (fm + 1.);
    if (m == 0) {
	return ret_val;
    }

    if (m == 1) {
	ret_val = -ret_val - 1. + 1. / *x;
    }
    if (m == 1) {
	return ret_val;
    }

    te = fm;
    t = 1.;
    s = t;
    mm1 = m - 1;
    i__1 = mm1;
    for (k = 1; k <= i__1; ++k) {
	fk = (doublereal) k;
	te = -(*x) * te / fk;
	t = te / (fm - fk);
	s += t;
	if (abs(t) < eps * abs(s)) {
	    goto L50;
	}
/* L40: */
    }

L50:
    i__1 = m;
    for (k = 1; k <= i__1; ++k) {
	ret_val += 1. / k;
/* L60: */
    }

    sgng = 1.;
    if (m % 2 == 1) {
	sgng = -1.;
    }
    d__1 = fm + 1.;
    alng = log(ret_val) - dlngam_(&d__1);

    ret_val = 0.;
    if (alng > bot) {
	ret_val = sgng * exp(alng);
    }
    if (s != 0.) {
	d__1 = exp(-fm * *alx + log(abs(s) / fm));
	ret_val += d_sign(&d__1, &s);
    }

    if (ret_val == 0. && s == 0.) {
	xermsg_("SLATEC", "D9GMIC", "RESULT UNDERFLOWS", &c__1, &c__1, (
		ftnlen)6, (ftnlen)6, (ftnlen)17);
    }
    return ret_val;

} /* d9gmic_ */

