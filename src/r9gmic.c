/* r9gmic.f -- translated by f2c (version 12.02.01).
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

/* DECK R9GMIC */
doublereal r9gmic_(real *a, real *x, real *alx)
{
    /* Initialized data */

    static real euler = .5772156649015329f;
    static real eps = 0.f;
    static real bot = 0.f;

    /* System generated locals */
    integer i__1;
    real ret_val, r__1;

    /* Local variables */
    static integer k, m;
    static real s, t;
    static integer ma;
    static real fk, fm, te;
    static integer mm1;
    static real fkp1, alng, sgng;
    extern doublereal r1mach_(integer *), alngam_(real *);
    extern /* Subroutine */ int xermsg_(char *, char *, char *, integer *, 
	    integer *, ftnlen, ftnlen, ftnlen);

/* ***BEGIN PROLOGUE  R9GMIC */
/* ***SUBSIDIARY */
/* ***PURPOSE  Compute the complementary incomplete Gamma function for A */
/*            near a negative integer and for small X. */
/* ***LIBRARY   SLATEC (FNLIB) */
/* ***CATEGORY  C7E */
/* ***TYPE      SINGLE PRECISION (R9GMIC-S, D9GMIC-D) */
/* ***KEYWORDS  COMPLEMENTARY INCOMPLETE GAMMA FUNCTION, FNLIB, SMALL X, */
/*             SPECIAL FUNCTIONS */
/* ***AUTHOR  Fullerton, W., (LANL) */
/* ***DESCRIPTION */

/* Compute the complementary incomplete gamma function for A near */
/* a negative integer and for small X. */

/* ***REFERENCES  (NONE) */
/* ***ROUTINES CALLED  ALNGAM, R1MACH, XERMSG */
/* ***REVISION HISTORY  (YYMMDD) */
/*   770701  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890531  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ) */
/*   900720  Routine changed from user-callable to subsidiary.  (WRB) */
/* ***END PROLOGUE  R9GMIC */
/* ***FIRST EXECUTABLE STATEMENT  R9GMIC */
    if (eps == 0.f) {
	eps = .5f * r1mach_(&c__3);
    }
    if (bot == 0.f) {
	bot = log(r1mach_(&c__1));
    }

    if (*a > 0.f) {
	xermsg_("SLATEC", "R9GMIC", "A MUST BE NEAR A NEGATIVE INTEGER", &
		c__2, &c__2, (ftnlen)6, (ftnlen)6, (ftnlen)33);
    }
    if (*x <= 0.f) {
	xermsg_("SLATEC", "R9GMIC", "X MUST BE GT ZERO", &c__3, &c__2, (
		ftnlen)6, (ftnlen)6, (ftnlen)17);
    }

    ma = *a - .5f;
    fm = (real) (-ma);
    m = -ma;

    te = 1.f;
    t = 1.f;
    s = t;
    for (k = 1; k <= 200; ++k) {
	fkp1 = (real) (k + 1);
	te = -(*x) * te / (fm + fkp1);
	t = te / fkp1;
	s += t;
	if (dabs(t) < eps * s) {
	    goto L30;
	}
/* L20: */
    }
    xermsg_("SLATEC", "R9GMIC", "NO CONVERGENCE IN 200 TERMS OF CONTINUED FR"
	    "ACTION", &c__4, &c__2, (ftnlen)6, (ftnlen)6, (ftnlen)49);

L30:
    ret_val = -(*alx) - euler + *x * s / (fm + 1.f);
    if (m == 0) {
	return ret_val;
    }

    if (m == 1) {
	ret_val = -ret_val - 1.f + 1.f / *x;
    }
    if (m == 1) {
	return ret_val;
    }

    te = fm;
    t = 1.f;
    s = t;
    mm1 = m - 1;
    i__1 = mm1;
    for (k = 1; k <= i__1; ++k) {
	fk = (real) k;
	te = -(*x) * te / fk;
	t = te / (fm - fk);
	s += t;
	if (dabs(t) < eps * dabs(s)) {
	    goto L50;
	}
/* L40: */
    }

L50:
    i__1 = m;
    for (k = 1; k <= i__1; ++k) {
	ret_val += 1.f / k;
/* L60: */
    }

    sgng = 1.f;
    if (m % 2 == 1) {
	sgng = -1.f;
    }
    r__1 = fm + 1.f;
    alng = log(ret_val) - alngam_(&r__1);

    ret_val = 0.f;
    if (alng > bot) {
	ret_val = sgng * exp(alng);
    }
    if (s != 0.f) {
	r__1 = exp(-fm * *alx + log(dabs(s) / fm));
	ret_val += r_sign(&r__1, &s);
    }

    if (ret_val == 0.f && s == 0.f) {
	xermsg_("SLATEC", "R9GMIC", "RESULT UNDERFLOWS", &c__1, &c__1, (
		ftnlen)6, (ftnlen)6, (ftnlen)17);
    }
    return ret_val;

} /* r9gmic_ */

