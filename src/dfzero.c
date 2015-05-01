/* dfzero.f -- translated by f2c (version 12.02.01).
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
static doublereal c_b3 = 1.;

/* DECK DFZERO */
/* Subroutine */ int dfzero_(D_fp f, doublereal *b, doublereal *c__, 
	doublereal *r__, doublereal *re, doublereal *ae, integer *iflag)
{
    /* System generated locals */
    doublereal d__1, d__2;

    /* Local variables */
    static doublereal a, p, q, t, z__, fa, fb, fc;
    static integer ic;
    static doublereal aw, er, fx, fz, rw, cmb, tol, acmb, acbs;
    extern doublereal d1mach_(integer *);
    static integer kount;

/* ***BEGIN PROLOGUE  DFZERO */
/* ***PURPOSE  Search for a zero of a function F(X) in a given interval */
/*            (B,C).  It is designed primarily for problems where F(B) */
/*            and F(C) have opposite signs. */
/* ***LIBRARY   SLATEC */
/* ***CATEGORY  F1B */
/* ***TYPE      DOUBLE PRECISION (FZERO-S, DFZERO-D) */
/* ***KEYWORDS  BISECTION, NONLINEAR, ROOTS, ZEROS */
/* ***AUTHOR  Shampine, L. F., (SNLA) */
/*           Watts, H. A., (SNLA) */
/* ***DESCRIPTION */

/*     DFZERO searches for a zero of a DOUBLE PRECISION function F(X) */
/*     between the given DOUBLE PRECISION values B and C until the width */
/*     of the interval (B,C) has collapsed to within a tolerance */
/*     specified by the stopping criterion, */
/*        ABS(B-C) .LE. 2.*(RW*ABS(B)+AE). */
/*     The method used is an efficient combination of bisection and the */
/*     secant rule and is due to T. J. Dekker. */

/*     Description Of Arguments */

/*   F     :EXT   - Name of the DOUBLE PRECISION external function.  This */
/*                  name must be in an EXTERNAL statement in the calling */
/*                  program.  F must be a function of one DOUBLE */
/*                  PRECISION argument. */

/*   B     :INOUT - One end of the DOUBLE PRECISION interval (B,C).  The */
/*                  value returned for B usually is the better */
/*                  approximation to a zero of F. */

/*   C     :INOUT - The other end of the DOUBLE PRECISION interval (B,C) */

/*   R     :IN    - A (better) DOUBLE PRECISION guess of a zero of F */
/*                  which could help in speeding up convergence.  If F(B) */
/*                  and F(R) have opposite signs, a root will be found in */
/*                  the interval (B,R);  if not, but F(R) and F(C) have */
/*                  opposite signs, a root will be found in the interval */
/*                  (R,C);  otherwise, the interval (B,C) will be */
/*                  searched for a possible root.  When no better guess */
/*                  is known, it is recommended that R be set to B or C, */
/*                  since if R is not interior to the interval (B,C), it */
/*                  will be ignored. */

/*   RE    :IN    - Relative error used for RW in the stopping criterion. */
/*                  If the requested RE is less than machine precision, */
/*                  then RW is set to approximately machine precision. */

/*   AE    :IN    - Absolute error used in the stopping criterion.  If */
/*                  the given interval (B,C) contains the origin, then a */
/*                  nonzero value should be chosen for AE. */

/*   IFLAG :OUT   - A status code.  User must check IFLAG after each */
/*                  call.  Control returns to the user from DFZERO in all */
/*                  cases. */

/*                1  B is within the requested tolerance of a zero. */
/*                   The interval (B,C) collapsed to the requested */
/*                   tolerance, the function changes sign in (B,C), and */
/*                   F(X) decreased in magnitude as (B,C) collapsed. */

/*                2  F(B) = 0.  However, the interval (B,C) may not have */
/*                   collapsed to the requested tolerance. */

/*                3  B may be near a singular point of F(X). */
/*                   The interval (B,C) collapsed to the requested tol- */
/*                   erance and the function changes sign in (B,C), but */
/*                   F(X) increased in magnitude as (B,C) collapsed, i.e. */
/*                     ABS(F(B out)) .GT. MAX(ABS(F(B in)),ABS(F(C in))) */

/*                4  No change in sign of F(X) was found although the */
/*                   interval (B,C) collapsed to the requested tolerance. */
/*                   The user must examine this case and decide whether */
/*                   B is near a local minimum of F(X), or B is near a */
/*                   zero of even multiplicity, or neither of these. */

/*                5  Too many (.GT. 500) function evaluations used. */

/* ***REFERENCES  L. F. Shampine and H. A. Watts, FZERO, a root-solving */
/*                 code, Report SC-TM-70-631, Sandia Laboratories, */
/*                 September 1970. */
/*               T. J. Dekker, Finding a zero by means of successive */
/*                 linear interpolation, Constructive Aspects of the */
/*                 Fundamental Theorem of Algebra, edited by B. Dejon */
/*                 and P. Henrici, Wiley-Interscience, 1969. */
/* ***ROUTINES CALLED  D1MACH */
/* ***REVISION HISTORY  (YYMMDD) */
/*   700901  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890531  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   920501  Reformatted the REFERENCES section.  (WRB) */
/* ***END PROLOGUE  DFZERO */

/* ***FIRST EXECUTABLE STATEMENT  DFZERO */

/*   ER is two times the computer unit roundoff value which is defined */
/*   here by the function D1MACH. */

    er = d1mach_(&c__4) * 2.;

/*   Initialize. */

    z__ = *r__;
    if (*r__ <= min(*b,*c__) || *r__ >= max(*b,*c__)) {
	z__ = *c__;
    }
    rw = max(*re,er);
    aw = max(*ae,0.);
    ic = 0;
    t = z__;
    fz = (*f)(&t);
    fc = fz;
    t = *b;
    fb = (*f)(&t);
    kount = 2;
    if (d_sign(&c_b3, &fz) == d_sign(&c_b3, &fb)) {
	goto L1;
    }
    *c__ = z__;
    goto L2;
L1:
    if (z__ == *c__) {
	goto L2;
    }
    t = *c__;
    fc = (*f)(&t);
    kount = 3;
    if (d_sign(&c_b3, &fz) == d_sign(&c_b3, &fc)) {
	goto L2;
    }
    *b = z__;
    fb = fz;
L2:
    a = *c__;
    fa = fc;
    acbs = (d__1 = *b - *c__, abs(d__1));
/* Computing MAX */
    d__1 = abs(fb), d__2 = abs(fc);
    fx = max(d__1,d__2);

L3:
    if (abs(fc) >= abs(fb)) {
	goto L4;
    }

/*   Perform interchange. */

    a = *b;
    fa = fb;
    *b = *c__;
    fb = fc;
    *c__ = a;
    fc = fa;

L4:
    cmb = (*c__ - *b) * .5;
    acmb = abs(cmb);
    tol = rw * abs(*b) + aw;

/*   Test stopping criterion and function count. */

    if (acmb <= tol) {
	goto L10;
    }
    if (fb == 0.) {
	goto L11;
    }
    if (kount >= 500) {
	goto L14;
    }

/*   Calculate new iterate implicitly as B+P/Q, where we arrange */
/*   P .GE. 0.  The implicit form is used to prevent overflow. */

    p = (*b - a) * fb;
    q = fa - fb;
    if (p >= 0.) {
	goto L5;
    }
    p = -p;
    q = -q;

/*   Update A and check for satisfactory reduction in the size of the */
/*   bracketing interval.  If not, perform bisection. */

L5:
    a = *b;
    fa = fb;
    ++ic;
    if (ic < 4) {
	goto L6;
    }
    if (acmb * 8. >= acbs) {
	goto L8;
    }
    ic = 0;
    acbs = acmb;

/*   Test for too small a change. */

L6:
    if (p > abs(q) * tol) {
	goto L7;
    }

/*   Increment by TOLerance. */

    *b += d_sign(&tol, &cmb);
    goto L9;

/*   Root ought to be between B and (C+B)/2. */

L7:
    if (p >= cmb * q) {
	goto L8;
    }

/*   Use secant rule. */

    *b += p / q;
    goto L9;

/*   Use bisection (C+B)/2. */

L8:
    *b += cmb;

/*   Have completed computation for new iterate B. */

L9:
    t = *b;
    fb = (*f)(&t);
    ++kount;

/*   Decide whether next step is interpolation or extrapolation. */

    if (d_sign(&c_b3, &fb) != d_sign(&c_b3, &fc)) {
	goto L3;
    }
    *c__ = a;
    fc = fa;
    goto L3;

/*   Finished.  Process results for proper setting of IFLAG. */

L10:
    if (d_sign(&c_b3, &fb) == d_sign(&c_b3, &fc)) {
	goto L13;
    }
    if (abs(fb) > fx) {
	goto L12;
    }
    *iflag = 1;
    return 0;
L11:
    *iflag = 2;
    return 0;
L12:
    *iflag = 3;
    return 0;
L13:
    *iflag = 4;
    return 0;
L14:
    *iflag = 5;
    return 0;
} /* dfzero_ */

