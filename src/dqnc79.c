/* dqnc79.f -- translated by f2c (version 12.02.01).
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

static integer c__5 = 5;
static integer c__14 = 14;
static doublereal c_b6 = 1.;
static doublereal c_b7 = 2.;
static integer c__4 = 4;
static integer c__2 = 2;
static integer c__1 = 1;
static integer c_n1 = -1;

/* DECK DQNC79 */
/* Subroutine */ int dqnc79_(D_fp fun, doublereal *a, doublereal *b, 
	doublereal *err, doublereal *ans, integer *ierr, integer *k)
{
    /* Initialized data */

    static integer kml = 7;
    static integer kmx = 5000;
    static integer nlmn = 2;
    static logical first = TRUE_;

    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1, d__2;

    /* Local variables */
    static doublereal c__, f[13];
    static integer i__, l;
    static doublereal f1[99], f2[99], f3[99], f4[99], f5[99], f6[99], f7[99], 
	    q7, w1, w2, w3, w4, aa[99], ae, ce, ee, ef, hh[99], q13;
    static integer lr[99];
    static doublereal vl[99], vr, q7l, sq2, q7r[99];
    static integer nib, lmn;
    static doublereal eps, tol;
    static integer lmx;
    static doublereal area, bank;
    static integer nlmx;
    static doublereal test;
    static integer nbits;
    extern doublereal d1mach_(integer *);
    extern integer i1mach_(integer *);
    static doublereal blocal;
    extern /* Subroutine */ int xermsg_(char *, char *, char *, integer *, 
	    integer *, ftnlen, ftnlen, ftnlen);

/* ***BEGIN PROLOGUE  DQNC79 */
/* ***PURPOSE  Integrate a function using a 7-point adaptive Newton-Cotes */
/*            quadrature rule. */
/* ***LIBRARY   SLATEC */
/* ***CATEGORY  H2A1A1 */
/* ***TYPE      DOUBLE PRECISION (QNC79-S, DQNC79-D) */
/* ***KEYWORDS  ADAPTIVE QUADRATURE, INTEGRATION, NEWTON-COTES */
/* ***AUTHOR  Kahaner, D. K., (NBS) */
/*           Jones, R. E., (SNLA) */
/* ***DESCRIPTION */

/*     Abstract  *** a DOUBLE PRECISION routine *** */
/*       DQNC79 is a general purpose program for evaluation of */
/*       one dimensional integrals of user defined functions. */
/*       DQNC79 will pick its own points for evaluation of the */
/*       integrand and these will vary from problem to problem. */
/*       Thus, DQNC79 is not designed to integrate over data sets. */
/*       Moderately smooth integrands will be integrated efficiently */
/*       and reliably.  For problems with strong singularities, */
/*       oscillations etc., the user may wish to use more sophis- */
/*       ticated routines such as those in QUADPACK.  One measure */
/*       of the reliability of DQNC79 is the output parameter K, */
/*       giving the number of integrand evaluations that were needed. */

/*     Description of Arguments */

/*     --Input--* FUN, A, B, ERR are DOUBLE PRECISION * */
/*       FUN  - name of external function to be integrated.  This name */
/*              must be in an EXTERNAL statement in your calling */
/*              program.  You must write a Fortran function to evaluate */
/*              FUN.  This should be of the form */
/*                    DOUBLE PRECISION FUNCTION FUN (X) */
/*              C */
/*              C     X can vary from A to B */
/*              C     FUN(X) should be finite for all X on interval. */
/*              C */
/*                    FUN = ... */
/*                    RETURN */
/*                    END */
/*       A    - lower limit of integration */
/*       B    - upper limit of integration (may be less than A) */
/*       ERR  - is a requested error tolerance.  Normally, pick a value */
/*              0 .LT. ERR .LT. 1.0D-8. */

/*     --Output-- */
/*       ANS  - computed value of the integral.  Hopefully, ANS is */
/*              accurate to within ERR * integral of ABS(FUN(X)). */
/*       IERR - a status code */
/*            - Normal codes */
/*               1  ANS most likely meets requested error tolerance. */
/*              -1  A equals B, or A and B are too nearly equal to */
/*                  allow normal integration.  ANS is set to zero. */
/*            - Abnormal code */
/*               2  ANS probably does not meet requested error tolerance. */
/*       K    - the number of function evaluations actually used to do */
/*              the integration.  A value of K .GT. 1000 indicates a */
/*              difficult problem; other programs may be more efficient. */
/*              DQNC79 will gracefully give up if K exceeds 2000. */

/* ***REFERENCES  (NONE) */
/* ***ROUTINES CALLED  D1MACH, I1MACH, XERMSG */
/* ***REVISION HISTORY  (YYMMDD) */
/*   790601  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890911  Removed unnecessary intrinsics.  (WRB) */
/*   890911  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ) */
/*   920218  Code redone to parallel QNC79.  (WRB) */
/*   930120  Increase array size 80->99, and KMX 2000->5000 for SUN -r8 */
/*           wordlength.  (RWC) */
/* ***END PROLOGUE  DQNC79 */
/*     .. Scalar Arguments .. */
/*     .. Function Arguments .. */
/*     .. Local Scalars .. */
/*     .. Local Arrays .. */
/*     .. External Functions .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */
/*     .. Save statement .. */
/*     .. Data statements .. */
/* ***FIRST EXECUTABLE STATEMENT  DQNC79 */
    if (first) {
	w1 = .29285714285714287;
	w2 = 1.5428571428571429;
	w3 = .19285714285714287;
	w4 = 1.9428571428571428;
	nbits = (integer) (d1mach_(&c__5) * i1mach_(&c__14) / .30102);
/* Computing MIN */
	i__1 = 99, i__2 = (nbits << 2) / 5;
	nlmx = min(i__1,i__2);
	sq2 = sqrt(2.);
    }
    first = FALSE_;
    *ans = 0.;
    *ierr = 1;
    ce = 0.;
    if (*a == *b) {
	goto L260;
    }
    lmx = nlmx;
    lmn = nlmn;
    if (*b == 0.) {
	goto L100;
    }
    if (d_sign(&c_b6, b) * *a <= 0.) {
	goto L100;
    }
    c__ = (d__1 = 1. - *a / *b, abs(d__1));
    if (c__ > .1) {
	goto L100;
    }
    if (c__ <= 0.) {
	goto L260;
    }
    nib = (integer) (.5 - log(c__) / log(2.));
/* Computing MIN */
    i__1 = nlmx, i__2 = nbits - nib - 4;
    lmx = min(i__1,i__2);
    if (lmx < 2) {
	goto L260;
    }
    lmn = min(lmn,lmx);
L100:
/* Computing MAX */
    i__1 = 5 - nbits;
    d__1 = abs(*err), d__2 = pow_di(&c_b7, &i__1);
    tol = max(d__1,d__2);
    if (*err == 0.) {
	tol = sqrt(d1mach_(&c__4));
    }
    eps = tol;
    hh[0] = (*b - *a) / 12.;
    aa[0] = *a;
    lr[0] = 1;
    for (i__ = 1; i__ <= 11; i__ += 2) {
	d__1 = *a + (i__ - 1) * hh[0];
	f[i__ - 1] = (*fun)(&d__1);
/* L110: */
    }
    blocal = *b;
    f[12] = (*fun)(&blocal);
    *k = 7;
    l = 1;
    area = 0.;
    q7 = 0.;
    ef = 1.003921568627451;
    bank = 0.;

/*     Compute refined estimates, estimate the error, etc. */

L120:
    for (i__ = 2; i__ <= 12; i__ += 2) {
	d__1 = aa[l - 1] + (i__ - 1) * hh[l - 1];
	f[i__ - 1] = (*fun)(&d__1);
/* L130: */
    }
    *k += 6;

/*     Compute left and right half estimates */

    q7l = hh[l - 1] * (w1 * (f[0] + f[6]) + w2 * (f[1] + f[5]) + (w3 * (f[2] 
	    + f[4]) + w4 * f[3]));
    q7r[l - 1] = hh[l - 1] * (w1 * (f[6] + f[12]) + w2 * (f[7] + f[11]) + (w3 
	    * (f[8] + f[10]) + w4 * f[9]));

/*     Update estimate of integral of absolute value */

    area += abs(q7l) + (d__1 = q7r[l - 1], abs(d__1)) - abs(q7);

/*     Do not bother to test convergence before minimum refinement level */

    if (l < lmn) {
	goto L180;
    }

/*     Estimate the error in new value for whole interval, Q13 */

    q13 = q7l + q7r[l - 1];
    ee = (d__1 = q7 - q13, abs(d__1)) * ef;

/*     Compute nominal allowed error */

    ae = eps * area;

/*     Borrow from bank account, but not too much */

/* Computing MIN */
    d__1 = ae + bank * .8, d__2 = ae * 10.;
    test = min(d__1,d__2);

/*     Don't ask for excessive accuracy */

/* Computing MAX */
    d__1 = test, d__2 = tol * abs(q13), d__1 = max(d__1,d__2), d__2 = tol * 
	    3e-5 * area;
    test = max(d__1,d__2);

/*     Now, did this interval pass or not? */

    if (ee - test <= 0.) {
	goto L150;
    } else {
	goto L170;
    }

/*     Have hit maximum refinement level -- penalize the cumulative error */

L140:
    ce += q7 - q13;
    goto L160;

/*     On good intervals accumulate the theoretical estimate */

L150:
    ce += (q7 - q13) / 255.;

/*     Update the bank account.  Don't go into debt. */

L160:
    bank += ae - ee;
    if (bank < 0.) {
	bank = 0.;
    }

/*     Did we just finish a left half or a right half? */

    if (lr[l - 1] <= 0) {
	goto L190;
    } else {
	goto L210;
    }

/*     Consider the left half of next deeper level */

L170:
    if (*k > kmx) {
	lmx = min(kml,lmx);
    }
    if (l >= lmx) {
	goto L140;
    }
L180:
    ++l;
    eps *= .5;
    if (l <= 17) {
	ef /= sq2;
    }
    hh[l - 1] = hh[l - 2] * .5;
    lr[l - 1] = -1;
    aa[l - 1] = aa[l - 2];
    q7 = q7l;
    f1[l - 1] = f[6];
    f2[l - 1] = f[7];
    f3[l - 1] = f[8];
    f4[l - 1] = f[9];
    f5[l - 1] = f[10];
    f6[l - 1] = f[11];
    f7[l - 1] = f[12];
    f[12] = f[6];
    f[10] = f[5];
    f[8] = f[4];
    f[6] = f[3];
    f[4] = f[2];
    f[2] = f[1];
    goto L120;

/*     Proceed to right half at this level */

L190:
    vl[l - 1] = q13;
L200:
    q7 = q7r[l - 2];
    lr[l - 1] = 1;
    aa[l - 1] += hh[l - 1] * 12.;
    f[0] = f1[l - 1];
    f[2] = f2[l - 1];
    f[4] = f3[l - 1];
    f[6] = f4[l - 1];
    f[8] = f5[l - 1];
    f[10] = f6[l - 1];
    f[12] = f7[l - 1];
    goto L120;

/*     Left and right halves are done, so go back up a level */

L210:
    vr = q13;
L220:
    if (l <= 1) {
	goto L250;
    }
    if (l <= 17) {
	ef *= sq2;
    }
    eps *= 2.;
    --l;
    if (lr[l - 1] <= 0) {
	goto L230;
    } else {
	goto L240;
    }
L230:
    vl[l - 1] = vl[l] + vr;
    goto L200;
L240:
    vr = vl[l] + vr;
    goto L220;

/*     Exit */

L250:
    *ans = vr;
    if (abs(ce) <= tol * 2. * area) {
	goto L270;
    }
    *ierr = 2;
    xermsg_("SLATEC", "DQNC79", "ANS is probably insufficiently accurate.", &
	    c__2, &c__1, (ftnlen)6, (ftnlen)6, (ftnlen)40);
    goto L270;
L260:
    *ierr = -1;
    xermsg_("SLATEC", "DQNC79", "A and B are too nearly equal to allow norma"
	    "l integration. $$ANS is set to zero and IERR to -1.", &c_n1, &
	    c_n1, (ftnlen)6, (ftnlen)6, (ftnlen)94);
L270:
    return 0;
} /* dqnc79_ */

