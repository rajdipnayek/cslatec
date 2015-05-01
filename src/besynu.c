/* besynu.f -- translated by f2c (version 12.02.01).
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
static integer c__2 = 2;
static integer c__1 = 1;

/* DECK BESYNU */
/* Subroutine */ int besynu_(real *x, real *fnu, integer *n, real *y)
{
    /* Initialized data */

    static real x1 = 3.f;
    static real x2 = 20.f;
    static real pi = 3.14159265358979f;
    static real rthpi = .797884560802865f;
    static real hpi = 1.5707963267949f;
    static real cc[8] = { .577215664901533f,-.0420026350340952f,
	    -.0421977345555443f,.007218943246663f,-2.152416741149e-4f,
	    -2.01348547807e-5f,1.133027232e-6f,6.116095e-9f };

    /* System generated locals */
    integer i__1;
    real r__1, r__2;

    /* Local variables */
    static real a[120], f, g;
    static integer i__, j, k;
    static real p, q, s, a1, a2, g1, g2, s1, s2, t1, t2, cb[120], fc, ak, bk, 
	    ck, fk, fn, rb[120];
    static integer kk;
    static real cs, sa, sb, cx;
    static integer nn;
    static real tb, fx, tm, pt, rs, ss, st, rx, cp1, cp2, cs1, cs2, rp1, rp2, 
	    rs1, rs2, cbk, cck, arg, rbk, rck, fhs, fks, cpt, dnu, fmu;
    static integer inu;
    static real tol, etx, smu, rpt, dnu2, coef, relb, flrx;
    extern doublereal gamma_(real *);
    static real etest;
    extern doublereal r1mach_(integer *);
    extern /* Subroutine */ int xermsg_(char *, char *, char *, integer *, 
	    integer *, ftnlen, ftnlen, ftnlen);

/* ***BEGIN PROLOGUE  BESYNU */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to BESY */
/* ***LIBRARY   SLATEC */
/* ***TYPE      SINGLE PRECISION (BESYNU-S, DBSYNU-D) */
/* ***AUTHOR  Amos, D. E., (SNLA) */
/* ***DESCRIPTION */

/*     Abstract */
/*         BESYNU computes N member sequences of Y Bessel functions */
/*         Y/SUB(FNU+I-1)/(X), I=1,N for non-negative orders FNU and */
/*         positive X. Equations of the references are implemented on */
/*         small orders DNU for Y/SUB(DNU)/(X) and Y/SUB(DNU+1)/(X). */
/*         Forward recursion with the three term recursion relation */
/*         generates higher orders FNU+I-1, I=1,...,N. */

/*         To start the recursion FNU is normalized to the interval */
/*         -0.5.LE.DNU.LT.0.5. A special form of the power series is */
/*         implemented on 0.LT.X.LE.X1 while the Miller algorithm for the */
/*         K Bessel function in terms of the confluent hypergeometric */
/*         function U(FNU+0.5,2*FNU+1,I*X) is implemented on X1.LT.X.LE.X */
/*         Here I is the complex number SQRT(-1.). */
/*         For X.GT.X2, the asymptotic expansion for large X is used. */
/*         When FNU is a half odd integer, a special formula for */
/*         DNU=-0.5 and DNU+1.0=0.5 is used to start the recursion. */

/*         BESYNU assumes that a significant digit SINH(X) function is */
/*         available. */

/*     Description of Arguments */

/*         Input */
/*           X      - X.GT.0.0E0 */
/*           FNU    - Order of initial Y function, FNU.GE.0.0E0 */
/*           N      - Number of members of the sequence, N.GE.1 */

/*         Output */
/*           Y      - A vector whose first N components contain values */
/*                    for the sequence Y(I)=Y/SUB(FNU+I-1), I=1,N. */

/*     Error Conditions */
/*         Improper input arguments - a fatal error */
/*         Overflow - a fatal error */

/* ***SEE ALSO  BESY */
/* ***REFERENCES  N. M. Temme, On the numerical evaluation of the ordinary */
/*                 Bessel function of the second kind, Journal of */
/*                 Computational Physics 21, (1976), pp. 343-350. */
/*               N. M. Temme, On the numerical evaluation of the modified */
/*                 Bessel function of the third kind, Journal of */
/*                 Computational Physics 19, (1975), pp. 324-337. */
/* ***ROUTINES CALLED  GAMMA, R1MACH, XERMSG */
/* ***REVISION HISTORY  (YYMMDD) */
/*   800501  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ) */
/*   900326  Removed duplicate information from DESCRIPTION section. */
/*           (WRB) */
/*   900328  Added TYPE section.  (WRB) */
/*   900727  Added EXTERNAL statement.  (WRB) */
/*   910408  Updated the AUTHOR and REFERENCES sections.  (WRB) */
/*   920501  Reformatted the REFERENCES section.  (WRB) */
/* ***END PROLOGUE  BESYNU */

    /* Parameter adjustments */
    --y;

    /* Function Body */
/* ***FIRST EXECUTABLE STATEMENT  BESYNU */
    ak = r1mach_(&c__3);
    tol = dmax(ak,1e-15f);
    if (*x <= 0.f) {
	goto L270;
    }
    if (*fnu < 0.f) {
	goto L280;
    }
    if (*n < 1) {
	goto L290;
    }
    rx = 2.f / *x;
    inu = (integer) (*fnu + .5f);
    dnu = *fnu - inu;
    if (dabs(dnu) == .5f) {
	goto L260;
    }
    dnu2 = 0.f;
    if (dabs(dnu) < tol) {
	goto L10;
    }
    dnu2 = dnu * dnu;
L10:
    if (*x > x1) {
	goto L120;
    }

/*     SERIES FOR X.LE.X1 */

    a1 = 1.f - dnu;
    a2 = dnu + 1.f;
    t1 = 1.f / gamma_(&a1);
    t2 = 1.f / gamma_(&a2);
    if (dabs(dnu) > .1f) {
	goto L40;
    }
/*     SERIES FOR F0 TO RESOLVE INDETERMINACY FOR SMALL ABS(DNU) */
    s = cc[0];
    ak = 1.f;
    for (k = 2; k <= 8; ++k) {
	ak *= dnu2;
	tm = cc[k - 1] * ak;
	s += tm;
	if (dabs(tm) < tol) {
	    goto L30;
	}
/* L20: */
    }
L30:
    g1 = -(s + s);
    goto L50;
L40:
    g1 = (t1 - t2) / dnu;
L50:
    g2 = t1 + t2;
    smu = 1.f;
    fc = 1.f / pi;
    flrx = log(rx);
    fmu = dnu * flrx;
    tm = 0.f;
    if (dnu == 0.f) {
	goto L60;
    }
    tm = sin(dnu * hpi) / dnu;
    tm = (dnu + dnu) * tm * tm;
    fc = dnu / sin(dnu * pi);
    if (fmu != 0.f) {
	smu = sinh(fmu) / fmu;
    }
L60:
    f = fc * (g1 * cosh(fmu) + g2 * flrx * smu);
    fx = exp(fmu);
    p = fc * t1 * fx;
    q = fc * t2 / fx;
    g = f + tm * q;
    ak = 1.f;
    ck = 1.f;
    bk = 1.f;
    s1 = g;
    s2 = p;
    if (inu > 0 || *n > 1) {
	goto L90;
    }
    if (*x < tol) {
	goto L80;
    }
    cx = *x * *x * .25f;
L70:
    f = (ak * f + p + q) / (bk - dnu2);
    p /= ak - dnu;
    q /= ak + dnu;
    g = f + tm * q;
    ck = -ck * cx / ak;
    t1 = ck * g;
    s1 += t1;
    bk = bk + ak + ak + 1.f;
    ak += 1.f;
    s = dabs(t1) / (dabs(s1) + 1.f);
    if (s > tol) {
	goto L70;
    }
L80:
    y[1] = -s1;
    return 0;
L90:
    if (*x < tol) {
	goto L110;
    }
    cx = *x * *x * .25f;
L100:
    f = (ak * f + p + q) / (bk - dnu2);
    p /= ak - dnu;
    q /= ak + dnu;
    g = f + tm * q;
    ck = -ck * cx / ak;
    t1 = ck * g;
    s1 += t1;
    t2 = ck * (p - ak * g);
    s2 += t2;
    bk = bk + ak + ak + 1.f;
    ak += 1.f;
    s = dabs(t1) / (dabs(s1) + 1.f) + dabs(t2) / (dabs(s2) + 1.f);
    if (s > tol) {
	goto L100;
    }
L110:
    s2 = -s2 * rx;
    s1 = -s1;
    goto L160;
L120:
    coef = rthpi / sqrt(*x);
    if (*x > x2) {
	goto L210;
    }

/*     MILLER ALGORITHM FOR X1.LT.X.LE.X2 */

    etest = cos(pi * dnu) / (pi * *x * tol);
    fks = 1.f;
    fhs = .25f;
    fk = 0.f;
    rck = 2.f;
    cck = *x + *x;
    rp1 = 0.f;
    cp1 = 0.f;
    rp2 = 1.f;
    cp2 = 0.f;
    k = 0;
L130:
    ++k;
    fk += 1.f;
    ak = (fhs - dnu2) / (fks + fk);
    pt = fk + 1.f;
    rbk = rck / pt;
    cbk = cck / pt;
    rpt = rp2;
    cpt = cp2;
    rp2 = rbk * rpt - cbk * cpt - ak * rp1;
    cp2 = cbk * rpt + rbk * cpt - ak * cp1;
    rp1 = rpt;
    cp1 = cpt;
    rb[k - 1] = rbk;
    cb[k - 1] = cbk;
    a[k - 1] = ak;
    rck += 2.f;
    fks = fks + fk + fk + 1.f;
    fhs = fhs + fk + fk;
/* Computing MAX */
    r__1 = dabs(rp1), r__2 = dabs(cp1);
    pt = dmax(r__1,r__2);
/* Computing 2nd power */
    r__1 = rp1 / pt;
/* Computing 2nd power */
    r__2 = cp1 / pt;
    fc = r__1 * r__1 + r__2 * r__2;
    pt = pt * sqrt(fc) * fk;
    if (etest > pt) {
	goto L130;
    }
    kk = k;
    rs = 1.f;
    cs = 0.f;
    rp1 = 0.f;
    cp1 = 0.f;
    rp2 = 1.f;
    cp2 = 0.f;
    i__1 = k;
    for (i__ = 1; i__ <= i__1; ++i__) {
	rpt = rp2;
	cpt = cp2;
	rp2 = (rb[kk - 1] * rpt - cb[kk - 1] * cpt - rp1) / a[kk - 1];
	cp2 = (cb[kk - 1] * rpt + rb[kk - 1] * cpt - cp1) / a[kk - 1];
	rp1 = rpt;
	cp1 = cpt;
	rs += rp2;
	cs += cp2;
	--kk;
/* L140: */
    }
/* Computing MAX */
    r__1 = dabs(rs), r__2 = dabs(cs);
    pt = dmax(r__1,r__2);
/* Computing 2nd power */
    r__1 = rs / pt;
/* Computing 2nd power */
    r__2 = cs / pt;
    fc = r__1 * r__1 + r__2 * r__2;
    pt *= sqrt(fc);
    rs1 = (rp2 * (rs / pt) + cp2 * (cs / pt)) / pt;
    cs1 = (cp2 * (rs / pt) - rp2 * (cs / pt)) / pt;
    fc = hpi * (dnu - .5f) - *x;
    p = cos(fc);
    q = sin(fc);
    s1 = (cs1 * q - rs1 * p) * coef;
    if (inu > 0 || *n > 1) {
	goto L150;
    }
    y[1] = s1;
    return 0;
L150:
/* Computing MAX */
    r__1 = dabs(rp2), r__2 = dabs(cp2);
    pt = dmax(r__1,r__2);
/* Computing 2nd power */
    r__1 = rp2 / pt;
/* Computing 2nd power */
    r__2 = cp2 / pt;
    fc = r__1 * r__1 + r__2 * r__2;
    pt *= sqrt(fc);
    rpt = dnu + .5f - (rp1 * (rp2 / pt) + cp1 * (cp2 / pt)) / pt;
    cpt = *x - (cp1 * (rp2 / pt) - rp1 * (cp2 / pt)) / pt;
    cs2 = cs1 * cpt - rs1 * rpt;
    rs2 = rpt * cs1 + rs1 * cpt;
    s2 = (rs2 * q + cs2 * p) * coef / *x;

/*     FORWARD RECURSION ON THE THREE TERM RECURSION RELATION */

L160:
    ck = (dnu + dnu + 2.f) / *x;
    if (*n == 1) {
	--inu;
    }
    if (inu > 0) {
	goto L170;
    }
    if (*n > 1) {
	goto L190;
    }
    s1 = s2;
    goto L190;
L170:
    i__1 = inu;
    for (i__ = 1; i__ <= i__1; ++i__) {
	st = s2;
	s2 = ck * s2 - s1;
	s1 = st;
	ck += rx;
/* L180: */
    }
    if (*n == 1) {
	s1 = s2;
    }
L190:
    y[1] = s1;
    if (*n == 1) {
	return 0;
    }
    y[2] = s2;
    if (*n == 2) {
	return 0;
    }
    i__1 = *n;
    for (i__ = 3; i__ <= i__1; ++i__) {
	y[i__] = ck * y[i__ - 1] - y[i__ - 2];
	ck += rx;
/* L200: */
    }
    return 0;

/*     ASYMPTOTIC EXPANSION FOR LARGE X, X.GT.X2 */

L210:
    nn = 2;
    if (inu == 0 && *n == 1) {
	nn = 1;
    }
    dnu2 = dnu + dnu;
    fmu = 0.f;
    if (dabs(dnu2) < tol) {
	goto L220;
    }
    fmu = dnu2 * dnu2;
L220:
    arg = *x - hpi * (dnu + .5f);
    sa = sin(arg);
    sb = cos(arg);
    etx = *x * 8.f;
    i__1 = nn;
    for (k = 1; k <= i__1; ++k) {
	s1 = s2;
	t2 = (fmu - 1.f) / etx;
	ss = t2;
	relb = tol * dabs(t2);
	t1 = etx;
	s = 1.f;
	fn = 1.f;
	ak = 0.f;
	for (j = 1; j <= 13; ++j) {
	    t1 += etx;
	    ak += 8.f;
	    fn += ak;
	    t2 = -t2 * (fmu - fn) / t1;
	    s += t2;
	    t1 += etx;
	    ak += 8.f;
	    fn += ak;
	    t2 = t2 * (fmu - fn) / t1;
	    ss += t2;
	    if (dabs(t2) <= relb) {
		goto L240;
	    }
/* L230: */
	}
L240:
	s2 = coef * (s * sa + ss * sb);
	fmu = fmu + dnu * 8.f + 4.f;
	tb = sa;
	sa = -sb;
	sb = tb;
/* L250: */
    }
    if (nn > 1) {
	goto L160;
    }
    s1 = s2;
    goto L190;

/*     FNU=HALF ODD INTEGER CASE */

L260:
    coef = rthpi / sqrt(*x);
    s1 = coef * sin(*x);
    s2 = -coef * cos(*x);
    goto L160;


L270:
    xermsg_("SLATEC", "BESYNU", "X NOT GREATER THAN ZERO", &c__2, &c__1, (
	    ftnlen)6, (ftnlen)6, (ftnlen)23);
    return 0;
L280:
    xermsg_("SLATEC", "BESYNU", "FNU NOT ZERO OR POSITIVE", &c__2, &c__1, (
	    ftnlen)6, (ftnlen)6, (ftnlen)24);
    return 0;
L290:
    xermsg_("SLATEC", "BESYNU", "N NOT GREATER THAN 0", &c__2, &c__1, (ftnlen)
	    6, (ftnlen)6, (ftnlen)20);
    return 0;
} /* besynu_ */

