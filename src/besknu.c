/* besknu.f -- translated by f2c (version 12.02.01).
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

static integer c__12 = 12;
static integer c__5 = 5;
static integer c__3 = 3;
static integer c__2 = 2;
static integer c__1 = 1;

/* DECK BESKNU */
/* Subroutine */ int besknu_(real *x, real *fnu, integer *kode, integer *n, 
	real *y, integer *nz)
{
    /* Initialized data */

    static real x1 = 2.f;
    static real x2 = 17.f;
    static real pi = 3.14159265358979f;
    static real rthpi = 1.2533141373155f;
    static real cc[8] = { .577215664901533f,-.0420026350340952f,
	    -.0421977345555443f,.007218943246663f,-2.152416741149e-4f,
	    -2.01348547807e-5f,1.133027232e-6f,6.116095e-9f };

    /* System generated locals */
    integer i__1;

    /* Local variables */
    static real a[160], b[160], f;
    static integer i__, j, k;
    static real p, q, s, a1, a2, g1, g2, p1, p2, s1, s2, t1, t2, fc, ak, bk, 
	    ck, dk, fk;
    static integer kk;
    static real cx;
    static integer nn;
    static real ex, tm, pt, st, rx, fhs, fks, dnu, fmu;
    static integer inu;
    static real sqk, tol, smu, dnu2, coef, elim, flrx;
    static integer iflag;
    extern doublereal gamma_(real *);
    static integer koded;
    static real etest;
    extern integer i1mach_(integer *);
    extern doublereal r1mach_(integer *);
    extern /* Subroutine */ int xermsg_(char *, char *, char *, integer *, 
	    integer *, ftnlen, ftnlen, ftnlen);

/* ***BEGIN PROLOGUE  BESKNU */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to BESK */
/* ***LIBRARY   SLATEC */
/* ***TYPE      SINGLE PRECISION (BESKNU-S, DBSKNU-D) */
/* ***AUTHOR  Amos, D. E., (SNLA) */
/* ***DESCRIPTION */

/*     Abstract */
/*         BESKNU computes N member sequences of K Bessel functions */
/*         K/SUB(FNU+I-1)/(X), I=1,N for non-negative orders FNU and */
/*         positive X. Equations of the references are implemented on */
/*         small orders DNU for K/SUB(DNU)/(X) and K/SUB(DNU+1)/(X). */
/*         Forward recursion with the three term recursion relation */
/*         generates higher orders FNU+I-1, I=1,...,N. The parameter */
/*         KODE permits K/SUB(FNU+I-1)/(X) values or scaled values */
/*         EXP(X)*K/SUB(FNU+I-1)/(X), I=1,N to be returned. */

/*         To start the recursion FNU is normalized to the interval */
/*         -0.5.LE.DNU.LT.0.5. A special form of the power series is */
/*         implemented on 0.LT.X.LE.X1 while the Miller algorithm for the */
/*         K Bessel function in terms of the confluent hypergeometric */
/*         function U(FNU+0.5,2*FNU+1,X) is implemented on X1.LT.X.LE.X2. */
/*         For X.GT.X2, the asymptotic expansion for large X is used. */
/*         When FNU is a half odd integer, a special formula for */
/*         DNU=-0.5 and DNU+1.0=0.5 is used to start the recursion. */

/*         BESKNU assumes that a significant digit SINH(X) function is */
/*         available. */

/*     Description of Arguments */

/*         Input */
/*           X      - X.GT.0.0E0 */
/*           FNU    - Order of initial K function, FNU.GE.0.0E0 */
/*           N      - Number of members of the sequence, N.GE.1 */
/*           KODE   - A parameter to indicate the scaling option */
/*                    KODE= 1  returns */
/*                             Y(I)=       K/SUB(FNU+I-1)/(X) */
/*                                  I=1,...,N */
/*                        = 2  returns */
/*                             Y(I)=EXP(X)*K/SUB(FNU+I-1)/(X) */
/*                                  I=1,...,N */

/*         Output */
/*           Y      - A vector whose first N components contain values */
/*                    for the sequence */
/*                    Y(I)=       K/SUB(FNU+I-1)/(X), I=1,...,N or */
/*                    Y(I)=EXP(X)*K/SUB(FNU+I-1)/(X), I=1,...,N */
/*                    depending on KODE */
/*           NZ     - Number of components set to zero due to */
/*                    underflow, */
/*                    NZ= 0   , Normal return */
/*                    NZ.NE.0 , First NZ components of Y set to zero */
/*                              due to underflow, Y(I)=0.0E0,I=1,...,NZ */

/*     Error Conditions */
/*         Improper input arguments - a fatal error */
/*         Overflow - a fatal error */
/*         Underflow with KODE=1 - a non-fatal error (NZ.NE.0) */

/* ***SEE ALSO  BESK */
/* ***REFERENCES  N. M. Temme, On the numerical evaluation of the modified */
/*                 Bessel function of the third kind, Journal of */
/*                 Computational Physics 19, (1975), pp. 324-337. */
/* ***ROUTINES CALLED  GAMMA, I1MACH, R1MACH, XERMSG */
/* ***REVISION HISTORY  (YYMMDD) */
/*   790201  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ) */
/*   900326  Removed duplicate information from DESCRIPTION section. */
/*           (WRB) */
/*   900328  Added TYPE section.  (WRB) */
/*   900727  Added EXTERNAL statement.  (WRB) */
/*   910408  Updated the AUTHOR and REFERENCES sections.  (WRB) */
/*   920501  Reformatted the REFERENCES section.  (WRB) */
/* ***END PROLOGUE  BESKNU */

    /* Parameter adjustments */
    --y;

    /* Function Body */
/* ***FIRST EXECUTABLE STATEMENT  BESKNU */
    kk = -i1mach_(&c__12);
    elim = (kk * r1mach_(&c__5) - 3.f) * 2.303f;
    ak = r1mach_(&c__3);
    tol = dmax(ak,1e-15f);
    if (*x <= 0.f) {
	goto L350;
    }
    if (*fnu < 0.f) {
	goto L360;
    }
    if (*kode < 1 || *kode > 2) {
	goto L370;
    }
    if (*n < 1) {
	goto L380;
    }
    *nz = 0;
    iflag = 0;
    koded = *kode;
    rx = 2.f / *x;
    inu = (integer) (*fnu + .5f);
    dnu = *fnu - inu;
    if (dabs(dnu) == .5f) {
	goto L120;
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
    g1 = -s;
    goto L50;
L40:
    g1 = (t1 - t2) / (dnu + dnu);
L50:
    g2 = (t1 + t2) * .5f;
    smu = 1.f;
    fc = 1.f;
    flrx = log(rx);
    fmu = dnu * flrx;
    if (dnu == 0.f) {
	goto L60;
    }
    fc = dnu * pi;
    fc /= sin(fc);
    if (fmu != 0.f) {
	smu = sinh(fmu) / fmu;
    }
L60:
    f = fc * (g1 * cosh(fmu) + g2 * flrx * smu);
    fc = exp(fmu);
    p = fc * .5f / t2;
    q = .5f / (fc * t1);
    ak = 1.f;
    ck = 1.f;
    bk = 1.f;
    s1 = f;
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
    ck = ck * cx / ak;
    t1 = ck * f;
    s1 += t1;
    bk = bk + ak + ak + 1.f;
    ak += 1.f;
    s = dabs(t1) / (dabs(s1) + 1.f);
    if (s > tol) {
	goto L70;
    }
L80:
    y[1] = s1;
    if (koded == 1) {
	return 0;
    }
    y[1] = s1 * exp(*x);
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
    ck = ck * cx / ak;
    t1 = ck * f;
    s1 += t1;
    t2 = ck * (p - ak * f);
    s2 += t2;
    bk = bk + ak + ak + 1.f;
    ak += 1.f;
    s = dabs(t1) / (dabs(s1) + 1.f) + dabs(t2) / (dabs(s2) + 1.f);
    if (s > tol) {
	goto L100;
    }
L110:
    s2 *= rx;
    if (koded == 1) {
	goto L170;
    }
    f = exp(*x);
    s1 *= f;
    s2 *= f;
    goto L170;
L120:
    coef = rthpi / sqrt(*x);
    if (koded == 2) {
	goto L130;
    }
    if (*x > elim) {
	goto L330;
    }
    coef *= exp(-(*x));
L130:
    if (dabs(dnu) == .5f) {
	goto L340;
    }
    if (*x > x2) {
	goto L280;
    }

/*     MILLER ALGORITHM FOR X1.LT.X.LE.X2 */

    etest = cos(pi * dnu) / (pi * *x * tol);
    fks = 1.f;
    fhs = .25f;
    fk = 0.f;
    ck = *x + *x + 2.f;
    p1 = 0.f;
    p2 = 1.f;
    k = 0;
L140:
    ++k;
    fk += 1.f;
    ak = (fhs - dnu2) / (fks + fk);
    bk = ck / (fk + 1.f);
    pt = p2;
    p2 = bk * p2 - ak * p1;
    p1 = pt;
    a[k - 1] = ak;
    b[k - 1] = bk;
    ck += 2.f;
    fks = fks + fk + fk + 1.f;
    fhs = fhs + fk + fk;
    if (etest > fk * p1) {
	goto L140;
    }
    kk = k;
    s = 1.f;
    p1 = 0.f;
    p2 = 1.f;
    i__1 = k;
    for (i__ = 1; i__ <= i__1; ++i__) {
	pt = p2;
	p2 = (b[kk - 1] * p2 - p1) / a[kk - 1];
	p1 = pt;
	s += p2;
	--kk;
/* L150: */
    }
    s1 = coef * (p2 / s);
    if (inu > 0 || *n > 1) {
	goto L160;
    }
    goto L200;
L160:
    s2 = s1 * (*x + dnu + .5f - p1 / p2) / *x;

/*     FORWARD RECURSION ON THE THREE TERM RECURSION RELATION */

L170:
    ck = (dnu + dnu + 2.f) / *x;
    if (*n == 1) {
	--inu;
    }
    if (inu > 0) {
	goto L180;
    }
    if (*n > 1) {
	goto L200;
    }
    s1 = s2;
    goto L200;
L180:
    i__1 = inu;
    for (i__ = 1; i__ <= i__1; ++i__) {
	st = s2;
	s2 = ck * s2 + s1;
	s1 = st;
	ck += rx;
/* L190: */
    }
    if (*n == 1) {
	s1 = s2;
    }
L200:
    if (iflag == 1) {
	goto L220;
    }
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
	y[i__] = ck * y[i__ - 1] + y[i__ - 2];
	ck += rx;
/* L210: */
    }
    return 0;
/*     IFLAG=1 CASES */
L220:
    s = -(*x) + log(s1);
    y[1] = 0.f;
    *nz = 1;
    if (s < -elim) {
	goto L230;
    }
    y[1] = exp(s);
    *nz = 0;
L230:
    if (*n == 1) {
	return 0;
    }
    s = -(*x) + log(s2);
    y[2] = 0.f;
    ++(*nz);
    if (s < -elim) {
	goto L240;
    }
    --(*nz);
    y[2] = exp(s);
L240:
    if (*n == 2) {
	return 0;
    }
    kk = 2;
    if (*nz < 2) {
	goto L260;
    }
    i__1 = *n;
    for (i__ = 3; i__ <= i__1; ++i__) {
	kk = i__;
	st = s2;
	s2 = ck * s2 + s1;
	s1 = st;
	ck += rx;
	s = -(*x) + log(s2);
	++(*nz);
	y[i__] = 0.f;
	if (s < -elim) {
	    goto L250;
	}
	y[i__] = exp(s);
	--(*nz);
	goto L260;
L250:
	;
    }
    return 0;
L260:
    if (kk == *n) {
	return 0;
    }
    s2 = s2 * ck + s1;
    ck += rx;
    ++kk;
    y[kk] = exp(-(*x) + log(s2));
    if (kk == *n) {
	return 0;
    }
    ++kk;
    i__1 = *n;
    for (i__ = kk; i__ <= i__1; ++i__) {
	y[i__] = ck * y[i__ - 1] + y[i__ - 2];
	ck += rx;
/* L270: */
    }
    return 0;

/*     ASYMPTOTIC EXPANSION FOR LARGE X, X.GT.X2 */

/*     IFLAG=0 MEANS NO UNDERFLOW OCCURRED */
/*     IFLAG=1 MEANS AN UNDERFLOW OCCURRED- COMPUTATION PROCEEDS WITH */
/*     KODED=2 AND A TEST FOR ON SCALE VALUES IS MADE DURING FORWARD */
/*     RECURSION */
L280:
    nn = 2;
    if (inu == 0 && *n == 1) {
	nn = 1;
    }
    dnu2 = dnu + dnu;
    fmu = 0.f;
    if (dabs(dnu2) < tol) {
	goto L290;
    }
    fmu = dnu2 * dnu2;
L290:
    ex = *x * 8.f;
    s2 = 0.f;
    i__1 = nn;
    for (k = 1; k <= i__1; ++k) {
	s1 = s2;
	s = 1.f;
	ak = 0.f;
	ck = 1.f;
	sqk = 1.f;
	dk = ex;
	for (j = 1; j <= 30; ++j) {
	    ck = ck * (fmu - sqk) / dk;
	    s += ck;
	    dk += ex;
	    ak += 8.f;
	    sqk += ak;
	    if (dabs(ck) < tol) {
		goto L310;
	    }
/* L300: */
	}
L310:
	s2 = s * coef;
	fmu = fmu + dnu * 8.f + 4.f;
/* L320: */
    }
    if (nn > 1) {
	goto L170;
    }
    s1 = s2;
    goto L200;
L330:
    koded = 2;
    iflag = 1;
    goto L120;

/*     FNU=HALF ODD INTEGER CASE */

L340:
    s1 = coef;
    s2 = coef;
    goto L170;


L350:
    xermsg_("SLATEC", "BESKNU", "X NOT GREATER THAN ZERO", &c__2, &c__1, (
	    ftnlen)6, (ftnlen)6, (ftnlen)23);
    return 0;
L360:
    xermsg_("SLATEC", "BESKNU", "FNU NOT ZERO OR POSITIVE", &c__2, &c__1, (
	    ftnlen)6, (ftnlen)6, (ftnlen)24);
    return 0;
L370:
    xermsg_("SLATEC", "BESKNU", "KODE NOT 1 OR 2", &c__2, &c__1, (ftnlen)6, (
	    ftnlen)6, (ftnlen)15);
    return 0;
L380:
    xermsg_("SLATEC", "BESKNU", "N NOT GREATER THAN 0", &c__2, &c__1, (ftnlen)
	    6, (ftnlen)6, (ftnlen)20);
    return 0;
} /* besknu_ */

