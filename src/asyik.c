/* asyik.f -- translated by f2c (version 12.02.01).
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

/* DECK ASYIK */
/* Subroutine */ int asyik_(real *x, real *fnu, integer *kode, real *flgik, 
	real *ra, real *arg, integer *in, real *y)
{
    /* Initialized data */

    static real con[2] = { .398942280401432678f,1.25331413731550025f };
    static real c__[65] = { -.208333333333333f,.125f,.334201388888889f,
	    -.401041666666667f,.0703125f,-1.02581259645062f,1.84646267361111f,
	    -.8912109375f,.0732421875f,4.66958442342625f,-11.207002616223f,
	    8.78912353515625f,-2.3640869140625f,.112152099609375f,
	    -28.2120725582002f,84.6362176746007f,-91.81824154324f,
	    42.5349987453885f,-7.36879435947963f,.227108001708984f,
	    212.570130039217f,-765.252468141182f,1059.990452528f,
	    -699.579627376133f,218.190511744212f,-26.4914304869516f,
	    .572501420974731f,-1919.45766231841f,8061.72218173731f,
	    -13586.5500064341f,11655.3933368645f,-5305.6469786134f,
	    1200.90291321635f,-108.090919788395f,1.72772750258446f,
	    20204.2913309661f,-96980.5983886375f,192547.001232532f,
	    -203400.177280416f,122200.464983017f,-41192.6549688976f,
	    7109.51430248936f,-493.915304773088f,6.07404200127348f,
	    -242919.187900551f,1311763.61466298f,-2998015.91853811f,
	    3763271.2976564f,-2813563.22658653f,1268365.27332162f,
	    -331645.172484564f,45218.7689813627f,-2499.83048181121f,
	    24.3805296995561f,3284469.85307204f,-19706819.1184322f,
	    50952602.4926646f,-74105148.2115327f,66344512.274729f,
	    -37567176.6607634f,13288767.1664218f,-2785618.12808645f,
	    308186.404612662f,-13886.089753717f,110.017140269247f };

    /* System generated locals */
    integer i__1, i__2;
    real r__1, r__2;

    /* Local variables */
    static integer j, k, l;
    static real t, z__, s1, s2, t2, ak, ap, fn;
    static integer kk, jn;
    static real gln, tol, etx, coef;
    extern doublereal r1mach_(integer *);

/* ***BEGIN PROLOGUE  ASYIK */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to BESI and BESK */
/* ***LIBRARY   SLATEC */
/* ***TYPE      SINGLE PRECISION (ASYIK-S, DASYIK-D) */
/* ***AUTHOR  Amos, D. E., (SNLA) */
/* ***DESCRIPTION */

/*                    ASYIK computes Bessel functions I and K */
/*                  for arguments X.GT.0.0 and orders FNU.GE.35 */
/*                  on FLGIK = 1 and FLGIK = -1 respectively. */

/*                                    INPUT */

/*      X    - argument, X.GT.0.0E0 */
/*      FNU  - order of first Bessel function */
/*      KODE - a parameter to indicate the scaling option */
/*             KODE=1 returns Y(I)=        I/SUB(FNU+I-1)/(X), I=1,IN */
/*                    or      Y(I)=        K/SUB(FNU+I-1)/(X), I=1,IN */
/*                    on FLGIK = 1.0E0 or FLGIK = -1.0E0 */
/*             KODE=2 returns Y(I)=EXP(-X)*I/SUB(FNU+I-1)/(X), I=1,IN */
/*                    or      Y(I)=EXP( X)*K/SUB(FNU+I-1)/(X), I=1,IN */
/*                    on FLGIK = 1.0E0 or FLGIK = -1.0E0 */
/*     FLGIK - selection parameter for I or K function */
/*             FLGIK =  1.0E0 gives the I function */
/*             FLGIK = -1.0E0 gives the K function */
/*        RA - SQRT(1.+Z*Z), Z=X/FNU */
/*       ARG - argument of the leading exponential */
/*        IN - number of functions desired, IN=1 or 2 */

/*                                    OUTPUT */

/*         Y - a vector whose first in components contain the sequence */

/*     Abstract */
/*         ASYIK implements the uniform asymptotic expansion of */
/*         the I and K Bessel functions for FNU.GE.35 and real */
/*         X.GT.0.0E0. The forms are identical except for a change */
/*         in sign of some of the terms. This change in sign is */
/*         accomplished by means of the flag FLGIK = 1 or -1. */

/* ***SEE ALSO  BESI, BESK */
/* ***ROUTINES CALLED  R1MACH */
/* ***REVISION HISTORY  (YYMMDD) */
/*   750101  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900328  Added TYPE section.  (WRB) */
/*   910408  Updated the AUTHOR section.  (WRB) */
/* ***END PROLOGUE  ASYIK */

    /* Parameter adjustments */
    --y;

    /* Function Body */
/* ***FIRST EXECUTABLE STATEMENT  ASYIK */
    tol = r1mach_(&c__3);
    tol = dmax(tol,1e-15f);
    fn = *fnu;
    z__ = (3.f - *flgik) / 2.f;
    kk = (integer) z__;
    i__1 = *in;
    for (jn = 1; jn <= i__1; ++jn) {
	if (jn == 1) {
	    goto L10;
	}
	fn -= *flgik;
	z__ = *x / fn;
	*ra = sqrt(z__ * z__ + 1.f);
	gln = log((*ra + 1.f) / z__);
	etx = (real) (*kode - 1);
	t = *ra * (1.f - etx) + etx / (z__ + *ra);
	*arg = fn * (t - gln) * *flgik;
L10:
	coef = exp(*arg);
	t = 1.f / *ra;
	t2 = t * t;
	t /= fn;
	t = r_sign(&t, flgik);
	s2 = 1.f;
	ap = 1.f;
	l = 0;
	for (k = 2; k <= 11; ++k) {
	    ++l;
	    s1 = c__[l - 1];
	    i__2 = k;
	    for (j = 2; j <= i__2; ++j) {
		++l;
		s1 = s1 * t2 + c__[l - 1];
/* L20: */
	    }
	    ap *= t;
	    ak = ap * s1;
	    s2 += ak;
/* Computing MAX */
	    r__1 = dabs(ak), r__2 = dabs(ap);
	    if (dmax(r__1,r__2) < tol) {
		goto L40;
	    }
/* L30: */
	}
L40:
	t = dabs(t);
	y[jn] = s2 * coef * sqrt(t) * con[kk - 1];
/* L50: */
    }
    return 0;
} /* asyik_ */

