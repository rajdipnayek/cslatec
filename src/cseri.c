/* cseri.f -- translated by f2c (version 12.02.01).
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

static integer c__1 = 1;

/* DECK CSERI */
/* Subroutine */ int cseri_(complex *z__, real *fnu, integer *kode, integer *
	n, complex *y, integer *nz, real *tol, real *elim, real *alim)
{
    /* Initialized data */

    static complex czero = {0.f,0.f};
    static complex cone = {1.f,0.f};

    /* System generated locals */
    integer i__1, i__2, i__3, i__4;
    real r__1, r__2;
    complex q__1, q__2, q__3, q__4;

    /* Local variables */
    static integer i__, k, l, m;
    static real s;
    static complex w[2];
    static real x;
    static complex s1, s2;
    static real aa;
    static integer ib;
    static real ak;
    static complex ck;
    static integer il;
    static real az;
    static integer nn;
    static complex cz, hz;
    static real rs, ss;
    static integer nw;
    static complex rz, ak1;
    static real acz, arm, rak1, rtr1;
    static complex coef, crsc;
    static real dfnu;
    static integer idum;
    static real atol, fnup;
    static integer iflag;
    static real ascle;
    extern /* Subroutine */ int cuchk_(complex *, integer *, real *, real *);
    extern doublereal gamln_(real *, integer *), r1mach_(integer *);

/* ***BEGIN PROLOGUE  CSERI */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to CBESI and CBESK */
/* ***LIBRARY   SLATEC */
/* ***TYPE      ALL (CSERI-A, ZSERI-A) */
/* ***AUTHOR  Amos, D. E., (SNL) */
/* ***DESCRIPTION */

/*     CSERI COMPUTES THE I BESSEL FUNCTION FOR REAL(Z).GE.0.0 BY */
/*     MEANS OF THE POWER SERIES FOR LARGE ABS(Z) IN THE */
/*     REGION ABS(Z).LE.2*SQRT(FNU+1). NZ=0 IS A NORMAL RETURN. */
/*     NZ.GT.0 MEANS THAT THE LAST NZ COMPONENTS WERE SET TO ZERO */
/*     DUE TO UNDERFLOW. NZ.LT.0 MEANS UNDERFLOW OCCURRED, BUT THE */
/*     CONDITION ABS(Z).LE.2*SQRT(FNU+1) WAS VIOLATED AND THE */
/*     COMPUTATION MUST BE COMPLETED IN ANOTHER ROUTINE WITH N=N-ABS(NZ). */

/* ***SEE ALSO  CBESI, CBESK */
/* ***ROUTINES CALLED  CUCHK, GAMLN, R1MACH */
/* ***REVISION HISTORY  (YYMMDD) */
/*   830501  DATE WRITTEN */
/*   910415  Prologue converted to Version 4.0 format.  (BAB) */
/* ***END PROLOGUE  CSERI */
    /* Parameter adjustments */
    --y;

    /* Function Body */
/* ***FIRST EXECUTABLE STATEMENT  CSERI */
    *nz = 0;
    az = c_abs(z__);
    if (az == 0.f) {
	goto L150;
    }
    x = z__->r;
    arm = r1mach_(&c__1) * 1e3f;
    rtr1 = sqrt(arm);
    crsc.r = 1.f, crsc.i = 0.f;
    iflag = 0;
    if (az < arm) {
	goto L140;
    }
    q__1.r = z__->r * .5f - z__->i * 0.f, q__1.i = z__->r * 0.f + z__->i * 
	    .5f;
    hz.r = q__1.r, hz.i = q__1.i;
    cz.r = czero.r, cz.i = czero.i;
    if (az > rtr1) {
	q__1.r = hz.r * hz.r - hz.i * hz.i, q__1.i = hz.r * hz.i + hz.i * 
		hz.r;
	cz.r = q__1.r, cz.i = q__1.i;
    }
    acz = c_abs(&cz);
    nn = *n;
    c_log(&q__1, &hz);
    ck.r = q__1.r, ck.i = q__1.i;
L10:
    dfnu = *fnu + (nn - 1);
    fnup = dfnu + 1.f;
/* ----------------------------------------------------------------------- */
/*     UNDERFLOW TEST */
/* ----------------------------------------------------------------------- */
    q__2.r = dfnu, q__2.i = 0.f;
    q__1.r = ck.r * q__2.r - ck.i * q__2.i, q__1.i = ck.r * q__2.i + ck.i * 
	    q__2.r;
    ak1.r = q__1.r, ak1.i = q__1.i;
    ak = gamln_(&fnup, &idum);
    q__2.r = ak, q__2.i = 0.f;
    q__1.r = ak1.r - q__2.r, q__1.i = ak1.i - q__2.i;
    ak1.r = q__1.r, ak1.i = q__1.i;
    if (*kode == 2) {
	q__2.r = x, q__2.i = 0.f;
	q__1.r = ak1.r - q__2.r, q__1.i = ak1.i - q__2.i;
	ak1.r = q__1.r, ak1.i = q__1.i;
    }
    rak1 = ak1.r;
    if (rak1 > -(*elim)) {
	goto L30;
    }
L20:
    ++(*nz);
    i__1 = nn;
    y[i__1].r = czero.r, y[i__1].i = czero.i;
    if (acz > dfnu) {
	goto L170;
    }
    --nn;
    if (nn == 0) {
	return 0;
    }
    goto L10;
L30:
    if (rak1 > -(*alim)) {
	goto L40;
    }
    iflag = 1;
    ss = 1.f / *tol;
    q__1.r = *tol, q__1.i = 0.f;
    crsc.r = q__1.r, crsc.i = q__1.i;
    ascle = arm * ss;
L40:
    ak = r_imag(&ak1);
    aa = exp(rak1);
    if (iflag == 1) {
	aa *= ss;
    }
    q__2.r = aa, q__2.i = 0.f;
    r__1 = cos(ak);
    r__2 = sin(ak);
    q__3.r = r__1, q__3.i = r__2;
    q__1.r = q__2.r * q__3.r - q__2.i * q__3.i, q__1.i = q__2.r * q__3.i + 
	    q__2.i * q__3.r;
    coef.r = q__1.r, coef.i = q__1.i;
    atol = *tol * acz / fnup;
    il = min(2,nn);
    i__1 = il;
    for (i__ = 1; i__ <= i__1; ++i__) {
	dfnu = *fnu + (nn - i__);
	fnup = dfnu + 1.f;
	s1.r = cone.r, s1.i = cone.i;
	if (acz < *tol * fnup) {
	    goto L60;
	}
	ak1.r = cone.r, ak1.i = cone.i;
	ak = fnup + 2.f;
	s = fnup;
	aa = 2.f;
L50:
	rs = 1.f / s;
	q__2.r = ak1.r * cz.r - ak1.i * cz.i, q__2.i = ak1.r * cz.i + ak1.i * 
		cz.r;
	q__3.r = rs, q__3.i = 0.f;
	q__1.r = q__2.r * q__3.r - q__2.i * q__3.i, q__1.i = q__2.r * q__3.i 
		+ q__2.i * q__3.r;
	ak1.r = q__1.r, ak1.i = q__1.i;
	q__1.r = s1.r + ak1.r, q__1.i = s1.i + ak1.i;
	s1.r = q__1.r, s1.i = q__1.i;
	s += ak;
	ak += 2.f;
	aa = aa * acz * rs;
	if (aa > atol) {
	    goto L50;
	}
L60:
	m = nn - i__ + 1;
	q__1.r = s1.r * coef.r - s1.i * coef.i, q__1.i = s1.r * coef.i + s1.i 
		* coef.r;
	s2.r = q__1.r, s2.i = q__1.i;
	i__2 = i__ - 1;
	w[i__2].r = s2.r, w[i__2].i = s2.i;
	if (iflag == 0) {
	    goto L70;
	}
	cuchk_(&s2, &nw, &ascle, tol);
	if (nw != 0) {
	    goto L20;
	}
L70:
	i__2 = m;
	q__1.r = s2.r * crsc.r - s2.i * crsc.i, q__1.i = s2.r * crsc.i + s2.i 
		* crsc.r;
	y[i__2].r = q__1.r, y[i__2].i = q__1.i;
	if (i__ != il) {
	    q__3.r = dfnu, q__3.i = 0.f;
	    q__2.r = coef.r * q__3.r - coef.i * q__3.i, q__2.i = coef.r * 
		    q__3.i + coef.i * q__3.r;
	    c_div(&q__1, &q__2, &hz);
	    coef.r = q__1.r, coef.i = q__1.i;
	}
/* L80: */
    }
    if (nn <= 2) {
	return 0;
    }
    k = nn - 2;
    ak = (real) k;
    q__2.r = cone.r + cone.r, q__2.i = cone.i + cone.i;
    c_div(&q__1, &q__2, z__);
    rz.r = q__1.r, rz.i = q__1.i;
    if (iflag == 1) {
	goto L110;
    }
    ib = 3;
L90:
    i__1 = nn;
    for (i__ = ib; i__ <= i__1; ++i__) {
	i__2 = k;
	r__1 = ak + *fnu;
	q__4.r = r__1, q__4.i = 0.f;
	q__3.r = q__4.r * rz.r - q__4.i * rz.i, q__3.i = q__4.r * rz.i + 
		q__4.i * rz.r;
	i__3 = k + 1;
	q__2.r = q__3.r * y[i__3].r - q__3.i * y[i__3].i, q__2.i = q__3.r * y[
		i__3].i + q__3.i * y[i__3].r;
	i__4 = k + 2;
	q__1.r = q__2.r + y[i__4].r, q__1.i = q__2.i + y[i__4].i;
	y[i__2].r = q__1.r, y[i__2].i = q__1.i;
	ak += -1.f;
	--k;
/* L100: */
    }
    return 0;
/* ----------------------------------------------------------------------- */
/*     RECUR BACKWARD WITH SCALED VALUES */
/* ----------------------------------------------------------------------- */
L110:
/* ----------------------------------------------------------------------- */
/*     EXP(-ALIM)=EXP(-ELIM)/TOL=APPROX. ONE PRECISION ABOVE THE */
/*     UNDERFLOW LIMIT = ASCLE = R1MACH(1)*CSCL*1.0E+3 */
/* ----------------------------------------------------------------------- */
    s1.r = w[0].r, s1.i = w[0].i;
    s2.r = w[1].r, s2.i = w[1].i;
    i__1 = nn;
    for (l = 3; l <= i__1; ++l) {
	ck.r = s2.r, ck.i = s2.i;
	r__1 = ak + *fnu;
	q__4.r = r__1, q__4.i = 0.f;
	q__3.r = q__4.r * rz.r - q__4.i * rz.i, q__3.i = q__4.r * rz.i + 
		q__4.i * rz.r;
	q__2.r = q__3.r * s2.r - q__3.i * s2.i, q__2.i = q__3.r * s2.i + 
		q__3.i * s2.r;
	q__1.r = s1.r + q__2.r, q__1.i = s1.i + q__2.i;
	s2.r = q__1.r, s2.i = q__1.i;
	s1.r = ck.r, s1.i = ck.i;
	q__1.r = s2.r * crsc.r - s2.i * crsc.i, q__1.i = s2.r * crsc.i + s2.i 
		* crsc.r;
	ck.r = q__1.r, ck.i = q__1.i;
	i__2 = k;
	y[i__2].r = ck.r, y[i__2].i = ck.i;
	ak += -1.f;
	--k;
	if (c_abs(&ck) > ascle) {
	    goto L130;
	}
/* L120: */
    }
    return 0;
L130:
    ib = l + 1;
    if (ib > nn) {
	return 0;
    }
    goto L90;
L140:
    *nz = *n;
    if (*fnu == 0.f) {
	--(*nz);
    }
L150:
    y[1].r = czero.r, y[1].i = czero.i;
    if (*fnu == 0.f) {
	y[1].r = cone.r, y[1].i = cone.i;
    }
    if (*n == 1) {
	return 0;
    }
    i__1 = *n;
    for (i__ = 2; i__ <= i__1; ++i__) {
	i__2 = i__;
	y[i__2].r = czero.r, y[i__2].i = czero.i;
/* L160: */
    }
    return 0;
/* ----------------------------------------------------------------------- */
/*     RETURN WITH NZ.LT.0 IF ABS(Z*Z/4).GT.FNU+N-NZ-1 COMPLETE */
/*     THE CALCULATION IN CBINU WITH N=N-ABS(NZ) */
/* ----------------------------------------------------------------------- */
L170:
    *nz = -(*nz);
    return 0;
} /* cseri_ */

