/* hkseq.f -- translated by f2c (version 12.02.01).
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
static integer c__5 = 5;
static integer c__11 = 11;

/* DECK HKSEQ */
/* Subroutine */ int hkseq_(real *x, integer *m, real *h__, integer *ierr)
{
    /* Initialized data */

    static real b[22] = { 1.f,-.5f,.25f,-.0625f,.046875f,-.06640625f,
	    .1513671875f,-.506103515625f,2.33319091796875f,
	    -14.1840972900390625f,109.941936492919922f,-1058.24747562408447f,
	    12384.2434241771698f,-173160.495905935764f,2851034.29084961116f,
	    -54596461.9322445132f,1203161746.68075304f,-30232631527.1452307f,
	    859229286072.319606f,-27423310409777.6039f,976664637943633.248f,
	    -38593158683845036.f };

    /* System generated locals */
    integer i__1, i__2;
    real r__1;

    /* Local variables */
    static integer i__, j, k;
    static real s, t, u[25], v[25], fk, fn, tk, xh;
    static integer mx, nx;
    static real xm, fln, fnp, r1m5, rln, hrx, trm[22], tst, xinc, trmh[25], 
	    xmin, xdmy, yint, trmr[25], rxsq, slope, wdtol;
    extern integer i1mach_(integer *);
    extern doublereal r1mach_(integer *);

/* ***BEGIN PROLOGUE  HKSEQ */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to BSKIN */
/* ***LIBRARY   SLATEC */
/* ***TYPE      SINGLE PRECISION (HKSEQ-S, DHKSEQ-D) */
/* ***AUTHOR  Amos, D. E., (SNLA) */
/* ***DESCRIPTION */

/*   HKSEQ is an adaptation of subroutine PSIFN described in the */
/*   reference below.  HKSEQ generates the sequence */
/*   H(K,X) = (-X)**(K+1)*(PSI(K,X) PSI(K,X+0.5))/GAMMA(K+1), for */
/*            K=0,...,M. */

/* ***SEE ALSO  BSKIN */
/* ***REFERENCES  D. E. Amos, A portable Fortran subroutine for */
/*                 derivatives of the Psi function, Algorithm 610, ACM */
/*                 Transactions on Mathematical Software 9, 4 (1983), */
/*                 pp. 494-502. */
/* ***ROUTINES CALLED  I1MACH, R1MACH */
/* ***REVISION HISTORY  (YYMMDD) */
/*   820601  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900328  Added TYPE section.  (WRB) */
/*   910722  Updated AUTHOR section.  (ALS) */
/*   920528  DESCRIPTION and REFERENCES sections revised.  (WRB) */
/* ***END PROLOGUE  HKSEQ */
/* ----------------------------------------------------------------------- */
/*             SCALED BERNOULLI NUMBERS 2.0*B(2K)*(1-2**(-2K)) */
/* ----------------------------------------------------------------------- */
    /* Parameter adjustments */
    --h__;

    /* Function Body */

/* ***FIRST EXECUTABLE STATEMENT  HKSEQ */
    *ierr = 0;
/* Computing MAX */
    r__1 = r1mach_(&c__4);
    wdtol = dmax(r__1,1e-18f);
    fn = (real) (*m - 1);
    fnp = fn + 1.f;
/* ----------------------------------------------------------------------- */
/*     COMPUTE XMIN */
/* ----------------------------------------------------------------------- */
    r1m5 = r1mach_(&c__5);
    rln = r1m5 * i1mach_(&c__11);
    rln = dmin(rln,18.06f);
    fln = dmax(rln,3.f) - 3.f;
    yint = fln * .4f + 3.5f;
    slope = fln * (fln * 6.038e-4f + .008677f) + .21f;
    xm = yint + slope * fn;
    mx = (integer) xm + 1;
    xmin = (real) mx;
/* ----------------------------------------------------------------------- */
/*     GENERATE H(M-1,XDMY)*XDMY**(M) BY THE ASYMPTOTIC EXPANSION */
/* ----------------------------------------------------------------------- */
    xdmy = *x;
    xinc = 0.f;
    if (*x >= xmin) {
	goto L10;
    }
    nx = (integer) (*x);
    xinc = xmin - nx;
    xdmy = *x + xinc;
L10:
    rxsq = 1.f / (xdmy * xdmy);
    hrx = .5f / xdmy;
    tst = wdtol * .5f;
    t = fnp * hrx;
/* ----------------------------------------------------------------------- */
/*     INITIALIZE COEFFICIENT ARRAY */
/* ----------------------------------------------------------------------- */
    s = t * b[2];
    if (dabs(s) < tst) {
	goto L30;
    }
    tk = 2.f;
    for (k = 4; k <= 22; ++k) {
	t = t * ((tk + fn + 1.f) / (tk + 1.f)) * ((tk + fn) / (tk + 2.f)) * 
		rxsq;
	trm[k - 1] = t * b[k - 1];
	if ((r__1 = trm[k - 1], dabs(r__1)) < tst) {
	    goto L30;
	}
	s += trm[k - 1];
	tk += 2.f;
/* L20: */
    }
    goto L110;
L30:
    h__[*m] = s + .5f;
    if (*m == 1) {
	goto L70;
    }
/* ----------------------------------------------------------------------- */
/*     GENERATE LOWER DERIVATIVES, I.LT.M-1 */
/* ----------------------------------------------------------------------- */
    i__1 = *m;
    for (i__ = 2; i__ <= i__1; ++i__) {
	fnp = fn;
	fn += -1.f;
	s = fnp * hrx * b[2];
	if (dabs(s) < tst) {
	    goto L50;
	}
	fk = fnp + 3.f;
	for (k = 4; k <= 22; ++k) {
	    trm[k - 1] = trm[k - 1] * fnp / fk;
	    if ((r__1 = trm[k - 1], dabs(r__1)) < tst) {
		goto L50;
	    }
	    s += trm[k - 1];
	    fk += 2.f;
/* L40: */
	}
	goto L110;
L50:
	mx = *m - i__ + 1;
	h__[mx] = s + .5f;
/* L60: */
    }
L70:
    if (xinc == 0.f) {
	return 0;
    }
/* ----------------------------------------------------------------------- */
/*     RECUR BACKWARD FROM XDMY TO X */
/* ----------------------------------------------------------------------- */
    xh = *x + .5f;
    s = 0.f;
    nx = (integer) xinc;
    i__1 = nx;
    for (i__ = 1; i__ <= i__1; ++i__) {
	trmr[i__ - 1] = *x / (*x + nx - i__);
	u[i__ - 1] = trmr[i__ - 1];
	trmh[i__ - 1] = *x / (xh + nx - i__);
	v[i__ - 1] = trmh[i__ - 1];
	s = s + u[i__ - 1] - v[i__ - 1];
/* L80: */
    }
    mx = nx + 1;
    trmr[mx - 1] = *x / xdmy;
    u[mx - 1] = trmr[mx - 1];
    h__[1] = h__[1] * trmr[mx - 1] + s;
    if (*m == 1) {
	return 0;
    }
    i__1 = *m;
    for (j = 2; j <= i__1; ++j) {
	s = 0.f;
	i__2 = nx;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    trmr[i__ - 1] *= u[i__ - 1];
	    trmh[i__ - 1] *= v[i__ - 1];
	    s = s + trmr[i__ - 1] - trmh[i__ - 1];
/* L90: */
	}
	trmr[mx - 1] *= u[mx - 1];
	h__[j] = h__[j] * trmr[mx - 1] + s;
/* L100: */
    }
    return 0;
L110:
    *ierr = 2;
    return 0;
} /* hkseq_ */

