/* cacai.f -- translated by f2c (version 12.02.01).
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

/* DECK CACAI */
/* Subroutine */ int cacai_(complex *z__, real *fnu, integer *kode, integer *
	mr, integer *n, complex *y, integer *nz, real *rl, real *tol, real *
	elim, real *alim)
{
    /* Initialized data */

    static real pi = 3.14159265358979324f;

    /* System generated locals */
    complex q__1, q__2, q__3;

    /* Local variables */
    static complex c1, c2;
    static real az;
    static complex cy[2];
    static integer nn, nw;
    static complex zn;
    static real yy, arg, cpn;
    static integer iuf;
    static real fmr, sgn;
    static integer inu;
    static real spn;
    static complex csgn;
    extern /* Subroutine */ int cs1s2_(complex *, complex *, complex *, 
	    integer *, real *, real *, integer *);
    static real dfnu;
    static complex cspn;
    static real ascle;
    extern /* Subroutine */ int cbknu_(complex *, real *, integer *, integer *
	    , complex *, integer *, real *, real *, real *), cseri_(complex *,
	     real *, integer *, integer *, complex *, integer *, real *, real 
	    *, real *), cmlri_(complex *, real *, integer *, integer *, 
	    complex *, integer *, real *), casyi_(complex *, real *, integer *
	    , integer *, complex *, integer *, real *, real *, real *, real *)
	    ;
    extern doublereal r1mach_(integer *);

/* ***BEGIN PROLOGUE  CACAI */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to CAIRY */
/* ***LIBRARY   SLATEC */
/* ***TYPE      ALL (CACAI-A, ZACAI-A) */
/* ***AUTHOR  Amos, D. E., (SNL) */
/* ***DESCRIPTION */

/*     CACAI APPLIES THE ANALYTIC CONTINUATION FORMULA */

/*         K(FNU,ZN*EXP(MP))=K(FNU,ZN)*EXP(-MP*FNU) - MP*I(FNU,ZN) */
/*                 MP=PI*MR*CMPLX(0.0,1.0) */

/*     TO CONTINUE THE K FUNCTION FROM THE RIGHT HALF TO THE LEFT */
/*     HALF Z PLANE FOR USE WITH CAIRY WHERE FNU=1/3 OR 2/3 AND N=1. */
/*     CACAI IS THE SAME AS CACON WITH THE PARTS FOR LARGER ORDERS AND */
/*     RECURRENCE REMOVED. A RECURSIVE CALL TO CACON CAN RESULT IF CACON */
/*     IS CALLED FROM CAIRY. */

/* ***SEE ALSO  CAIRY */
/* ***ROUTINES CALLED  CASYI, CBKNU, CMLRI, CS1S2, CSERI, R1MACH */
/* ***REVISION HISTORY  (YYMMDD) */
/*   830501  DATE WRITTEN */
/*   910415  Prologue converted to Version 4.0 format.  (BAB) */
/* ***END PROLOGUE  CACAI */
    /* Parameter adjustments */
    --y;

    /* Function Body */
/* ***FIRST EXECUTABLE STATEMENT  CACAI */
    *nz = 0;
    q__1.r = -z__->r, q__1.i = -z__->i;
    zn.r = q__1.r, zn.i = q__1.i;
    az = c_abs(z__);
    nn = *n;
    dfnu = *fnu + (*n - 1);
    if (az <= 2.f) {
	goto L10;
    }
    if (az * az * .25f > dfnu + 1.f) {
	goto L20;
    }
L10:
/* ----------------------------------------------------------------------- */
/*     POWER SERIES FOR THE I FUNCTION */
/* ----------------------------------------------------------------------- */
    cseri_(&zn, fnu, kode, &nn, &y[1], &nw, tol, elim, alim);
    goto L40;
L20:
    if (az < *rl) {
	goto L30;
    }
/* ----------------------------------------------------------------------- */
/*     ASYMPTOTIC EXPANSION FOR LARGE Z FOR THE I FUNCTION */
/* ----------------------------------------------------------------------- */
    casyi_(&zn, fnu, kode, &nn, &y[1], &nw, rl, tol, elim, alim);
    if (nw < 0) {
	goto L70;
    }
    goto L40;
L30:
/* ----------------------------------------------------------------------- */
/*     MILLER ALGORITHM NORMALIZED BY THE SERIES FOR THE I FUNCTION */
/* ----------------------------------------------------------------------- */
    cmlri_(&zn, fnu, kode, &nn, &y[1], &nw, tol);
    if (nw < 0) {
	goto L70;
    }
L40:
/* ----------------------------------------------------------------------- */
/*     ANALYTIC CONTINUATION TO THE LEFT HALF PLANE FOR THE K FUNCTION */
/* ----------------------------------------------------------------------- */
    cbknu_(&zn, fnu, kode, &c__1, cy, &nw, tol, elim, alim);
    if (nw != 0) {
	goto L70;
    }
    fmr = (real) (*mr);
    sgn = -r_sign(&pi, &fmr);
    q__1.r = 0.f, q__1.i = sgn;
    csgn.r = q__1.r, csgn.i = q__1.i;
    if (*kode == 1) {
	goto L50;
    }
    yy = -r_imag(&zn);
    cpn = cos(yy);
    spn = sin(yy);
    q__2.r = cpn, q__2.i = spn;
    q__1.r = csgn.r * q__2.r - csgn.i * q__2.i, q__1.i = csgn.r * q__2.i + 
	    csgn.i * q__2.r;
    csgn.r = q__1.r, csgn.i = q__1.i;
L50:
/* ----------------------------------------------------------------------- */
/*     CALCULATE CSPN=EXP(FNU*PI*I) TO MINIMIZE LOSSES OF SIGNIFICANCE */
/*     WHEN FNU IS LARGE */
/* ----------------------------------------------------------------------- */
    inu = *fnu;
    arg = (*fnu - inu) * sgn;
    cpn = cos(arg);
    spn = sin(arg);
    q__1.r = cpn, q__1.i = spn;
    cspn.r = q__1.r, cspn.i = q__1.i;
    if (inu % 2 == 1) {
	q__1.r = -cspn.r, q__1.i = -cspn.i;
	cspn.r = q__1.r, cspn.i = q__1.i;
    }
    c1.r = cy[0].r, c1.i = cy[0].i;
    c2.r = y[1].r, c2.i = y[1].i;
    if (*kode == 1) {
	goto L60;
    }
    iuf = 0;
    ascle = r1mach_(&c__1) * 1e3f / *tol;
    cs1s2_(&zn, &c1, &c2, &nw, &ascle, alim, &iuf);
    *nz += nw;
L60:
    q__2.r = cspn.r * c1.r - cspn.i * c1.i, q__2.i = cspn.r * c1.i + cspn.i * 
	    c1.r;
    q__3.r = csgn.r * c2.r - csgn.i * c2.i, q__3.i = csgn.r * c2.i + csgn.i * 
	    c2.r;
    q__1.r = q__2.r + q__3.r, q__1.i = q__2.i + q__3.i;
    y[1].r = q__1.r, y[1].i = q__1.i;
    return 0;
L70:
    *nz = -1;
    if (nw == -2) {
	*nz = -2;
    }
    return 0;
} /* cacai_ */

