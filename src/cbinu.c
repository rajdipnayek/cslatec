/* cbinu.f -- translated by f2c (version 12.02.01).
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
static integer c__2 = 2;

/* DECK CBINU */
/* Subroutine */ int cbinu_(complex *z__, real *fnu, integer *kode, integer *
	n, complex *cy, integer *nz, real *rl, real *fnul, real *tol, real *
	elim, real *alim)
{
    /* Initialized data */

    static complex czero = {0.f,0.f};

    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer i__;
    static complex cw[2];
    static real az;
    static integer nn, nw, nui, inw;
    static real dfnu;
    extern /* Subroutine */ int cbuni_(complex *, real *, integer *, integer *
	    , complex *, integer *, integer *, integer *, real *, real *, 
	    real *, real *), cseri_(complex *, real *, integer *, integer *, 
	    complex *, integer *, real *, real *, real *), cmlri_(complex *, 
	    real *, integer *, integer *, complex *, integer *, real *), 
	    casyi_(complex *, real *, integer *, integer *, complex *, 
	    integer *, real *, real *, real *, real *), cuoik_(complex *, 
	    real *, integer *, integer *, integer *, complex *, integer *, 
	    real *, real *, real *);
    static integer nlast;
    extern /* Subroutine */ int cwrsk_(complex *, real *, integer *, integer *
	    , complex *, integer *, complex *, real *, real *, real *);

/* ***BEGIN PROLOGUE  CBINU */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to CAIRY, CBESH, CBESI, CBESJ, CBESK and CBIRY */
/* ***LIBRARY   SLATEC */
/* ***TYPE      ALL (CBINU-A, ZBINU-A) */
/* ***AUTHOR  Amos, D. E., (SNL) */
/* ***DESCRIPTION */

/*     CBINU COMPUTES THE I FUNCTION IN THE RIGHT HALF Z PLANE */

/* ***SEE ALSO  CAIRY, CBESH, CBESI, CBESJ, CBESK, CBIRY */
/* ***ROUTINES CALLED  CASYI, CBUNI, CMLRI, CSERI, CUOIK, CWRSK */
/* ***REVISION HISTORY  (YYMMDD) */
/*   830501  DATE WRITTEN */
/*   910415  Prologue converted to Version 4.0 format.  (BAB) */
/* ***END PROLOGUE  CBINU */
    /* Parameter adjustments */
    --cy;

    /* Function Body */
/* ***FIRST EXECUTABLE STATEMENT  CBINU */
    *nz = 0;
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
/*     POWER SERIES */
/* ----------------------------------------------------------------------- */
    cseri_(z__, fnu, kode, &nn, &cy[1], &nw, tol, elim, alim);
    inw = abs(nw);
    *nz += inw;
    nn -= inw;
    if (nn == 0) {
	return 0;
    }
    if (nw >= 0) {
	goto L120;
    }
    dfnu = *fnu + (nn - 1);
L20:
    if (az < *rl) {
	goto L40;
    }
    if (dfnu <= 1.f) {
	goto L30;
    }
    if (az + az < dfnu * dfnu) {
	goto L50;
    }
/* ----------------------------------------------------------------------- */
/*     ASYMPTOTIC EXPANSION FOR LARGE Z */
/* ----------------------------------------------------------------------- */
L30:
    casyi_(z__, fnu, kode, &nn, &cy[1], &nw, rl, tol, elim, alim);
    if (nw < 0) {
	goto L130;
    }
    goto L120;
L40:
    if (dfnu <= 1.f) {
	goto L70;
    }
L50:
/* ----------------------------------------------------------------------- */
/*     OVERFLOW AND UNDERFLOW TEST ON I SEQUENCE FOR MILLER ALGORITHM */
/* ----------------------------------------------------------------------- */
    cuoik_(z__, fnu, kode, &c__1, &nn, &cy[1], &nw, tol, elim, alim);
    if (nw < 0) {
	goto L130;
    }
    *nz += nw;
    nn -= nw;
    if (nn == 0) {
	return 0;
    }
    dfnu = *fnu + (nn - 1);
    if (dfnu > *fnul) {
	goto L110;
    }
    if (az > *fnul) {
	goto L110;
    }
L60:
    if (az > *rl) {
	goto L80;
    }
L70:
/* ----------------------------------------------------------------------- */
/*     MILLER ALGORITHM NORMALIZED BY THE SERIES */
/* ----------------------------------------------------------------------- */
    cmlri_(z__, fnu, kode, &nn, &cy[1], &nw, tol);
    if (nw < 0) {
	goto L130;
    }
    goto L120;
L80:
/* ----------------------------------------------------------------------- */
/*     MILLER ALGORITHM NORMALIZED BY THE WRONSKIAN */
/* ----------------------------------------------------------------------- */
/* ----------------------------------------------------------------------- */
/*     OVERFLOW TEST ON K FUNCTIONS USED IN WRONSKIAN */
/* ----------------------------------------------------------------------- */
    cuoik_(z__, fnu, kode, &c__2, &c__2, cw, &nw, tol, elim, alim);
    if (nw >= 0) {
	goto L100;
    }
    *nz = nn;
    i__1 = nn;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = i__;
	cy[i__2].r = czero.r, cy[i__2].i = czero.i;
/* L90: */
    }
    return 0;
L100:
    if (nw > 0) {
	goto L130;
    }
    cwrsk_(z__, fnu, kode, &nn, &cy[1], &nw, cw, tol, elim, alim);
    if (nw < 0) {
	goto L130;
    }
    goto L120;
L110:
/* ----------------------------------------------------------------------- */
/*     INCREMENT FNU+NN-1 UP TO FNUL, COMPUTE AND RECUR BACKWARD */
/* ----------------------------------------------------------------------- */
    nui = *fnul - dfnu + 1;
    nui = max(nui,0);
    cbuni_(z__, fnu, kode, &nn, &cy[1], &nw, &nui, &nlast, fnul, tol, elim, 
	    alim);
    if (nw < 0) {
	goto L130;
    }
    *nz += nw;
    if (nlast == 0) {
	goto L120;
    }
    nn = nlast;
    goto L60;
L120:
    return 0;
L130:
    *nz = -1;
    if (nw == -2) {
	*nz = -2;
    }
    return 0;
} /* cbinu_ */

