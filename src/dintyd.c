/* dintyd.f -- translated by f2c (version 12.02.01).
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

/* Common Block Declarations */

struct {
    doublereal rownd, rowns[210], el0, h__, hmin, hmxi, hu, tn, uround;
    integer iownd[14], iowns[6], ier, jstart, kflag, l, meth, miter, maxord, 
	    n, nq, nst, nfe, nje, nqu;
} ddebd1_;

#define ddebd1_1 ddebd1_

/* DECK DINTYD */
/* Subroutine */ int dintyd_(doublereal *t, integer *k, doublereal *yh, 
	integer *nyh, doublereal *dky, integer *iflag)
{
    /* System generated locals */
    integer yh_dim1, yh_offset, i__1, i__2;

    /* Local variables */
    static doublereal c__;
    static integer i__, j;
    static doublereal r__, s;
    static integer ic, jb, jj;
    static doublereal tp;
    static integer jb2, jj1, jp1;

/* ***BEGIN PROLOGUE  DINTYD */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to DDEBDF */
/* ***LIBRARY   SLATEC */
/* ***TYPE      DOUBLE PRECISION (INTYD-S, DINTYD-D) */
/* ***AUTHOR  Watts, H. A., (SNLA) */
/* ***DESCRIPTION */

/*   DINTYD approximates the solution and derivatives at T by polynomial */
/*   interpolation. Must be used in conjunction with the integrator */
/*   package DDEBDF. */
/* ---------------------------------------------------------------------- */
/* DINTYD computes interpolated values of the K-th derivative of the */
/* dependent variable vector Y, and stores it in DKY. */
/* This routine is called by DDEBDF with K = 0,1 and T = TOUT, but may */
/* also be called by the user for any K up to the current order. */
/* (see detailed instructions in LSODE usage documentation.) */
/* ---------------------------------------------------------------------- */
/* The computed values in DKY are gotten by interpolation using the */
/* Nordsieck history array YH.  This array corresponds uniquely to a */
/* vector-valued polynomial of degree NQCUR or less, and DKY is set */
/* to the K-th derivative of this polynomial at T. */
/* The formula for DKY is.. */
/*              Q */
/*  DKY(I)  =  Sum  C(J,K) * (T - TN)**(J-K) * H**(-J) * YH(I,J+1) */
/*             J=K */
/* where  C(J,K) = J*(J-1)*...*(J-K+1), Q = NQCUR, TN = TCUR, H = HCUR. */
/* The quantities  NQ = NQCUR, L = NQ+1, N = NEQ, TN, and H are */
/* communicated by common.  The above sum is done in reverse order. */
/* IFLAG is returned negative if either K or T is out of bounds. */
/* ---------------------------------------------------------------------- */

/* ***SEE ALSO  DDEBDF */
/* ***ROUTINES CALLED  (NONE) */
/* ***COMMON BLOCKS    DDEBD1 */
/* ***REVISION HISTORY  (YYMMDD) */
/*   820301  DATE WRITTEN */
/*   890911  Removed unnecessary intrinsics.  (WRB) */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900328  Added TYPE section.  (WRB) */
/*   910722  Updated AUTHOR section.  (ALS) */
/* ***END PROLOGUE  DINTYD */


/*     BEGIN BLOCK PERMITTING ...EXITS TO 130 */
/* ***FIRST EXECUTABLE STATEMENT  DINTYD */
    /* Parameter adjustments */
    yh_dim1 = *nyh;
    yh_offset = 1 + yh_dim1;
    yh -= yh_offset;
    --dky;

    /* Function Body */
    *iflag = 0;
    if (*k < 0 || *k > ddebd1_1.nq) {
	goto L110;
    }
    tp = ddebd1_1.tn - ddebd1_1.hu * (ddebd1_1.uround * 100. + 1.);
    if ((*t - tp) * (*t - ddebd1_1.tn) <= 0.) {
	goto L10;
    }
    *iflag = -2;
/*     .........EXIT */
    goto L130;
L10:

    s = (*t - ddebd1_1.tn) / ddebd1_1.h__;
    ic = 1;
    if (*k == 0) {
	goto L30;
    }
    jj1 = ddebd1_1.l - *k;
    i__1 = ddebd1_1.nq;
    for (jj = jj1; jj <= i__1; ++jj) {
	ic *= jj;
/* L20: */
    }
L30:
    c__ = (doublereal) ic;
    i__1 = ddebd1_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	dky[i__] = c__ * yh[i__ + ddebd1_1.l * yh_dim1];
/* L40: */
    }
    if (*k == ddebd1_1.nq) {
	goto L90;
    }
    jb2 = ddebd1_1.nq - *k;
    i__1 = jb2;
    for (jb = 1; jb <= i__1; ++jb) {
	j = ddebd1_1.nq - jb;
	jp1 = j + 1;
	ic = 1;
	if (*k == 0) {
	    goto L60;
	}
	jj1 = jp1 - *k;
	i__2 = j;
	for (jj = jj1; jj <= i__2; ++jj) {
	    ic *= jj;
/* L50: */
	}
L60:
	c__ = (doublereal) ic;
	i__2 = ddebd1_1.n;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    dky[i__] = c__ * yh[i__ + jp1 * yh_dim1] + s * dky[i__];
/* L70: */
	}
/* L80: */
    }
/*     .........EXIT */
    if (*k == 0) {
	goto L130;
    }
L90:
    i__1 = -(*k);
    r__ = pow_di(&ddebd1_1.h__, &i__1);
    i__1 = ddebd1_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	dky[i__] = r__ * dky[i__];
/* L100: */
    }
    goto L120;
L110:

    *iflag = -1;
L120:
L130:
    return 0;
/*     ----------------------- END OF SUBROUTINE DINTYD */
/*     ----------------------- */
} /* dintyd_ */

