/* intyd.f -- translated by f2c (version 12.02.01).
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
    real rownd, rowns[210], el0, h__, hmin, hmxi, hu, tn, uround;
    integer iownd[14], iowns[6], ier, jstart, kflag, l, meth, miter, maxord, 
	    n, nq, nst, nfe, nje, nqu;
} debdf1_;

#define debdf1_1 debdf1_

/* DECK INTYD */
/* Subroutine */ int intyd_(real *t, integer *k, real *yh, integer *nyh, real 
	*dky, integer *iflag)
{
    /* System generated locals */
    integer yh_dim1, yh_offset, i__1, i__2;

    /* Local variables */
    static real c__;
    static integer i__, j;
    static real r__, s;
    static integer ic, jb, jj;
    static real tp;
    static integer jb2, jj1, jp1;

/* ***BEGIN PROLOGUE  INTYD */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to DEBDF */
/* ***LIBRARY   SLATEC */
/* ***TYPE      SINGLE PRECISION (INTYD-S, DINTYD-D) */
/* ***AUTHOR  Watts, H. A., (SNLA) */
/* ***DESCRIPTION */

/*   INTYD approximates the solution and derivatives at T by polynomial */
/*   interpolation. Must be used in conjunction with the integrator */
/*   package DEBDF. */
/* ---------------------------------------------------------------------- */
/* INTYD computes interpolated values of the K-th derivative of the */
/* dependent variable vector Y, and stores it in DKY. */
/* This routine is called by DEBDF with K = 0,1 and T = TOUT, but may */
/* also be called by the user for any K up to the current order. */
/* (see detailed instructions in LSODE usage documentation.) */
/* ---------------------------------------------------------------------- */
/* The computed values in DKY are gotten by interpolation using the */
/* Nordsieck history array YH.  This array corresponds uniquely to a */
/* vector-valued polynomial of degree NQCUR or less, and DKY is set */
/* to the K-th derivative of this polynomial at T. */
/* The formula for DKY is.. */
/*              Q */
/*  DKY(I)  =  sum  C(J,K) * (T - TN)**(J-K) * H**(-J) * YH(I,J+1) */
/*             J=K */
/* where  C(J,K) = J*(J-1)*...*(J-K+1), Q = NQCUR, TN = TCUR, H = HCUR. */
/* The quantities  NQ = NQCUR, L = NQ+1, N = NEQ, TN, and H are */
/* communicated by common.  The above sum is done in reverse order. */
/* IFLAG is returned negative if either K or T is out of bounds. */
/* ---------------------------------------------------------------------- */

/* ***SEE ALSO  DEBDF */
/* ***ROUTINES CALLED  (NONE) */
/* ***COMMON BLOCKS    DEBDF1 */
/* ***REVISION HISTORY  (YYMMDD) */
/*   800901  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900328  Added TYPE section.  (WRB) */
/*   910722  Updated AUTHOR section.  (ALS) */
/* ***END PROLOGUE  INTYD */

/* LLL. OPTIMIZE */

/* ***FIRST EXECUTABLE STATEMENT  INTYD */
    /* Parameter adjustments */
    yh_dim1 = *nyh;
    yh_offset = 1 + yh_dim1;
    yh -= yh_offset;
    --dky;

    /* Function Body */
    *iflag = 0;
    if (*k < 0 || *k > debdf1_1.nq) {
	goto L80;
    }
    tp = debdf1_1.tn - debdf1_1.hu * (debdf1_1.uround * 100.f + 1.f);
    if ((*t - tp) * (*t - debdf1_1.tn) > 0.f) {
	goto L90;
    }

    s = (*t - debdf1_1.tn) / debdf1_1.h__;
    ic = 1;
    if (*k == 0) {
	goto L15;
    }
    jj1 = debdf1_1.l - *k;
    i__1 = debdf1_1.nq;
    for (jj = jj1; jj <= i__1; ++jj) {
/* L10: */
	ic *= jj;
    }
L15:
    c__ = (real) ic;
    i__1 = debdf1_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L20: */
	dky[i__] = c__ * yh[i__ + debdf1_1.l * yh_dim1];
    }
    if (*k == debdf1_1.nq) {
	goto L55;
    }
    jb2 = debdf1_1.nq - *k;
    i__1 = jb2;
    for (jb = 1; jb <= i__1; ++jb) {
	j = debdf1_1.nq - jb;
	jp1 = j + 1;
	ic = 1;
	if (*k == 0) {
	    goto L35;
	}
	jj1 = jp1 - *k;
	i__2 = j;
	for (jj = jj1; jj <= i__2; ++jj) {
/* L30: */
	    ic *= jj;
	}
L35:
	c__ = (real) ic;
	i__2 = debdf1_1.n;
	for (i__ = 1; i__ <= i__2; ++i__) {
/* L40: */
	    dky[i__] = c__ * yh[i__ + jp1 * yh_dim1] + s * dky[i__];
	}
/* L50: */
    }
    if (*k == 0) {
	return 0;
    }
L55:
    i__1 = -(*k);
    r__ = pow_ri(&debdf1_1.h__, &i__1);
    i__1 = debdf1_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L60: */
	dky[i__] = r__ * dky[i__];
    }
    return 0;

L80:
    *iflag = -1;
    return 0;
L90:
    *iflag = -2;
    return 0;
/* ----------------------- END OF SUBROUTINE INTYD ----------------------- */
} /* intyd_ */

