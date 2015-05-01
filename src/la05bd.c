/* la05bd.f -- translated by f2c (version 12.02.01).
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
    doublereal small;
    integer lp, lenl, lenu, ncp, lrow, lcol;
} la05dd_;

#define la05dd_1 la05dd_

/* Table of constant values */

static integer c_n8 = -8;
static integer c__2 = 2;

/* DECK LA05BD */
/* Subroutine */ int la05bd_(doublereal *a, integer *ind, integer *ia, 
	integer *n, integer *ip, integer *iw, doublereal *w, doublereal *g, 
	doublereal *b, logical *trans)
{
    /* System generated locals */
    integer ind_dim1, ind_offset, iw_dim1, iw_offset, ip_dim1, ip_offset, 
	    i__1, i__2, i__3;

    /* Local variables */
    static integer i__, j, k, l1, k2, n1;
    static doublereal am;
    static integer ii, kk, kl, kp, nz, kpc, kll;
    extern /* Subroutine */ int xermsg_(char *, char *, char *, integer *, 
	    integer *, ftnlen, ftnlen, ftnlen), xsetun_(integer *);

/* ***BEGIN PROLOGUE  LA05BD */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to DSPLP */
/* ***LIBRARY   SLATEC */
/* ***TYPE      DOUBLE PRECISION (LA05BS-S, LA05BD-D) */
/* ***AUTHOR  (UNKNOWN) */
/* ***DESCRIPTION */

/*     THIS SUBPROGRAM IS A SLIGHT MODIFICATION OF A SUBPROGRAM */
/*     FROM THE C. 1979 AERE HARWELL LIBRARY.  THE NAME OF THE */
/*     CORRESPONDING HARWELL CODE CAN BE OBTAINED BY DELETING */
/*     THE FINAL LETTER =D= IN THE NAMES USED HERE. */
/*     REVISED SEP. 13, 1979. */

/*     ROYALTIES HAVE BEEN PAID TO AERE-UK FOR USE OF THEIR CODES */
/*     IN THE PACKAGE GIVEN HERE.  ANY PRIMARY USAGE OF THE HARWELL */
/*     SUBROUTINES REQUIRES A ROYALTY AGREEMENT AND PAYMENT BETWEEN */
/*     THE USER AND AERE-UK.  ANY USAGE OF THE SANDIA WRITTEN CODES */
/*     DSPLP( ) (WHICH USES THE HARWELL SUBROUTINES) IS PERMITTED. */

/* IP(I,1),IP(I,2) POINT TO START OF ROW/COLUMN I OF U. */
/* IW(I,1),IW(I,2) ARE LENGTHS OF ROW/COL I OF U. */
/* IW(.,3),IW(.,4) HOLD ROW/COL NUMBERS IN PIVOTAL ORDER. */

/* ***SEE ALSO  DSPLP */
/* ***ROUTINES CALLED  XERMSG, XSETUN */
/* ***COMMON BLOCKS    LA05DD */
/* ***REVISION HISTORY  (YYMMDD) */
/*   811215  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890831  Modified array declarations.  (WRB) */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ) */
/*   900402  Added TYPE section.  (WRB) */
/*   920410  Corrected second dimension on IW declaration.  (WRB) */
/* ***END PROLOGUE  LA05BD */
/* ***FIRST EXECUTABLE STATEMENT  LA05BD */
    /* Parameter adjustments */
    --a;
    ind_dim1 = *ia;
    ind_offset = 1 + ind_dim1;
    ind -= ind_offset;
    iw_dim1 = *n;
    iw_offset = 1 + iw_dim1;
    iw -= iw_offset;
    ip_dim1 = *n;
    ip_offset = 1 + ip_dim1;
    ip -= ip_offset;
    --w;
    --b;

    /* Function Body */
    if (*g < 0.) {
	goto L130;
    }
    kll = *ia - la05dd_1.lenl + 1;
    if (*trans) {
	goto L80;
    }

/*     MULTIPLY VECTOR BY INVERSE OF L */
    if (la05dd_1.lenl <= 0) {
	goto L20;
    }
    l1 = *ia + 1;
    i__1 = la05dd_1.lenl;
    for (kk = 1; kk <= i__1; ++kk) {
	k = l1 - kk;
	i__ = ind[k + ind_dim1];
	if (b[i__] == 0.) {
	    goto L10;
	}
	j = ind[k + (ind_dim1 << 1)];
	b[j] += a[k] * b[i__];
L10:
	;
    }
L20:
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	w[i__] = b[i__];
	b[i__] = 0.;
/* L30: */
    }

/*     MULTIPLY VECTOR BY INVERSE OF U */
    n1 = *n + 1;
    i__1 = *n;
    for (ii = 1; ii <= i__1; ++ii) {
	i__ = n1 - ii;
	i__ = iw[i__ + iw_dim1 * 3];
	am = w[i__];
	kp = ip[i__ + ip_dim1];
	if (kp > 0) {
	    goto L50;
	}
	kp = -kp;
	ip[i__ + ip_dim1] = kp;
	nz = iw[i__ + iw_dim1];
	kl = kp - 1 + nz;
	k2 = kp + 1;
	i__2 = kl;
	for (k = k2; k <= i__2; ++k) {
	    j = ind[k + (ind_dim1 << 1)];
	    am -= a[k] * b[j];
/* L40: */
	}
L50:
	if (am == 0.f) {
	    goto L70;
	}
	j = ind[kp + (ind_dim1 << 1)];
	b[j] = am / a[kp];
	kpc = ip[j + (ip_dim1 << 1)];
	kl = iw[j + (iw_dim1 << 1)] + kpc - 1;
	if (kl == kpc) {
	    goto L70;
	}
	k2 = kpc + 1;
	i__2 = kl;
	for (k = k2; k <= i__2; ++k) {
	    i__ = ind[k + ind_dim1];
	    ip[i__ + ip_dim1] = -(i__3 = ip[i__ + ip_dim1], abs(i__3));
/* L60: */
	}
L70:
	;
    }
    goto L140;

/*     MULTIPLY VECTOR BY INVERSE OF TRANSPOSE OF U */
L80:
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	w[i__] = b[i__];
	b[i__] = 0.;
/* L90: */
    }
    i__1 = *n;
    for (ii = 1; ii <= i__1; ++ii) {
	i__ = iw[ii + (iw_dim1 << 2)];
	am = w[i__];
	if (am == 0.) {
	    goto L110;
	}
	j = iw[ii + iw_dim1 * 3];
	kp = ip[j + ip_dim1];
	am /= a[kp];
	b[j] = am;
	kl = iw[j + iw_dim1] + kp - 1;
	if (kp == kl) {
	    goto L110;
	}
	k2 = kp + 1;
	i__2 = kl;
	for (k = k2; k <= i__2; ++k) {
	    i__ = ind[k + (ind_dim1 << 1)];
	    w[i__] -= am * a[k];
/* L100: */
	}
L110:
	;
    }

/*     MULTIPLY VECTOR BY INVERSE OF TRANSPOSE OF L */
    if (kll > *ia) {
	return 0;
    }
    i__1 = *ia;
    for (k = kll; k <= i__1; ++k) {
	j = ind[k + (ind_dim1 << 1)];
	if (b[j] == 0.) {
	    goto L120;
	}
	i__ = ind[k + ind_dim1];
	b[i__] += a[k] * b[j];
L120:
	;
    }
    goto L140;

L130:
    xsetun_(&la05dd_1.lp);
    if (la05dd_1.lp > 0) {
	xermsg_("SLATEC", "LA05BD", "EARLIER ENTRY GAVE ERROR RETURN.", &c_n8,
		 &c__2, (ftnlen)6, (ftnlen)6, (ftnlen)32);
    }
L140:
    return 0;
} /* la05bd_ */

