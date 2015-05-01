/* bsplvd.f -- translated by f2c (version 12.02.01).
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
static integer c__0 = 0;
static integer c__2 = 2;

/* DECK BSPLVD */
/* Subroutine */ int bsplvd_(real *t, integer *k, real *x, integer *ileft, 
	real *vnikx, integer *nderiv)
{
    /* System generated locals */
    integer vnikx_dim1, vnikx_offset, i__1, i__2, i__3;

    /* Local variables */
    static real a[400]	/* was [20][20] */;
    static integer i__, j, l, m;
    static real v;
    static integer jm1, kmd;
    static real diff, fkmd;
    static integer jlow, ipkmd, ideriv, idervm;
    extern /* Subroutine */ int bsplvn_(real *, integer *, integer *, real *, 
	    integer *, real *);

/* ***BEGIN PROLOGUE  BSPLVD */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to FC */
/* ***LIBRARY   SLATEC */
/* ***TYPE      SINGLE PRECISION (BSPLVD-S, DFSPVD-D) */
/* ***AUTHOR  (UNKNOWN) */
/* ***DESCRIPTION */

/* Calculates value and deriv.s of all B-splines which do not vanish at X */

/*  Fill VNIKX(J,IDERIV), J=IDERIV, ... ,K  with nonzero values of */
/*  B-splines of order K+1-IDERIV , IDERIV=NDERIV, ... ,1, by repeated */
/*  calls to BSPLVN */

/* ***SEE ALSO  FC */
/* ***ROUTINES CALLED  BSPLVN */
/* ***REVISION HISTORY  (YYMMDD) */
/*   780801  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890831  Modified array declarations.  (WRB) */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900328  Added TYPE section.  (WRB) */
/* ***END PROLOGUE  BSPLVD */
/* ***FIRST EXECUTABLE STATEMENT  BSPLVD */
    /* Parameter adjustments */
    --t;
    vnikx_dim1 = *k;
    vnikx_offset = 1 + vnikx_dim1;
    vnikx -= vnikx_offset;

    /* Function Body */
    i__1 = *k + 1 - *nderiv;
    bsplvn_(&t[1], &i__1, &c__1, x, ileft, &vnikx[*nderiv + *nderiv * 
	    vnikx_dim1]);
    if (*nderiv <= 1) {
	goto L99;
    }
    ideriv = *nderiv;
    i__1 = *nderiv;
    for (i__ = 2; i__ <= i__1; ++i__) {
	idervm = ideriv - 1;
	i__2 = *k;
	for (j = ideriv; j <= i__2; ++j) {
/* L11: */
	    vnikx[j - 1 + idervm * vnikx_dim1] = vnikx[j + ideriv * 
		    vnikx_dim1];
	}
	ideriv = idervm;
	bsplvn_(&t[1], &c__0, &c__2, x, ileft, &vnikx[ideriv + ideriv * 
		vnikx_dim1]);
/* L15: */
    }

    i__1 = *k;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = *k;
	for (j = 1; j <= i__2; ++j) {
/* L19: */
	    a[i__ + j * 20 - 21] = 0.f;
	}
/* L20: */
	a[i__ + i__ * 20 - 21] = 1.f;
    }
    kmd = *k;
    i__1 = *nderiv;
    for (m = 2; m <= i__1; ++m) {
	--kmd;
	fkmd = (real) kmd;
	i__ = *ileft;
	j = *k;
L21:
	jm1 = j - 1;
	ipkmd = i__ + kmd;
	diff = t[ipkmd] - t[i__];
	if (jm1 == 0) {
	    goto L26;
	}
	if (diff == 0.f) {
	    goto L25;
	}
	i__2 = j;
	for (l = 1; l <= i__2; ++l) {
/* L24: */
	    a[l + j * 20 - 21] = (a[l + j * 20 - 21] - a[l + (j - 1) * 20 - 
		    21]) / diff * fkmd;
	}
L25:
	j = jm1;
	--i__;
	goto L21;
L26:
	if (diff == 0.f) {
	    goto L30;
	}
	a[0] = a[0] / diff * fkmd;

L30:
	i__2 = *k;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    v = 0.f;
	    jlow = max(i__,m);
	    i__3 = *k;
	    for (j = jlow; j <= i__3; ++j) {
/* L35: */
		v = a[i__ + j * 20 - 21] * vnikx[j + m * vnikx_dim1] + v;
	    }
/* L40: */
	    vnikx[i__ + m * vnikx_dim1] = v;
	}
    }
L99:
    return 0;
} /* bsplvd_ */

