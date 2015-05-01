/* reort.f -- translated by f2c (version 12.02.01).
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
    real c__, xsav;
    integer igofx, inhomo, ivp, ncompd, nfc;
} ml8sz_;

#define ml8sz_1 ml8sz_

struct {
    real px, pwcnd, tnd, x, xbeg, xend, xot, xop;
    integer info[15], istkop, knswot, kop, lotjp, mnswot, nswot;
} ml15to_;

#define ml15to_1 ml15to_

struct {
    real ae, re, tol;
    integer nxpts, nic, nopg, mxnon, ndisk, ntape, neq, indpvt, integ, nps, 
	    ntp, neqivp, numort, nfcc, icoco;
} ml18jr_;

#define ml18jr_1 ml18jr_

/* Table of constant values */

static integer c__1 = 1;
static integer c__0 = 0;

/* DECK REORT */
/* Subroutine */ int reort_(integer *ncomp, real *y, real *yp, real *yhp, 
	integer *niv, real *w, real *s, real *p, integer *ip, real *stowa, 
	integer *iflag)
{
    /* System generated locals */
    integer y_dim1, y_offset, yhp_dim1, yhp_offset, i__1, i__2;
    real r__1, r__2;

    /* Local variables */
    static integer j, k, l, kk;
    static real dx, dnd;
    static integer ijk;
    static real srp;
    static integer nfcp;
    static real dndt, wcnd;
    extern doublereal sdot_(integer *, real *, integer *, real *, integer *);
    static real ypnm;
    extern /* Subroutine */ int stor1_(real *, real *, real *, real *, 
	    integer *, integer *, integer *);
    static integer mflag;
    extern /* Subroutine */ int mgsbv_(integer *, integer *, real *, integer *
	    , integer *, integer *, real *, real *, integer *, integer *, 
	    real *, real *, real *);
    static real vnorm;
    extern /* Subroutine */ int stway_(real *, real *, real *, integer *, 
	    real *);

/* ***BEGIN PROLOGUE  REORT */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to BVSUP */
/* ***LIBRARY   SLATEC */
/* ***TYPE      SINGLE PRECISION (REORT-S, DREORT-D) */
/* ***AUTHOR  Watts, H. A., (SNLA) */
/* ***DESCRIPTION */

/* ********************************************************************** */
/*   INPUT */
/* ********* */
/*     Y, YP and YHP = homogeneous solution matrix and particular */
/*                     solution vector to be orthonormalized. */
/*     IFLAG = 1 --  store YHP into Y and YP, test for */
/*                   reorthonormalization, orthonormalize if needed, */
/*                   save restart data. */
/*             2 --  store YHP into Y and YP, reorthonormalization, */
/*                   no restarts. */
/*                   (preset orthonormalization mode) */
/*             3 --  store YHP into Y and YP, reorthonormalization */
/*                   (when INHOMO=3 and X=XEND). */
/* ********************************************************************** */
/*   OUTPUT */
/* ********* */
/*     Y, YP = orthonormalized solutions. */
/*     NIV = number of independent vectors returned from DMGSBV. */
/*     IFLAG = 0 --  reorthonormalization was performed. */
/*            10 --  solution process must be restarted at the last */
/*                   orthonormalization point. */
/*            30 --  solutions are linearly dependent, problem must */
/*                   be restarted from the beginning. */
/*     W, P, IP = orthonormalization information. */
/* ********************************************************************** */

/* ***SEE ALSO  BVSUP */
/* ***ROUTINES CALLED  MGSBV, SDOT, STOR1, STWAY */
/* ***COMMON BLOCKS    ML15TO, ML18JR, ML8SZ */
/* ***REVISION HISTORY  (YYMMDD) */
/*   750601  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890831  Modified array declarations.  (WRB) */
/*   890921  Realigned order of variables in certain COMMON blocks. */
/*           (WRB) */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900328  Added TYPE section.  (WRB) */
/*   910722  Updated AUTHOR section.  (ALS) */
/* ***END PROLOGUE  REORT */


/* ********************************************************************** */


/* ********************************************************************** */
/* ***FIRST EXECUTABLE STATEMENT  REORT */
    /* Parameter adjustments */
    yhp_dim1 = *ncomp;
    yhp_offset = 1 + yhp_dim1;
    yhp -= yhp_offset;
    y_dim1 = *ncomp;
    y_offset = 1 + y_dim1;
    y -= y_offset;
    --yp;
    --w;
    --s;
    --p;
    --ip;
    --stowa;

    /* Function Body */
    nfcp = ml8sz_1.nfc + 1;

/*     CHECK TO SEE IF ORTHONORMALIZATION TEST IS TO BE PERFORMED */

    if (*iflag != 1) {
	goto L5;
    }
    ++ml15to_1.knswot;
    if (ml15to_1.knswot >= ml15to_1.nswot) {
	goto L5;
    }
    if ((ml15to_1.xend - ml15to_1.x) * (ml15to_1.x - ml15to_1.xot) < 0.f) {
	return 0;
    }
L5:
    stor1_(&y[y_offset], &yhp[yhp_offset], &yp[1], &yhp[nfcp * yhp_dim1 + 1], 
	    &c__1, &c__0, &c__0);

/*     **************************************** */

/*     ORTHOGONALIZE THE HOMOGENEOUS SOLUTIONS Y */
/*     AND PARTICULAR SOLUTION YP. */

    *niv = ml8sz_1.nfc;
    mgsbv_(ncomp, &ml8sz_1.nfc, &y[y_offset], ncomp, niv, &mflag, &s[1], &p[1]
	    , &ip[1], &ml8sz_1.inhomo, &yp[1], &w[1], &wcnd);

/*     **************************************** */

/*  CHECK FOR LINEAR DEPENDENCE OF THE SOLUTIONS. */

    if (mflag == 0) {
	goto L25;
    }
    if (*iflag == 2) {
	goto L15;
    }
    if (ml15to_1.nswot > 1 || ml15to_1.lotjp == 0) {
	goto L20;
    }
L15:
    *iflag = 30;
    return 0;

/*     RETRIEVE DATA FOR A RESTART AT LAST ORTHONORMALIZATION POINT */

L20:
    stway_(&y[y_offset], &yp[1], &yhp[yhp_offset], &c__1, &stowa[1]);
    ml15to_1.lotjp = 1;
    ml15to_1.nswot = 1;
    ml15to_1.knswot = 0;
    ml15to_1.mnswot /= 2;
    ml15to_1.tnd += 1.f;
    *iflag = 10;
    return 0;

/*     **************************************** */

L25:
    if (*iflag != 1) {
	goto L60;
    }

/*     TEST FOR ORTHONORMALIZATION */

    if (wcnd < ml18jr_1.tol * 50.f) {
	goto L60;
    }
    i__1 = nfcp;
    for (ijk = 1; ijk <= i__1; ++ijk) {
	if (s[ijk] > 1e20f) {
	    goto L60;
	}
/* L30: */
    }

/*     USE LINEAR EXTRAPOLATION ON LOGARITHMIC VALUES OF THE NORM */
/*     DECREMENTS TO DETERMINE NEXT ORTHONORMALIZATION CHECKPOINT. */
/*     OTHER CONTROLS ON THE NUMBER OF STEPS TO THE NEXT CHECKPOINT */
/*     ARE ADDED FOR SAFETY PURPOSES. */

    ml15to_1.nswot = ml15to_1.knswot;
    ml15to_1.knswot = 0;
    ml15to_1.lotjp = 0;
    wcnd = r_lg10(&wcnd);
    if (wcnd > ml15to_1.tnd + 3.f) {
	ml15to_1.nswot <<= 1;
    }
    if (wcnd >= ml15to_1.pwcnd) {
	goto L40;
    }
    dx = ml15to_1.x - ml15to_1.px;
    dnd = ml15to_1.pwcnd - wcnd;
    if (dnd >= 4.f) {
	ml15to_1.nswot /= 2;
    }
    dndt = wcnd - ml15to_1.tnd;
    if ((r__2 = dx * dndt, dabs(r__2)) > dnd * (r__1 = ml15to_1.xend - 
	    ml15to_1.x, dabs(r__1))) {
	goto L40;
    }
    ml15to_1.xot = ml15to_1.x + dx * dndt / dnd;
    goto L50;
L40:
    ml15to_1.xot = ml15to_1.xend;
L50:
    ml15to_1.nswot = min(ml15to_1.mnswot,ml15to_1.nswot);
    ml15to_1.pwcnd = wcnd;
    ml15to_1.px = ml15to_1.x;
    return 0;

/*     **************************************** */

/*     ORTHONORMALIZATION NECESSARY SO WE NORMALIZE THE HOMOGENEOUS */
/*     SOLUTION VECTORS AND CHANGE W ACCORDINGLY. */

L60:
    ml15to_1.nswot = 1;
    ml15to_1.knswot = 0;
    ml15to_1.lotjp = 1;
    kk = 1;
    l = 1;
    i__1 = ml18jr_1.nfcc;
    for (k = 1; k <= i__1; ++k) {
	srp = sqrt(p[kk]);
	if (ml8sz_1.inhomo == 1) {
	    w[k] = srp * w[k];
	}
	vnorm = 1.f / srp;
	p[kk] = vnorm;
	kk = kk + ml18jr_1.nfcc + 1 - k;
	if (ml8sz_1.nfc == ml18jr_1.nfcc) {
	    goto L63;
	}
	if (l != k / 2) {
	    goto L70;
	}
L63:
	i__2 = *ncomp;
	for (j = 1; j <= i__2; ++j) {
/* L65: */
	    y[j + l * y_dim1] *= vnorm;
	}
	++l;
L70:
	;
    }

    if (ml8sz_1.inhomo != 1 || ml18jr_1.nps == 1) {
	goto L100;
    }

/*     NORMALIZE THE PARTICULAR SOLUTION */

    ypnm = sdot_(ncomp, &yp[1], &c__1, &yp[1], &c__1);
    if (ypnm == 0.f) {
	ypnm = 1.f;
    }
    ypnm = sqrt(ypnm);
    s[nfcp] = ypnm;
    i__1 = *ncomp;
    for (j = 1; j <= i__1; ++j) {
/* L80: */
	yp[j] /= ypnm;
    }
    i__1 = ml18jr_1.nfcc;
    for (j = 1; j <= i__1; ++j) {
/* L90: */
	w[j] = ml8sz_1.c__ * w[j];
    }

L100:
    if (*iflag == 1) {
	stway_(&y[y_offset], &yp[1], &yhp[yhp_offset], &c__0, &stowa[1]);
    }
    *iflag = 0;
    return 0;
} /* reort_ */

