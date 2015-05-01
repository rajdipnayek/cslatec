/* dreort.f -- translated by f2c (version 12.02.01).
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
    doublereal c__, xsav;
    integer igofx, inhomo, ivp, ncompd, nfc;
} dml8sz_;

#define dml8sz_1 dml8sz_

struct {
    doublereal px, pwcnd, tnd, x, xbeg, xend, xot, xop;
    integer info[15], istkop, knswot, kop, lotjp, mnswot, nswot;
} dml15t_;

#define dml15t_1 dml15t_

struct {
    doublereal ae, re, tol;
    integer nxpts, nic, nopg, mxnon, ndisk, ntape, neq, indpvt, integ, nps, 
	    ntp, neqivp, numort, nfcc, icoco;
} dml18j_;

#define dml18j_1 dml18j_

/* Table of constant values */

static integer c__1 = 1;
static integer c__0 = 0;

/* DECK DREORT */
/* Subroutine */ int dreort_(integer *ncomp, doublereal *y, doublereal *yp, 
	doublereal *yhp, integer *niv, doublereal *w, doublereal *s, 
	doublereal *p, integer *ip, doublereal *stowa, integer *iflag)
{
    /* System generated locals */
    integer y_dim1, y_offset, yhp_dim1, yhp_offset, i__1, i__2;
    doublereal d__1, d__2;

    /* Local variables */
    static integer j, k, l, kk;
    static doublereal dx, dnd;
    static integer ijk;
    static doublereal srp;
    static integer nfcp;
    static doublereal dndt;
    extern doublereal ddot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    static doublereal wcnd, ypnm;
    static integer mflag;
    static doublereal vnorm;
    extern /* Subroutine */ int dstor1_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, integer *, integer *, integer *), 
	    dmgsbv_(integer *, integer *, doublereal *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *), dstway_(doublereal *, 
	    doublereal *, doublereal *, integer *, doublereal *);

/* ***BEGIN PROLOGUE  DREORT */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to DBVSUP */
/* ***LIBRARY   SLATEC */
/* ***TYPE      DOUBLE PRECISION (REORT-S, DREORT-D) */
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

/* ***SEE ALSO  DBVSUP */
/* ***ROUTINES CALLED  DDOT, DMGSBV, DSTOR1, DSTWAY */
/* ***COMMON BLOCKS    DML15T, DML18J, DML8SZ */
/* ***REVISION HISTORY  (YYMMDD) */
/*   750601  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890831  Modified array declarations.  (WRB) */
/*   890921  Realigned order of variables in certain COMMON blocks. */
/*           (WRB) */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900328  Added TYPE section.  (WRB) */
/*   910722  Updated AUTHOR section.  (ALS) */
/* ***END PROLOGUE  DREORT */


/*     ****************************************************************** */


/* ********************************************************************** */
/*     BEGIN BLOCK PERMITTING ...EXITS TO 210 */
/*        BEGIN BLOCK PERMITTING ...EXITS TO 10 */
/* ***FIRST EXECUTABLE STATEMENT  DREORT */
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
    nfcp = dml8sz_1.nfc + 1;

/*           CHECK TO SEE IF ORTHONORMALIZATION TEST IS TO BE PERFORMED */

/*        ...EXIT */
    if (*iflag != 1) {
	goto L10;
    }
    ++dml15t_1.knswot;
/*        ...EXIT */
    if (dml15t_1.knswot >= dml15t_1.nswot) {
	goto L10;
    }
/*     ......EXIT */
    if ((dml15t_1.xend - dml15t_1.x) * (dml15t_1.x - dml15t_1.xot) < 0.) {
	goto L210;
    }
L10:
    dstor1_(&y[y_offset], &yhp[yhp_offset], &yp[1], &yhp[nfcp * yhp_dim1 + 1],
	     &c__1, &c__0, &c__0);

/*        *************************************************************** */

/*        ORTHOGONALIZE THE HOMOGENEOUS SOLUTIONS Y */
/*        AND PARTICULAR SOLUTION YP. */

    *niv = dml8sz_1.nfc;
    dmgsbv_(ncomp, &dml8sz_1.nfc, &y[y_offset], ncomp, niv, &mflag, &s[1], &p[
	    1], &ip[1], &dml8sz_1.inhomo, &yp[1], &w[1], &wcnd);

/*           ************************************************************ */

/*        CHECK FOR LINEAR DEPENDENCE OF THE SOLUTIONS. */

    if (mflag == 0) {
	goto L50;
    }
/*           BEGIN BLOCK PERMITTING ...EXITS TO 40 */
    if (*iflag == 2) {
	goto L30;
    }
    if (dml15t_1.nswot <= 1 && dml15t_1.lotjp != 0) {
	goto L20;
    }

/*                    RETRIEVE DATA FOR A RESTART AT LAST */
/*                    ORTHONORMALIZATION POINT */

    dstway_(&y[y_offset], &yp[1], &yhp[yhp_offset], &c__1, &stowa[1]);
    dml15t_1.lotjp = 1;
    dml15t_1.nswot = 1;
    dml15t_1.knswot = 0;
    dml15t_1.mnswot /= 2;
    dml15t_1.tnd += 1.;
    *iflag = 10;
/*           .........EXIT */
    goto L40;
L20:
L30:
    *iflag = 30;
L40:
    goto L200;
L50:
/*           BEGIN BLOCK PERMITTING ...EXITS TO 190 */
/*              BEGIN BLOCK PERMITTING ...EXITS TO 110 */

/*                 ****************************************************** */

/*              ...EXIT */
    if (*iflag != 1) {
	goto L110;
    }

/*                 TEST FOR ORTHONORMALIZATION */

/*              ...EXIT */
    if (wcnd < dml18j_1.tol * 50.) {
	goto L110;
    }
    i__1 = nfcp;
    for (ijk = 1; ijk <= i__1; ++ijk) {
/*              ......EXIT */
	if (s[ijk] > 1e20) {
	    goto L110;
	}
/* L60: */
    }

/*                 USE LINEAR EXTRAPOLATION ON LOGARITHMIC VALUES OF THE */
/*                 NORM DECREMENTS TO DETERMINE NEXT ORTHONORMALIZATION */
/*                 CHECKPOINT.  OTHER CONTROLS ON THE NUMBER OF STEPS TO */
/*                 THE NEXT CHECKPOINT ARE ADDED FOR SAFETY PURPOSES. */

    dml15t_1.nswot = dml15t_1.knswot;
    dml15t_1.knswot = 0;
    dml15t_1.lotjp = 0;
    wcnd = d_lg10(&wcnd);
    if (wcnd > dml15t_1.tnd + 3.) {
	dml15t_1.nswot <<= 1;
    }
    if (wcnd < dml15t_1.pwcnd) {
	goto L70;
    }
    dml15t_1.xot = dml15t_1.xend;
    dml15t_1.nswot = min(dml15t_1.mnswot,dml15t_1.nswot);
    dml15t_1.pwcnd = wcnd;
    dml15t_1.px = dml15t_1.x;
    goto L100;
L70:
    dx = dml15t_1.x - dml15t_1.px;
    dnd = dml15t_1.pwcnd - wcnd;
    if (dnd >= 4.) {
	dml15t_1.nswot /= 2;
    }
    dndt = wcnd - dml15t_1.tnd;
    if ((d__2 = dx * dndt, abs(d__2)) <= dnd * (d__1 = dml15t_1.xend - 
	    dml15t_1.x, abs(d__1))) {
	goto L80;
    }
    dml15t_1.xot = dml15t_1.xend;
    dml15t_1.nswot = min(dml15t_1.mnswot,dml15t_1.nswot);
    dml15t_1.pwcnd = wcnd;
    dml15t_1.px = dml15t_1.x;
    goto L90;
L80:
    dml15t_1.xot = dml15t_1.x + dx * dndt / dnd;
    dml15t_1.nswot = min(dml15t_1.mnswot,dml15t_1.nswot);
    dml15t_1.pwcnd = wcnd;
    dml15t_1.px = dml15t_1.x;
L90:
L100:
/*           ......EXIT */
    goto L190;
L110:

/*              ********************************************************* */

/*              ORTHONORMALIZATION NECESSARY SO WE NORMALIZE THE */
/*              HOMOGENEOUS SOLUTION VECTORS AND CHANGE W ACCORDINGLY. */

    dml15t_1.nswot = 1;
    dml15t_1.knswot = 0;
    dml15t_1.lotjp = 1;
    kk = 1;
    l = 1;
    i__1 = dml18j_1.nfcc;
    for (k = 1; k <= i__1; ++k) {
/*                 BEGIN BLOCK PERMITTING ...EXITS TO 140 */
	srp = sqrt(p[kk]);
	if (dml8sz_1.inhomo == 1) {
	    w[k] = srp * w[k];
	}
	vnorm = 1. / srp;
	p[kk] = vnorm;
	kk = kk + dml18j_1.nfcc + 1 - k;
	if (dml8sz_1.nfc == dml18j_1.nfcc) {
	    goto L120;
	}
/*                 ......EXIT */
	if (l != k / 2) {
	    goto L140;
	}
L120:
	i__2 = *ncomp;
	for (j = 1; j <= i__2; ++j) {
	    y[j + l * y_dim1] *= vnorm;
/* L130: */
	}
	++l;
L140:
/* L150: */
	;
    }

    if (dml8sz_1.inhomo != 1 || dml18j_1.nps == 1) {
	goto L180;
    }

/*                 NORMALIZE THE PARTICULAR SOLUTION */

    ypnm = ddot_(ncomp, &yp[1], &c__1, &yp[1], &c__1);
    if (ypnm == 0.) {
	ypnm = 1.;
    }
    ypnm = sqrt(ypnm);
    s[nfcp] = ypnm;
    i__1 = *ncomp;
    for (j = 1; j <= i__1; ++j) {
	yp[j] /= ypnm;
/* L160: */
    }
    i__1 = dml18j_1.nfcc;
    for (j = 1; j <= i__1; ++j) {
	w[j] = dml8sz_1.c__ * w[j];
/* L170: */
    }
L180:

    if (*iflag == 1) {
	dstway_(&y[y_offset], &yp[1], &yhp[yhp_offset], &c__0, &stowa[1]);
    }
    *iflag = 0;
L190:
L200:
L210:
    return 0;
} /* dreort_ */

