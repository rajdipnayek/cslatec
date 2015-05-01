/* rkfab.f -- translated by f2c (version 12.02.01).
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
    integer igofx, inhomo, ivp, ncompd, nfcd;
} ml8sz_;

#define ml8sz_1 ml8sz_

struct {
    real px, pwcnd, tnd, x, xbeg, xend, xot, xop;
    integer info[15], istkop, knswot, kop, lotjp, mnswot, nswot;
} ml15to_;

#define ml15to_1 ml15to_

struct {
    real ae, re, tol;
    integer nxptsd, nic, nopg, mxnond, ndisk, ntape, neq, indpvt, integ, nps, 
	    ntpd, neqivp, numort, nfccd, icoco;
} ml18jr_;

#define ml18jr_1 ml18jr_

struct {
    integer kkkzpw, needw, neediw, k1, k2, k3, k4, k5, k6, k7, k8, k9, k10, 
	    k11, l1, l2, kkkint, lllint;
} ml17bw_;

#define ml17bw_1 ml17bw_

/* Table of constant values */

static integer c__1 = 1;
static integer c__0 = 0;

/* DECK RKFAB */
/* Subroutine */ int rkfab_(integer *ncomp, real *xpts, integer *nxpts, 
	integer *nfc, integer *iflag, real *z__, integer *mxnon, real *p, 
	integer *ntp, integer *ip, real *yhp, integer *niv, real *u, real *v, 
	real *w, real *s, real *stowa, real *g, real *work, integer *iwork, 
	integer *nfcc)
{
    /* System generated locals */
    integer p_dim1, p_offset, ip_dim1, ip_offset, u_dim1, u_dim2, u_offset, 
	    v_dim1, v_offset, w_dim1, w_offset, yhp_dim1, yhp_offset, i__1, 
	    i__2, i__3;

    /* Local variables */
    static integer j, kod, jon, non, idid, ipar, kopp;
    static real xxop;
    static integer nfcp1;
    extern /* Subroutine */ int deabm_(U_fp, integer *, real *, real *, real *
	    , integer *, real *, real *, integer *, real *, integer *, 
	    integer *, integer *, real *, integer *), stor1_(real *, real *, 
	    real *, real *, integer *, integer *, integer *);
    static integer jflag;
    extern /* Subroutine */ int derkf_(U_fp, integer *, real *, real *, real *
	    , integer *, real *, real *, integer *, real *, integer *, 
	    integer *, integer *, real *, integer *);
    extern /* Subroutine */ int bvder_();
    extern /* Subroutine */ int reort_(integer *, real *, real *, real *, 
	    integer *, real *, real *, real *, integer *, real *, integer *);

    /* Fortran I/O blocks */
    static cilist io___10 = { 0, 0, 0, 0, 0 };
    static cilist io___12 = { 0, 0, 0, 0, 0 };


/* ***BEGIN PROLOGUE  RKFAB */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to BVSUP */
/* ***LIBRARY   SLATEC */
/* ***TYPE      SINGLE PRECISION (RKFAB-S, DRKFAB-D) */
/* ***AUTHOR  Watts, H. A., (SNLA) */
/* ***DESCRIPTION */

/* ********************************************************************** */

/*     Subroutine RKFAB integrates the initial value equations using */
/*     the variable-step RUNGE-KUTTA-FEHLBERG integration scheme or */
/*     the variable-order ADAMS method and orthonormalization */
/*     determined by a linear dependence test. */

/* ********************************************************************** */

/* ***SEE ALSO  BVSUP */
/* ***ROUTINES CALLED  BVDER, DEABM, DERKF, REORT, STOR1 */
/* ***COMMON BLOCKS    ML15TO, ML17BW, ML18JR, ML8SZ */
/* ***REVISION HISTORY  (YYMMDD) */
/*   750601  DATE WRITTEN */
/*   890831  Modified array declarations.  (WRB) */
/*   890921  Realigned order of variables in certain COMMON blocks. */
/*           (WRB) */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900328  Added TYPE section.  (WRB) */
/*   910722  Updated AUTHOR section.  (ALS) */
/* ***END PROLOGUE  RKFAB */


/* ********************************************************************** */



/* ********************************************************************** */
/*  INITIALIZATION OF COUNTERS AND VARIABLES. */

/* ***FIRST EXECUTABLE STATEMENT  RKFAB */
    /* Parameter adjustments */
    v_dim1 = *ncomp;
    v_offset = 1 + v_dim1;
    v -= v_offset;
    yhp_dim1 = *ncomp;
    yhp_offset = 1 + yhp_dim1;
    yhp -= yhp_offset;
    --xpts;
    u_dim1 = *ncomp;
    u_dim2 = *nfc;
    u_offset = 1 + u_dim1 * (1 + u_dim2);
    u -= u_offset;
    --z__;
    p_dim1 = *ntp;
    p_offset = 1 + p_dim1;
    p -= p_offset;
    --s;
    --stowa;
    --g;
    --work;
    --iwork;
    w_dim1 = *nfcc;
    w_offset = 1 + w_dim1;
    w -= w_offset;
    ip_dim1 = *nfcc;
    ip_offset = 1 + ip_dim1;
    ip -= ip_offset;

    /* Function Body */
    kod = 1;
    non = 1;
    ml15to_1.x = ml15to_1.xbeg;
    jon = 1;
    ml15to_1.info[0] = 0;
    ml15to_1.info[1] = 0;
    ml15to_1.info[2] = 1;
    ml15to_1.info[3] = 1;
    work[1] = ml15to_1.xend;
    if (ml18jr_1.nopg == 0) {
	goto L1;
    }
    ml15to_1.info[2] = 0;
    if (ml15to_1.x == z__[1]) {
	jon = 2;
    }
L1:
    nfcp1 = *nfc + 1;

/* ********************************************************************** */
/* *****BEGINNING OF INTEGRATION LOOP AT OUTPUT POINTS.****************** */
/* ********************************************************************** */

    i__1 = *nxpts;
    for (kopp = 2; kopp <= i__1; ++kopp) {
	ml15to_1.kop = kopp;

L5:
	ml15to_1.xop = xpts[ml15to_1.kop];
	if (ml18jr_1.ndisk == 0) {
	    kod = ml15to_1.kop;
	}

/*     STEP BY STEP INTEGRATION LOOP BETWEEN OUTPUT POINTS. */

L10:
	xxop = ml15to_1.xop;
	if (ml18jr_1.nopg == 0) {
	    goto L15;
	}
	if (ml15to_1.xend > ml15to_1.xbeg && ml15to_1.xop > z__[jon]) {
	    xxop = z__[jon];
	}
	if (ml15to_1.xend < ml15to_1.xbeg && ml15to_1.xop < z__[jon]) {
	    xxop = z__[jon];
	}

/* ********************************************************************** */
L15:
	switch (ml18jr_1.integ) {
	    case 1:  goto L20;
	    case 2:  goto L25;
	}
/*     DERKF INTEGRATOR */

L20:
	derkf_((U_fp)bvder_, &ml18jr_1.neq, &ml15to_1.x, &yhp[yhp_offset], &
		xxop, ml15to_1.info, &ml18jr_1.re, &ml18jr_1.ae, &idid, &work[
		1], &ml17bw_1.kkkint, &iwork[1], &ml17bw_1.lllint, &g[1], &
		ipar);
	goto L28;
/*     DEABM INTEGRATOR */

L25:
	deabm_((U_fp)bvder_, &ml18jr_1.neq, &ml15to_1.x, &yhp[yhp_offset], &
		xxop, ml15to_1.info, &ml18jr_1.re, &ml18jr_1.ae, &idid, &work[
		1], &ml17bw_1.kkkint, &iwork[1], &ml17bw_1.lllint, &g[1], &
		ipar);
L28:
	if (idid >= 1) {
	    goto L30;
	}
	ml15to_1.info[0] = 1;
	if (idid == -1) {
	    goto L15;
	}
	*iflag = 20 - idid;
	return 0;

/* ********************************************************************** */
/*     GRAM-SCHMIDT ORTHOGONALIZATION TEST FOR ORTHONORMALIZATION */
/*     (TEMPORARILY USING U AND V IN THE TEST) */

L30:
	if (ml18jr_1.nopg == 0) {
	    goto L35;
	}
	if (xxop != z__[jon]) {
	    goto L100;
	}
	jflag = 2;
	goto L40;
L35:
	jflag = 1;
	if (ml8sz_1.inhomo == 3 && ml15to_1.x == ml15to_1.xend) {
	    jflag = 3;
	}

L40:
	if (ml18jr_1.ndisk == 0) {
	    non = ml18jr_1.numort + 1;
	}
	reort_(ncomp, &u[(kod * u_dim2 + 1) * u_dim1 + 1], &v[kod * v_dim1 + 
		1], &yhp[yhp_offset], niv, &w[non * w_dim1 + 1], &s[1], &p[
		non * p_dim1 + 1], &ip[non * ip_dim1 + 1], &stowa[1], &jflag);

	if (jflag != 30) {
	    goto L45;
	}
	*iflag = 30;
	return 0;

L45:
	if (jflag == 10) {
	    goto L5;
	}

	if (jflag != 0) {
	    goto L100;
	}

/* ********************************************************************** */
/*     STORE ORTHONORMALIZED VECTORS INTO SOLUTION VECTORS. */

	if (ml18jr_1.numort < *mxnon) {
	    goto L65;
	}
	if (ml15to_1.x == ml15to_1.xend) {
	    goto L65;
	}
	*iflag = 13;
	return 0;

L65:
	++ml18jr_1.numort;
	stor1_(&yhp[yhp_offset], &u[(kod * u_dim2 + 1) * u_dim1 + 1], &yhp[
		nfcp1 * yhp_dim1 + 1], &v[kod * v_dim1 + 1], &c__1, &
		ml18jr_1.ndisk, &ml18jr_1.ntape);

/* ********************************************************************** */
/*     STORE ORTHONORMALIZATION INFORMATION, INITIALIZE */
/*     INTEGRATION FLAG, AND CONTINUE INTEGRATION TO THE NEXT */
/*     ORTHONORMALIZATION POINT OR OUTPUT POINT. */

	z__[ml18jr_1.numort] = ml15to_1.x;
	if (ml8sz_1.inhomo == 1 && ml18jr_1.nps == 0) {
	    ml8sz_1.c__ = s[nfcp1] * ml8sz_1.c__;
	}
	if (ml18jr_1.ndisk == 0) {
	    goto L90;
	}
	if (ml8sz_1.inhomo == 1) {
	    io___10.ciunit = ml18jr_1.ntape;
	    s_wsue(&io___10);
	    i__2 = *nfcc;
	    for (j = 1; j <= i__2; ++j) {
		do_uio(&c__1, (char *)&w[j + w_dim1], (ftnlen)sizeof(real));
	    }
	    e_wsue();
	}
	io___12.ciunit = ml18jr_1.ntape;
	s_wsue(&io___12);
	i__2 = *nfcc;
	for (j = 1; j <= i__2; ++j) {
	    do_uio(&c__1, (char *)&ip[j + ip_dim1], (ftnlen)sizeof(integer));
	}
	i__3 = *ntp;
	for (j = 1; j <= i__3; ++j) {
	    do_uio(&c__1, (char *)&p[j + p_dim1], (ftnlen)sizeof(real));
	}
	e_wsue();
L90:
	ml15to_1.info[0] = 0;
	++jon;
	if (ml18jr_1.nopg == 1 && ml15to_1.x != ml15to_1.xop) {
	    goto L10;
	}

/* ********************************************************************** */
/*     CONTINUE INTEGRATION IF WE ARE NOT AT AN OUTPUT POINT. */

L100:
	if (idid == 1) {
	    goto L15;
	}

/*     STORAGE OF HOMOGENEOUS SOLUTIONS IN U AND THE PARTICULAR */
/*     SOLUTION IN V AT THE OUTPUT POINTS. */

	stor1_(&u[(kod * u_dim2 + 1) * u_dim1 + 1], &yhp[yhp_offset], &v[kod *
		 v_dim1 + 1], &yhp[nfcp1 * yhp_dim1 + 1], &c__0, &
		ml18jr_1.ndisk, &ml18jr_1.ntape);
/* L110: */
    }
/* ********************************************************************** */
/* ********************************************************************** */

    *iflag = 0;
    return 0;
} /* rkfab_ */

