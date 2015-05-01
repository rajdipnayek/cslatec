/* drkfab.f -- translated by f2c (version 12.02.01).
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
    integer igofx, inhomo, ivp, ncompd, nfcd;
} dml8sz_;

#define dml8sz_1 dml8sz_

struct {
    doublereal px, pwcnd, tnd, x, xbeg, xend, xot, xop;
    integer info[15], istkop, knswot, kop, lotjp, mnswot, nswot;
} dml15t_;

#define dml15t_1 dml15t_

struct {
    doublereal ae, re, tol;
    integer nxptsd, nic, nopg, mxnond, ndisk, ntape, neq, indpvt, integ, nps, 
	    ntpd, neqivp, numort, nfccd, icoco;
} dml18j_;

#define dml18j_1 dml18j_

struct {
    integer kkkzpw, needw, neediw, k1, k2, k3, k4, k5, k6, k7, k8, k9, k10, 
	    k11, l1, l2, kkkint, lllint;
} dml17b_;

#define dml17b_1 dml17b_

/* Table of constant values */

static integer c__1 = 1;
static integer c__0 = 0;

/* DECK DRKFAB */
/* Subroutine */ int drkfab_(integer *ncomp, doublereal *xpts, integer *nxpts,
	 integer *nfc, integer *iflag, doublereal *z__, integer *mxnon, 
	doublereal *p, integer *ntp, integer *ip, doublereal *yhp, integer *
	niv, doublereal *u, doublereal *v, doublereal *w, doublereal *s, 
	doublereal *stowa, doublereal *g, doublereal *work, integer *iwork, 
	integer *nfcc)
{
    /* System generated locals */
    integer ip_dim1, ip_offset, p_dim1, p_offset, u_dim1, u_dim2, u_offset, 
	    v_dim1, v_offset, w_dim1, w_offset, yhp_dim1, yhp_offset, i__1, 
	    i__2, i__3;

    /* Local variables */
    static integer j, kod, jon, non, idid, ipar, kopp;
    static doublereal xxop;
    static integer nfcp1, jflag;
    extern /* Subroutine */ int ddeabm_(U_fp, integer *, doublereal *, 
	    doublereal *, doublereal *, integer *, doublereal *, doublereal *,
	     integer *, doublereal *, integer *, integer *, integer *, 
	    doublereal *, integer *), dstor1_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, integer *, integer *, integer *), 
	    dderkf_(U_fp, integer *, doublereal *, doublereal *, doublereal *,
	     integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *, integer *, integer *, doublereal *, integer *);
    extern /* Subroutine */ int dbvder_();
    extern /* Subroutine */ int dreort_(integer *, doublereal *, doublereal *,
	     doublereal *, integer *, doublereal *, doublereal *, doublereal *
	    , integer *, doublereal *, integer *);

    /* Fortran I/O blocks */
    static cilist io___10 = { 0, 0, 0, 0, 0 };
    static cilist io___12 = { 0, 0, 0, 0, 0 };


/* ***BEGIN PROLOGUE  DRKFAB */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to DBVSUP */
/* ***LIBRARY   SLATEC */
/* ***TYPE      DOUBLE PRECISION (RKFAB-S, DRKFAB-D) */
/* ***AUTHOR  Watts, H. A., (SNLA) */
/* ***DESCRIPTION */

/* ********************************************************************** */

/*     Subroutine DRKFAB integrates the initial value equations using */
/*     the variable-step Runge-Kutta-Fehlberg integration scheme or */
/*     the variable-order Adams method and orthonormalization */
/*     determined by a linear dependence test. */

/* ********************************************************************** */

/* ***SEE ALSO  DBVSUP */
/* ***ROUTINES CALLED  DBVDER, DDEABM, DDERKF, DREORT, DSTOR1 */
/* ***COMMON BLOCKS    DML15T, DML17B, DML18J, DML8SZ */
/* ***REVISION HISTORY  (YYMMDD) */
/*   750601  DATE WRITTEN */
/*   890921  Realigned order of variables in certain COMMON blocks. */
/*           (WRB) */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900328  Added TYPE section.  (WRB) */
/*   910722  Updated AUTHOR section.  (ALS) */
/* ***END PROLOGUE  DRKFAB */


/*     ****************************************************************** */



/*      ***************************************************************** */
/*       INITIALIZATION OF COUNTERS AND VARIABLES. */

/*     BEGIN BLOCK PERMITTING ...EXITS TO 220 */
/*        BEGIN BLOCK PERMITTING ...EXITS TO 10 */
/* ***FIRST EXECUTABLE STATEMENT  DRKFAB */
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
    dml15t_1.x = dml15t_1.xbeg;
    jon = 1;
    dml15t_1.info[0] = 0;
    dml15t_1.info[1] = 0;
    dml15t_1.info[2] = 1;
    dml15t_1.info[3] = 1;
    work[1] = dml15t_1.xend;
/*        ...EXIT */
    if (dml18j_1.nopg == 0) {
	goto L10;
    }
    dml15t_1.info[2] = 0;
    if (dml15t_1.x == z__[1]) {
	jon = 2;
    }
L10:
    nfcp1 = *nfc + 1;

/*        *************************************************************** */
/*        *****BEGINNING OF INTEGRATION LOOP AT OUTPUT */
/*        POINTS.****************** */
/*        *************************************************************** */

    i__1 = *nxpts;
    for (kopp = 2; kopp <= i__1; ++kopp) {
	dml15t_1.kop = kopp;
	dml15t_1.xop = xpts[dml15t_1.kop];
	if (dml18j_1.ndisk == 0) {
	    kod = dml15t_1.kop;
	}

L20:

/*              STEP BY STEP INTEGRATION LOOP BETWEEN OUTPUT POINTS. */

/*              BEGIN BLOCK PERMITTING ...EXITS TO 190 */
/*                 BEGIN BLOCK PERMITTING ...EXITS TO 30 */
	xxop = dml15t_1.xop;
/*                 ...EXIT */
	if (dml18j_1.nopg == 0) {
	    goto L30;
	}
	if (dml15t_1.xend > dml15t_1.xbeg && dml15t_1.xop > z__[jon]) {
	    xxop = z__[jon];
	}
	if (dml15t_1.xend < dml15t_1.xbeg && dml15t_1.xop < z__[jon]) {
	    xxop = z__[jon];
	}
L30:

/*                 ****************************************************** */
L40:
/*                    BEGIN BLOCK PERMITTING ...EXITS TO 170 */
	switch (dml18j_1.integ) {
	    case 1:  goto L50;
	    case 2:  goto L60;
	}
/*                       DDERKF INTEGRATOR */

L50:
	dderkf_((U_fp)dbvder_, &dml18j_1.neq, &dml15t_1.x, &yhp[yhp_offset], &
		xxop, dml15t_1.info, &dml18j_1.re, &dml18j_1.ae, &idid, &work[
		1], &dml17b_1.kkkint, &iwork[1], &dml17b_1.lllint, &g[1], &
		ipar);
	goto L70;
/*                       DDEABM INTEGRATOR */

L60:
	ddeabm_((U_fp)dbvder_, &dml18j_1.neq, &dml15t_1.x, &yhp[yhp_offset], &
		xxop, dml15t_1.info, &dml18j_1.re, &dml18j_1.ae, &idid, &work[
		1], &dml17b_1.kkkint, &iwork[1], &dml17b_1.lllint, &g[1], &
		ipar);
L70:
	if (idid >= 1) {
	    goto L80;
	}
	dml15t_1.info[0] = 1;
/*                    ......EXIT */
	if (idid == -1) {
	    goto L170;
	}
	*iflag = 20 - idid;
/*     .....................EXIT */
	goto L220;
L80:

/*                       ************************************************ */
/*                           GRAM-SCHMIDT ORTHOGONALIZATION TEST FOR */
/*                           ORTHONORMALIZATION (TEMPORARILY USING U AND */
/*                           V IN THE TEST) */

	if (dml18j_1.nopg == 0) {
	    goto L100;
	}
	if (xxop == z__[jon]) {
	    goto L90;
	}

/*                             ****************************************** */
/*                                 CONTINUE INTEGRATION IF WE ARE NOT AT */
/*                                 AN OUTPUT POINT. */

/*           ..................EXIT */
	if (idid != 1) {
	    goto L200;
	}
/*                    .........EXIT */
	goto L170;
L90:
	jflag = 2;
	goto L110;
L100:
	jflag = 1;
	if (dml8sz_1.inhomo == 3 && dml15t_1.x == dml15t_1.xend) {
	    jflag = 3;
	}
L110:

	if (dml18j_1.ndisk == 0) {
	    non = dml18j_1.numort + 1;
	}
	dreort_(ncomp, &u[(kod * u_dim2 + 1) * u_dim1 + 1], &v[kod * v_dim1 + 
		1], &yhp[yhp_offset], niv, &w[non * w_dim1 + 1], &s[1], &p[
		non * p_dim1 + 1], &ip[non * ip_dim1 + 1], &stowa[1], &jflag);

	if (jflag != 30) {
	    goto L120;
	}
	*iflag = 30;
/*     .....................EXIT */
	goto L220;
L120:

	if (jflag != 10) {
	    goto L130;
	}
	dml15t_1.xop = xpts[dml15t_1.kop];
	if (dml18j_1.ndisk == 0) {
	    kod = dml15t_1.kop;
	}
/*              ............EXIT */
	goto L190;
L130:

	if (jflag == 0) {
	    goto L140;
	}

/*                          ********************************************* */
/*                              CONTINUE INTEGRATION IF WE ARE NOT AT AN */
/*                              OUTPUT POINT. */

/*           ...............EXIT */
	if (idid != 1) {
	    goto L200;
	}
/*                    ......EXIT */
	goto L170;
L140:

/*                       ************************************************ */
/*                           STORE ORTHONORMALIZED VECTORS INTO SOLUTION */
/*                           VECTORS. */

	if (dml18j_1.numort < *mxnon) {
	    goto L150;
	}
	if (dml15t_1.x == dml15t_1.xend) {
	    goto L150;
	}
	*iflag = 13;
/*     .....................EXIT */
	goto L220;
L150:

	++dml18j_1.numort;
	dstor1_(&yhp[yhp_offset], &u[(kod * u_dim2 + 1) * u_dim1 + 1], &yhp[
		nfcp1 * yhp_dim1 + 1], &v[kod * v_dim1 + 1], &c__1, &
		dml18j_1.ndisk, &dml18j_1.ntape);

/*                       ************************************************ */
/*                           STORE ORTHONORMALIZATION INFORMATION, */
/*                           INITIALIZE INTEGRATION FLAG, AND CONTINUE */
/*                           INTEGRATION TO THE NEXT ORTHONORMALIZATION */
/*                           POINT OR OUTPUT POINT. */

	z__[dml18j_1.numort] = dml15t_1.x;
	if (dml8sz_1.inhomo == 1 && dml18j_1.nps == 0) {
	    dml8sz_1.c__ = s[nfcp1] * dml8sz_1.c__;
	}
	if (dml18j_1.ndisk == 0) {
	    goto L160;
	}
	if (dml8sz_1.inhomo == 1) {
	    io___10.ciunit = dml18j_1.ntape;
	    s_wsue(&io___10);
	    i__2 = *nfcc;
	    for (j = 1; j <= i__2; ++j) {
		do_uio(&c__1, (char *)&w[j + w_dim1], (ftnlen)sizeof(
			doublereal));
	    }
	    e_wsue();
	}
	io___12.ciunit = dml18j_1.ntape;
	s_wsue(&io___12);
	i__2 = *nfcc;
	for (j = 1; j <= i__2; ++j) {
	    do_uio(&c__1, (char *)&ip[j + ip_dim1], (ftnlen)sizeof(integer));
	}
	i__3 = *ntp;
	for (j = 1; j <= i__3; ++j) {
	    do_uio(&c__1, (char *)&p[j + p_dim1], (ftnlen)sizeof(doublereal));
	}
	e_wsue();
L160:
	dml15t_1.info[0] = 0;
	++jon;
/*                 ......EXIT */
	if (dml18j_1.nopg == 1 && dml15t_1.x != dml15t_1.xop) {
	    goto L180;
	}

/*                       ************************************************ */
/*                           CONTINUE INTEGRATION IF WE ARE NOT AT AN */
/*                           OUTPUT POINT. */

/*           ............EXIT */
	if (idid != 1) {
	    goto L200;
	}
L170:
	goto L40;
L180:
L190:
	goto L20;
L200:

/*           STORAGE OF HOMOGENEOUS SOLUTIONS IN U AND THE PARTICULAR */
/*           SOLUTION IN V AT THE OUTPUT POINTS. */

	dstor1_(&u[(kod * u_dim2 + 1) * u_dim1 + 1], &yhp[yhp_offset], &v[kod 
		* v_dim1 + 1], &yhp[nfcp1 * yhp_dim1 + 1], &c__0, &
		dml18j_1.ndisk, &dml18j_1.ntape);
/* L210: */
    }
/*        *************************************************************** */
/*        *************************************************************** */

    *iflag = 0;
L220:
    return 0;
} /* drkfab_ */

