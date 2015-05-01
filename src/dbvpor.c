/* dbvpor.f -- translated by f2c (version 12.02.01).
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
    integer nxptsd, nicd, nopg, mxnond, ndisk, ntape, neq, indpvt, integ, nps,
	     ntpd, neqivp, numort, nfccd, icoco;
} dml18j_;

#define dml18j_1 dml18j_

/* Table of constant values */

static integer c__1 = 1;
static integer c__0 = 0;
static integer c__2 = 2;

/* DECK DBVPOR */
/* Subroutine */ int dbvpor_(doublereal *y, integer *nrowy, integer *ncomp, 
	doublereal *xpts, integer *nxpts, doublereal *a, integer *nrowa, 
	doublereal *alpha, integer *nic, doublereal *b, integer *nrowb, 
	doublereal *beta, integer *nfc, integer *iflag, doublereal *z__, 
	integer *mxnon, doublereal *p, integer *ntp, integer *ip, doublereal *
	w, integer *niv, doublereal *yhp, doublereal *u, doublereal *v, 
	doublereal *coef, doublereal *s, doublereal *stowa, doublereal *g, 
	doublereal *work, integer *iwork, integer *nfcc)
{
    /* System generated locals */
    integer ip_dim1, ip_offset, a_dim1, a_offset, b_dim1, b_offset, p_dim1, 
	    p_offset, u_dim1, u_dim2, u_offset, v_dim1, v_offset, w_dim1, 
	    w_offset, y_dim1, y_offset, yhp_dim1, yhp_offset, i__1, i__2, 
	    i__3, i__4;
    doublereal d__1;
    alist al__1;

    /* Local variables */
    static integer i__, j, k, l, m, n, i1, i2, ic, nn, ira, kod, kwc, kwd, 
	    ndw, non, kws, kwt;
    extern doublereal ddot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    static integer kpts, nfcp1, nfcp2;
    extern /* Subroutine */ int dcoef_(doublereal *, doublereal *, integer *, 
	    integer *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, doublereal *, doublereal *, doublereal *,
	     doublereal *, doublereal *, integer *, integer *, integer *);
    static integer isflg;
    extern /* Subroutine */ int dvecs_(integer *, integer *, doublereal *, 
	    doublereal *, integer *, integer *, integer *);
    static integer ncomp2;
    extern /* Subroutine */ int dstor1_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, integer *, integer *, integer *), 
	    drkfab_(integer *, doublereal *, integer *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, doublereal *,
	     doublereal *, doublereal *, doublereal *, doublereal *, integer *
	    , integer *), dbksol_(integer *, doublereal *, doublereal *), 
	    dlssud_(doublereal *, doublereal *, doublereal *, integer *, 
	    integer *, integer *, doublereal *, integer *, integer *, integer 
	    *, integer *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, doublereal *)
	    , dstway_(doublereal *, doublereal *, doublereal *, integer *, 
	    doublereal *);

    /* Fortran I/O blocks */
    static cilist io___19 = { 0, 0, 0, 0, 0 };
    static cilist io___21 = { 0, 0, 0, 0, 0 };
    static cilist io___25 = { 0, 0, 0, 0, 0 };


/* ***BEGIN PROLOGUE  DBVPOR */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to DBVSUP */
/* ***LIBRARY   SLATEC */
/* ***TYPE      DOUBLE PRECISION (BVPOR-S, DBVPOR-D) */
/* ***AUTHOR  Watts, H. A., (SNLA) */
/* ***DESCRIPTION */

/* ********************************************************************** */
/*     INPUT to DBVPOR    (items not defined in DBVSUP comments) */
/* ********************************************************************** */

/*     NOPG = 0 -- orthonormalization points not pre-assigned */
/*          = 1 -- orthonormalization points pre-assigned */

/*     MXNON = maximum number of orthogonalizations allowed. */

/*     NDISK = 0 -- in-core storage */
/*           = 1 -- disk storage.  Value of NTAPE in data statement */
/*                  is set to 13.  If another value is desired, */
/*                  the data statement must be changed. */

/*     INTEG = type of integrator and associated test to be used */
/*             to determine when to orthonormalize. */

/*             1 -- use GRAM-SCHMIDT test and DDERKF */
/*             2 -- use GRAM-SCHMIDT test and DDEABM */

/*     TOL = tolerance for allowable error in orthogonalization test. */

/*     NPS = 0 normalize particular solution to unit length at each */
/*             point of orthonormalization. */
/*         = 1 do not normalize particular solution. */

/*     NTP = must be .GE. NFC*(NFC+1)/2. */

/*     NFCC = 2*NFC for special treatment of a COMPLEX*16 valued problem */

/*     ICOCO = 0 skip final computations (superposition coefficients */
/*               and, hence, boundary problem solution) */
/*           = 1 calculate superposition coefficients and obtain */
/*               solution to the boundary value problem */

/* ********************************************************************** */
/*     OUTPUT from DBVPOR */
/* ********************************************************************** */

/*     Y(NROWY,NXPTS) = solution at specified output points. */

/*     MXNON = number of orthonormalizations performed by DBVPOR. */

/*     Z(MXNON+1) = locations of orthonormalizations performed by DBVPOR. */

/*     NIV = number of independent vectors returned from DMGSBV. Normally */
/*           this parameter will be meaningful only when DMGSBV returns */
/*           with MFLAG = 2. */

/* ********************************************************************** */

/*     The following variables are in the argument list because of */
/*     variable dimensioning.  In general, they contain no information of */
/*     use to the user.  The amount of storage set aside by the user must */
/*     be greater than or equal to that indicated by the dimension */
/*     statements.  For the disk storage mode, NON = 0 and KPTS = 1, */
/*     while for the in-core storage mode, NON = MXNON and KPTS = NXPTS. */

/*     P(NTP,NON+1) */
/*     IP(NFCC,NON+1) */
/*     YHP(NCOMP,NFC+1)  plus an additional column of the length  NEQIVP */
/*     U(NCOMP,NFC,KPTS) */
/*     V(NCOMP,KPTS) */
/*     W(NFCC,NON+1) */
/*     COEF(NFCC) */
/*     S(NFC+1) */
/*     STOWA(NCOMP*(NFC+1)+NEQIVP+1) */
/*     G(NCOMP) */
/*     WORK(KKKWS) */
/*     IWORK(LLLIWS) */

/* ********************************************************************** */
/*     SUBROUTINES used by DBVPOR */
/*         DLSSUD -- solves an underdetermined system of linear */
/*                   equations.  This routine is used to get a full */
/*                   set of initial conditions for integration. */
/*                   Called by DBVPOR. */

/*         DVECS -- obtains starting vectors for special treatment */
/*                   of COMPLEX*16 valued problems, called by DBVPOR. */

/*         DRKFAB -- routine which conducts integration using DDERKF or */
/*                   DDEABM. */

/*         DSTWAY -- storage for backup capability, called by */
/*                   DBVPOR and DREORT. */

/*         DSTOR1 -- storage at output points, called by DBVPOR, */
/*                   DRKFAB, DREORT and DSTWAY. */

/*         DDOT -- single precision vector inner product routine, */
/*                   called by DBVPOR, DCOEF, DLSSUD, DMGSBV, */
/*                   DBKSOL, DREORT and DPRVEC. */
/*         ** NOTE ** */
/*         a considerable improvement in speed can be achieved if a */
/*         machine language version is used for DDOT. */

/*         DCOEF -- computes the superposition constants from the */
/*                   boundary conditions at XFINAL. */

/*         DBKSOL -- solves an upper triangular set of linear equations. */

/* ********************************************************************** */

/* ***SEE ALSO  DBVSUP */
/* ***ROUTINES CALLED  DBKSOL, DCOEF, DDOT, DLSSUD, DRKFAB, DSTOR1, */
/*                    DSTWAY, DVECS */
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
/* ***END PROLOGUE  DBVPOR */


/*     ****************************************************************** */


/*      ***************************************************************** */

/* ***FIRST EXECUTABLE STATEMENT  DBVPOR */
    /* Parameter adjustments */
    y_dim1 = *nrowy;
    y_offset = 1 + y_dim1;
    y -= y_offset;
    v_dim1 = *ncomp;
    v_offset = 1 + v_dim1;
    v -= v_offset;
    yhp_dim1 = *ncomp;
    yhp_offset = 1 + yhp_dim1;
    yhp -= yhp_offset;
    --xpts;
    a_dim1 = *nrowa;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --alpha;
    b_dim1 = *nrowb;
    b_offset = 1 + b_dim1;
    b -= b_offset;
    --beta;
    u_dim1 = *ncomp;
    u_dim2 = *nfc;
    u_offset = 1 + u_dim1 * (1 + u_dim2);
    u -= u_offset;
    --z__;
    p_dim1 = *ntp;
    p_offset = 1 + p_dim1;
    p -= p_offset;
    --coef;
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
    nfcp1 = *nfc + 1;
    dml18j_1.numort = 0;
    dml8sz_1.c__ = 1.;

/*     ****************************************************************** */
/*         CALCULATE INITIAL CONDITIONS WHICH SATISFY */
/*                       A*YH(XINITIAL)=0  AND  A*YP(XINITIAL)=ALPHA. */
/*         WHEN NFC .NE. NFCC DLSSUD DEFINES VALUES YHP IN A MATRIX OF */
/*         SIZE (NFCC+1)*NCOMP AND ,HENCE, OVERFLOWS THE STORAGE */
/*         ALLOCATION INTO THE U ARRAY. HOWEVER, THIS IS OKAY SINCE */
/*         PLENTY OF SPACE IS AVAILABLE IN U AND IT HAS NOT YET BEEN */
/*         USED. */

    ndw = *nrowa * *ncomp;
    kws = ndw + *nic + 1;
    kwd = kws + *nic;
    kwt = kwd + *nic;
    kwc = kwt + *nic;
    *iflag = 0;
    dlssud_(&a[a_offset], &yhp[(*nfcc + 1) * yhp_dim1 + 1], &alpha[1], nic, 
	    ncomp, nrowa, &yhp[yhp_offset], ncomp, iflag, &c__1, &ira, &c__0, 
	    &work[1], &work[ndw + 1], &iwork[1], &work[kws], &work[kwd], &
	    work[kwt], &isflg, &work[kwc]);
    if (*iflag == 1) {
	goto L10;
    }
    *iflag = -4;
    goto L200;
L10:
    if (*nfc != *nfcc) {
	dvecs_(ncomp, nfc, &yhp[yhp_offset], &work[1], &iwork[1], &
		dml8sz_1.inhomo, iflag);
    }
    if (*iflag == 1) {
	goto L20;
    }
    *iflag = -5;
    goto L190;
L20:

/*           ************************************************************ */
/*               DETERMINE THE NUMBER OF DIFFERENTIAL EQUATIONS TO BE */
/*               INTEGRATED, INITIALIZE VARIABLES FOR AUXILIARY INITIAL */
/*               VALUE PROBLEM AND STORE INITIAL CONDITIONS. */

    dml18j_1.neq = *ncomp * *nfc;
    if (dml8sz_1.inhomo == 1) {
	dml18j_1.neq += *ncomp;
    }
    dml8sz_1.ivp = 0;
    if (dml18j_1.neqivp == 0) {
	goto L40;
    }
    dml8sz_1.ivp = dml18j_1.neq;
    dml18j_1.neq += dml18j_1.neqivp;
    nfcp2 = nfcp1;
    if (dml8sz_1.inhomo == 1) {
	nfcp2 = nfcp1 + 1;
    }
    i__1 = dml18j_1.neqivp;
    for (k = 1; k <= i__1; ++k) {
	yhp[k + nfcp2 * yhp_dim1] = alpha[*nic + k];
/* L30: */
    }
L40:
    dstor1_(&u[u_offset], &yhp[yhp_offset], &v[v_offset], &yhp[nfcp1 * 
	    yhp_dim1 + 1], &c__0, &dml18j_1.ndisk, &dml18j_1.ntape);

/*           ************************************************************ */
/*               SET UP DATA FOR THE ORTHONORMALIZATION TESTING PROCEDURE */
/*               AND SAVE INITIAL CONDITIONS IN CASE A RESTART IS */
/*               NECESSARY. */

    dml15t_1.nswot = 1;
    dml15t_1.knswot = 0;
    dml15t_1.lotjp = 1;
    d__1 = dml18j_1.tol * 10.;
    dml15t_1.tnd = d_lg10(&d__1);
    d__1 = sqrt(dml18j_1.tol);
    dml15t_1.pwcnd = d_lg10(&d__1);
    dml15t_1.x = dml15t_1.xbeg;
    dml15t_1.px = dml15t_1.x;
    dml15t_1.xot = dml15t_1.xend;
    dml15t_1.xop = dml15t_1.x;
    dml15t_1.kop = 1;
    dstway_(&u[u_offset], &v[v_offset], &yhp[yhp_offset], &c__0, &stowa[1]);

/*           ************************************************************ */
/*           ******** FORWARD INTEGRATION OF ALL INITIAL VALUE EQUATIONS */
/*           ********** */
/*           ************************************************************ */

    drkfab_(ncomp, &xpts[1], nxpts, nfc, iflag, &z__[1], mxnon, &p[p_offset], 
	    ntp, &ip[ip_offset], &yhp[yhp_offset], niv, &u[u_offset], &v[
	    v_offset], &w[w_offset], &s[1], &stowa[1], &g[1], &work[1], &
	    iwork[1], nfcc);
    if (*iflag != 0 || dml18j_1.icoco == 0) {
	goto L180;
    }

/*              ********************************************************* */
/*              **************** BACKWARD SWEEP TO OBTAIN SOLUTION */
/*              ******************* */
/*              ********************************************************* */

/*                  CALCULATE SUPERPOSITION COEFFICIENTS AT XFINAL. */

/*                FOR THE DISK STORAGE VERSION, IT IS NOT NECESSARY TO */
/*                READ  U  AND  V AT THE LAST OUTPUT POINT, SINCE THE */
/*                LOCAL COPY OF EACH STILL EXISTS. */

    kod = 1;
    if (dml18j_1.ndisk == 0) {
	kod = *nxpts;
    }
    i1 = *nfcc * *nfcc + 1;
    i2 = i1 + *nfcc;
    dcoef_(&u[(kod * u_dim2 + 1) * u_dim1 + 1], &v[kod * v_dim1 + 1], ncomp, 
	    nrowb, nfc, nic, &b[b_offset], &beta[1], &coef[1], &
	    dml8sz_1.inhomo, &dml18j_1.re, &dml18j_1.ae, &work[1], &work[i1], 
	    &work[i2], &iwork[1], iflag, nfcc);

/*              ********************************************************* */
/*                  CALCULATE SOLUTION AT OUTPUT POINTS BY RECURRING */
/*                  BACKWARDS.  AS WE RECUR BACKWARDS FROM XFINAL TO */
/*                  XINITIAL WE MUST CALCULATE NEW SUPERPOSITION */
/*                  COEFFICIENTS EACH TIME WE CROSS A POINT OF */
/*                  ORTHONORMALIZATION. */

    k = dml18j_1.numort;
    ncomp2 = *ncomp / 2;
    ic = 1;
    if (*nfc != *nfcc) {
	ic = 2;
    }
    i__1 = *nxpts;
    for (j = 1; j <= i__1; ++j) {
	kpts = *nxpts - j + 1;
	kod = kpts;
	if (dml18j_1.ndisk == 1) {
	    kod = 1;
	}
L50:
/*                 ...EXIT */
	if (k == 0) {
	    goto L120;
	}
/*                 ...EXIT */
	if (dml15t_1.xend > dml15t_1.xbeg && xpts[kpts] >= z__[k]) {
	    goto L120;
	}
/*                 ...EXIT */
	if (dml15t_1.xend < dml15t_1.xbeg && xpts[kpts] <= z__[k]) {
	    goto L120;
	}
	non = k;
	if (dml18j_1.ndisk == 0) {
	    goto L60;
	}
	non = 1;
	al__1.aerr = 0;
	al__1.aunit = dml18j_1.ntape;
	f_back(&al__1);
	io___19.ciunit = dml18j_1.ntape;
	s_rsue(&io___19);
	i__2 = *nfcc;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    do_uio(&c__1, (char *)&ip[i__ + ip_dim1], (ftnlen)sizeof(integer))
		    ;
	}
	i__3 = *ntp;
	for (i__ = 1; i__ <= i__3; ++i__) {
	    do_uio(&c__1, (char *)&p[i__ + p_dim1], (ftnlen)sizeof(doublereal)
		    );
	}
	e_rsue();
	al__1.aerr = 0;
	al__1.aunit = dml18j_1.ntape;
	f_back(&al__1);
L60:
	if (dml8sz_1.inhomo != 1) {
	    goto L90;
	}
	if (dml18j_1.ndisk == 0) {
	    goto L70;
	}
	al__1.aerr = 0;
	al__1.aunit = dml18j_1.ntape;
	f_back(&al__1);
	io___21.ciunit = dml18j_1.ntape;
	s_rsue(&io___21);
	i__2 = *nfcc;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    do_uio(&c__1, (char *)&w[i__ + w_dim1], (ftnlen)sizeof(doublereal)
		    );
	}
	e_rsue();
	al__1.aerr = 0;
	al__1.aunit = dml18j_1.ntape;
	f_back(&al__1);
L70:
	i__2 = *nfcc;
	for (n = 1; n <= i__2; ++n) {
	    coef[n] -= w[n + non * w_dim1];
/* L80: */
	}
L90:
	dbksol_(nfcc, &p[non * p_dim1 + 1], &coef[1]);
	i__2 = *nfcc;
	for (m = 1; m <= i__2; ++m) {
	    work[m] = coef[m];
/* L100: */
	}
	i__2 = *nfcc;
	for (m = 1; m <= i__2; ++m) {
	    l = ip[m + non * ip_dim1];
	    coef[l] = work[m];
/* L110: */
	}
	--k;
	goto L50;
L120:
	if (dml18j_1.ndisk == 0) {
	    goto L130;
	}
	al__1.aerr = 0;
	al__1.aunit = dml18j_1.ntape;
	f_back(&al__1);
	io___25.ciunit = dml18j_1.ntape;
	s_rsue(&io___25);
	i__2 = *ncomp;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    do_uio(&c__1, (char *)&v[i__ + v_dim1], (ftnlen)sizeof(doublereal)
		    );
	}
	i__3 = *nfc;
	for (m = 1; m <= i__3; ++m) {
	    i__4 = *ncomp;
	    for (i__ = 1; i__ <= i__4; ++i__) {
		do_uio(&c__1, (char *)&u[i__ + (m + u_dim2) * u_dim1], (
			ftnlen)sizeof(doublereal));
	    }
	}
	e_rsue();
	al__1.aerr = 0;
	al__1.aunit = dml18j_1.ntape;
	f_back(&al__1);
L130:
	i__2 = *ncomp;
	for (n = 1; n <= i__2; ++n) {
	    y[n + kpts * y_dim1] = v[n + kod * v_dim1] + ddot_(nfc, &u[n + (
		    kod * u_dim2 + 1) * u_dim1], ncomp, &coef[1], &ic);
/* L140: */
	}
	if (*nfc == *nfcc) {
	    goto L160;
	}
	i__2 = ncomp2;
	for (n = 1; n <= i__2; ++n) {
	    nn = ncomp2 + n;
	    y[n + kpts * y_dim1] -= ddot_(nfc, &u[nn + (kod * u_dim2 + 1) * 
		    u_dim1], ncomp, &coef[2], &c__2);
	    y[nn + kpts * y_dim1] += ddot_(nfc, &u[n + (kod * u_dim2 + 1) * 
		    u_dim1], ncomp, &coef[2], &c__2);
/* L150: */
	}
L160:
/* L170: */
	;
    }
L180:
L190:
L200:

/*     ****************************************************************** */

    *mxnon = dml18j_1.numort;
    return 0;
} /* dbvpor_ */

