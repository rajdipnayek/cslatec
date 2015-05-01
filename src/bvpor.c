/* bvpor.f -- translated by f2c (version 12.02.01).
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
    integer nxptsd, nicd, nopg, mxnond, ndisk, ntape, neq, indpvt, integ, nps,
	     ntpd, neqivp, numort, nfccd, icoco;
} ml18jr_;

#define ml18jr_1 ml18jr_

/* Table of constant values */

static integer c__1 = 1;
static integer c__0 = 0;
static integer c__2 = 2;

/* DECK BVPOR */
/* Subroutine */ int bvpor_(real *y, integer *nrowy, integer *ncomp, real *
	xpts, integer *nxpts, real *a, integer *nrowa, real *alpha, integer *
	nic, real *b, integer *nrowb, real *beta, integer *nfc, integer *
	iflag, real *z__, integer *mxnon, real *p, integer *ntp, integer *ip, 
	real *w, integer *niv, real *yhp, real *u, real *v, real *coef, real *
	s, real *stowa, real *g, real *work, integer *iwork, integer *nfcc)
{
    /* System generated locals */
    integer y_dim1, y_offset, a_dim1, a_offset, b_dim1, b_offset, p_dim1, 
	    p_offset, ip_dim1, ip_offset, u_dim1, u_dim2, u_offset, v_dim1, 
	    v_offset, w_dim1, w_offset, yhp_dim1, yhp_offset, i__1, i__2, 
	    i__3, i__4;
    real r__1;
    alist al__1;

    /* Local variables */
    static integer i__, j, k, l, m, n, i1, i2, ic, nn, ira, kod, kwc, kwd, 
	    ndw, non, kws, kwt;
    extern doublereal sdot_(integer *, real *, integer *, real *, integer *);
    static integer kpts, nfcp1, nfcp2;
    extern /* Subroutine */ int stor1_(real *, real *, real *, real *, 
	    integer *, integer *, integer *), rkfab_(integer *, real *, 
	    integer *, integer *, integer *, real *, integer *, real *, 
	    integer *, integer *, real *, integer *, real *, real *, real *, 
	    real *, real *, real *, real *, integer *, integer *), scoef_(
	    real *, real *, integer *, integer *, integer *, integer *, real *
	    , real *, real *, integer *, real *, real *, real *, real *, real 
	    *, integer *, integer *, integer *);
    static integer isflg;
    extern /* Subroutine */ int bksol_(integer *, real *, real *), svecs_(
	    integer *, integer *, real *, real *, integer *, integer *, 
	    integer *), stway_(real *, real *, real *, integer *, real *);
    static integer ncomp2;
    extern /* Subroutine */ int lssuds_(real *, real *, real *, integer *, 
	    integer *, integer *, real *, integer *, integer *, integer *, 
	    integer *, integer *, real *, real *, integer *, real *, real *, 
	    real *, integer *, real *);

    /* Fortran I/O blocks */
    static cilist io___19 = { 0, 0, 0, 0, 0 };
    static cilist io___21 = { 0, 0, 0, 0, 0 };
    static cilist io___25 = { 0, 0, 0, 0, 0 };


/* ***BEGIN PROLOGUE  BVPOR */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to BVSUP */
/* ***LIBRARY   SLATEC */
/* ***TYPE      SINGLE PRECISION (BVPOR-S, DBVPOR-D) */
/* ***AUTHOR  Watts, H. A., (SNLA) */
/* ***DESCRIPTION */

/* ********************************************************************** */
/*     INPUT to BVPOR    (items not defined in BVSUP comments) */
/* ********************************************************************** */

/*     NOPG = 0 -- Orthonormalization points not pre-assigned */
/*          = 1 -- Orthonormalization points pre-assigned */

/*     MXNON = Maximum number of orthogonalizations allowed. */

/*     NDISK = 0 -- IN-CORE storage */
/*           = 1 -- DISK storage.  Value of NTAPE in data statement */
/*                  is set to 13.  If another value is desired, */
/*                  the data statement must be changed. */

/*     INTEG = Type of integrator and associated test to be used */
/*             to determine when to orthonormalize. */

/*             1 -- Use GRAM-SCHMIDT test and DERKF */
/*             2 -- Use GRAM-SCHMIDT test and DEABM */

/*     TOL = Tolerance for allowable error in orthogonalization test. */

/*     NPS = 0 Normalize particular solution to unit length at each */
/*             point of orthonormalization. */
/*         = 1 Do not normalize particular solution. */

/*     NTP = Must be .GE. NFC*(NFC+1)/2. */


/*     NFCC = 2*NFC for special treatment of a complex valued problem */

/*     ICOCO = 0 Skip final computations (superposition coefficients */
/*               and ,hence, boundary problem solution) */
/*           = 1 Calculate superposition coefficients and obtain */
/*               solution to the boundary value problem */

/* ********************************************************************** */
/*     OUTPUT from BVPOR */
/* ********************************************************************** */

/*     Y(NROWY,NXPTS) = Solution at specified output points. */

/*     MXNON = Number of orthonormalizations performed by BVPOR. */

/*     Z(MXNON+1) = Locations of orthonormalizations performed by BVPOR. */

/*     NIV = Number of independent vectors returned from MGSBV. Normally */
/*        this parameter will be meaningful only when MGSBV returns with */
/*           MFLAG = 2. */

/* ********************************************************************** */

/*     The following variables are in the argument list because of */
/*     variable dimensioning. In general, they contain no information of */
/*     use to the user.  The amount of storage set aside by the user must */
/*     be greater than or equal to that indicated by the dimension */
/*     statements.   For the DISK storage mode, NON = 0 and KPTS = 1, */
/*     while for the IN-CORE storage mode, NON = MXNON and KPTS = NXPTS. */

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
/*     Subroutines used by BVPOR */
/*         LSSUDS -- Solves an underdetermined system of linear */
/*                   equations.  This routine is used to get a full */
/*                   set of initial conditions for integration. */
/*                   Called by BVPOR */

/*         SVECS -- Obtains starting vectors for special treatment */
/*                  of complex valued problems , called by BVPOR */

/*         RKFAB -- Routine which conducts integration using DERKF or */
/*                   DEABM */

/*         STWAY -- Storage for backup capability, called by */
/*                   BVPOR and REORT */

/*         STOR1 -- Storage at output points, called by BVPOR, */
/*                  RKFAB, REORT and STWAY. */

/*         SDOT -- Single precision vector inner product routine, */
/*                   called by BVPOR, SCOEF, LSSUDS, MGSBV, */
/*                   BKSOL, REORT and PRVEC. */
/*         ** NOTE ** */
/*         A considerable improvement in speed can be achieved if a */
/*         machine language version is used for SDOT. */

/*         SCOEF -- Computes the superposition constants from the */
/*                  boundary conditions at Xfinal. */

/*         BKSOL -- Solves an upper triangular set of linear equations. */

/* ********************************************************************** */

/* ***SEE ALSO  BVSUP */
/* ***ROUTINES CALLED  BKSOL, LSSUDS, RKFAB, SCOEF, SDOT, STOR1, STWAY, */
/*                    SVECS */
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
/* ***END PROLOGUE  BVPOR */


/* ********************************************************************** */


/* ********************************************************************** */

/* ***FIRST EXECUTABLE STATEMENT  BVPOR */
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
    ml18jr_1.numort = 0;
    ml8sz_1.c__ = 1.f;

/* ********************************************************************** */
/*     CALCULATE INITIAL CONDITIONS WHICH SATISFY */
/*                   A*YH(XINITIAL)=0  AND  A*YP(XINITIAL)=ALPHA. */
/*     WHEN NFC .NE. NFCC LSSUDS DEFINES VALUES YHP IN A MATRIX OF SIZE */
/*     (NFCC+1)*NCOMP AND ,HENCE, OVERFLOWS THE STORAGE ALLOCATION INTO */
/*     THE U ARRAY. HOWEVER, THIS IS OKAY SINCE PLENTY OF SPACE IS */
/*     AVAILABLE IN U AND IT HAS NOT YET BEEN USED. */

    ndw = *nrowa * *ncomp;
    kws = ndw + *nic + 1;
    kwd = kws + *nic;
    kwt = kwd + *nic;
    kwc = kwt + *nic;
    *iflag = 0;
    lssuds_(&a[a_offset], &yhp[(*nfcc + 1) * yhp_dim1 + 1], &alpha[1], nic, 
	    ncomp, nrowa, &yhp[yhp_offset], ncomp, iflag, &c__1, &ira, &c__0, 
	    &work[1], &work[ndw + 1], &iwork[1], &work[kws], &work[kwd], &
	    work[kwt], &isflg, &work[kwc]);
    if (*iflag == 1) {
	goto L3;
    }
    *iflag = -4;
    goto L250;
L3:
    if (*nfc != *nfcc) {
	svecs_(ncomp, nfc, &yhp[yhp_offset], &work[1], &iwork[1], &
		ml8sz_1.inhomo, iflag);
    }
    if (*iflag == 1) {
	goto L5;
    }
    *iflag = -5;
    goto L250;

/* ********************************************************************** */
/*     DETERMINE THE NUMBER OF DIFFERENTIAL EQUATIONS TO BE INTEGRATED, */
/*     INITIALIZE VARIABLES FOR AUXILIARY INITIAL VALUE PROBLEM AND */
/*     STORE INITIAL CONDITIONS. */

L5:
    ml18jr_1.neq = *ncomp * *nfc;
    if (ml8sz_1.inhomo == 1) {
	ml18jr_1.neq += *ncomp;
    }
    ml8sz_1.ivp = 0;
    if (ml18jr_1.neqivp == 0) {
	goto L10;
    }
    ml8sz_1.ivp = ml18jr_1.neq;
    ml18jr_1.neq += ml18jr_1.neqivp;
    nfcp2 = nfcp1;
    if (ml8sz_1.inhomo == 1) {
	nfcp2 = nfcp1 + 1;
    }
    i__1 = ml18jr_1.neqivp;
    for (k = 1; k <= i__1; ++k) {
/* L7: */
	yhp[k + nfcp2 * yhp_dim1] = alpha[*nic + k];
    }
L10:
    stor1_(&u[u_offset], &yhp[yhp_offset], &v[v_offset], &yhp[nfcp1 * 
	    yhp_dim1 + 1], &c__0, &ml18jr_1.ndisk, &ml18jr_1.ntape);

/* ********************************************************************** */
/*     SET UP DATA FOR THE ORTHONORMALIZATION TESTING PROCEDURE AND */
/*     SAVE INITIAL CONDITIONS IN CASE A RESTART IS NECESSARY. */

    ml15to_1.nswot = 1;
    ml15to_1.knswot = 0;
    ml15to_1.lotjp = 1;
    r__1 = ml18jr_1.tol * 10.f;
    ml15to_1.tnd = r_lg10(&r__1);
    r__1 = sqrt(ml18jr_1.tol);
    ml15to_1.pwcnd = r_lg10(&r__1);
    ml15to_1.x = ml15to_1.xbeg;
    ml15to_1.px = ml15to_1.x;
    ml15to_1.xot = ml15to_1.xend;
    ml15to_1.xop = ml15to_1.x;
    ml15to_1.kop = 1;
    stway_(&u[u_offset], &v[v_offset], &yhp[yhp_offset], &c__0, &stowa[1]);

/* ********************************************************************** */
/* ******** FORWARD INTEGRATION OF ALL INITIAL VALUE EQUATIONS ********** */
/* ********************************************************************** */

    rkfab_(ncomp, &xpts[1], nxpts, nfc, iflag, &z__[1], mxnon, &p[p_offset], 
	    ntp, &ip[ip_offset], &yhp[yhp_offset], niv, &u[u_offset], &v[
	    v_offset], &w[w_offset], &s[1], &stowa[1], &g[1], &work[1], &
	    iwork[1], nfcc);
    if (*iflag != 0 || ml18jr_1.icoco == 0) {
	goto L250;
    }

/* ********************************************************************** */
/* **************** BACKWARD SWEEP TO OBTAIN SOLUTION ******************* */
/* ********************************************************************** */

/*     CALCULATE SUPERPOSITION COEFFICIENTS AT XFINAL. */

/*   FOR THE DISK STORAGE VERSION, IT IS NOT NECESSARY TO READ  U  AND  V */
/*   AT THE LAST OUTPUT POINT, SINCE THE LOCAL COPY OF EACH STILL EXISTS. */

    kod = 1;
    if (ml18jr_1.ndisk == 0) {
	kod = *nxpts;
    }
    i1 = *nfcc * *nfcc + 1;
    i2 = i1 + *nfcc;
    scoef_(&u[(kod * u_dim2 + 1) * u_dim1 + 1], &v[kod * v_dim1 + 1], ncomp, 
	    nrowb, nfc, nic, &b[b_offset], &beta[1], &coef[1], &
	    ml8sz_1.inhomo, &ml18jr_1.re, &ml18jr_1.ae, &work[1], &work[i1], &
	    work[i2], &iwork[1], iflag, nfcc);

/* ********************************************************************** */
/*     CALCULATE SOLUTION AT OUTPUT POINTS BY RECURRING BACKWARDS. */
/*     AS WE RECUR BACKWARDS FROM XFINAL TO XINITIAL WE MUST CALCULATE */
/*     NEW SUPERPOSITION COEFFICIENTS EACH TIME WE CROSS A POINT OF */
/*     ORTHONORMALIZATION. */

    k = ml18jr_1.numort;
    ncomp2 = *ncomp / 2;
    ic = 1;
    if (*nfc != *nfcc) {
	ic = 2;
    }
    i__1 = *nxpts;
    for (j = 1; j <= i__1; ++j) {
	kpts = *nxpts - j + 1;
	kod = kpts;
	if (ml18jr_1.ndisk == 1) {
	    kod = 1;
	}
L135:
	if (k == 0) {
	    goto L170;
	}
	if (ml15to_1.xend > ml15to_1.xbeg && xpts[kpts] >= z__[k]) {
	    goto L170;
	}
	if (ml15to_1.xend < ml15to_1.xbeg && xpts[kpts] <= z__[k]) {
	    goto L170;
	}
	non = k;
	if (ml18jr_1.ndisk == 0) {
	    goto L136;
	}
	non = 1;
	al__1.aerr = 0;
	al__1.aunit = ml18jr_1.ntape;
	f_back(&al__1);
	io___19.ciunit = ml18jr_1.ntape;
	s_rsue(&io___19);
	i__2 = *nfcc;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    do_uio(&c__1, (char *)&ip[i__ + ip_dim1], (ftnlen)sizeof(integer))
		    ;
	}
	i__3 = *ntp;
	for (i__ = 1; i__ <= i__3; ++i__) {
	    do_uio(&c__1, (char *)&p[i__ + p_dim1], (ftnlen)sizeof(real));
	}
	e_rsue();
	al__1.aerr = 0;
	al__1.aunit = ml18jr_1.ntape;
	f_back(&al__1);
L136:
	if (ml8sz_1.inhomo != 1) {
	    goto L150;
	}
	if (ml18jr_1.ndisk == 0) {
	    goto L138;
	}
	al__1.aerr = 0;
	al__1.aunit = ml18jr_1.ntape;
	f_back(&al__1);
	io___21.ciunit = ml18jr_1.ntape;
	s_rsue(&io___21);
	i__2 = *nfcc;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    do_uio(&c__1, (char *)&w[i__ + w_dim1], (ftnlen)sizeof(real));
	}
	e_rsue();
	al__1.aerr = 0;
	al__1.aunit = ml18jr_1.ntape;
	f_back(&al__1);
L138:
	i__2 = *nfcc;
	for (n = 1; n <= i__2; ++n) {
/* L140: */
	    coef[n] -= w[n + non * w_dim1];
	}
L150:
	bksol_(nfcc, &p[non * p_dim1 + 1], &coef[1]);
	i__2 = *nfcc;
	for (m = 1; m <= i__2; ++m) {
/* L155: */
	    work[m] = coef[m];
	}
	i__2 = *nfcc;
	for (m = 1; m <= i__2; ++m) {
	    l = ip[m + non * ip_dim1];
/* L160: */
	    coef[l] = work[m];
	}
	--k;
	goto L135;
L170:
	if (ml18jr_1.ndisk == 0) {
	    goto L175;
	}
	al__1.aerr = 0;
	al__1.aunit = ml18jr_1.ntape;
	f_back(&al__1);
	io___25.ciunit = ml18jr_1.ntape;
	s_rsue(&io___25);
	i__2 = *ncomp;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    do_uio(&c__1, (char *)&v[i__ + v_dim1], (ftnlen)sizeof(real));
	}
	i__3 = *nfc;
	for (m = 1; m <= i__3; ++m) {
	    i__4 = *ncomp;
	    for (i__ = 1; i__ <= i__4; ++i__) {
		do_uio(&c__1, (char *)&u[i__ + (m + u_dim2) * u_dim1], (
			ftnlen)sizeof(real));
	    }
	}
	e_rsue();
	al__1.aerr = 0;
	al__1.aunit = ml18jr_1.ntape;
	f_back(&al__1);
L175:
	i__2 = *ncomp;
	for (n = 1; n <= i__2; ++n) {
/* L180: */
	    y[n + kpts * y_dim1] = v[n + kod * v_dim1] + sdot_(nfc, &u[n + (
		    kod * u_dim2 + 1) * u_dim1], ncomp, &coef[1], &ic);
	}
	if (*nfc == *nfcc) {
	    goto L200;
	}
	i__2 = ncomp2;
	for (n = 1; n <= i__2; ++n) {
	    nn = ncomp2 + n;
	    y[n + kpts * y_dim1] -= sdot_(nfc, &u[nn + (kod * u_dim2 + 1) * 
		    u_dim1], ncomp, &coef[2], &c__2);
/* L190: */
	    y[nn + kpts * y_dim1] += sdot_(nfc, &u[n + (kod * u_dim2 + 1) * 
		    u_dim1], ncomp, &coef[2], &c__2);
	}
L200:
	;
    }

/* ********************************************************************** */

L250:
    *mxnon = ml18jr_1.numort;
    return 0;
} /* bvpor_ */

