/* dexbvp.f -- translated by f2c (version 12.02.01).
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
    integer igofx, inhomo, ivp, ncomp, nfc;
} dml8sz_;

#define dml8sz_1 dml8sz_

struct {
    doublereal ae, re, tol;
    integer nxpts, nic, nopg, mxnon, ndisk, ntape, neq, indpvt, integ, nps, 
	    ntp, neqivp, numort, nfcc, icoco;
} dml18j_;

#define dml18j_1 dml18j_

struct {
    doublereal px, pwcnd, tnd, x, xbeg, xend, xot, xop;
    integer info[15], istkop, knswot, kop, lotjp, mnswot, nswot;
} dml15t_;

#define dml15t_1 dml15t_

struct {
    integer kkkzpw, needw, neediw, k1, k2, k3, k4, k5, k6, k7, k8, k9, k10, 
	    k11, l1, l2, kkkint, lllint;
} dml17b_;

#define dml17b_1 dml17b_

struct {
    doublereal uro, sru, eps, sqovfl, twou, fouru;
    integer lpar;
} dml5mc_;

#define dml5mc_1 dml5mc_

/* Table of constant values */

static doublereal c_b3 = 10.;
static integer c__1 = 1;
static integer c__0 = 0;
static integer c__4 = 4;

/* DECK DEXBVP */
/* Subroutine */ int dexbvp_(doublereal *y, integer *nrowy, doublereal *xpts, 
	doublereal *a, integer *nrowa, doublereal *alpha, doublereal *b, 
	integer *nrowb, doublereal *beta, integer *iflag, doublereal *work, 
	integer *iwork)
{
    /* System generated locals */
    address a__1[4];
    integer a_dim1, a_offset, b_dim1, b_offset, y_dim1, y_offset, i__1, i__2[
	    4];
    doublereal d__1;
    char ch__1[124];

    /* Local variables */
    static doublereal xl;
    static integer inc, kotc, iexp;
    static char xern1[8], xern2[8];
    static integer nsafw;
    static doublereal zquit;
    static integer nsafiw;
    extern /* Subroutine */ int dbvpor_(doublereal *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, integer *, integer *), xermsg_(char *,
	     char *, char *, integer *, integer *, ftnlen, ftnlen, ftnlen);

    /* Fortran I/O blocks */
    static icilist io___9 = { 0, xern1, 0, "(I8)", 8, 1 };
    static icilist io___11 = { 0, xern2, 0, "(I8)", 8, 1 };


/* ***BEGIN PROLOGUE  DEXBVP */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to DBVSUP */
/* ***LIBRARY   SLATEC */
/* ***TYPE      DOUBLE PRECISION (EXBVP-S, DEXBVP-D) */
/* ***AUTHOR  Watts, H. A., (SNLA) */
/* ***DESCRIPTION */

/*  This subroutine is used to execute the basic technique for solving */
/*  the two-point boundary value problem. */

/* ***SEE ALSO  DBVSUP */
/* ***ROUTINES CALLED  DBVPOR, XERMSG */
/* ***COMMON BLOCKS    DML15T, DML17B, DML18J, DML5MC, DML8SZ */
/* ***REVISION HISTORY  (YYMMDD) */
/*   750601  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890831  Modified array declarations.  (WRB) */
/*   890911  Removed unnecessary intrinsics.  (WRB) */
/*   890921  Realigned order of variables in certain COMMON blocks. */
/*           (WRB) */
/*   890921  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900328  Added TYPE section.  (WRB) */
/*   900510  Convert XERRWV calls to XERMSG calls.  (RWC) */
/*   910722  Updated AUTHOR section.  (ALS) */
/* ***END PROLOGUE  DEXBVP */


/*     ****************************************************************** */



/* ***FIRST EXECUTABLE STATEMENT  DEXBVP */
    /* Parameter adjustments */
    y_dim1 = *nrowy;
    y_offset = 1 + y_dim1;
    y -= y_offset;
    --xpts;
    a_dim1 = *nrowa;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --alpha;
    b_dim1 = *nrowb;
    b_offset = 1 + b_dim1;
    b -= b_offset;
    --beta;
    --work;
    --iwork;

    /* Function Body */
    kotc = 1;
    iexp = 0;
    if (iwork[7] == -1) {
	iexp = iwork[8];
    }

/*     COMPUTE ORTHONORMALIZATION TOLERANCES. */

L10:
    i__1 = -dml5mc_1.lpar - iexp << 1;
    dml18j_1.tol = pow_di(&c_b3, &i__1);

    iwork[8] = iexp;
    dml18j_1.mxnon = iwork[2];

/* ********************************************************************** */
/* ********************************************************************** */

    dbvpor_(&y[y_offset], nrowy, &dml8sz_1.ncomp, &xpts[1], &dml18j_1.nxpts, &
	    a[a_offset], nrowa, &alpha[1], &dml18j_1.nic, &b[b_offset], nrowb,
	     &beta[1], &dml8sz_1.nfc, iflag, &work[1], &dml18j_1.mxnon, &work[
	    dml17b_1.k1], &dml18j_1.ntp, &iwork[18], &work[dml17b_1.k2], &
	    iwork[16], &work[dml17b_1.k3], &work[dml17b_1.k4], &work[
	    dml17b_1.k5], &work[dml17b_1.k6], &work[dml17b_1.k7], &work[
	    dml17b_1.k8], &work[dml17b_1.k9], &work[dml17b_1.k10], &iwork[
	    dml17b_1.l1], &dml18j_1.nfcc);

/* ********************************************************************** */
/* ********************************************************************** */
/*     IF DMGSBV RETURNS WITH MESSAGE OF DEPENDENT VECTORS, WE REDUCE */
/*     ORTHONORMALIZATION TOLERANCE AND TRY AGAIN. THIS IS DONE */
/*     A MAXIMUM OF 2 TIMES. */

    if (*iflag != 30) {
	goto L20;
    }
    if (kotc == 3 || dml18j_1.nopg == 1) {
	goto L30;
    }
    ++kotc;
    iexp += -2;
    goto L10;

/* ********************************************************************** */
/*     IF DBVPOR RETURNS MESSAGE THAT THE MAXIMUM NUMBER OF */
/*     ORTHONORMALIZATIONS HAS BEEN ATTAINED AND WE CANNOT CONTINUE, THEN */
/*     WE ESTIMATE THE NEW STORAGE REQUIREMENTS IN ORDER TO SOLVE PROBLEM */

L20:
    if (*iflag != 13) {
	goto L30;
    }
    xl = (d__1 = dml15t_1.xend - dml15t_1.xbeg, abs(d__1));
    zquit = (d__1 = dml15t_1.x - dml15t_1.xbeg, abs(d__1));
    inc = (integer) (xl * 1.5 / zquit * (dml18j_1.mxnon + 1));
    if (dml18j_1.ndisk != 1) {
	nsafw = inc * dml17b_1.kkkzpw + dml17b_1.needw;
	nsafiw = inc * dml18j_1.nfcc + dml17b_1.neediw;
    } else {
	nsafw = dml17b_1.needw + inc;
	nsafiw = dml17b_1.neediw;
    }

    s_wsfi(&io___9);
    do_fio(&c__1, (char *)&nsafw, (ftnlen)sizeof(integer));
    e_wsfi();
    s_wsfi(&io___11);
    do_fio(&c__1, (char *)&nsafiw, (ftnlen)sizeof(integer));
    e_wsfi();
/* Writing concatenation */
    i__2[0] = 58, a__1[0] = "IN DBVSUP, PREDICTED STORAGE ALLOCATION FOR WOR"
	    "K ARRAY IS ";
    i__2[1] = 8, a__1[1] = xern1;
    i__2[2] = 50, a__1[2] = ", PREDICTED STORAGE ALLOCATION FOR IWORK ARRAY "
	    "IS ";
    i__2[3] = 8, a__1[3] = xern2;
    s_cat(ch__1, a__1, i__2, &c__4, (ftnlen)124);
    xermsg_("SLATEC", "DEXBVP", ch__1, &c__1, &c__0, (ftnlen)6, (ftnlen)6, (
	    ftnlen)124);

L30:
    iwork[1] = dml18j_1.mxnon;
    return 0;
} /* dexbvp_ */

