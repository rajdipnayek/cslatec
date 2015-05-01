/* exbvp.f -- translated by f2c (version 12.02.01).
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
    integer igofx, inhomo, ivp, ncomp, nfc;
} ml8sz_;

#define ml8sz_1 ml8sz_

struct {
    real ae, re, tol;
    integer nxpts, nic, nopg, mxnon, ndisk, ntape, neq, indpvt, integ, nps, 
	    ntp, neqivp, numort, nfcc, icoco;
} ml18jr_;

#define ml18jr_1 ml18jr_

struct {
    real px, pwcnd, tnd, x, xbeg, xend, xot, xop;
    integer info[15], istkop, knswot, kop, lotjp, mnswot, nswot;
} ml15to_;

#define ml15to_1 ml15to_

struct {
    integer kkkzpw, needw, neediw, k1, k2, k3, k4, k5, k6, k7, k8, k9, k10, 
	    k11, l1, l2, kkkint, lllint;
} ml17bw_;

#define ml17bw_1 ml17bw_

struct {
    real uro, sru, eps, sqovfl, twou, fouru;
    integer lpar;
} ml5mco_;

#define ml5mco_1 ml5mco_

/* Table of constant values */

static real c_b3 = 10.f;
static integer c__1 = 1;
static integer c__0 = 0;
static integer c__4 = 4;

/* DECK EXBVP */
/* Subroutine */ int exbvp_(real *y, integer *nrowy, real *xpts, real *a, 
	integer *nrowa, real *alpha, real *b, integer *nrowb, real *beta, 
	integer *iflag, real *work, integer *iwork)
{
    /* System generated locals */
    address a__1[4];
    integer y_dim1, y_offset, a_dim1, a_offset, b_dim1, b_offset, i__1, i__2[
	    4];
    real r__1;
    char ch__1[123];

    /* Local variables */
    static real xl;
    static integer inc, kotc, iexp;
    static char xern1[8], xern2[8];
    static integer nsafw;
    extern /* Subroutine */ int bvpor_(real *, integer *, integer *, real *, 
	    integer *, real *, integer *, real *, integer *, real *, integer *
	    , real *, integer *, integer *, real *, integer *, real *, 
	    integer *, integer *, real *, integer *, real *, real *, real *, 
	    real *, real *, real *, real *, real *, integer *, integer *);
    static real zquit;
    static integer nsafiw;
    extern /* Subroutine */ int xermsg_(char *, char *, char *, integer *, 
	    integer *, ftnlen, ftnlen, ftnlen);

    /* Fortran I/O blocks */
    static icilist io___9 = { 0, xern1, 0, "(I8)", 8, 1 };
    static icilist io___11 = { 0, xern2, 0, "(I8)", 8, 1 };


/* ***BEGIN PROLOGUE  EXBVP */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to BVSUP */
/* ***LIBRARY   SLATEC */
/* ***TYPE      SINGLE PRECISION (EXBVP-S, DEXBVP-D) */
/* ***AUTHOR  Watts, H. A., (SNLA) */
/* ***DESCRIPTION */

/*  This subroutine is used to execute the basic technique for solving */
/*  the two-point boundary value problem */

/* ***SEE ALSO  BVSUP */
/* ***ROUTINES CALLED  BVPOR, XERMSG */
/* ***COMMON BLOCKS    ML15TO, ML17BW, ML18JR, ML5MCO, ML8SZ */
/* ***REVISION HISTORY  (YYMMDD) */
/*   750601  DATE WRITTEN */
/*   890921  Realigned order of variables in certain COMMON blocks. */
/*           (WRB) */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900328  Added TYPE section.  (WRB) */
/*   900510  Convert XERRWV calls to XERMSG calls.  (RWC) */
/*   910722  Updated AUTHOR section.  (ALS) */
/* ***END PROLOGUE  EXBVP */


/*     **************************************************************** */



/* ***FIRST EXECUTABLE STATEMENT  EXBVP */
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
    i__1 = -ml5mco_1.lpar - iexp << 1;
    ml18jr_1.tol = pow_ri(&c_b3, &i__1);

    iwork[8] = iexp;
    ml18jr_1.mxnon = iwork[2];

/* ********************************************************************** */
/* ********************************************************************** */

    bvpor_(&y[y_offset], nrowy, &ml8sz_1.ncomp, &xpts[1], &ml18jr_1.nxpts, &a[
	    a_offset], nrowa, &alpha[1], &ml18jr_1.nic, &b[b_offset], nrowb, &
	    beta[1], &ml8sz_1.nfc, iflag, &work[1], &ml18jr_1.mxnon, &work[
	    ml17bw_1.k1], &ml18jr_1.ntp, &iwork[18], &work[ml17bw_1.k2], &
	    iwork[16], &work[ml17bw_1.k3], &work[ml17bw_1.k4], &work[
	    ml17bw_1.k5], &work[ml17bw_1.k6], &work[ml17bw_1.k7], &work[
	    ml17bw_1.k8], &work[ml17bw_1.k9], &work[ml17bw_1.k10], &iwork[
	    ml17bw_1.l1], &ml18jr_1.nfcc);

/* ********************************************************************** */
/* ********************************************************************** */
/*     IF MGSBV RETURNS WITH MESSAGE OF DEPENDENT VECTORS, WE REDUCE */
/*     ORTHONORMALIZATION TOLERANCE AND TRY AGAIN. THIS IS DONE */
/*     A MAXIMUM OF 2 TIMES. */

    if (*iflag != 30) {
	goto L20;
    }
    if (kotc == 3 || ml18jr_1.nopg == 1) {
	goto L30;
    }
    ++kotc;
    iexp += -2;
    goto L10;

/* ********************************************************************** */
/*     IF BVPOR RETURNS MESSAGE THAT THE MAXIMUM NUMBER OF */
/*     ORTHONORMALIZATIONS HAS BEEN ATTAINED AND WE CANNOT CONTINUE, THEN */
/*     WE ESTIMATE THE NEW STORAGE REQUIREMENTS IN ORDER TO SOLVE PROBLEM */

L20:
    if (*iflag != 13) {
	goto L30;
    }
    xl = (r__1 = ml15to_1.xend - ml15to_1.xbeg, dabs(r__1));
    zquit = (r__1 = ml15to_1.x - ml15to_1.xbeg, dabs(r__1));
    inc = xl * 1.5f / zquit * (ml18jr_1.mxnon + 1);
    if (ml18jr_1.ndisk != 1) {
	nsafw = inc * ml17bw_1.kkkzpw + ml17bw_1.needw;
	nsafiw = inc * ml18jr_1.nfcc + ml17bw_1.neediw;
    } else {
	nsafw = ml17bw_1.needw + inc;
	nsafiw = ml17bw_1.neediw;
    }

    s_wsfi(&io___9);
    do_fio(&c__1, (char *)&nsafw, (ftnlen)sizeof(integer));
    e_wsfi();
    s_wsfi(&io___11);
    do_fio(&c__1, (char *)&nsafiw, (ftnlen)sizeof(integer));
    e_wsfi();
/* Writing concatenation */
    i__2[0] = 57, a__1[0] = "IN BVSUP, PREDICTED STORAGE ALLOCATION FOR WORK"
	    " ARRAY IS ";
    i__2[1] = 8, a__1[1] = xern1;
    i__2[2] = 50, a__1[2] = ", PREDICTED STORAGE ALLOCATION FOR IWORK ARRAY "
	    "IS ";
    i__2[3] = 8, a__1[3] = xern2;
    s_cat(ch__1, a__1, i__2, &c__4, (ftnlen)123);
    xermsg_("SLATEC", "EXBVP", ch__1, &c__1, &c__0, (ftnlen)6, (ftnlen)5, (
	    ftnlen)123);

L30:
    iwork[1] = ml18jr_1.mxnon;
    return 0;
} /* exbvp_ */

