/* svecs.f -- translated by f2c (version 12.02.01).
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
    real ae, re, tol;
    integer nxpts, nic, nopg, mxnon, ndisk, ntape, neq, indpvt, integ, nps, 
	    ntp, neqivp, numort, lnfcc, icoco;
} ml18jr_;

#define ml18jr_1 ml18jr_

/* DECK SVECS */
/* Subroutine */ int svecs_(integer *ncomp, integer *lnfc, real *yhp, real *
	work, integer *iwork, integer *inhomo, integer *iflag)
{
    /* System generated locals */
    integer yhp_dim1, yhp_offset, i__1;

    /* Local variables */
    static integer k, kp, idp;
    static real dum;
    static integer niv;
    extern /* Subroutine */ int mgsbv_(integer *, integer *, real *, integer *
	    , integer *, integer *, real *, real *, integer *, integer *, 
	    real *, real *, real *);

/* ***BEGIN PROLOGUE  SVECS */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to BVSUP */
/* ***LIBRARY   SLATEC */
/* ***TYPE      SINGLE PRECISION (SVECS-S, DVECS-D) */
/* ***AUTHOR  Watts, H. A., (SNLA) */
/* ***DESCRIPTION */

/*  This subroutine is used for the special structure of complex valued */
/*  problems. MGSBV is called upon to obtain LNFC vectors from an */
/*  original set of 2*LNFC independent vectors so that the resulting */
/*  LNFC vectors together with their imaginary product or mate vectors */
/*  form an independent set. */

/* ***SEE ALSO  BVSUP */
/* ***ROUTINES CALLED  MGSBV */
/* ***COMMON BLOCKS    ML18JR */
/* ***REVISION HISTORY  (YYMMDD) */
/*   750601  DATE WRITTEN */
/*   890921  Realigned order of variables in certain COMMON blocks. */
/*           (WRB) */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900328  Added TYPE section.  (WRB) */
/*   910722  Updated AUTHOR section.  (ALS) */
/* ***END PROLOGUE  SVECS */

/* ***FIRST EXECUTABLE STATEMENT  SVECS */
    /* Parameter adjustments */
    yhp_dim1 = *ncomp;
    yhp_offset = 1 + yhp_dim1;
    yhp -= yhp_offset;
    --work;
    --iwork;

    /* Function Body */
    if (*lnfc == 1) {
	goto L5;
    }
    niv = *lnfc;
    *lnfc <<= 1;
    ml18jr_1.lnfcc <<= 1;
    kp = *lnfc + 2 + ml18jr_1.lnfcc;
    idp = ml18jr_1.indpvt;
    ml18jr_1.indpvt = 0;
    mgsbv_(ncomp, lnfc, &yhp[yhp_offset], ncomp, &niv, iflag, &work[1], &work[
	    kp], &iwork[1], inhomo, &yhp[(*lnfc + 1) * yhp_dim1 + 1], &work[*
	    lnfc + 2], &dum);
    *lnfc /= 2;
    ml18jr_1.lnfcc /= 2;
    ml18jr_1.indpvt = idp;
    if (*iflag == 0 && niv == *lnfc) {
	goto L5;
    }
    *iflag = 99;
    return 0;
L5:
    i__1 = *ncomp;
    for (k = 1; k <= i__1; ++k) {
/* L6: */
	yhp[k + (*lnfc + 1) * yhp_dim1] = yhp[k + (ml18jr_1.lnfcc + 1) * 
		yhp_dim1];
    }
    *iflag = 1;
    return 0;
} /* svecs_ */

