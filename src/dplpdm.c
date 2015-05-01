/* dplpdm.f -- translated by f2c (version 12.02.01).
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
    doublereal small;
    integer lp, lenl, lenu, ncp, lrow, lcol;
} la05dd_;

#define la05dd_1 la05dd_

/* Table of constant values */

static integer c__1 = 1;
static integer c__28 = 28;
static integer c__27 = 27;
static integer c__2 = 2;

/* DECK DPLPDM */
/* Subroutine */ int dplpdm_(integer *mrelas, integer *nvars__, integer *lmx, 
	integer *lbm, integer *nredc, integer *info, integer *iopt, integer *
	ibasis, integer *imat, integer *ibrc, integer *ipr, integer *iwr, 
	integer *ind, integer *ibb, doublereal *anorm, doublereal *eps, 
	doublereal *uu, doublereal *gg, doublereal *amat, doublereal *basmat, 
	doublereal *csc, doublereal *wr, logical *singlr, logical *redbas)
{
    /* System generated locals */
    address a__1[2];
    integer ibrc_dim1, ibrc_offset, i__1, i__2[2];
    char ch__1[55];

    /* Local variables */
    static integer i__, j, k;
    static doublereal aij, one;
    static integer nzbm;
    static doublereal zero;
    static char xern3[16];
    extern /* Subroutine */ int la05ad_(doublereal *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, doublereal *, 
	    doublereal *, doublereal *);
    extern doublereal dasum_(integer *, doublereal *, integer *);
    static integer iplace;
    extern /* Subroutine */ int xermsg_(char *, char *, char *, integer *, 
	    integer *, ftnlen, ftnlen, ftnlen), dpnnzr_(integer *, doublereal 
	    *, integer *, doublereal *, integer *, integer *);

    /* Fortran I/O blocks */
    static icilist io___10 = { 0, xern3, 0, "(1PE15.6)", 16, 1 };


/* ***BEGIN PROLOGUE  DPLPDM */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to DSPLP */
/* ***LIBRARY   SLATEC */
/* ***TYPE      DOUBLE PRECISION (SPLPDM-S, DPLPDM-D) */
/* ***AUTHOR  (UNKNOWN) */
/* ***DESCRIPTION */

/*     THIS SUBPROGRAM IS FROM THE DSPLP( ) PACKAGE.  IT PERFORMS THE */
/*     TASK OF DEFINING THE ENTRIES OF THE BASIS MATRIX AND */
/*     DECOMPOSING IT USING THE LA05 PACKAGE. */
/*     IT IS THE MAIN PART OF THE PROCEDURE (DECOMPOSE BASIS MATRIX). */

/* ***SEE ALSO  DSPLP */
/* ***ROUTINES CALLED  DASUM, DPNNZR, LA05AD, XERMSG */
/* ***COMMON BLOCKS    LA05DD */
/* ***REVISION HISTORY  (YYMMDD) */
/*   811215  DATE WRITTEN */
/*   890605  Added DASUM to list of DOUBLE PRECISION variables. */
/*   890605  Removed unreferenced labels.  (WRB) */
/*   891009  Removed unreferenced variable.  (WRB) */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ) */
/*   900328  Added TYPE section.  (WRB) */
/*   900510  Convert XERRWV calls to XERMSG calls, convert do-it-yourself */
/*           DO loops to DO loops.  (RWC) */
/* ***END PROLOGUE  DPLPDM */

/*     COMMON BLOCK USED BY LA05 () PACKAGE.. */

/* ***FIRST EXECUTABLE STATEMENT  DPLPDM */
    /* Parameter adjustments */
    ibrc_dim1 = *lbm;
    ibrc_offset = 1 + ibrc_dim1;
    ibrc -= ibrc_offset;
    --ibasis;
    --imat;
    --ipr;
    --iwr;
    --ind;
    --ibb;
    --amat;
    --basmat;
    --csc;
    --wr;

    /* Function Body */
    zero = 0.;
    one = 1.;

/*     DEFINE BASIS MATRIX BY COLUMNS FOR SPARSE MATRIX EQUATION SOLVER. */
/*     THE LA05AD() SUBPROGRAM REQUIRES THE NONZERO ENTRIES OF THE MATRIX */
/*     TOGETHER WITH THE ROW AND COLUMN INDICES. */

    nzbm = 0;

/*     DEFINE DEPENDENT VARIABLE COLUMNS. THESE ARE */
/*     COLS. OF THE IDENTITY MATRIX AND IMPLICITLY GENERATED. */

    i__1 = *mrelas;
    for (k = 1; k <= i__1; ++k) {
	j = ibasis[k];
	if (j > *nvars__) {
	    ++nzbm;
	    if (ind[j] == 2) {
		basmat[nzbm] = one;
	    } else {
		basmat[nzbm] = -one;
	    }
	    ibrc[nzbm + ibrc_dim1] = j - *nvars__;
	    ibrc[nzbm + (ibrc_dim1 << 1)] = k;
	} else {

/*           DEFINE THE INDEP. VARIABLE COLS.  THIS REQUIRES RETRIEVING */
/*           THE COLS. FROM THE SPARSE MATRIX DATA STRUCTURE. */

	    i__ = 0;
L10:
	    dpnnzr_(&i__, &aij, &iplace, &amat[1], &imat[1], &j);
	    if (i__ > 0) {
		++nzbm;
		basmat[nzbm] = aij * csc[j];
		ibrc[nzbm + ibrc_dim1] = i__;
		ibrc[nzbm + (ibrc_dim1 << 1)] = k;
		goto L10;
	    }
	}
/* L20: */
    }

    *singlr = FALSE_;

/*     RECOMPUTE MATRIX NORM USING CRUDE NORM  =  SUM OF MAGNITUDES. */

    *anorm = dasum_(&nzbm, &basmat[1], &c__1);
    la05dd_1.small = *eps * *anorm;

/*     GET AN L-U FACTORIZATION OF THE BASIS MATRIX. */

    ++(*nredc);
    *redbas = TRUE_;
    la05ad_(&basmat[1], &ibrc[ibrc_offset], &nzbm, lbm, mrelas, &ipr[1], &iwr[
	    1], &wr[1], gg, uu);

/*     CHECK RETURN VALUE OF ERROR FLAG, GG. */

    if (*gg >= zero) {
	return 0;
    }
    if (*gg == -7.f) {
	xermsg_("SLATEC", "DPLPDM", "IN DSPLP, SHORT ON STORAGE FOR LA05AD. "
		" USE PRGOPT(*) TO GIVE MORE.", &c__28, iopt, (ftnlen)6, (
		ftnlen)6, (ftnlen)67);
	*info = -28;
    } else if (*gg == -5.f) {
	*singlr = TRUE_;
    } else {
	s_wsfi(&io___10);
	do_fio(&c__1, (char *)&(*gg), (ftnlen)sizeof(doublereal));
	e_wsfi();
/* Writing concatenation */
	i__2[0] = 39, a__1[0] = "IN DSPLP, LA05AD RETURNED ERROR FLAG = ";
	i__2[1] = 16, a__1[1] = xern3;
	s_cat(ch__1, a__1, i__2, &c__2, (ftnlen)55);
	xermsg_("SLATEC", "DPLPDM", ch__1, &c__27, iopt, (ftnlen)6, (ftnlen)6,
		 (ftnlen)55);
	*info = -27;
    }
    return 0;
} /* dplpdm_ */

