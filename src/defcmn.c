/* defcmn.f -- translated by f2c (version 12.02.01).
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

/* Table of constant values */

static doublereal c_b2 = 0.;
static integer c__0 = 0;
static integer c__1 = 1;
static integer c__3 = 3;
static integer c__4 = 4;
static integer c__5 = 5;
static integer c__6 = 6;
static integer c__7 = 7;
static integer c__2 = 2;

/* DECK DEFCMN */
/* Subroutine */ int defcmn_(integer *ndata, doublereal *xdata, doublereal *
	ydata, doublereal *sddata, integer *nord, integer *nbkpt, doublereal *
	bkptin, integer *mdein, integer *mdeout, doublereal *coeff, 
	doublereal *bf, doublereal *xtemp, doublereal *ptemp, doublereal *
	bkpt, doublereal *g, integer *mdg, doublereal *w, integer *mdw, 
	integer *lw)
{
    /* System generated locals */
    address a__1[4];
    integer bf_dim1, bf_offset, g_dim1, g_offset, w_dim1, w_offset, i__1, 
	    i__2[4], i__3;
    doublereal d__1, d__2;
    char ch__1[112];

    /* Local variables */
    static integer i__, l, n, nb, ip, ir, mt, np1;
    static doublereal xval, xmin, xmax;
    static integer irow;
    static char xern1[8], xern2[8];
    static integer idata;
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *);
    static integer ileft;
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), dsort_(doublereal *, doublereal *, 
	    integer *, integer *);
    static doublereal dummy, rnorm;
    static integer nordm1, nordp1;
    extern /* Subroutine */ int dbndac_(doublereal *, integer *, integer *, 
	    integer *, integer *, integer *, integer *), dbndsl_(integer *, 
	    doublereal *, integer *, integer *, integer *, integer *, 
	    doublereal *, integer *, doublereal *), dfspvn_(doublereal *, 
	    integer *, integer *, doublereal *, integer *, doublereal *);
    static integer intseq;
    extern /* Subroutine */ int xermsg_(char *, char *, char *, integer *, 
	    integer *, ftnlen, ftnlen, ftnlen);

    /* Fortran I/O blocks */
    static icilist io___5 = { 0, xern1, 0, "(I8)", 8, 1 };
    static icilist io___7 = { 0, xern2, 0, "(I8)", 8, 1 };


/* ***BEGIN PROLOGUE  DEFCMN */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to DEFC */
/* ***LIBRARY   SLATEC */
/* ***TYPE      DOUBLE PRECISION (EFCMN-S, DEFCMN-D) */
/* ***AUTHOR  Hanson, R. J., (SNLA) */
/* ***DESCRIPTION */

/*     This is a companion subprogram to DEFC( ). */
/*     This subprogram does weighted least squares fitting of data by */
/*     B-spline curves. */
/*     The documentation for DEFC( ) has complete usage instructions. */

/* ***SEE ALSO  DEFC */
/* ***ROUTINES CALLED  DBNDAC, DBNDSL, DCOPY, DFSPVN, DSCAL, DSORT, XERMSG */
/* ***REVISION HISTORY  (YYMMDD) */
/*   800801  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890618  Completely restructured and extensively revised (WRB & RWC) */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ) */
/*   900328  Added TYPE section.  (WRB) */
/*   900510  Convert XERRWV calls to XERMSG calls.  (RWC) */
/*   900604  DP version created from SP version.  (RWC) */
/* ***END PROLOGUE  DEFCMN */



/* ***FIRST EXECUTABLE STATEMENT  DEFCMN */

/*     Initialize variables and analyze input. */

    /* Parameter adjustments */
    --xdata;
    --ydata;
    --sddata;
    bf_dim1 = *nord;
    bf_offset = 1 + bf_dim1;
    bf -= bf_offset;
    --bkptin;
    --coeff;
    --xtemp;
    --ptemp;
    --bkpt;
    g_dim1 = *mdg;
    g_offset = 1 + g_dim1;
    g -= g_offset;
    w_dim1 = *mdw;
    w_offset = 1 + w_dim1;
    w -= w_offset;

    /* Function Body */
    n = *nbkpt - *nord;
    np1 = n + 1;

/*     Initially set all output coefficients to zero. */

    dcopy_(&n, &c_b2, &c__0, &coeff[1], &c__1);
    *mdeout = -1;
    if (*nord < 1 || *nord > 20) {
	xermsg_("SLATEC", "DEFCMN", "IN DEFC, THE ORDER OF THE B-SPLINE MUST"
		" BE 1 THRU 20.", &c__3, &c__1, (ftnlen)6, (ftnlen)6, (ftnlen)
		53);
	return 0;
    }

    if (*nbkpt < *nord << 1) {
	xermsg_("SLATEC", "DEFCMN", "IN DEFC, THE NUMBER OF KNOTS MUST BE AT"
		" LEAST TWICE THE B-SPLINE ORDER.", &c__4, &c__1, (ftnlen)6, (
		ftnlen)6, (ftnlen)71);
	return 0;
    }

    if (*ndata < 0) {
	xermsg_("SLATEC", "DEFCMN", "IN DEFC, THE NUMBER OF DATA POINTS MUST"
		" BE NONNEGATIVE.", &c__5, &c__1, (ftnlen)6, (ftnlen)6, (
		ftnlen)55);
	return 0;
    }

/* Computing 2nd power */
    i__1 = *nord;
    nb = (*nbkpt - *nord + 3) * (*nord + 1) + (*nbkpt + 1) * (*nord + 1) + (
	    max(*nbkpt,*ndata) << 1) + *nbkpt + i__1 * i__1;
    if (*lw < nb) {
	s_wsfi(&io___5);
	do_fio(&c__1, (char *)&nb, (ftnlen)sizeof(integer));
	e_wsfi();
	s_wsfi(&io___7);
	do_fio(&c__1, (char *)&(*lw), (ftnlen)sizeof(integer));
	e_wsfi();
/* Writing concatenation */
	i__2[0] = 87, a__1[0] = "IN DEFC, INSUFFICIENT STORAGE FOR W(*).  CH"
		"ECK FORMULA THAT READS LW.GE. ... .  NEED = ";
	i__2[1] = 8, a__1[1] = xern1;
	i__2[2] = 9, a__1[2] = " GIVEN = ";
	i__2[3] = 8, a__1[3] = xern2;
	s_cat(ch__1, a__1, i__2, &c__4, (ftnlen)112);
	xermsg_("SLATEC", "DEFCMN", ch__1, &c__6, &c__1, (ftnlen)6, (ftnlen)6,
		 (ftnlen)112);
	*mdeout = -1;
	return 0;
    }

    if (*mdein != 1 && *mdein != 2) {
	xermsg_("SLATEC", "DEFCMN", "IN DEFC, INPUT VALUE OF MDEIN MUST BE 1"
		"-2.", &c__7, &c__1, (ftnlen)6, (ftnlen)6, (ftnlen)42);
	return 0;
    }

/*     Sort the breakpoints. */

    dcopy_(nbkpt, &bkptin[1], &c__1, &bkpt[1], &c__1);
    dsort_(&bkpt[1], &dummy, nbkpt, &c__1);

/*     Save interval containing knots. */

    xmin = bkpt[*nord];
    xmax = bkpt[np1];
    nordm1 = *nord - 1;
    nordp1 = *nord + 1;

/*     Process least squares equations. */

/*     Sort data and an array of pointers. */

    dcopy_(ndata, &xdata[1], &c__1, &xtemp[1], &c__1);
    i__1 = *ndata;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ptemp[i__] = (doublereal) i__;
/* L100: */
    }

    if (*ndata > 0) {
	dsort_(&xtemp[1], &ptemp[1], ndata, &c__2);
	xmin = min(xmin,xtemp[1]);
/* Computing MAX */
	d__1 = xmax, d__2 = xtemp[*ndata];
	xmax = max(d__1,d__2);
    }

/*     Fix breakpoint array if needed. This should only involve very */
/*     minor differences with the input array of breakpoints. */

    i__1 = *nord;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing MIN */
	d__1 = bkpt[i__];
	bkpt[i__] = min(d__1,xmin);
/* L110: */
    }

    i__1 = *nbkpt;
    for (i__ = np1; i__ <= i__1; ++i__) {
/* Computing MAX */
	d__1 = bkpt[i__];
	bkpt[i__] = max(d__1,xmax);
/* L120: */
    }

/*     Initialize parameters of banded matrix processor, DBNDAC( ). */

    mt = 0;
    ip = 1;
    ir = 1;
    ileft = *nord;
    intseq = 1;
    i__1 = *ndata;
    for (idata = 1; idata <= i__1; ++idata) {

/*        Sorted indices are in PTEMP(*). */

	l = (integer) ptemp[idata];
	xval = xdata[l];

/*        When interval changes, process equations in the last block. */

	if (xval >= bkpt[ileft + 1]) {
	    i__3 = ileft - nordm1;
	    dbndac_(&g[g_offset], mdg, nord, &ip, &ir, &mt, &i__3);
	    mt = 0;

/*           Move pointer up to have BKPT(ILEFT).LE.XVAL, ILEFT.LE.N. */

	    i__3 = n;
	    for (ileft = ileft; ileft <= i__3; ++ileft) {
		if (xval < bkpt[ileft + 1]) {
		    goto L140;
		}
		if (*mdein == 2) {

/*                 Data is being sequentially accumulated. */
/*                 Transfer previously accumulated rows from W(*,*) to */
/*                 G(*,*) and process them. */

		    dcopy_(&nordp1, &w[intseq + w_dim1], mdw, &g[ir + g_dim1],
			     mdg);
		    dbndac_(&g[g_offset], mdg, nord, &ip, &ir, &c__1, &intseq)
			    ;
		    ++intseq;
		}
/* L130: */
	    }
	}

/*        Obtain B-spline function value. */

L140:
	dfspvn_(&bkpt[1], nord, &c__1, &xval, &ileft, &bf[bf_offset]);

/*        Move row into place. */

	irow = ir + mt;
	++mt;
	dcopy_(nord, &bf[bf_offset], &c__1, &g[irow + g_dim1], mdg);
	g[irow + nordp1 * g_dim1] = ydata[l];

/*        Scale data if uncertainty is nonzero. */

	if (sddata[l] != 0.) {
	    d__1 = 1. / sddata[l];
	    dscal_(&nordp1, &d__1, &g[irow + g_dim1], mdg);
	}

/*        When staging work area is exhausted, process rows. */

	if (irow == *mdg - 1) {
	    i__3 = ileft - nordm1;
	    dbndac_(&g[g_offset], mdg, nord, &ip, &ir, &mt, &i__3);
	    mt = 0;
	}
/* L150: */
    }

/*     Process last block of equations. */

    i__1 = ileft - nordm1;
    dbndac_(&g[g_offset], mdg, nord, &ip, &ir, &mt, &i__1);

/*     Finish processing any previously accumulated rows from W(*,*) */
/*     to G(*,*). */

    if (*mdein == 2) {
	i__1 = np1;
	for (i__ = intseq; i__ <= i__1; ++i__) {
	    dcopy_(&nordp1, &w[i__ + w_dim1], mdw, &g[ir + g_dim1], mdg);
	    i__3 = min(n,i__);
	    dbndac_(&g[g_offset], mdg, nord, &ip, &ir, &c__1, &i__3);
/* L160: */
	}
    }

/*     Last call to adjust block positioning. */

    dcopy_(&nordp1, &c_b2, &c__0, &g[ir + g_dim1], mdg);
    dbndac_(&g[g_offset], mdg, nord, &ip, &ir, &c__1, &np1);

/*     Transfer accumulated rows from G(*,*) to W(*,*) for */
/*     possible later sequential accumulation. */

    i__1 = np1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	dcopy_(&nordp1, &g[i__ + g_dim1], mdg, &w[i__ + w_dim1], mdw);
/* L170: */
    }

/*     Solve for coefficients when possible. */

    i__1 = n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (g[i__ + g_dim1] == 0.) {
	    *mdeout = 2;
	    return 0;
	}
/* L180: */
    }

/*     All the diagonal terms in the accumulated triangular */
/*     matrix are nonzero.  The solution can be computed but */
/*     it may be unsuitable for further use due to poor */
/*     conditioning or the lack of constraints.  No checking */
/*     for either of these is done here. */

    dbndsl_(&c__1, &g[g_offset], mdg, nord, &ip, &ir, &coeff[1], &n, &rnorm);
    *mdeout = 1;
    return 0;
} /* defcmn_ */

