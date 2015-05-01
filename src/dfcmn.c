/* dfcmn.f -- translated by f2c (version 12.02.01).
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

static integer c__2 = 2;
static integer c__1 = 1;
static doublereal c_b45 = 0.;
static integer c__0 = 0;
static doublereal c_b74 = -1.;

/* DECK DFCMN */
/* Subroutine */ int dfcmn_(integer *ndata, doublereal *xdata, doublereal *
	ydata, doublereal *sddata, integer *nord, integer *nbkpt, doublereal *
	bkptin, integer *nconst, doublereal *xconst, doublereal *yconst, 
	integer *nderiv, integer *mode, doublereal *coeff, doublereal *bf, 
	doublereal *xtemp, doublereal *ptemp, doublereal *bkpt, doublereal *g,
	 integer *mdg, doublereal *w, integer *mdw, doublereal *work, integer 
	*iwork)
{
    /* System generated locals */
    address a__1[2];
    integer bf_dim1, bf_offset, g_dim1, g_offset, w_dim1, w_offset, i__1, 
	    i__2[2], i__3, i__4;
    doublereal d__1, d__2;
    char ch__1[59], ch__2[61];

    /* Local variables */
    static integer i__, l, n, nb, ip, ir, mt, lw, np1, iw1, iw2;
    static logical var, new__, band;
    static doublereal xval, xmin, yval, xmax;
    static integer irow;
    static char xern1[8];
    static integer intw1, idata;
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *), dlsei_(doublereal *, integer *, integer *, integer *, 
	    integer *, integer *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *);
    static integer ileft;
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), daxpy_(integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *);
    static integer itype;
    extern /* Subroutine */ int dsort_(doublereal *, doublereal *, integer *, 
	    integer *);
    static doublereal dummy, rnorm;
    static integer nordm1, nordp1;
    extern /* Subroutine */ int dbndac_(doublereal *, integer *, integer *, 
	    integer *, integer *, integer *, integer *), dbndsl_(integer *, 
	    doublereal *, integer *, integer *, integer *, integer *, 
	    doublereal *, integer *, doublereal *);
    static integer ideriv, neqcon, nincon;
    extern /* Subroutine */ int dfspvd_(doublereal *, integer *, doublereal *,
	     integer *, doublereal *, integer *), dfspvn_(doublereal *, 
	    integer *, integer *, doublereal *, integer *, doublereal *);
    static doublereal rnorme;
    extern /* Subroutine */ int xermsg_(char *, char *, char *, integer *, 
	    integer *, ftnlen, ftnlen, ftnlen);
    static doublereal rnorml, prgopt[10];
    static integer intrvl;

    /* Fortran I/O blocks */
    static icilist io___5 = { 0, xern1, 0, "(I8)", 8, 1 };
    static icilist io___32 = { 0, xern1, 0, "(I8)", 8, 1 };
    static icilist io___33 = { 0, xern1, 0, "(I8)", 8, 1 };


/* ***BEGIN PROLOGUE  DFCMN */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to FC */
/* ***LIBRARY   SLATEC */
/* ***TYPE      DOUBLE PRECISION (FCMN-S, DFCMN-D) */
/* ***AUTHOR  (UNKNOWN) */
/* ***DESCRIPTION */

/*     This is a companion subprogram to DFC( ). */
/*     The documentation for DFC( ) has complete usage instructions. */

/* ***SEE ALSO  DFC */
/* ***ROUTINES CALLED  DAXPY, DBNDAC, DBNDSL, DCOPY, DFSPVD, DFSPVN, */
/*                    DLSEI, DSCAL, DSORT, XERMSG */
/* ***REVISION HISTORY  (YYMMDD) */
/*   780801  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890618  Completely restructured and extensively revised (WRB & RWC) */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ) */
/*   900328  Added TYPE section.  (WRB) */
/*   900510  Convert XERRWV calls to XERMSG calls.  (RWC) */
/*   900604  DP version created from SP version.  (RWC) */
/* ***END PROLOGUE  DFCMN */



/* ***FIRST EXECUTABLE STATEMENT  DFCMN */

/*     Analyze input. */

    /* Parameter adjustments */
    --xdata;
    --ydata;
    --sddata;
    bf_dim1 = *nord;
    bf_offset = 1 + bf_dim1;
    bf -= bf_offset;
    --bkptin;
    --xconst;
    --yconst;
    --nderiv;
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
    --work;
    --iwork;

    /* Function Body */
    if (*nord < 1 || *nord > 20) {
	xermsg_("SLATEC", "DFCMN", "IN DFC, THE ORDER OF THE B-SPLINE MUST B"
		"E 1 THRU 20.", &c__2, &c__1, (ftnlen)6, (ftnlen)5, (ftnlen)52)
		;
	*mode = -1;
	return 0;

    } else if (*nbkpt < *nord << 1) {
	xermsg_("SLATEC", "DFCMN", "IN DFC, THE NUMBER OF KNOTS MUST BE AT L"
		"EAST TWICE THE B-SPLINE ORDER.", &c__2, &c__1, (ftnlen)6, (
		ftnlen)5, (ftnlen)70);
	*mode = -1;
	return 0;
    }

    if (*ndata < 0) {
	xermsg_("SLATEC", "DFCMN", "IN DFC, THE NUMBER OF DATA POINTS MUST B"
		"E NONNEGATIVE.", &c__2, &c__1, (ftnlen)6, (ftnlen)5, (ftnlen)
		54);
	*mode = -1;
	return 0;
    }

/*     Amount of storage allocated for W(*), IW(*). */

    iw1 = iwork[1];
    iw2 = iwork[2];
/* Computing 2nd power */
    i__1 = *nord;
    nb = (*nbkpt - *nord + 3) * (*nord + 1) + (max(*ndata,*nbkpt) << 1) + *
	    nbkpt + i__1 * i__1;

/*     See if sufficient storage has been allocated. */

    if (iw1 < nb) {
	s_wsfi(&io___5);
	do_fio(&c__1, (char *)&nb, (ftnlen)sizeof(integer));
	e_wsfi();
/* Writing concatenation */
	i__2[0] = 51, a__1[0] = "IN DFC, INSUFFICIENT STORAGE FOR W(*).  CHE"
		"CK NB = ";
	i__2[1] = 8, a__1[1] = xern1;
	s_cat(ch__1, a__1, i__2, &c__2, (ftnlen)59);
	xermsg_("SLATEC", "DFCMN", ch__1, &c__2, &c__1, (ftnlen)6, (ftnlen)5, 
		(ftnlen)59);
	*mode = -1;
	return 0;
    }

    if (*mode == 1) {
	band = TRUE_;
	var = FALSE_;
	new__ = TRUE_;
    } else if (*mode == 2) {
	band = FALSE_;
	var = TRUE_;
	new__ = TRUE_;
    } else if (*mode == 3) {
	band = TRUE_;
	var = FALSE_;
	new__ = FALSE_;
    } else if (*mode == 4) {
	band = FALSE_;
	var = TRUE_;
	new__ = FALSE_;
    } else {
	xermsg_("SLATEC", "DFCMN", "IN DFC, INPUT VALUE OF MODE MUST BE 1-4.",
		 &c__2, &c__1, (ftnlen)6, (ftnlen)5, (ftnlen)40);
	*mode = -1;
	return 0;
    }
    *mode = 0;

/*     Sort the breakpoints. */

    dcopy_(nbkpt, &bkptin[1], &c__1, &bkpt[1], &c__1);
    dsort_(&bkpt[1], &dummy, nbkpt, &c__1);

/*     Initialize variables. */

    neqcon = 0;
    nincon = 0;
    i__1 = *nconst;
    for (i__ = 1; i__ <= i__1; ++i__) {
	l = nderiv[i__];
	itype = l % 4;
	if (itype < 2) {
	    ++nincon;
	} else {
	    ++neqcon;
	}
/* L100: */
    }

/*     Compute the number of variables. */

    n = *nbkpt - *nord;
    np1 = n + 1;
    lw = nb + (np1 + *nconst) * np1 + (neqcon + np1 << 1) + (nincon + np1) + (
	    nincon + 2) * (np1 + 6);
    intw1 = nincon + (np1 << 1);

/*     Save interval containing knots. */

    xmin = bkpt[*nord];
    xmax = bkpt[np1];

/*     Find the smallest referenced independent variable value in any */
/*     constraint. */

    i__1 = *nconst;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing MIN */
	d__1 = xmin, d__2 = xconst[i__];
	xmin = min(d__1,d__2);
/* Computing MAX */
	d__1 = xmax, d__2 = xconst[i__];
	xmax = max(d__1,d__2);
/* L110: */
    }
    nordm1 = *nord - 1;
    nordp1 = *nord + 1;

/*     Define the option vector PRGOPT(1-10) for use in DLSEI( ). */

    prgopt[0] = 4.;

/*     Set the covariance matrix computation flag. */

    prgopt[1] = 1.;
    if (var) {
	prgopt[2] = 1.;
    } else {
	prgopt[2] = 0.;
    }

/*     Increase the rank determination tolerances for both equality */
/*     constraint equations and least squares equations. */

    prgopt[3] = 7.;
    prgopt[4] = 4.;
    prgopt[5] = 1e-4;

    prgopt[6] = 10.;
    prgopt[7] = 5.;
    prgopt[8] = 1e-4;

    prgopt[9] = 1.;

/*     Turn off work array length checking in DLSEI( ). */

    iwork[1] = 0;
    iwork[2] = 0;

/*     Initialize variables and analyze input. */

    if (new__) {

/*        To process least squares equations sort data and an array of */
/*        pointers. */

	dcopy_(ndata, &xdata[1], &c__1, &xtemp[1], &c__1);
	i__1 = *ndata;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    ptemp[i__] = (doublereal) i__;
/* L120: */
	}

	if (*ndata > 0) {
	    dsort_(&xtemp[1], &ptemp[1], ndata, &c__2);
	    xmin = min(xmin,xtemp[1]);
/* Computing MAX */
	    d__1 = xmax, d__2 = xtemp[*ndata];
	    xmax = max(d__1,d__2);
	}

/*        Fix breakpoint array if needed. */

	i__1 = *nord;
	for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing MIN */
	    d__1 = bkpt[i__];
	    bkpt[i__] = min(d__1,xmin);
/* L130: */
	}

	i__1 = *nbkpt;
	for (i__ = np1; i__ <= i__1; ++i__) {
/* Computing MAX */
	    d__1 = bkpt[i__];
	    bkpt[i__] = max(d__1,xmax);
/* L140: */
	}

/*        Initialize parameters of banded matrix processor, DBNDAC( ). */

	mt = 0;
	ip = 1;
	ir = 1;
	ileft = *nord;
	i__1 = *ndata;
	for (idata = 1; idata <= i__1; ++idata) {

/*           Sorted indices are in PTEMP(*). */

	    l = (integer) ptemp[idata];
	    xval = xdata[l];

/*           When interval changes, process equations in the last block. */

	    if (xval >= bkpt[ileft + 1]) {
		i__3 = ileft - nordm1;
		dbndac_(&g[g_offset], mdg, nord, &ip, &ir, &mt, &i__3);
		mt = 0;

/*              Move pointer up to have BKPT(ILEFT).LE.XVAL, */
/*                 ILEFT.LT.NP1. */

L150:
		if (xval >= bkpt[ileft + 1] && ileft < n) {
		    ++ileft;
		    goto L150;
		}
	    }

/*           Obtain B-spline function value. */

	    dfspvn_(&bkpt[1], nord, &c__1, &xval, &ileft, &bf[bf_offset]);

/*           Move row into place. */

	    irow = ir + mt;
	    ++mt;
	    dcopy_(nord, &bf[bf_offset], &c__1, &g[irow + g_dim1], mdg);
	    g[irow + nordp1 * g_dim1] = ydata[l];

/*           Scale data if uncertainty is nonzero. */

	    if (sddata[l] != 0.) {
		d__1 = 1. / sddata[l];
		dscal_(&nordp1, &d__1, &g[irow + g_dim1], mdg);
	    }

/*           When staging work area is exhausted, process rows. */

	    if (irow == *mdg - 1) {
		i__3 = ileft - nordm1;
		dbndac_(&g[g_offset], mdg, nord, &ip, &ir, &mt, &i__3);
		mt = 0;
	    }
/* L160: */
	}

/*        Process last block of equations. */

	i__1 = ileft - nordm1;
	dbndac_(&g[g_offset], mdg, nord, &ip, &ir, &mt, &i__1);

/*        Last call to adjust block positioning. */

	dcopy_(&nordp1, &c_b45, &c__0, &g[ir + g_dim1], mdg);
	dbndac_(&g[g_offset], mdg, nord, &ip, &ir, &c__1, &np1);
    }

    band = band && *nconst == 0;
    i__1 = n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	band = band && g[i__ + g_dim1] != 0.;
/* L170: */
    }

/*     Process banded least squares equations. */

    if (band) {
	dbndsl_(&c__1, &g[g_offset], mdg, nord, &ip, &ir, &coeff[1], &n, &
		rnorm);
	return 0;
    }

/*     Check further for sufficient storage in working arrays. */

    if (iw1 < lw) {
	s_wsfi(&io___32);
	do_fio(&c__1, (char *)&lw, (ftnlen)sizeof(integer));
	e_wsfi();
/* Writing concatenation */
	i__2[0] = 51, a__1[0] = "IN DFC, INSUFFICIENT STORAGE FOR W(*).  CHE"
		"CK LW = ";
	i__2[1] = 8, a__1[1] = xern1;
	s_cat(ch__1, a__1, i__2, &c__2, (ftnlen)59);
	xermsg_("SLATEC", "DFCMN", ch__1, &c__2, &c__1, (ftnlen)6, (ftnlen)5, 
		(ftnlen)59);
	*mode = -1;
	return 0;
    }

    if (iw2 < intw1) {
	s_wsfi(&io___33);
	do_fio(&c__1, (char *)&intw1, (ftnlen)sizeof(integer));
	e_wsfi();
/* Writing concatenation */
	i__2[0] = 53, a__1[0] = "IN DFC, INSUFFICIENT STORAGE FOR IW(*).  CH"
		"ECK IW1 = ";
	i__2[1] = 8, a__1[1] = xern1;
	s_cat(ch__2, a__1, i__2, &c__2, (ftnlen)61);
	xermsg_("SLATEC", "DFCMN", ch__2, &c__2, &c__1, (ftnlen)6, (ftnlen)5, 
		(ftnlen)61);
	*mode = -1;
	return 0;
    }

/*     Write equality constraints. */
/*     Analyze constraint indicators for an equality constraint. */

    neqcon = 0;
    i__1 = *nconst;
    for (idata = 1; idata <= i__1; ++idata) {
	l = nderiv[idata];
	itype = l % 4;
	if (itype > 1) {
	    ideriv = l / 4;
	    ++neqcon;
	    ileft = *nord;
	    xval = xconst[idata];

L180:
	    if (xval < bkpt[ileft + 1] || ileft >= n) {
		goto L190;
	    }
	    ++ileft;
	    goto L180;

L190:
	    i__3 = ideriv + 1;
	    dfspvd_(&bkpt[1], nord, &xval, &ileft, &bf[bf_offset], &i__3);
	    dcopy_(&np1, &c_b45, &c__0, &w[neqcon + w_dim1], mdw);
	    dcopy_(nord, &bf[(ideriv + 1) * bf_dim1 + 1], &c__1, &w[neqcon + (
		    ileft - nordm1) * w_dim1], mdw);

	    if (itype == 2) {
		w[neqcon + np1 * w_dim1] = yconst[idata];
	    } else {
		ileft = *nord;
		yval = yconst[idata];

L200:
		if (yval < bkpt[ileft + 1] || ileft >= n) {
		    goto L210;
		}
		++ileft;
		goto L200;

L210:
		i__3 = ideriv + 1;
		dfspvd_(&bkpt[1], nord, &yval, &ileft, &bf[bf_offset], &i__3);
		daxpy_(nord, &c_b74, &bf[(ideriv + 1) * bf_dim1 + 1], &c__1, &
			w[neqcon + (ileft - nordm1) * w_dim1], mdw);
	    }
	}
/* L220: */
    }

/*     Transfer least squares data. */

    i__1 = np1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	irow = i__ + neqcon;
	dcopy_(&n, &c_b45, &c__0, &w[irow + w_dim1], mdw);
/* Computing MIN */
	i__4 = np1 - i__;
	i__3 = min(i__4,*nord);
	dcopy_(&i__3, &g[i__ + g_dim1], mdg, &w[irow + i__ * w_dim1], mdw);
	w[irow + np1 * w_dim1] = g[i__ + nordp1 * g_dim1];
/* L230: */
    }

/*     Write inequality constraints. */
/*     Analyze constraint indicators for inequality constraints. */

    nincon = 0;
    i__1 = *nconst;
    for (idata = 1; idata <= i__1; ++idata) {
	l = nderiv[idata];
	itype = l % 4;
	if (itype < 2) {
	    ideriv = l / 4;
	    ++nincon;
	    ileft = *nord;
	    xval = xconst[idata];

L240:
	    if (xval < bkpt[ileft + 1] || ileft >= n) {
		goto L250;
	    }
	    ++ileft;
	    goto L240;

L250:
	    i__3 = ideriv + 1;
	    dfspvd_(&bkpt[1], nord, &xval, &ileft, &bf[bf_offset], &i__3);
	    irow = neqcon + np1 + nincon;
	    dcopy_(&n, &c_b45, &c__0, &w[irow + w_dim1], mdw);
	    intrvl = ileft - nordm1;
	    dcopy_(nord, &bf[(ideriv + 1) * bf_dim1 + 1], &c__1, &w[irow + 
		    intrvl * w_dim1], mdw);

	    if (itype == 1) {
		w[irow + np1 * w_dim1] = yconst[idata];
	    } else {
		w[irow + np1 * w_dim1] = -yconst[idata];
		dscal_(nord, &c_b74, &w[irow + intrvl * w_dim1], mdw);
	    }
	}
/* L260: */
    }

/*     Solve constrained least squares equations. */

    dlsei_(&w[w_offset], mdw, &neqcon, &np1, &nincon, &n, prgopt, &coeff[1], &
	    rnorme, &rnorml, mode, &work[1], &iwork[1]);
    return 0;
} /* dfcmn_ */

