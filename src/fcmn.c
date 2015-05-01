/* fcmn.f -- translated by f2c (version 12.02.01).
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
static real c_b45 = 0.f;
static integer c__0 = 0;
static real c_b74 = -1.f;

/* DECK FCMN */
/* Subroutine */ int fcmn_(integer *ndata, real *xdata, real *ydata, real *
	sddata, integer *nord, integer *nbkpt, real *bkptin, integer *nconst, 
	real *xconst, real *yconst, integer *nderiv, integer *mode, real *
	coeff, real *bf, real *xtemp, real *ptemp, real *bkpt, real *g, 
	integer *mdg, real *w, integer *mdw, real *work, integer *iwork)
{
    /* System generated locals */
    address a__1[2];
    integer bf_dim1, bf_offset, g_dim1, g_offset, w_dim1, w_offset, i__1, 
	    i__2[2], i__3, i__4;
    real r__1, r__2;
    char ch__1[58], ch__2[60];

    /* Local variables */
    static integer i__, l, n, nb, ip, ir, mt, lw, np1, iw1, iw2;
    static logical var, new__, band;
    extern /* Subroutine */ int lsei_(real *, integer *, integer *, integer *,
	     integer *, integer *, real *, real *, real *, real *, integer *, 
	    real *, integer *);
    static real xval, xmin, yval, xmax;
    static integer irow;
    static char xern1[8];
    static integer intw1, idata, ileft;
    extern /* Subroutine */ int sscal_(integer *, real *, real *, integer *);
    static integer itype;
    static real dummy;
    extern /* Subroutine */ int scopy_(integer *, real *, integer *, real *, 
	    integer *);
    static real rnorm;
    extern /* Subroutine */ int saxpy_(integer *, real *, real *, integer *, 
	    real *, integer *), ssort_(real *, real *, integer *, integer *);
    static integer nordm1, nordp1;
    extern /* Subroutine */ int bndacc_(real *, integer *, integer *, integer 
	    *, integer *, integer *, integer *), bndsol_(integer *, real *, 
	    integer *, integer *, integer *, integer *, real *, integer *, 
	    real *);
    static integer ideriv, neqcon, nincon;
    extern /* Subroutine */ int bsplvd_(real *, integer *, real *, integer *, 
	    real *, integer *);
    static real rnorme;
    extern /* Subroutine */ int bsplvn_(real *, integer *, integer *, real *, 
	    integer *, real *), xermsg_(char *, char *, char *, integer *, 
	    integer *, ftnlen, ftnlen, ftnlen);
    static real rnorml, prgopt[10];
    static integer intrvl;

    /* Fortran I/O blocks */
    static icilist io___5 = { 0, xern1, 0, "(I8)", 8, 1 };
    static icilist io___32 = { 0, xern1, 0, "(I8)", 8, 1 };
    static icilist io___33 = { 0, xern1, 0, "(I8)", 8, 1 };


/* ***BEGIN PROLOGUE  FCMN */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to FC */
/* ***LIBRARY   SLATEC */
/* ***TYPE      SINGLE PRECISION (FCMN-S, DFCMN-D) */
/* ***AUTHOR  (UNKNOWN) */
/* ***DESCRIPTION */

/*     This is a companion subprogram to FC( ). */
/*     The documentation for FC( ) has complete usage instructions. */

/* ***SEE ALSO  FC */
/* ***ROUTINES CALLED  BNDACC, BNDSOL, BSPLVD, BSPLVN, LSEI, SAXPY, SCOPY, */
/*                    SSCAL, SSORT, XERMSG */
/* ***REVISION HISTORY  (YYMMDD) */
/*   780801  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890618  Completely restructured and extensively revised (WRB & RWC) */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ) */
/*   900328  Added TYPE section.  (WRB) */
/*   900510  Convert XERRWV calls to XERMSG calls.  (RWC) */
/* ***END PROLOGUE  FCMN */



/* ***FIRST EXECUTABLE STATEMENT  FCMN */

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
	xermsg_("SLATEC", "FCMN", "IN FC, THE ORDER OF THE B-SPLINE MUST BE "
		"1 THRU 20.", &c__2, &c__1, (ftnlen)6, (ftnlen)4, (ftnlen)51);
	*mode = -1;
	return 0;

    } else if (*nbkpt < *nord << 1) {
	xermsg_("SLATEC", "FCMN", "IN FC, THE NUMBER OF KNOTS MUST BE AT LEA"
		"ST TWICE THE B-SPLINE ORDER.", &c__2, &c__1, (ftnlen)6, (
		ftnlen)4, (ftnlen)69);
	*mode = -1;
	return 0;
    }

    if (*ndata < 0) {
	xermsg_("SLATEC", "FCMN", "IN FC, THE NUMBER OF DATA POINTS MUST BE "
		"NONNEGATIVE.", &c__2, &c__1, (ftnlen)6, (ftnlen)4, (ftnlen)53)
		;
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
	i__2[0] = 50, a__1[0] = "IN FC, INSUFFICIENT STORAGE FOR W(*).  CHEC"
		"K NB = ";
	i__2[1] = 8, a__1[1] = xern1;
	s_cat(ch__1, a__1, i__2, &c__2, (ftnlen)58);
	xermsg_("SLATEC", "FCMN", ch__1, &c__2, &c__1, (ftnlen)6, (ftnlen)4, (
		ftnlen)58);
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
	xermsg_("SLATEC", "FCMN", "IN FC, INPUT VALUE OF MODE MUST BE 1-4.", &
		c__2, &c__1, (ftnlen)6, (ftnlen)4, (ftnlen)39);
	*mode = -1;
	return 0;
    }
    *mode = 0;

/*     Sort the breakpoints. */

    scopy_(nbkpt, &bkptin[1], &c__1, &bkpt[1], &c__1);
    ssort_(&bkpt[1], &dummy, nbkpt, &c__1);

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
	r__1 = xmin, r__2 = xconst[i__];
	xmin = dmin(r__1,r__2);
/* Computing MAX */
	r__1 = xmax, r__2 = xconst[i__];
	xmax = dmax(r__1,r__2);
/* L110: */
    }
    nordm1 = *nord - 1;
    nordp1 = *nord + 1;

/*     Define the option vector PRGOPT(1-10) for use in LSEI( ). */

    prgopt[0] = 4.f;

/*     Set the covariance matrix computation flag. */

    prgopt[1] = 1.f;
    if (var) {
	prgopt[2] = 1.f;
    } else {
	prgopt[2] = 0.f;
    }

/*     Increase the rank determination tolerances for both equality */
/*     constraint equations and least squares equations. */

    prgopt[3] = 7.f;
    prgopt[4] = 4.f;
    prgopt[5] = 1e-4f;

    prgopt[6] = 10.f;
    prgopt[7] = 5.f;
    prgopt[8] = 1e-4f;

    prgopt[9] = 1.f;

/*     Turn off work array length checking in LSEI( ). */

    iwork[1] = 0;
    iwork[2] = 0;

/*     Initialize variables and analyze input. */

    if (new__) {

/*        To process least squares equations sort data and an array of */
/*        pointers. */

	scopy_(ndata, &xdata[1], &c__1, &xtemp[1], &c__1);
	i__1 = *ndata;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    ptemp[i__] = (real) i__;
/* L120: */
	}

	if (*ndata > 0) {
	    ssort_(&xtemp[1], &ptemp[1], ndata, &c__2);
	    xmin = dmin(xmin,xtemp[1]);
/* Computing MAX */
	    r__1 = xmax, r__2 = xtemp[*ndata];
	    xmax = dmax(r__1,r__2);
	}

/*        Fix breakpoint array if needed. */

	i__1 = *nord;
	for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing MIN */
	    r__1 = bkpt[i__];
	    bkpt[i__] = dmin(r__1,xmin);
/* L130: */
	}

	i__1 = *nbkpt;
	for (i__ = np1; i__ <= i__1; ++i__) {
/* Computing MAX */
	    r__1 = bkpt[i__];
	    bkpt[i__] = dmax(r__1,xmax);
/* L140: */
	}

/*        Initialize parameters of banded matrix processor, BNDACC( ). */

	mt = 0;
	ip = 1;
	ir = 1;
	ileft = *nord;
	i__1 = *ndata;
	for (idata = 1; idata <= i__1; ++idata) {

/*           Sorted indices are in PTEMP(*). */

	    l = ptemp[idata];
	    xval = xdata[l];

/*           When interval changes, process equations in the last block. */

	    if (xval >= bkpt[ileft + 1]) {
		i__3 = ileft - nordm1;
		bndacc_(&g[g_offset], mdg, nord, &ip, &ir, &mt, &i__3);
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

	    bsplvn_(&bkpt[1], nord, &c__1, &xval, &ileft, &bf[bf_offset]);

/*           Move row into place. */

	    irow = ir + mt;
	    ++mt;
	    scopy_(nord, &bf[bf_offset], &c__1, &g[irow + g_dim1], mdg);
	    g[irow + nordp1 * g_dim1] = ydata[l];

/*           Scale data if uncertainty is nonzero. */

	    if (sddata[l] != 0.f) {
		r__1 = 1.f / sddata[l];
		sscal_(&nordp1, &r__1, &g[irow + g_dim1], mdg);
	    }

/*           When staging work area is exhausted, process rows. */

	    if (irow == *mdg - 1) {
		i__3 = ileft - nordm1;
		bndacc_(&g[g_offset], mdg, nord, &ip, &ir, &mt, &i__3);
		mt = 0;
	    }
/* L160: */
	}

/*        Process last block of equations. */

	i__1 = ileft - nordm1;
	bndacc_(&g[g_offset], mdg, nord, &ip, &ir, &mt, &i__1);

/*        Last call to adjust block positioning. */

	scopy_(&nordp1, &c_b45, &c__0, &g[ir + g_dim1], mdg);
	bndacc_(&g[g_offset], mdg, nord, &ip, &ir, &c__1, &np1);
    }

    band = band && *nconst == 0;
    i__1 = n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	band = band && g[i__ + g_dim1] != 0.f;
/* L170: */
    }

/*     Process banded least squares equations. */

    if (band) {
	bndsol_(&c__1, &g[g_offset], mdg, nord, &ip, &ir, &coeff[1], &n, &
		rnorm);
	return 0;
    }

/*     Check further for sufficient storage in working arrays. */

    if (iw1 < lw) {
	s_wsfi(&io___32);
	do_fio(&c__1, (char *)&lw, (ftnlen)sizeof(integer));
	e_wsfi();
/* Writing concatenation */
	i__2[0] = 50, a__1[0] = "IN FC, INSUFFICIENT STORAGE FOR W(*).  CHEC"
		"K LW = ";
	i__2[1] = 8, a__1[1] = xern1;
	s_cat(ch__1, a__1, i__2, &c__2, (ftnlen)58);
	xermsg_("SLATEC", "FCMN", ch__1, &c__2, &c__1, (ftnlen)6, (ftnlen)4, (
		ftnlen)58);
	*mode = -1;
	return 0;
    }

    if (iw2 < intw1) {
	s_wsfi(&io___33);
	do_fio(&c__1, (char *)&intw1, (ftnlen)sizeof(integer));
	e_wsfi();
/* Writing concatenation */
	i__2[0] = 52, a__1[0] = "IN FC, INSUFFICIENT STORAGE FOR IW(*).  CHE"
		"CK IW1 = ";
	i__2[1] = 8, a__1[1] = xern1;
	s_cat(ch__2, a__1, i__2, &c__2, (ftnlen)60);
	xermsg_("SLATEC", "FCMN", ch__2, &c__2, &c__1, (ftnlen)6, (ftnlen)4, (
		ftnlen)60);
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
	    bsplvd_(&bkpt[1], nord, &xval, &ileft, &bf[bf_offset], &i__3);
	    scopy_(&np1, &c_b45, &c__0, &w[neqcon + w_dim1], mdw);
	    scopy_(nord, &bf[(ideriv + 1) * bf_dim1 + 1], &c__1, &w[neqcon + (
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
		bsplvd_(&bkpt[1], nord, &yval, &ileft, &bf[bf_offset], &i__3);
		saxpy_(nord, &c_b74, &bf[(ideriv + 1) * bf_dim1 + 1], &c__1, &
			w[neqcon + (ileft - nordm1) * w_dim1], mdw);
	    }
	}
/* L220: */
    }

/*     Transfer least squares data. */

    i__1 = np1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	irow = i__ + neqcon;
	scopy_(&n, &c_b45, &c__0, &w[irow + w_dim1], mdw);
/* Computing MIN */
	i__4 = np1 - i__;
	i__3 = min(i__4,*nord);
	scopy_(&i__3, &g[i__ + g_dim1], mdg, &w[irow + i__ * w_dim1], mdw);
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
	    bsplvd_(&bkpt[1], nord, &xval, &ileft, &bf[bf_offset], &i__3);
	    irow = neqcon + np1 + nincon;
	    scopy_(&n, &c_b45, &c__0, &w[irow + w_dim1], mdw);
	    intrvl = ileft - nordm1;
	    scopy_(nord, &bf[(ideriv + 1) * bf_dim1 + 1], &c__1, &w[irow + 
		    intrvl * w_dim1], mdw);

	    if (itype == 1) {
		w[irow + np1 * w_dim1] = yconst[idata];
	    } else {
		w[irow + np1 * w_dim1] = -yconst[idata];
		sscal_(nord, &c_b74, &w[irow + intrvl * w_dim1], mdw);
	    }
	}
/* L260: */
    }

/*     Solve constrained least squares equations. */

    lsei_(&w[w_offset], mdw, &neqcon, &np1, &nincon, &n, prgopt, &coeff[1], &
	    rnorme, &rnorml, mode, &work[1], &iwork[1]);
    return 0;
} /* fcmn_ */

