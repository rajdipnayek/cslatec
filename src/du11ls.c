/* du11ls.f -- translated by f2c (version 12.02.01).
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

static integer c__1 = 1;
static integer c__8 = 8;
static integer c__0 = 0;

/* DECK DU11LS */
/* Subroutine */ int du11ls_(doublereal *a, integer *mda, integer *m, integer 
	*n, doublereal *ub, doublereal *db, integer *mode, integer *np, 
	integer *krank, integer *ksure, doublereal *h__, doublereal *w, 
	doublereal *eb, integer *ic, integer *ir)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;
    doublereal d__1, d__2;

    /* Local variables */
    static integer i__, j, l;
    static doublereal t, r2, bb;
    static integer ii, kk, mm, is;
    static doublereal tn;
    static integer kz;
    static doublereal tt;
    static integer im1, jm1, km1, lm1, jp1, kp1, kmi, nmk;
    static doublereal sum;
    extern doublereal ddot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    static integer imin, jmax;
    static doublereal rmin, temp;
    extern doublereal dnrm2_(integer *, doublereal *, integer *);
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *), dswap_(integer *, doublereal *, integer *, doublereal 
	    *, integer *), iswap_(integer *, integer *, integer *, integer *, 
	    integer *), daxpy_(integer *, doublereal *, doublereal *, integer 
	    *, doublereal *, integer *);
    extern integer idamax_(integer *, doublereal *, integer *);
    extern /* Subroutine */ int xermsg_(char *, char *, char *, integer *, 
	    integer *, ftnlen, ftnlen, ftnlen);

/* ***BEGIN PROLOGUE  DU11LS */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to DLLSIA */
/* ***LIBRARY   SLATEC */
/* ***TYPE      DOUBLE PRECISION (U11LS-S, DU11LS-D) */
/* ***AUTHOR  (UNKNOWN) */
/* ***DESCRIPTION */

/*   **** Double Precision version of U11LS **** */

/*       This routine performs a QR factorization of A */
/*       using Householder transformations. Row and */
/*       column pivots are chosen to reduce the growth */
/*       of round-off and to help detect possible rank */
/*       deficiency. */

/* ***SEE ALSO  DLLSIA */
/* ***ROUTINES CALLED  DAXPY, DDOT, DNRM2, DSCAL, DSWAP, IDAMAX, ISWAP, */
/*                    XERMSG */
/* ***REVISION HISTORY  (YYMMDD) */
/*   810801  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890831  Modified array declarations.  (WRB) */
/*   891009  Removed unreferenced variable.  (WRB) */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ) */
/*   900328  Added TYPE section.  (WRB) */
/* ***END PROLOGUE  DU11LS */

/*        INITIALIZATION */

/* ***FIRST EXECUTABLE STATEMENT  DU11LS */
    /* Parameter adjustments */
    a_dim1 = *mda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --ub;
    --db;
    --h__;
    --w;
    --eb;
    --ic;
    --ir;

    /* Function Body */
    j = 0;
    *krank = *n;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ic[i__] = i__;
/* L10: */
    }
    i__1 = *m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ir[i__] = i__;
/* L12: */
    }

/*        DETERMINE REL AND ABS ERROR VECTORS */



/*        CALCULATE COL LENGTH */

    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	h__[i__] = dnrm2_(m, &a[i__ * a_dim1 + 1], &c__1);
	w[i__] = h__[i__];
/* L30: */
    }

/*         INITIALIZE ERROR BOUNDS */

    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing MAX */
	d__1 = db[i__], d__2 = ub[i__] * h__[i__];
	eb[i__] = max(d__1,d__2);
	ub[i__] = eb[i__];
	db[i__] = 0.;
/* L40: */
    }

/*          DISCARD SELF DEPENDENT COLUMNS */

    i__ = 1;
L50:
    if (eb[i__] >= h__[i__]) {
	goto L60;
    }
    if (i__ == *krank) {
	goto L70;
    }
    ++i__;
    goto L50;

/*          MATRIX REDUCTION */

L60:
    kk = *krank;
    --(*krank);
    if (*mode == 0) {
	return 0;
    }
    if (i__ > *np) {
	goto L64;
    }
    xermsg_("SLATEC", "DU11LS", "FIRST NP COLUMNS ARE LINEARLY DEPENDENT", &
	    c__8, &c__0, (ftnlen)6, (ftnlen)6, (ftnlen)39);
    *krank = i__ - 1;
    return 0;
L64:
    if (i__ > *krank) {
	goto L70;
    }
    dswap_(&c__1, &eb[i__], &c__1, &eb[kk], &c__1);
    dswap_(&c__1, &ub[i__], &c__1, &ub[kk], &c__1);
    dswap_(&c__1, &w[i__], &c__1, &w[kk], &c__1);
    dswap_(&c__1, &h__[i__], &c__1, &h__[kk], &c__1);
    iswap_(&c__1, &ic[i__], &c__1, &ic[kk], &c__1);
    dswap_(m, &a[i__ * a_dim1 + 1], &c__1, &a[kk * a_dim1 + 1], &c__1);
    goto L50;

/*           TEST FOR ZERO RANK */

L70:
    if (*krank > 0) {
	goto L80;
    }
    *krank = 0;
    *ksure = 0;
    return 0;
L80:

/*        M A I N    L O O P */

L110:
    ++j;
    jp1 = j + 1;
    jm1 = j - 1;
    kz = *krank;
    if (j <= *np) {
	kz = j;
    }

/*        EACH COL HAS MM=M-J+1 COMPONENTS */

    mm = *m - j + 1;

/*         UB DETERMINES COLUMN PIVOT */

L115:
    imin = j;
    if (h__[j] == 0.) {
	goto L170;
    }
    rmin = ub[j] / h__[j];
    i__1 = kz;
    for (i__ = j; i__ <= i__1; ++i__) {
	if (ub[i__] >= h__[i__] * rmin) {
	    goto L120;
	}
	rmin = ub[i__] / h__[i__];
	imin = i__;
L120:
	;
    }

/*     TEST FOR RANK DEFICIENCY */

    if (rmin < 1.) {
	goto L200;
    }
    tt = (eb[imin] + (d__1 = db[imin], abs(d__1))) / h__[imin];
    if (tt >= 1.) {
	goto L170;
    }
/*     COMPUTE EXACT UB */
    i__1 = jm1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	w[i__] = a[i__ + imin * a_dim1];
/* L125: */
    }
    l = jm1;
L130:
    w[l] /= a[l + l * a_dim1];
    if (l == 1) {
	goto L150;
    }
    lm1 = l - 1;
    i__1 = jm1;
    for (i__ = l; i__ <= i__1; ++i__) {
	w[lm1] -= a[lm1 + i__ * a_dim1] * w[i__];
/* L140: */
    }
    l = lm1;
    goto L130;
L150:
    tt = eb[imin];
    i__1 = jm1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	tt += (d__1 = w[i__], abs(d__1)) * eb[i__];
/* L160: */
    }
    ub[imin] = tt;
    if (ub[imin] / h__[imin] >= 1.) {
	goto L170;
    }
    goto L200;

/*        MATRIX REDUCTION */

L170:
    kk = *krank;
    --(*krank);
    kz = *krank;
    if (*mode == 0) {
	return 0;
    }
    if (j > *np) {
	goto L172;
    }
    xermsg_("SLATEC", "DU11LS", "FIRST NP COLUMNS ARE LINEARLY DEPENDENT", &
	    c__8, &c__0, (ftnlen)6, (ftnlen)6, (ftnlen)39);
    *krank = j - 1;
    return 0;
L172:
    if (imin > *krank) {
	goto L180;
    }
    iswap_(&c__1, &ic[imin], &c__1, &ic[kk], &c__1);
    dswap_(m, &a[imin * a_dim1 + 1], &c__1, &a[kk * a_dim1 + 1], &c__1);
    dswap_(&c__1, &eb[imin], &c__1, &eb[kk], &c__1);
    dswap_(&c__1, &ub[imin], &c__1, &ub[kk], &c__1);
    dswap_(&c__1, &db[imin], &c__1, &db[kk], &c__1);
    dswap_(&c__1, &w[imin], &c__1, &w[kk], &c__1);
    dswap_(&c__1, &h__[imin], &c__1, &h__[kk], &c__1);
L180:
    if (j > *krank) {
	goto L300;
    }
    goto L115;

/*        COLUMN PIVOT */

L200:
    if (imin == j) {
	goto L230;
    }
    dswap_(&c__1, &h__[j], &c__1, &h__[imin], &c__1);
    dswap_(m, &a[j * a_dim1 + 1], &c__1, &a[imin * a_dim1 + 1], &c__1);
    dswap_(&c__1, &eb[j], &c__1, &eb[imin], &c__1);
    dswap_(&c__1, &ub[j], &c__1, &ub[imin], &c__1);
    dswap_(&c__1, &db[j], &c__1, &db[imin], &c__1);
    dswap_(&c__1, &w[j], &c__1, &w[imin], &c__1);
    iswap_(&c__1, &ic[j], &c__1, &ic[imin], &c__1);

/*        ROW PIVOT */

L230:
    jmax = idamax_(&mm, &a[j + j * a_dim1], &c__1);
    jmax = jmax + j - 1;
    if (jmax == j) {
	goto L240;
    }
    dswap_(n, &a[j + a_dim1], mda, &a[jmax + a_dim1], mda);
    iswap_(&c__1, &ir[j], &c__1, &ir[jmax], &c__1);
L240:

/*     APPLY HOUSEHOLDER TRANSFORMATION */

    tn = dnrm2_(&mm, &a[j + j * a_dim1], &c__1);
    if (tn == 0.) {
	goto L170;
    }
    if (a[j + j * a_dim1] != 0.) {
	tn = d_sign(&tn, &a[j + j * a_dim1]);
    }
    d__1 = 1. / tn;
    dscal_(&mm, &d__1, &a[j + j * a_dim1], &c__1);
    a[j + j * a_dim1] += 1.;
    if (j == *n) {
	goto L250;
    }
    i__1 = *n;
    for (i__ = jp1; i__ <= i__1; ++i__) {
	bb = -ddot_(&mm, &a[j + j * a_dim1], &c__1, &a[j + i__ * a_dim1], &
		c__1) / a[j + j * a_dim1];
	daxpy_(&mm, &bb, &a[j + j * a_dim1], &c__1, &a[j + i__ * a_dim1], &
		c__1);
	if (i__ <= *np) {
	    goto L248;
	}
	if (h__[i__] == 0.) {
	    goto L248;
	}
/* Computing 2nd power */
	d__2 = (d__1 = a[j + i__ * a_dim1], abs(d__1)) / h__[i__];
	tt = 1. - d__2 * d__2;
	tt = max(tt,0.);
	t = tt;
/* Computing 2nd power */
	d__1 = h__[i__] / w[i__];
	tt = tt * .05 * (d__1 * d__1) + 1.;
	if (tt == 1.) {
	    goto L244;
	}
	h__[i__] *= sqrt(t);
	goto L246;
L244:
	i__2 = *m - j;
	h__[i__] = dnrm2_(&i__2, &a[j + 1 + i__ * a_dim1], &c__1);
	w[i__] = h__[i__];
L246:
L248:
	;
    }
L250:
    h__[j] = a[j + j * a_dim1];
    a[j + j * a_dim1] = -tn;


/*          UPDATE UB, DB */

    ub[j] /= (d__1 = a[j + j * a_dim1], abs(d__1));
    db[j] = (d_sign(&eb[j], &db[j]) + db[j]) / a[j + j * a_dim1];
    if (j == *krank) {
	goto L300;
    }
    i__1 = *krank;
    for (i__ = jp1; i__ <= i__1; ++i__) {
	ub[i__] += (d__1 = a[j + i__ * a_dim1], abs(d__1)) * ub[j];
	db[i__] -= a[j + i__ * a_dim1] * db[j];
/* L260: */
    }
    goto L110;

/*        E N D    M A I N    L O O P */

L300:

/*        COMPUTE KSURE */

    km1 = *krank - 1;
    i__1 = km1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	is = 0;
	kmi = *krank - i__;
	i__2 = kmi;
	for (ii = 1; ii <= i__2; ++ii) {
	    if (ub[ii] <= ub[ii + 1]) {
		goto L315;
	    }
	    is = 1;
	    temp = ub[ii];
	    ub[ii] = ub[ii + 1];
	    ub[ii + 1] = temp;
L315:
	    ;
	}
	if (is == 0) {
	    goto L320;
	}
/* L318: */
    }
L320:
    *ksure = 0;
    sum = 0.;
    i__1 = *krank;
    for (i__ = 1; i__ <= i__1; ++i__) {
	r2 = ub[i__] * ub[i__];
	if (r2 + sum >= 1.) {
	    goto L330;
	}
	sum += r2;
	++(*ksure);
/* L328: */
    }
L330:

/*     IF SYSTEM IS OF REDUCED RANK AND MODE = 2 */
/*     COMPLETE THE DECOMPOSITION FOR SHORTEST LEAST SQUARES SOLUTION */

    if (*krank == *n || *mode < 2) {
	goto L360;
    }
    nmk = *n - *krank;
    kp1 = *krank + 1;
    i__ = *krank;
L340:
    tn = dnrm2_(&nmk, &a[i__ + kp1 * a_dim1], mda) / a[i__ + i__ * a_dim1];
    tn = a[i__ + i__ * a_dim1] * sqrt(tn * tn + 1.);
    d__1 = 1. / tn;
    dscal_(&nmk, &d__1, &a[i__ + kp1 * a_dim1], mda);
    w[i__] = a[i__ + i__ * a_dim1] / tn + 1.;
    a[i__ + i__ * a_dim1] = -tn;
    if (i__ == 1) {
	goto L350;
    }
    im1 = i__ - 1;
    i__1 = im1;
    for (ii = 1; ii <= i__1; ++ii) {
	tt = -ddot_(&nmk, &a[ii + kp1 * a_dim1], mda, &a[i__ + kp1 * a_dim1], 
		mda) / w[i__];
	tt -= a[ii + i__ * a_dim1];
	daxpy_(&nmk, &tt, &a[i__ + kp1 * a_dim1], mda, &a[ii + kp1 * a_dim1], 
		mda);
	a[ii + i__ * a_dim1] += tt * w[i__];
/* L345: */
    }
    --i__;
    goto L340;
L350:
L360:
    return 0;
} /* du11ls_ */

