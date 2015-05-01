/* hstcs1.f -- translated by f2c (version 12.02.01).
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

static integer c__0 = 0;
static integer c__1 = 1;

/* DECK HSTCS1 */
/* Subroutine */ int hstcs1_(integer *intl, real *a, real *b, integer *m, 
	integer *mbdcnd, real *bda, real *bdb, real *c__, real *d__, integer *
	n, integer *nbdcnd, real *bdc, real *bdd, real *elmbda, real *f, 
	integer *idimf, real *pertrb, integer *ierr1, real *am, real *bm, 
	real *cm, real *an, real *bn, real *cn, real *snth, real *rsq, real *
	wrk)
{
    /* System generated locals */
    integer f_dim1, f_offset, i__1, i__2;
    real r__1;

    /* Local variables */
    static integer i__, j;
    static real x, y, a1, a2, a3;
    static integer nb;
    static real dr, dth;
    static integer isw;
    static real dthsq;
    extern /* Subroutine */ int blktri_(integer *, integer *, integer *, real 
	    *, real *, real *, integer *, integer *, real *, real *, real *, 
	    integer *, real *, integer *, real *);

/* ***BEGIN PROLOGUE  HSTCS1 */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to HSTCSP */
/* ***LIBRARY   SLATEC */
/* ***TYPE      SINGLE PRECISION (HSTCS1-S) */
/* ***AUTHOR  (UNKNOWN) */
/* ***SEE ALSO  HSTCSP */
/* ***ROUTINES CALLED  BLKTRI */
/* ***REVISION HISTORY  (YYMMDD) */
/*   801001  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900402  Added TYPE section.  (WRB) */
/* ***END PROLOGUE  HSTCS1 */
/* ***FIRST EXECUTABLE STATEMENT  HSTCS1 */
    /* Parameter adjustments */
    --bda;
    --bdb;
    --bdc;
    --bdd;
    f_dim1 = *idimf;
    f_offset = 1 + f_dim1;
    f -= f_offset;
    --am;
    --bm;
    --cm;
    --an;
    --bn;
    --cn;
    --snth;
    --rsq;
    --wrk;

    /* Function Body */
    dth = (*b - *a) / *m;
    dthsq = dth * dth;
    i__1 = *m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	snth[i__] = sin(*a + (i__ - .5f) * dth);
/* L101: */
    }
    dr = (*d__ - *c__) / *n;
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
/* Computing 2nd power */
	r__1 = *c__ + (j - .5f) * dr;
	rsq[j] = r__1 * r__1;
/* L102: */
    }

/*     MULTIPLY RIGHT SIDE BY R(J)**2 */

    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	x = rsq[j];
	i__2 = *m;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    f[i__ + j * f_dim1] = x * f[i__ + j * f_dim1];
/* L103: */
	}
/* L104: */
    }

/*      DEFINE COEFFICIENTS AM,BM,CM */

    x = 1.f / (cos(dth / 2.f) * 2.f);
    i__1 = *m;
    for (i__ = 2; i__ <= i__1; ++i__) {
	am[i__] = (snth[i__ - 1] + snth[i__]) * x;
	cm[i__ - 1] = am[i__];
/* L105: */
    }
    am[1] = sin(*a);
    cm[*m] = sin(*b);
    i__1 = *m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	x = 1.f / snth[i__];
	y = x / dthsq;
	am[i__] *= y;
	cm[i__] *= y;
	bm[i__] = *elmbda * x * x - am[i__] - cm[i__];
/* L106: */
    }

/*     DEFINE COEFFICIENTS AN,BN,CN */

    x = *c__ / dr;
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
/* Computing 2nd power */
	r__1 = x + j - 1;
	an[j] = r__1 * r__1;
/* Computing 2nd power */
	r__1 = x + j;
	cn[j] = r__1 * r__1;
	bn[j] = -(an[j] + cn[j]);
/* L107: */
    }
    isw = 1;
    nb = *nbdcnd;
    if (*c__ == 0.f && nb == 2) {
	nb = 6;
    }

/*     ENTER DATA ON THETA BOUNDARIES */

    switch (*mbdcnd) {
	case 1:  goto L108;
	case 2:  goto L108;
	case 3:  goto L110;
	case 4:  goto L110;
	case 5:  goto L112;
	case 6:  goto L112;
	case 7:  goto L108;
	case 8:  goto L110;
	case 9:  goto L112;
    }
L108:
    bm[1] -= am[1];
    x = am[1] * 2.f;
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	f[j * f_dim1 + 1] -= x * bda[j];
/* L109: */
    }
    goto L112;
L110:
    bm[1] += am[1];
    x = dth * am[1];
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	f[j * f_dim1 + 1] += x * bda[j];
/* L111: */
    }
L112:
    switch (*mbdcnd) {
	case 1:  goto L113;
	case 2:  goto L115;
	case 3:  goto L115;
	case 4:  goto L113;
	case 5:  goto L113;
	case 6:  goto L115;
	case 7:  goto L117;
	case 8:  goto L117;
	case 9:  goto L117;
    }
L113:
    bm[*m] -= cm[*m];
    x = cm[*m] * 2.f;
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	f[*m + j * f_dim1] -= x * bdb[j];
/* L114: */
    }
    goto L117;
L115:
    bm[*m] += cm[*m];
    x = dth * cm[*m];
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	f[*m + j * f_dim1] -= x * bdb[j];
/* L116: */
    }
L117:

/*     ENTER DATA ON R BOUNDARIES */

    switch (nb) {
	case 1:  goto L118;
	case 2:  goto L118;
	case 3:  goto L120;
	case 4:  goto L120;
	case 5:  goto L122;
	case 6:  goto L122;
    }
L118:
    bn[1] -= an[1];
    x = an[1] * 2.f;
    i__1 = *m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	f[i__ + f_dim1] -= x * bdc[i__];
/* L119: */
    }
    goto L122;
L120:
    bn[1] += an[1];
    x = dr * an[1];
    i__1 = *m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	f[i__ + f_dim1] += x * bdc[i__];
/* L121: */
    }
L122:
    switch (nb) {
	case 1:  goto L123;
	case 2:  goto L125;
	case 3:  goto L125;
	case 4:  goto L123;
	case 5:  goto L123;
	case 6:  goto L125;
    }
L123:
    bn[*n] -= cn[*n];
    x = cn[*n] * 2.f;
    i__1 = *m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	f[i__ + *n * f_dim1] -= x * bdd[i__];
/* L124: */
    }
    goto L127;
L125:
    bn[*n] += cn[*n];
    x = dr * cn[*n];
    i__1 = *m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	f[i__ + *n * f_dim1] -= x * bdd[i__];
/* L126: */
    }
L127:

/*     CHECK FOR SINGULAR PROBLEM.  IF SINGULAR, PERTURB F. */

    *pertrb = 0.f;
    switch (*mbdcnd) {
	case 1:  goto L137;
	case 2:  goto L137;
	case 3:  goto L128;
	case 4:  goto L137;
	case 5:  goto L137;
	case 6:  goto L128;
	case 7:  goto L137;
	case 8:  goto L128;
	case 9:  goto L128;
    }
L128:
    switch (nb) {
	case 1:  goto L137;
	case 2:  goto L137;
	case 3:  goto L129;
	case 4:  goto L137;
	case 5:  goto L137;
	case 6:  goto L129;
    }
L129:
    if (*elmbda < 0.f) {
	goto L137;
    } else if (*elmbda == 0) {
	goto L131;
    } else {
	goto L130;
    }
L130:
    *ierr1 = 10;
    goto L137;
L131:
    isw = 2;
    i__1 = *m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	x = 0.f;
	i__2 = *n;
	for (j = 1; j <= i__2; ++j) {
	    x += f[i__ + j * f_dim1];
/* L132: */
	}
	*pertrb += x * snth[i__];
/* L133: */
    }
    x = 0.f;
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	x += rsq[j];
/* L134: */
    }
    *pertrb = *pertrb * sin(dth / 2.f) * 2.f / (x * (cos(*a) - cos(*b)));
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	x = rsq[j] * *pertrb;
	i__2 = *m;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    f[i__ + j * f_dim1] -= x;
/* L135: */
	}
/* L136: */
    }
L137:
    a2 = 0.f;
    i__1 = *m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	a2 += f[i__ + f_dim1];
/* L138: */
    }
    a2 /= rsq[1];

/*     INITIALIZE BLKTRI */

    if (*intl != 0) {
	goto L139;
    }
    blktri_(&c__0, &c__1, n, &an[1], &bn[1], &cn[1], &c__1, m, &am[1], &bm[1],
	     &cm[1], idimf, &f[f_offset], ierr1, &wrk[1]);
L139:

/*     CALL BLKTRI TO SOLVE SYSTEM OF EQUATIONS. */

    blktri_(&c__1, &c__1, n, &an[1], &bn[1], &cn[1], &c__1, m, &am[1], &bm[1],
	     &cm[1], idimf, &f[f_offset], ierr1, &wrk[1]);
    if (isw != 2 || *c__ != 0.f || *nbdcnd != 2) {
	goto L143;
    }
    a1 = 0.f;
    a3 = 0.f;
    i__1 = *m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	a1 += snth[i__] * f[i__ + f_dim1];
	a3 += snth[i__];
/* L140: */
    }
    a1 += rsq[1] * a2 / 2.f;
    if (*mbdcnd == 3) {
	a1 += (sin(*b) * bdb[1] - sin(*a) * bda[1]) / ((*b - *a) * 2.f);
    }
    a1 /= a3;
    a1 = bdc[1] - a1;
    i__1 = *m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = *n;
	for (j = 1; j <= i__2; ++j) {
	    f[i__ + j * f_dim1] += a1;
/* L141: */
	}
/* L142: */
    }
L143:
    return 0;
} /* hstcs1_ */

