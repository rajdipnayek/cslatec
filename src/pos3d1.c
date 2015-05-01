/* pos3d1.f -- translated by f2c (version 12.02.01).
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

/* DECK POS3D1 */
/* Subroutine */ int pos3d1_(integer *lp, integer *l, integer *mp, integer *m,
	 integer *n, real *a, real *b, real *c__, integer *ldimf, integer *
	mdimf, real *f, real *xrt, real *yrt, real *t, real *d__, real *wx, 
	real *wy, real *c1, real *c2, real *bb)
{
    /* System generated locals */
    integer f_dim1, f_dim2, f_offset, i__1, i__2, i__3;
    real r__1;

    /* Local variables */
    static integer i__, j, k;
    static real di, dj, pi;
    static integer lr, mr, nr;
    static real dx, dy, dum;
    extern /* Subroutine */ int cost_(integer *, real *, real *), sint_(
	    integer *, real *, real *);
    static integer lrdel, mrdel;
    extern /* Subroutine */ int rfftb_(integer *, real *, real *), rfftf_(
	    integer *, real *, real *), cosqb_(integer *, real *, real *);
    static real scalx;
    extern /* Subroutine */ int rffti_(integer *, real *);
    static real scaly;
    static integer ifwrd;
    extern /* Subroutine */ int cosqi_(integer *, real *), sinqb_(integer *, 
	    real *, real *), sinqf_(integer *, real *, real *), costi_(
	    integer *, real *), cosqf_(integer *, real *, real *), sinqi_(
	    integer *, real *), tridq_(integer *, real *, real *, real *, 
	    real *, real *), sinti_(integer *, real *);
    extern doublereal pimach_(real *);

/* ***BEGIN PROLOGUE  POS3D1 */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to POIS3D */
/* ***LIBRARY   SLATEC */
/* ***TYPE      SINGLE PRECISION (POS3D1-S) */
/* ***AUTHOR  (UNKNOWN) */
/* ***SEE ALSO  POIS3D */
/* ***ROUTINES CALLED  COSQB, COSQF, COSQI, COST, COSTI, PIMACH, RFFTB, */
/*                    RFFTF, RFFTI, SINQB, SINQF, SINQI, SINT, SINTI, */
/*                    TRIDQ */
/* ***REVISION HISTORY  (YYMMDD) */
/*   801001  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   891009  Removed unreferenced variable.  (WRB) */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900308  Changed call to TRID to call to TRIDQ.  (WRB) */
/*   900402  Added TYPE section.  (WRB) */
/* ***END PROLOGUE  POS3D1 */
/* ***FIRST EXECUTABLE STATEMENT  POS3D1 */
    /* Parameter adjustments */
    --a;
    --b;
    --c__;
    f_dim1 = *ldimf;
    f_dim2 = *mdimf;
    f_offset = 1 + f_dim1 * (1 + f_dim2);
    f -= f_offset;
    --xrt;
    --yrt;
    --t;
    --d__;
    --wx;
    --wy;
    --bb;

    /* Function Body */
    pi = pimach_(&dum);
    lr = *l;
    mr = *m;
    nr = *n;

/*     GENERATE TRANSFORM ROOTS */

    lrdel = (*lp - 1) * (*lp - 3) * (*lp - 5) / 3;
    scalx = (real) (lr + lrdel);
    dx = pi / (scalx * 2.f);
    switch (*lp) {
	case 1:  goto L108;
	case 2:  goto L103;
	case 3:  goto L101;
	case 4:  goto L102;
	case 5:  goto L101;
    }
L101:
    di = .5f;
    scalx *= 2.f;
    goto L104;
L102:
    di = 1.f;
    goto L104;
L103:
    di = 0.f;
L104:
    i__1 = lr;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing 2nd power */
	r__1 = sin((i__ - di) * dx);
	xrt[i__] = *c1 * -4.f * (r__1 * r__1);
/* L105: */
    }
    scalx *= 2.f;
    switch (*lp) {
	case 1:  goto L112;
	case 2:  goto L106;
	case 3:  goto L110;
	case 4:  goto L107;
	case 5:  goto L111;
    }
L106:
    sinti_(&lr, &wx[1]);
    goto L112;
L107:
    costi_(&lr, &wx[1]);
    goto L112;
L108:
    xrt[1] = 0.f;
    xrt[lr] = *c1 * -4.f;
    i__1 = lr;
    for (i__ = 3; i__ <= i__1; i__ += 2) {
/* Computing 2nd power */
	r__1 = sin((i__ - 1) * dx);
	xrt[i__ - 1] = *c1 * -4.f * (r__1 * r__1);
	xrt[i__] = xrt[i__ - 1];
/* L109: */
    }
    rffti_(&lr, &wx[1]);
    goto L112;
L110:
    sinqi_(&lr, &wx[1]);
    goto L112;
L111:
    cosqi_(&lr, &wx[1]);
L112:
    mrdel = (*mp - 1) * (*mp - 3) * (*mp - 5) / 3;
    scaly = (real) (mr + mrdel);
    dy = pi / (scaly * 2.f);
    switch (*mp) {
	case 1:  goto L120;
	case 2:  goto L115;
	case 3:  goto L113;
	case 4:  goto L114;
	case 5:  goto L113;
    }
L113:
    dj = .5f;
    scaly *= 2.f;
    goto L116;
L114:
    dj = 1.f;
    goto L116;
L115:
    dj = 0.f;
L116:
    i__1 = mr;
    for (j = 1; j <= i__1; ++j) {
/* Computing 2nd power */
	r__1 = sin((j - dj) * dy);
	yrt[j] = *c2 * -4.f * (r__1 * r__1);
/* L117: */
    }
    scaly *= 2.f;
    switch (*mp) {
	case 1:  goto L124;
	case 2:  goto L118;
	case 3:  goto L122;
	case 4:  goto L119;
	case 5:  goto L123;
    }
L118:
    sinti_(&mr, &wy[1]);
    goto L124;
L119:
    costi_(&mr, &wy[1]);
    goto L124;
L120:
    yrt[1] = 0.f;
    yrt[mr] = *c2 * -4.f;
    i__1 = mr;
    for (j = 3; j <= i__1; j += 2) {
/* Computing 2nd power */
	r__1 = sin((j - 1) * dy);
	yrt[j - 1] = *c2 * -4.f * (r__1 * r__1);
	yrt[j] = yrt[j - 1];
/* L121: */
    }
    rffti_(&mr, &wy[1]);
    goto L124;
L122:
    sinqi_(&mr, &wy[1]);
    goto L124;
L123:
    cosqi_(&mr, &wy[1]);
L124:
    ifwrd = 1;
L125:

/*     TRANSFORM X */

    i__1 = mr;
    for (j = 1; j <= i__1; ++j) {
	i__2 = nr;
	for (k = 1; k <= i__2; ++k) {
	    i__3 = lr;
	    for (i__ = 1; i__ <= i__3; ++i__) {
		t[i__] = f[i__ + (j + k * f_dim2) * f_dim1];
/* L126: */
	    }
	    switch (*lp) {
		case 1:  goto L127;
		case 2:  goto L130;
		case 3:  goto L131;
		case 4:  goto L134;
		case 5:  goto L135;
	    }
L127:
	    switch (ifwrd) {
		case 1:  goto L128;
		case 2:  goto L129;
	    }
L128:
	    rfftf_(&lr, &t[1], &wx[1]);
	    goto L138;
L129:
	    rfftb_(&lr, &t[1], &wx[1]);
	    goto L138;
L130:
	    sint_(&lr, &t[1], &wx[1]);
	    goto L138;
L131:
	    switch (ifwrd) {
		case 1:  goto L132;
		case 2:  goto L133;
	    }
L132:
	    sinqf_(&lr, &t[1], &wx[1]);
	    goto L138;
L133:
	    sinqb_(&lr, &t[1], &wx[1]);
	    goto L138;
L134:
	    cost_(&lr, &t[1], &wx[1]);
	    goto L138;
L135:
	    switch (ifwrd) {
		case 1:  goto L136;
		case 2:  goto L137;
	    }
L136:
	    cosqf_(&lr, &t[1], &wx[1]);
	    goto L138;
L137:
	    cosqb_(&lr, &t[1], &wx[1]);
L138:
	    i__3 = lr;
	    for (i__ = 1; i__ <= i__3; ++i__) {
		f[i__ + (j + k * f_dim2) * f_dim1] = t[i__];
/* L139: */
	    }
/* L140: */
	}
/* L141: */
    }
    switch (ifwrd) {
	case 1:  goto L142;
	case 2:  goto L164;
    }

/*     TRANSFORM Y */

L142:
    i__1 = lr;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = nr;
	for (k = 1; k <= i__2; ++k) {
	    i__3 = mr;
	    for (j = 1; j <= i__3; ++j) {
		t[j] = f[i__ + (j + k * f_dim2) * f_dim1];
/* L143: */
	    }
	    switch (*mp) {
		case 1:  goto L144;
		case 2:  goto L147;
		case 3:  goto L148;
		case 4:  goto L151;
		case 5:  goto L152;
	    }
L144:
	    switch (ifwrd) {
		case 1:  goto L145;
		case 2:  goto L146;
	    }
L145:
	    rfftf_(&mr, &t[1], &wy[1]);
	    goto L155;
L146:
	    rfftb_(&mr, &t[1], &wy[1]);
	    goto L155;
L147:
	    sint_(&mr, &t[1], &wy[1]);
	    goto L155;
L148:
	    switch (ifwrd) {
		case 1:  goto L149;
		case 2:  goto L150;
	    }
L149:
	    sinqf_(&mr, &t[1], &wy[1]);
	    goto L155;
L150:
	    sinqb_(&mr, &t[1], &wy[1]);
	    goto L155;
L151:
	    cost_(&mr, &t[1], &wy[1]);
	    goto L155;
L152:
	    switch (ifwrd) {
		case 1:  goto L153;
		case 2:  goto L154;
	    }
L153:
	    cosqf_(&mr, &t[1], &wy[1]);
	    goto L155;
L154:
	    cosqb_(&mr, &t[1], &wy[1]);
L155:
	    i__3 = mr;
	    for (j = 1; j <= i__3; ++j) {
		f[i__ + (j + k * f_dim2) * f_dim1] = t[j];
/* L156: */
	    }
/* L157: */
	}
/* L158: */
    }
    switch (ifwrd) {
	case 1:  goto L159;
	case 2:  goto L125;
    }
L159:

/*     SOLVE TRIDIAGONAL SYSTEMS IN Z */

    i__1 = lr;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = mr;
	for (j = 1; j <= i__2; ++j) {
	    i__3 = nr;
	    for (k = 1; k <= i__3; ++k) {
		bb[k] = b[k] + xrt[i__] + yrt[j];
		t[k] = f[i__ + (j + k * f_dim2) * f_dim1];
/* L160: */
	    }
	    tridq_(&nr, &a[1], &bb[1], &c__[1], &t[1], &d__[1]);
	    i__3 = nr;
	    for (k = 1; k <= i__3; ++k) {
		f[i__ + (j + k * f_dim2) * f_dim1] = t[k];
/* L161: */
	    }
/* L162: */
	}
/* L163: */
    }
    ifwrd = 2;
    goto L142;
L164:
    i__1 = lr;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = mr;
	for (j = 1; j <= i__2; ++j) {
	    i__3 = nr;
	    for (k = 1; k <= i__3; ++k) {
		f[i__ + (j + k * f_dim2) * f_dim1] /= scalx * scaly;
/* L165: */
	    }
/* L166: */
	}
/* L167: */
    }
    return 0;
} /* pos3d1_ */

