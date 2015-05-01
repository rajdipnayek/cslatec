/* poisn2.f -- translated by f2c (version 12.02.01).
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
static real c_b11 = .5f;
static real c_b12 = 0.f;
static integer c__0 = 0;
static real c_b120 = 1.f;
static integer c__2 = 2;

/* DECK POISN2 */
/* Subroutine */ int poisn2_(integer *m, integer *n, integer *istag, integer *
	mixbnd, real *a, real *bb, real *c__, real *q, integer *idimq, real *
	b, real *b2, real *b3, real *w, real *w2, real *w3, real *d__, real *
	tcos, real *p)
{
    /* System generated locals */
    integer q_dim1, q_offset, i__1, i__2, i__3;
    static integer equiv_3[4];

    /* Local variables */
    static integer i__, j;
#define k (equiv_3)
    static real t;
    static integer i1, i2;
#define k1 (equiv_3)
#define k2 (equiv_3 + 1)
#define k3 (equiv_3 + 2)
#define k4 (equiv_3 + 3)
    static real fi;
    static integer ii, ip, jr, kr, lr, mr, nr, jm1, jm2, jm3, jp1, jp2, i2r, 
	    jp3, jr2;
    extern /* Subroutine */ int tri3_(integer *, real *, real *, real *, 
	    integer *, real *, real *, real *, real *, real *, real *, real *,
	     real *);
    static real fden;
    static integer nrod;
    static real fnum;
    extern /* Subroutine */ int trix_(integer *, integer *, integer *, real *,
	     real *, real *, real *, real *, real *, real *);
    static integer i2rby2, nlast, jstep, jstop;
    extern /* Subroutine */ int s1merg_(real *, integer *, integer *, integer 
	    *, integer *, integer *);
    static real fistag;
    extern /* Subroutine */ int cosgen_(integer *, integer *, real *, real *, 
	    real *);
    static integer nlastp, nrodpr, jstart, ipstor;

/* ***BEGIN PROLOGUE  POISN2 */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to GENBUN */
/* ***LIBRARY   SLATEC */
/* ***TYPE      SINGLE PRECISION (POISN2-S, CMPOSN-C) */
/* ***AUTHOR  (UNKNOWN) */
/* ***DESCRIPTION */

/*     Subroutine to solve Poisson's equation with Neumann boundary */
/*     conditions. */

/*     ISTAG = 1 if the last diagonal block is A. */
/*     ISTAG = 2 if the last diagonal block is A-I. */
/*     MIXBND = 1 if have Neumann boundary conditions at both boundaries. */
/*     MIXBND = 2 if have Neumann boundary conditions at bottom and */
/*     Dirichlet condition at top.  (for this case, must have ISTAG = 1.) */

/* ***SEE ALSO  GENBUN */
/* ***ROUTINES CALLED  COSGEN, S1MERG, TRI3, TRIX */
/* ***REVISION HISTORY  (YYMMDD) */
/*   801001  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900402  Added TYPE section.  (WRB) */
/*   920130  Modified to use merge routine S1MERG rather than deleted */
/*           routine MERGE.  (WRB) */
/* ***END PROLOGUE  POISN2 */

/* ***FIRST EXECUTABLE STATEMENT  POISN2 */
    /* Parameter adjustments */
    --a;
    --bb;
    --c__;
    q_dim1 = *idimq;
    q_offset = 1 + q_dim1;
    q -= q_offset;
    --b;
    --b2;
    --b3;
    --w;
    --w2;
    --w3;
    --d__;
    --tcos;
    --p;

    /* Function Body */
    fistag = (real) (3 - *istag);
    fnum = 1.f / *istag;
    fden = (*istag - 1) * .5f;
    mr = *m;
    ip = -mr;
    ipstor = 0;
    i2r = 1;
    jr = 2;
    nr = *n;
    nlast = *n;
    kr = 1;
    lr = 0;
    switch (*istag) {
	case 1:  goto L101;
	case 2:  goto L103;
    }
L101:
    i__1 = mr;
    for (i__ = 1; i__ <= i__1; ++i__) {
	q[i__ + *n * q_dim1] *= .5f;
/* L102: */
    }
    switch (*mixbnd) {
	case 1:  goto L103;
	case 2:  goto L104;
    }
L103:
    if (*n <= 3) {
	goto L155;
    }
L104:
    jr = i2r << 1;
    nrod = 1;
    if (nr / 2 << 1 == nr) {
	nrod = 0;
    }
    switch (*mixbnd) {
	case 1:  goto L105;
	case 2:  goto L106;
    }
L105:
    jstart = 1;
    goto L107;
L106:
    jstart = jr;
    nrod = 1 - nrod;
L107:
    jstop = nlast - jr;
    if (nrod == 0) {
	jstop -= i2r;
    }
    cosgen_(&i2r, &c__1, &c_b11, &c_b12, &tcos[1]);
    i2rby2 = i2r / 2;
    if (jstop >= jstart) {
	goto L108;
    }
    j = jr;
    goto L116;
L108:

/*     REGULAR REDUCTION. */

    i__1 = jstop;
    i__2 = jr;
    for (j = jstart; i__2 < 0 ? j >= i__1 : j <= i__1; j += i__2) {
	jp1 = j + i2rby2;
	jp2 = j + i2r;
	jp3 = jp2 + i2rby2;
	jm1 = j - i2rby2;
	jm2 = j - i2r;
	jm3 = jm2 - i2rby2;
	if (j != 1) {
	    goto L109;
	}
	jm1 = jp1;
	jm2 = jp2;
	jm3 = jp3;
L109:
	if (i2r != 1) {
	    goto L111;
	}
	if (j == 1) {
	    jm2 = jp2;
	}
	i__3 = mr;
	for (i__ = 1; i__ <= i__3; ++i__) {
	    b[i__] = q[i__ + j * q_dim1] * 2.f;
	    q[i__ + j * q_dim1] = q[i__ + jm2 * q_dim1] + q[i__ + jp2 * 
		    q_dim1];
/* L110: */
	}
	goto L113;
L111:
	i__3 = mr;
	for (i__ = 1; i__ <= i__3; ++i__) {
	    fi = q[i__ + j * q_dim1];
	    q[i__ + j * q_dim1] = q[i__ + j * q_dim1] - q[i__ + jm1 * q_dim1] 
		    - q[i__ + jp1 * q_dim1] + q[i__ + jm2 * q_dim1] + q[i__ + 
		    jp2 * q_dim1];
	    b[i__] = fi + q[i__ + j * q_dim1] - q[i__ + jm3 * q_dim1] - q[i__ 
		    + jp3 * q_dim1];
/* L112: */
	}
L113:
	trix_(&i2r, &c__0, &mr, &a[1], &bb[1], &c__[1], &b[1], &tcos[1], &d__[
		1], &w[1]);
	i__3 = mr;
	for (i__ = 1; i__ <= i__3; ++i__) {
	    q[i__ + j * q_dim1] += b[i__];
/* L114: */
	}

/*     END OF REDUCTION FOR REGULAR UNKNOWNS. */

/* L115: */
    }

/*     BEGIN SPECIAL REDUCTION FOR LAST UNKNOWN. */

    j = jstop + jr;
L116:
    nlast = j;
    jm1 = j - i2rby2;
    jm2 = j - i2r;
    jm3 = jm2 - i2rby2;
    if (nrod == 0) {
	goto L128;
    }

/*     ODD NUMBER OF UNKNOWNS */

    if (i2r != 1) {
	goto L118;
    }
    i__2 = mr;
    for (i__ = 1; i__ <= i__2; ++i__) {
	b[i__] = fistag * q[i__ + j * q_dim1];
	q[i__ + j * q_dim1] = q[i__ + jm2 * q_dim1];
/* L117: */
    }
    goto L126;
L118:
    i__2 = mr;
    for (i__ = 1; i__ <= i__2; ++i__) {
	b[i__] = q[i__ + j * q_dim1] + (q[i__ + jm2 * q_dim1] - q[i__ + jm1 * 
		q_dim1] - q[i__ + jm3 * q_dim1]) * .5f;
/* L119: */
    }
    if (nrodpr != 0) {
	goto L121;
    }
    i__2 = mr;
    for (i__ = 1; i__ <= i__2; ++i__) {
	ii = ip + i__;
	q[i__ + j * q_dim1] = q[i__ + jm2 * q_dim1] + p[ii];
/* L120: */
    }
    ip -= mr;
    goto L123;
L121:
    i__2 = mr;
    for (i__ = 1; i__ <= i__2; ++i__) {
	q[i__ + j * q_dim1] = q[i__ + j * q_dim1] - q[i__ + jm1 * q_dim1] + q[
		i__ + jm2 * q_dim1];
/* L122: */
    }
L123:
    if (lr == 0) {
	goto L124;
    }
    cosgen_(&lr, &c__1, &c_b11, &fden, &tcos[kr + 1]);
    goto L126;
L124:
    i__2 = mr;
    for (i__ = 1; i__ <= i__2; ++i__) {
	b[i__] = fistag * b[i__];
/* L125: */
    }
L126:
    cosgen_(&kr, &c__1, &c_b11, &fden, &tcos[1]);
    trix_(&kr, &lr, &mr, &a[1], &bb[1], &c__[1], &b[1], &tcos[1], &d__[1], &w[
	    1]);
    i__2 = mr;
    for (i__ = 1; i__ <= i__2; ++i__) {
	q[i__ + j * q_dim1] += b[i__];
/* L127: */
    }
    kr += i2r;
    goto L151;
L128:

/*     EVEN NUMBER OF UNKNOWNS */

    jp1 = j + i2rby2;
    jp2 = j + i2r;
    if (i2r != 1) {
	goto L135;
    }
    i__2 = mr;
    for (i__ = 1; i__ <= i__2; ++i__) {
	b[i__] = q[i__ + j * q_dim1];
/* L129: */
    }
    trix_(&c__1, &c__0, &mr, &a[1], &bb[1], &c__[1], &b[1], &tcos[1], &d__[1],
	     &w[1]);
    ip = 0;
    ipstor = mr;
    switch (*istag) {
	case 1:  goto L133;
	case 2:  goto L130;
    }
L130:
    i__2 = mr;
    for (i__ = 1; i__ <= i__2; ++i__) {
	p[i__] = b[i__];
	b[i__] += q[i__ + *n * q_dim1];
/* L131: */
    }
    tcos[1] = 1.f;
    tcos[2] = 0.f;
    trix_(&c__1, &c__1, &mr, &a[1], &bb[1], &c__[1], &b[1], &tcos[1], &d__[1],
	     &w[1]);
    i__2 = mr;
    for (i__ = 1; i__ <= i__2; ++i__) {
	q[i__ + j * q_dim1] = q[i__ + jm2 * q_dim1] + p[i__] + b[i__];
/* L132: */
    }
    goto L150;
L133:
    i__2 = mr;
    for (i__ = 1; i__ <= i__2; ++i__) {
	p[i__] = b[i__];
	q[i__ + j * q_dim1] = q[i__ + jm2 * q_dim1] + q[i__ + jp2 * q_dim1] * 
		2.f + b[i__] * 3.f;
/* L134: */
    }
    goto L150;
L135:
    i__2 = mr;
    for (i__ = 1; i__ <= i__2; ++i__) {
	b[i__] = q[i__ + j * q_dim1] + (q[i__ + jm2 * q_dim1] - q[i__ + jm1 * 
		q_dim1] - q[i__ + jm3 * q_dim1]) * .5f;
/* L136: */
    }
    if (nrodpr != 0) {
	goto L138;
    }
    i__2 = mr;
    for (i__ = 1; i__ <= i__2; ++i__) {
	ii = ip + i__;
	b[i__] += p[ii];
/* L137: */
    }
    goto L140;
L138:
    i__2 = mr;
    for (i__ = 1; i__ <= i__2; ++i__) {
	b[i__] = b[i__] + q[i__ + jp2 * q_dim1] - q[i__ + jp1 * q_dim1];
/* L139: */
    }
L140:
    trix_(&i2r, &c__0, &mr, &a[1], &bb[1], &c__[1], &b[1], &tcos[1], &d__[1], 
	    &w[1]);
    ip += mr;
/* Computing MAX */
    i__2 = ipstor, i__1 = ip + mr;
    ipstor = max(i__2,i__1);
    i__2 = mr;
    for (i__ = 1; i__ <= i__2; ++i__) {
	ii = ip + i__;
	p[ii] = b[i__] + (q[i__ + j * q_dim1] - q[i__ + jm1 * q_dim1] - q[i__ 
		+ jp1 * q_dim1]) * .5f;
	b[i__] = p[ii] + q[i__ + jp2 * q_dim1];
/* L141: */
    }
    if (lr == 0) {
	goto L142;
    }
    cosgen_(&lr, &c__1, &c_b11, &fden, &tcos[i2r + 1]);
    s1merg_(&tcos[1], &c__0, &i2r, &i2r, &lr, &kr);
    goto L144;
L142:
    i__2 = i2r;
    for (i__ = 1; i__ <= i__2; ++i__) {
	ii = kr + i__;
	tcos[ii] = tcos[i__];
/* L143: */
    }
L144:
    cosgen_(&kr, &c__1, &c_b11, &fden, &tcos[1]);
    if (lr != 0) {
	goto L145;
    }
    switch (*istag) {
	case 1:  goto L146;
	case 2:  goto L145;
    }
L145:
    trix_(&kr, &kr, &mr, &a[1], &bb[1], &c__[1], &b[1], &tcos[1], &d__[1], &w[
	    1]);
    goto L148;
L146:
    i__2 = mr;
    for (i__ = 1; i__ <= i__2; ++i__) {
	b[i__] = fistag * b[i__];
/* L147: */
    }
L148:
    i__2 = mr;
    for (i__ = 1; i__ <= i__2; ++i__) {
	ii = ip + i__;
	q[i__ + j * q_dim1] = q[i__ + jm2 * q_dim1] + p[ii] + b[i__];
/* L149: */
    }
L150:
    lr = kr;
    kr += jr;
L151:
    switch (*mixbnd) {
	case 1:  goto L152;
	case 2:  goto L153;
    }
L152:
    nr = (nlast - 1) / jr + 1;
    if (nr <= 3) {
	goto L155;
    }
    goto L154;
L153:
    nr = nlast / jr;
    if (nr <= 1) {
	goto L192;
    }
L154:
    i2r = jr;
    nrodpr = nrod;
    goto L104;
L155:

/*      BEGIN SOLUTION */

    j = jr + 1;
    jm1 = j - i2r;
    jp1 = j + i2r;
    jm2 = nlast - i2r;
    if (nr == 2) {
	goto L184;
    }
    if (lr != 0) {
	goto L170;
    }
    if (*n != 3) {
	goto L161;
    }

/*     CASE N = 3. */

    switch (*istag) {
	case 1:  goto L156;
	case 2:  goto L168;
    }
L156:
    i__2 = mr;
    for (i__ = 1; i__ <= i__2; ++i__) {
	b[i__] = q[i__ + (q_dim1 << 1)];
/* L157: */
    }
    tcos[1] = 0.f;
    trix_(&c__1, &c__0, &mr, &a[1], &bb[1], &c__[1], &b[1], &tcos[1], &d__[1],
	     &w[1]);
    i__2 = mr;
    for (i__ = 1; i__ <= i__2; ++i__) {
	q[i__ + (q_dim1 << 1)] = b[i__];
	b[i__] = b[i__] * 4.f + q[i__ + q_dim1] + q[i__ + q_dim1 * 3] * 2.f;
/* L158: */
    }
    tcos[1] = -2.f;
    tcos[2] = 2.f;
    i1 = 2;
    i2 = 0;
    trix_(&i1, &i2, &mr, &a[1], &bb[1], &c__[1], &b[1], &tcos[1], &d__[1], &w[
	    1]);
    i__2 = mr;
    for (i__ = 1; i__ <= i__2; ++i__) {
	q[i__ + (q_dim1 << 1)] += b[i__];
	b[i__] = q[i__ + q_dim1] + q[i__ + (q_dim1 << 1)] * 2.f;
/* L159: */
    }
    tcos[1] = 0.f;
    trix_(&c__1, &c__0, &mr, &a[1], &bb[1], &c__[1], &b[1], &tcos[1], &d__[1],
	     &w[1]);
    i__2 = mr;
    for (i__ = 1; i__ <= i__2; ++i__) {
	q[i__ + q_dim1] = b[i__];
/* L160: */
    }
    jr = 1;
    i2r = 0;
    goto L194;

/*     CASE N = 2**P+1 */

L161:
    switch (*istag) {
	case 1:  goto L162;
	case 2:  goto L170;
    }
L162:
    i__2 = mr;
    for (i__ = 1; i__ <= i__2; ++i__) {
	b[i__] = q[i__ + j * q_dim1] + q[i__ + q_dim1] * .5f - q[i__ + jm1 * 
		q_dim1] + q[i__ + nlast * q_dim1] - q[i__ + jm2 * q_dim1];
/* L163: */
    }
    cosgen_(&jr, &c__1, &c_b11, &c_b12, &tcos[1]);
    trix_(&jr, &c__0, &mr, &a[1], &bb[1], &c__[1], &b[1], &tcos[1], &d__[1], &
	    w[1]);
    i__2 = mr;
    for (i__ = 1; i__ <= i__2; ++i__) {
	q[i__ + j * q_dim1] = (q[i__ + j * q_dim1] - q[i__ + jm1 * q_dim1] - 
		q[i__ + jp1 * q_dim1]) * .5f + b[i__];
	b[i__] = q[i__ + q_dim1] + q[i__ + nlast * q_dim1] * 2.f + q[i__ + j *
		 q_dim1] * 4.f;
/* L164: */
    }
    jr2 = jr << 1;
    cosgen_(&jr, &c__1, &c_b12, &c_b12, &tcos[1]);
    i__2 = jr;
    for (i__ = 1; i__ <= i__2; ++i__) {
	i1 = jr + i__;
	i2 = jr + 1 - i__;
	tcos[i1] = -tcos[i2];
/* L165: */
    }
    trix_(&jr2, &c__0, &mr, &a[1], &bb[1], &c__[1], &b[1], &tcos[1], &d__[1], 
	    &w[1]);
    i__2 = mr;
    for (i__ = 1; i__ <= i__2; ++i__) {
	q[i__ + j * q_dim1] += b[i__];
	b[i__] = q[i__ + q_dim1] + q[i__ + j * q_dim1] * 2.f;
/* L166: */
    }
    cosgen_(&jr, &c__1, &c_b11, &c_b12, &tcos[1]);
    trix_(&jr, &c__0, &mr, &a[1], &bb[1], &c__[1], &b[1], &tcos[1], &d__[1], &
	    w[1]);
    i__2 = mr;
    for (i__ = 1; i__ <= i__2; ++i__) {
	q[i__ + q_dim1] = q[i__ + q_dim1] * .5f - q[i__ + jm1 * q_dim1] + b[
		i__];
/* L167: */
    }
    goto L194;

/*     CASE OF GENERAL N WITH NR = 3 . */

L168:
    i__2 = mr;
    for (i__ = 1; i__ <= i__2; ++i__) {
	b[i__] = q[i__ + (q_dim1 << 1)];
	q[i__ + (q_dim1 << 1)] = 0.f;
	b2[i__] = q[i__ + q_dim1 * 3];
	b3[i__] = q[i__ + q_dim1];
/* L169: */
    }
    jr = 1;
    i2r = 0;
    j = 2;
    goto L177;
L170:
    i__2 = mr;
    for (i__ = 1; i__ <= i__2; ++i__) {
	b[i__] = q[i__ + q_dim1] * .5f - q[i__ + jm1 * q_dim1] + q[i__ + j * 
		q_dim1];
/* L171: */
    }
    if (nrod != 0) {
	goto L173;
    }
    i__2 = mr;
    for (i__ = 1; i__ <= i__2; ++i__) {
	ii = ip + i__;
	b[i__] += p[ii];
/* L172: */
    }
    goto L175;
L173:
    i__2 = mr;
    for (i__ = 1; i__ <= i__2; ++i__) {
	b[i__] = b[i__] + q[i__ + nlast * q_dim1] - q[i__ + jm2 * q_dim1];
/* L174: */
    }
L175:
    i__2 = mr;
    for (i__ = 1; i__ <= i__2; ++i__) {
	t = (q[i__ + j * q_dim1] - q[i__ + jm1 * q_dim1] - q[i__ + jp1 * 
		q_dim1]) * .5f;
	q[i__ + j * q_dim1] = t;
	b2[i__] = q[i__ + nlast * q_dim1] + t;
	b3[i__] = q[i__ + q_dim1] + t * 2.f;
/* L176: */
    }
L177:
    *k1 = kr + (jr << 1) - 1;
    *k2 = kr + jr;
    tcos[*k1 + 1] = -2.f;
    *k4 = *k1 + 3 - *istag;
    i__2 = *k2 + *istag - 2;
    cosgen_(&i__2, &c__1, &c_b12, &fnum, &tcos[*k4]);
    *k4 = *k1 + *k2 + 1;
    i__2 = jr - 1;
    cosgen_(&i__2, &c__1, &c_b12, &c_b120, &tcos[*k4]);
    i__2 = *k1 + *k2;
    i__1 = jr - 1;
    s1merg_(&tcos[1], k1, k2, &i__2, &i__1, &c__0);
    *k3 = *k1 + *k2 + lr;
    cosgen_(&jr, &c__1, &c_b11, &c_b12, &tcos[*k3 + 1]);
    *k4 = *k3 + jr + 1;
    cosgen_(&kr, &c__1, &c_b11, &fden, &tcos[*k4]);
    i__2 = *k3 + jr;
    s1merg_(&tcos[1], k3, &jr, &i__2, &kr, k1);
    if (lr == 0) {
	goto L178;
    }
    cosgen_(&lr, &c__1, &c_b11, &fden, &tcos[*k4]);
    i__2 = *k3 + jr;
    i__1 = *k3 - lr;
    s1merg_(&tcos[1], k3, &jr, &i__2, &lr, &i__1);
    cosgen_(&kr, &c__1, &c_b11, &fden, &tcos[*k4]);
L178:
    *k3 = kr;
    *k4 = kr;
    tri3_(&mr, &a[1], &bb[1], &c__[1], k, &b[1], &b2[1], &b3[1], &tcos[1], &
	    d__[1], &w[1], &w2[1], &w3[1]);
    i__2 = mr;
    for (i__ = 1; i__ <= i__2; ++i__) {
	b[i__] = b[i__] + b2[i__] + b3[i__];
/* L179: */
    }
    tcos[1] = 2.f;
    trix_(&c__1, &c__0, &mr, &a[1], &bb[1], &c__[1], &b[1], &tcos[1], &d__[1],
	     &w[1]);
    i__2 = mr;
    for (i__ = 1; i__ <= i__2; ++i__) {
	q[i__ + j * q_dim1] += b[i__];
	b[i__] = q[i__ + q_dim1] + q[i__ + j * q_dim1] * 2.f;
/* L180: */
    }
    cosgen_(&jr, &c__1, &c_b11, &c_b12, &tcos[1]);
    trix_(&jr, &c__0, &mr, &a[1], &bb[1], &c__[1], &b[1], &tcos[1], &d__[1], &
	    w[1]);
    if (jr != 1) {
	goto L182;
    }
    i__2 = mr;
    for (i__ = 1; i__ <= i__2; ++i__) {
	q[i__ + q_dim1] = b[i__];
/* L181: */
    }
    goto L194;
L182:
    i__2 = mr;
    for (i__ = 1; i__ <= i__2; ++i__) {
	q[i__ + q_dim1] = q[i__ + q_dim1] * .5f - q[i__ + jm1 * q_dim1] + b[
		i__];
/* L183: */
    }
    goto L194;
L184:
    if (*n != 2) {
	goto L188;
    }

/*     CASE  N = 2 */

    i__2 = mr;
    for (i__ = 1; i__ <= i__2; ++i__) {
	b[i__] = q[i__ + q_dim1];
/* L185: */
    }
    tcos[1] = 0.f;
    trix_(&c__1, &c__0, &mr, &a[1], &bb[1], &c__[1], &b[1], &tcos[1], &d__[1],
	     &w[1]);
    i__2 = mr;
    for (i__ = 1; i__ <= i__2; ++i__) {
	q[i__ + q_dim1] = b[i__];
	b[i__] = (q[i__ + (q_dim1 << 1)] + b[i__]) * 2.f * fistag;
/* L186: */
    }
    tcos[1] = -fistag;
    tcos[2] = 2.f;
    trix_(&c__2, &c__0, &mr, &a[1], &bb[1], &c__[1], &b[1], &tcos[1], &d__[1],
	     &w[1]);
    i__2 = mr;
    for (i__ = 1; i__ <= i__2; ++i__) {
	q[i__ + q_dim1] += b[i__];
/* L187: */
    }
    jr = 1;
    i2r = 0;
    goto L194;
L188:

/*     CASE OF GENERAL N AND NR = 2 . */

    i__2 = mr;
    for (i__ = 1; i__ <= i__2; ++i__) {
	ii = ip + i__;
	b3[i__] = 0.f;
	b[i__] = q[i__ + q_dim1] + p[ii] * 2.f;
	q[i__ + q_dim1] = q[i__ + q_dim1] * .5f - q[i__ + jm1 * q_dim1];
	b2[i__] = (q[i__ + q_dim1] + q[i__ + nlast * q_dim1]) * 2.f;
/* L189: */
    }
    *k1 = kr + jr - 1;
    tcos[*k1 + 1] = -2.f;
    *k4 = *k1 + 3 - *istag;
    i__2 = kr + *istag - 2;
    cosgen_(&i__2, &c__1, &c_b12, &fnum, &tcos[*k4]);
    *k4 = *k1 + kr + 1;
    i__2 = jr - 1;
    cosgen_(&i__2, &c__1, &c_b12, &c_b120, &tcos[*k4]);
    i__2 = *k1 + kr;
    i__1 = jr - 1;
    s1merg_(&tcos[1], k1, &kr, &i__2, &i__1, &c__0);
    cosgen_(&kr, &c__1, &c_b11, &fden, &tcos[*k1 + 1]);
    *k2 = kr;
    *k4 = *k1 + *k2 + 1;
    cosgen_(&lr, &c__1, &c_b11, &fden, &tcos[*k4]);
    *k3 = lr;
    *k4 = 0;
    tri3_(&mr, &a[1], &bb[1], &c__[1], k, &b[1], &b2[1], &b3[1], &tcos[1], &
	    d__[1], &w[1], &w2[1], &w3[1]);
    i__2 = mr;
    for (i__ = 1; i__ <= i__2; ++i__) {
	b[i__] += b2[i__];
/* L190: */
    }
    tcos[1] = 2.f;
    trix_(&c__1, &c__0, &mr, &a[1], &bb[1], &c__[1], &b[1], &tcos[1], &d__[1],
	     &w[1]);
    i__2 = mr;
    for (i__ = 1; i__ <= i__2; ++i__) {
	q[i__ + q_dim1] += b[i__];
/* L191: */
    }
    goto L194;
L192:
    i__2 = mr;
    for (i__ = 1; i__ <= i__2; ++i__) {
	b[i__] = q[i__ + nlast * q_dim1];
/* L193: */
    }
    goto L196;
L194:

/*     START BACK SUBSTITUTION. */

    j = nlast - jr;
    i__2 = mr;
    for (i__ = 1; i__ <= i__2; ++i__) {
	b[i__] = q[i__ + nlast * q_dim1] + q[i__ + j * q_dim1];
/* L195: */
    }
L196:
    jm2 = nlast - i2r;
    if (jr != 1) {
	goto L198;
    }
    i__2 = mr;
    for (i__ = 1; i__ <= i__2; ++i__) {
	q[i__ + nlast * q_dim1] = 0.f;
/* L197: */
    }
    goto L202;
L198:
    if (nrod != 0) {
	goto L200;
    }
    i__2 = mr;
    for (i__ = 1; i__ <= i__2; ++i__) {
	ii = ip + i__;
	q[i__ + nlast * q_dim1] = p[ii];
/* L199: */
    }
    ip -= mr;
    goto L202;
L200:
    i__2 = mr;
    for (i__ = 1; i__ <= i__2; ++i__) {
	q[i__ + nlast * q_dim1] -= q[i__ + jm2 * q_dim1];
/* L201: */
    }
L202:
    cosgen_(&kr, &c__1, &c_b11, &fden, &tcos[1]);
    cosgen_(&lr, &c__1, &c_b11, &fden, &tcos[kr + 1]);
    if (lr != 0) {
	goto L204;
    }
    i__2 = mr;
    for (i__ = 1; i__ <= i__2; ++i__) {
	b[i__] = fistag * b[i__];
/* L203: */
    }
L204:
    trix_(&kr, &lr, &mr, &a[1], &bb[1], &c__[1], &b[1], &tcos[1], &d__[1], &w[
	    1]);
    i__2 = mr;
    for (i__ = 1; i__ <= i__2; ++i__) {
	q[i__ + nlast * q_dim1] += b[i__];
/* L205: */
    }
    nlastp = nlast;
L206:
    jstep = jr;
    jr = i2r;
    i2r /= 2;
    if (jr == 0) {
	goto L222;
    }
    switch (*mixbnd) {
	case 1:  goto L207;
	case 2:  goto L208;
    }
L207:
    jstart = jr + 1;
    goto L209;
L208:
    jstart = jr;
L209:
    kr -= jr;
    if (nlast + jr > *n) {
	goto L210;
    }
    kr -= jr;
    nlast += jr;
    jstop = nlast - jstep;
    goto L211;
L210:
    jstop = nlast - jr;
L211:
    lr = kr - jr;
    cosgen_(&jr, &c__1, &c_b11, &c_b12, &tcos[1]);
    i__2 = jstop;
    i__1 = jstep;
    for (j = jstart; i__1 < 0 ? j >= i__2 : j <= i__2; j += i__1) {
	jm2 = j - jr;
	jp2 = j + jr;
	if (j != jr) {
	    goto L213;
	}
	i__3 = mr;
	for (i__ = 1; i__ <= i__3; ++i__) {
	    b[i__] = q[i__ + j * q_dim1] + q[i__ + jp2 * q_dim1];
/* L212: */
	}
	goto L215;
L213:
	i__3 = mr;
	for (i__ = 1; i__ <= i__3; ++i__) {
	    b[i__] = q[i__ + j * q_dim1] + q[i__ + jm2 * q_dim1] + q[i__ + 
		    jp2 * q_dim1];
/* L214: */
	}
L215:
	if (jr != 1) {
	    goto L217;
	}
	i__3 = mr;
	for (i__ = 1; i__ <= i__3; ++i__) {
	    q[i__ + j * q_dim1] = 0.f;
/* L216: */
	}
	goto L219;
L217:
	jm1 = j - i2r;
	jp1 = j + i2r;
	i__3 = mr;
	for (i__ = 1; i__ <= i__3; ++i__) {
	    q[i__ + j * q_dim1] = (q[i__ + j * q_dim1] - q[i__ + jm1 * q_dim1]
		     - q[i__ + jp1 * q_dim1]) * .5f;
/* L218: */
	}
L219:
	trix_(&jr, &c__0, &mr, &a[1], &bb[1], &c__[1], &b[1], &tcos[1], &d__[
		1], &w[1]);
	i__3 = mr;
	for (i__ = 1; i__ <= i__3; ++i__) {
	    q[i__ + j * q_dim1] += b[i__];
/* L220: */
	}
/* L221: */
    }
    nrod = 1;
    if (nlast + i2r <= *n) {
	nrod = 0;
    }
    if (nlastp != nlast) {
	goto L194;
    }
    goto L206;
L222:

/*     RETURN STORAGE REQUIREMENTS FOR P VECTORS. */

    w[1] = (real) ipstor;
    return 0;
} /* poisn2_ */

#undef k4
#undef k3
#undef k2
#undef k1
#undef k


