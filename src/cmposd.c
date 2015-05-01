/* cmposd.f -- translated by f2c (version 12.02.01).
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
static integer c__0 = 0;
static real c_b20 = .5f;
static real c_b21 = 0.f;

/* DECK CMPOSD */
/* Subroutine */ int cmposd_(integer *mr, integer *nr, integer *istag, 
	complex *ba, complex *bb, complex *bc, complex *q, integer *idimq, 
	complex *b, complex *w, complex *d__, complex *tcos, complex *p)
{
    /* System generated locals */
    integer q_dim1, q_offset, i__1, i__2, i__3, i__4, i__5, i__6, i__7, i__8;
    real r__1;
    complex q__1, q__2, q__3, q__4, q__5, q__6;

    /* Local variables */
    static integer i__, j, l, m, n;
    static complex t;
    static real fi;
    static integer ip, kr, lr, jm1, jm2, jm3, jp1, jp2, jp3, ip1, jsh, jsp, 
	    nun, jst, ideg, jdeg, nodd, krpi, irreg;
    extern /* Subroutine */ int c1merg_(complex *, integer *, integer *, 
	    integer *, integer *, integer *), cmpcsg_(integer *, integer *, 
	    real *, real *, complex *);
    static integer noddpr, jstsav;
    extern /* Subroutine */ int cmptrx_(integer *, integer *, integer *, 
	    complex *, complex *, complex *, complex *, complex *, complex *, 
	    complex *);
    static integer ipstor;

/* ***BEGIN PROLOGUE  CMPOSD */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to CMGNBN */
/* ***LIBRARY   SLATEC */
/* ***TYPE      COMPLEX (POISD2-S, CMPOSD-C) */
/* ***AUTHOR  (UNKNOWN) */
/* ***DESCRIPTION */

/*     Subroutine to solve Poisson's equation for Dirichlet boundary */
/*     conditions. */

/*     ISTAG = 1 if the last diagonal block is the matrix A. */
/*     ISTAG = 2 if the last diagonal block is the matrix A+I. */

/* ***SEE ALSO  CMGNBN */
/* ***ROUTINES CALLED  C1MERG, CMPCSG, CMPTRX */
/* ***REVISION HISTORY  (YYMMDD) */
/*   801001  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900402  Added TYPE section.  (WRB) */
/*   920130  Modified to use merge routine C1MERG rather than deleted */
/*           routine CMPMRG.  (WRB) */
/* ***END PROLOGUE  CMPOSD */

/* ***FIRST EXECUTABLE STATEMENT  CMPOSD */
    /* Parameter adjustments */
    --ba;
    --bb;
    --bc;
    q_dim1 = *idimq;
    q_offset = 1 + q_dim1;
    q -= q_offset;
    --b;
    --w;
    --d__;
    --tcos;
    --p;

    /* Function Body */
    m = *mr;
    n = *nr;
    fi = 1.f / *istag;
    ip = -m;
    ipstor = 0;
    jsh = 0;
    switch (*istag) {
	case 1:  goto L101;
	case 2:  goto L102;
    }
L101:
    kr = 0;
    irreg = 1;
    if (n > 1) {
	goto L106;
    }
    tcos[1].r = 0.f, tcos[1].i = 0.f;
    goto L103;
L102:
    kr = 1;
    jstsav = 1;
    irreg = 2;
    if (n > 1) {
	goto L106;
    }
    tcos[1].r = -1.f, tcos[1].i = 0.f;
L103:
    i__1 = m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = i__;
	i__3 = i__ + q_dim1;
	b[i__2].r = q[i__3].r, b[i__2].i = q[i__3].i;
/* L104: */
    }
    cmptrx_(&c__1, &c__0, &m, &ba[1], &bb[1], &bc[1], &b[1], &tcos[1], &d__[1]
	    , &w[1]);
    i__1 = m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = i__ + q_dim1;
	i__3 = i__;
	q[i__2].r = b[i__3].r, q[i__2].i = b[i__3].i;
/* L105: */
    }
    goto L183;
L106:
    lr = 0;
    i__1 = m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = i__;
	p[i__2].r = 0.f, p[i__2].i = 0.f;
/* L107: */
    }
    nun = n;
    jst = 1;
    jsp = n;

/*     IRREG = 1 WHEN NO IRREGULARITIES HAVE OCCURRED, OTHERWISE IT IS 2. */

L108:
    l = jst << 1;
    nodd = 2 - ((nun + 1) / 2 << 1) + nun;

/*     NODD = 1 WHEN NUN IS ODD, OTHERWISE IT IS 2. */

    switch (nodd) {
	case 1:  goto L110;
	case 2:  goto L109;
    }
L109:
    jsp -= l;
    goto L111;
L110:
    jsp -= jst;
    if (irreg != 1) {
	jsp -= l;
    }
L111:

/*     REGULAR REDUCTION */

    cmpcsg_(&jst, &c__1, &c_b20, &c_b21, &tcos[1]);
    if (l > jsp) {
	goto L118;
    }
    i__1 = jsp;
    i__2 = l;
    for (j = l; i__2 < 0 ? j >= i__1 : j <= i__1; j += i__2) {
	jm1 = j - jsh;
	jp1 = j + jsh;
	jm2 = j - jst;
	jp2 = j + jst;
	jm3 = jm2 - jsh;
	jp3 = jp2 + jsh;
	if (jst != 1) {
	    goto L113;
	}
	i__3 = m;
	for (i__ = 1; i__ <= i__3; ++i__) {
	    i__4 = i__;
	    i__5 = i__ + j * q_dim1;
	    q__1.r = q[i__5].r * 2.f, q__1.i = q[i__5].i * 2.f;
	    b[i__4].r = q__1.r, b[i__4].i = q__1.i;
	    i__4 = i__ + j * q_dim1;
	    i__5 = i__ + jm2 * q_dim1;
	    i__6 = i__ + jp2 * q_dim1;
	    q__1.r = q[i__5].r + q[i__6].r, q__1.i = q[i__5].i + q[i__6].i;
	    q[i__4].r = q__1.r, q[i__4].i = q__1.i;
/* L112: */
	}
	goto L115;
L113:
	i__3 = m;
	for (i__ = 1; i__ <= i__3; ++i__) {
	    i__4 = i__ + j * q_dim1;
	    i__5 = i__ + jm1 * q_dim1;
	    q__4.r = q[i__4].r - q[i__5].r, q__4.i = q[i__4].i - q[i__5].i;
	    i__6 = i__ + jp1 * q_dim1;
	    q__3.r = q__4.r - q[i__6].r, q__3.i = q__4.i - q[i__6].i;
	    i__7 = i__ + jm2 * q_dim1;
	    q__2.r = q__3.r + q[i__7].r, q__2.i = q__3.i + q[i__7].i;
	    i__8 = i__ + jp2 * q_dim1;
	    q__1.r = q__2.r + q[i__8].r, q__1.i = q__2.i + q[i__8].i;
	    t.r = q__1.r, t.i = q__1.i;
	    i__4 = i__;
	    i__5 = i__ + j * q_dim1;
	    q__3.r = t.r + q[i__5].r, q__3.i = t.i + q[i__5].i;
	    i__6 = i__ + jm3 * q_dim1;
	    q__2.r = q__3.r - q[i__6].r, q__2.i = q__3.i - q[i__6].i;
	    i__7 = i__ + jp3 * q_dim1;
	    q__1.r = q__2.r - q[i__7].r, q__1.i = q__2.i - q[i__7].i;
	    b[i__4].r = q__1.r, b[i__4].i = q__1.i;
	    i__4 = i__ + j * q_dim1;
	    q[i__4].r = t.r, q[i__4].i = t.i;
/* L114: */
	}
L115:
	cmptrx_(&jst, &c__0, &m, &ba[1], &bb[1], &bc[1], &b[1], &tcos[1], &
		d__[1], &w[1]);
	i__3 = m;
	for (i__ = 1; i__ <= i__3; ++i__) {
	    i__4 = i__ + j * q_dim1;
	    i__5 = i__ + j * q_dim1;
	    i__6 = i__;
	    q__1.r = q[i__5].r + b[i__6].r, q__1.i = q[i__5].i + b[i__6].i;
	    q[i__4].r = q__1.r, q[i__4].i = q__1.i;
/* L116: */
	}
/* L117: */
    }

/*     REDUCTION FOR LAST UNKNOWN */

L118:
    switch (nodd) {
	case 1:  goto L119;
	case 2:  goto L136;
    }
L119:
    switch (irreg) {
	case 1:  goto L152;
	case 2:  goto L120;
    }

/*     ODD NUMBER OF UNKNOWNS */

L120:
    jsp += l;
    j = jsp;
    jm1 = j - jsh;
    jp1 = j + jsh;
    jm2 = j - jst;
    jp2 = j + jst;
    jm3 = jm2 - jsh;
    switch (*istag) {
	case 1:  goto L123;
	case 2:  goto L121;
    }
L121:
    if (jst != 1) {
	goto L123;
    }
    i__2 = m;
    for (i__ = 1; i__ <= i__2; ++i__) {
	i__1 = i__;
	i__3 = i__ + j * q_dim1;
	b[i__1].r = q[i__3].r, b[i__1].i = q[i__3].i;
	i__1 = i__ + j * q_dim1;
	q[i__1].r = 0.f, q[i__1].i = 0.f;
/* L122: */
    }
    goto L130;
L123:
    switch (noddpr) {
	case 1:  goto L124;
	case 2:  goto L126;
    }
L124:
    i__2 = m;
    for (i__ = 1; i__ <= i__2; ++i__) {
	ip1 = ip + i__;
	i__1 = i__;
	i__3 = i__ + jm2 * q_dim1;
	i__4 = i__ + jm1 * q_dim1;
	q__5.r = q[i__3].r - q[i__4].r, q__5.i = q[i__3].i - q[i__4].i;
	i__5 = i__ + jm3 * q_dim1;
	q__4.r = q__5.r - q[i__5].r, q__4.i = q__5.i - q[i__5].i;
	q__3.r = q__4.r * .5f, q__3.i = q__4.i * .5f;
	i__6 = ip1;
	q__2.r = q__3.r + p[i__6].r, q__2.i = q__3.i + p[i__6].i;
	i__7 = i__ + j * q_dim1;
	q__1.r = q__2.r + q[i__7].r, q__1.i = q__2.i + q[i__7].i;
	b[i__1].r = q__1.r, b[i__1].i = q__1.i;
/* L125: */
    }
    goto L128;
L126:
    i__2 = m;
    for (i__ = 1; i__ <= i__2; ++i__) {
	i__1 = i__;
	i__3 = i__ + jm2 * q_dim1;
	i__4 = i__ + jm1 * q_dim1;
	q__6.r = q[i__3].r - q[i__4].r, q__6.i = q[i__3].i - q[i__4].i;
	i__5 = i__ + jm3 * q_dim1;
	q__5.r = q__6.r - q[i__5].r, q__5.i = q__6.i - q[i__5].i;
	q__4.r = q__5.r * .5f, q__4.i = q__5.i * .5f;
	i__6 = i__ + jp2 * q_dim1;
	q__3.r = q__4.r + q[i__6].r, q__3.i = q__4.i + q[i__6].i;
	i__7 = i__ + jp1 * q_dim1;
	q__2.r = q__3.r - q[i__7].r, q__2.i = q__3.i - q[i__7].i;
	i__8 = i__ + j * q_dim1;
	q__1.r = q__2.r + q[i__8].r, q__1.i = q__2.i + q[i__8].i;
	b[i__1].r = q__1.r, b[i__1].i = q__1.i;
/* L127: */
    }
L128:
    i__2 = m;
    for (i__ = 1; i__ <= i__2; ++i__) {
	i__1 = i__ + j * q_dim1;
	i__3 = i__ + j * q_dim1;
	i__4 = i__ + jm1 * q_dim1;
	q__3.r = q[i__3].r - q[i__4].r, q__3.i = q[i__3].i - q[i__4].i;
	i__5 = i__ + jp1 * q_dim1;
	q__2.r = q__3.r - q[i__5].r, q__2.i = q__3.i - q[i__5].i;
	q__1.r = q__2.r * .5f, q__1.i = q__2.i * .5f;
	q[i__1].r = q__1.r, q[i__1].i = q__1.i;
/* L129: */
    }
L130:
    cmptrx_(&jst, &c__0, &m, &ba[1], &bb[1], &bc[1], &b[1], &tcos[1], &d__[1],
	     &w[1]);
    ip += m;
/* Computing MAX */
    i__2 = ipstor, i__1 = ip + m;
    ipstor = max(i__2,i__1);
    i__2 = m;
    for (i__ = 1; i__ <= i__2; ++i__) {
	ip1 = ip + i__;
	i__1 = ip1;
	i__3 = i__ + j * q_dim1;
	i__4 = i__;
	q__1.r = q[i__3].r + b[i__4].r, q__1.i = q[i__3].i + b[i__4].i;
	p[i__1].r = q__1.r, p[i__1].i = q__1.i;
	i__1 = i__;
	i__3 = i__ + jp2 * q_dim1;
	i__4 = ip1;
	q__1.r = q[i__3].r + p[i__4].r, q__1.i = q[i__3].i + p[i__4].i;
	b[i__1].r = q__1.r, b[i__1].i = q__1.i;
/* L131: */
    }
    if (lr != 0) {
	goto L133;
    }
    i__2 = jst;
    for (i__ = 1; i__ <= i__2; ++i__) {
	krpi = kr + i__;
	i__1 = krpi;
	i__3 = i__;
	tcos[i__1].r = tcos[i__3].r, tcos[i__1].i = tcos[i__3].i;
/* L132: */
    }
    goto L134;
L133:
    cmpcsg_(&lr, &jstsav, &c_b21, &fi, &tcos[jst + 1]);
    c1merg_(&tcos[1], &c__0, &jst, &jst, &lr, &kr);
L134:
    cmpcsg_(&kr, &jstsav, &c_b21, &fi, &tcos[1]);
    cmptrx_(&kr, &kr, &m, &ba[1], &bb[1], &bc[1], &b[1], &tcos[1], &d__[1], &
	    w[1]);
    i__2 = m;
    for (i__ = 1; i__ <= i__2; ++i__) {
	ip1 = ip + i__;
	i__1 = i__ + j * q_dim1;
	i__3 = i__ + jm2 * q_dim1;
	i__4 = i__;
	q__2.r = q[i__3].r + b[i__4].r, q__2.i = q[i__3].i + b[i__4].i;
	i__5 = ip1;
	q__1.r = q__2.r + p[i__5].r, q__1.i = q__2.i + p[i__5].i;
	q[i__1].r = q__1.r, q[i__1].i = q__1.i;
/* L135: */
    }
    lr = kr;
    kr += l;
    goto L152;

/*     EVEN NUMBER OF UNKNOWNS */

L136:
    jsp += l;
    j = jsp;
    jm1 = j - jsh;
    jp1 = j + jsh;
    jm2 = j - jst;
    jp2 = j + jst;
    jm3 = jm2 - jsh;
    switch (irreg) {
	case 1:  goto L137;
	case 2:  goto L138;
    }
L137:
    jstsav = jst;
    ideg = jst;
    kr = l;
    goto L139;
L138:
    cmpcsg_(&kr, &jstsav, &c_b21, &fi, &tcos[1]);
    cmpcsg_(&lr, &jstsav, &c_b21, &fi, &tcos[kr + 1]);
    ideg = kr;
    kr += jst;
L139:
    if (jst != 1) {
	goto L141;
    }
    irreg = 2;
    i__2 = m;
    for (i__ = 1; i__ <= i__2; ++i__) {
	i__1 = i__;
	i__3 = i__ + j * q_dim1;
	b[i__1].r = q[i__3].r, b[i__1].i = q[i__3].i;
	i__1 = i__ + j * q_dim1;
	i__3 = i__ + jm2 * q_dim1;
	q[i__1].r = q[i__3].r, q[i__1].i = q[i__3].i;
/* L140: */
    }
    goto L150;
L141:
    i__2 = m;
    for (i__ = 1; i__ <= i__2; ++i__) {
	i__1 = i__;
	i__3 = i__ + j * q_dim1;
	i__4 = i__ + jm2 * q_dim1;
	i__5 = i__ + jm1 * q_dim1;
	q__4.r = q[i__4].r - q[i__5].r, q__4.i = q[i__4].i - q[i__5].i;
	i__6 = i__ + jm3 * q_dim1;
	q__3.r = q__4.r - q[i__6].r, q__3.i = q__4.i - q[i__6].i;
	q__2.r = q__3.r * .5f, q__2.i = q__3.i * .5f;
	q__1.r = q[i__3].r + q__2.r, q__1.i = q[i__3].i + q__2.i;
	b[i__1].r = q__1.r, b[i__1].i = q__1.i;
/* L142: */
    }
    switch (irreg) {
	case 1:  goto L143;
	case 2:  goto L145;
    }
L143:
    i__2 = m;
    for (i__ = 1; i__ <= i__2; ++i__) {
	i__1 = i__ + j * q_dim1;
	i__3 = i__ + jm2 * q_dim1;
	i__4 = i__ + j * q_dim1;
	i__5 = i__ + jm1 * q_dim1;
	q__4.r = q[i__4].r - q[i__5].r, q__4.i = q[i__4].i - q[i__5].i;
	i__6 = i__ + jp1 * q_dim1;
	q__3.r = q__4.r - q[i__6].r, q__3.i = q__4.i - q[i__6].i;
	q__2.r = q__3.r * .5f, q__2.i = q__3.i * .5f;
	q__1.r = q[i__3].r + q__2.r, q__1.i = q[i__3].i + q__2.i;
	q[i__1].r = q__1.r, q[i__1].i = q__1.i;
/* L144: */
    }
    irreg = 2;
    goto L150;
L145:
    switch (noddpr) {
	case 1:  goto L146;
	case 2:  goto L148;
    }
L146:
    i__2 = m;
    for (i__ = 1; i__ <= i__2; ++i__) {
	ip1 = ip + i__;
	i__1 = i__ + j * q_dim1;
	i__3 = i__ + jm2 * q_dim1;
	i__4 = ip1;
	q__1.r = q[i__3].r + p[i__4].r, q__1.i = q[i__3].i + p[i__4].i;
	q[i__1].r = q__1.r, q[i__1].i = q__1.i;
/* L147: */
    }
    ip -= m;
    goto L150;
L148:
    i__2 = m;
    for (i__ = 1; i__ <= i__2; ++i__) {
	i__1 = i__ + j * q_dim1;
	i__3 = i__ + jm2 * q_dim1;
	i__4 = i__ + j * q_dim1;
	q__2.r = q[i__3].r + q[i__4].r, q__2.i = q[i__3].i + q[i__4].i;
	i__5 = i__ + jm1 * q_dim1;
	q__1.r = q__2.r - q[i__5].r, q__1.i = q__2.i - q[i__5].i;
	q[i__1].r = q__1.r, q[i__1].i = q__1.i;
/* L149: */
    }
L150:
    cmptrx_(&ideg, &lr, &m, &ba[1], &bb[1], &bc[1], &b[1], &tcos[1], &d__[1], 
	    &w[1]);
    i__2 = m;
    for (i__ = 1; i__ <= i__2; ++i__) {
	i__1 = i__ + j * q_dim1;
	i__3 = i__ + j * q_dim1;
	i__4 = i__;
	q__1.r = q[i__3].r + b[i__4].r, q__1.i = q[i__3].i + b[i__4].i;
	q[i__1].r = q__1.r, q[i__1].i = q__1.i;
/* L151: */
    }
L152:
    nun /= 2;
    noddpr = nodd;
    jsh = jst;
    jst <<= 1;
    if (nun >= 2) {
	goto L108;
    }

/*     START SOLUTION. */

    j = jsp;
    i__2 = m;
    for (i__ = 1; i__ <= i__2; ++i__) {
	i__1 = i__;
	i__3 = i__ + j * q_dim1;
	b[i__1].r = q[i__3].r, b[i__1].i = q[i__3].i;
/* L153: */
    }
    switch (irreg) {
	case 1:  goto L154;
	case 2:  goto L155;
    }
L154:
    cmpcsg_(&jst, &c__1, &c_b20, &c_b21, &tcos[1]);
    ideg = jst;
    goto L156;
L155:
    kr = lr + jst;
    cmpcsg_(&kr, &jstsav, &c_b21, &fi, &tcos[1]);
    cmpcsg_(&lr, &jstsav, &c_b21, &fi, &tcos[kr + 1]);
    ideg = kr;
L156:
    cmptrx_(&ideg, &lr, &m, &ba[1], &bb[1], &bc[1], &b[1], &tcos[1], &d__[1], 
	    &w[1]);
    jm1 = j - jsh;
    jp1 = j + jsh;
    switch (irreg) {
	case 1:  goto L157;
	case 2:  goto L159;
    }
L157:
    i__2 = m;
    for (i__ = 1; i__ <= i__2; ++i__) {
	i__1 = i__ + j * q_dim1;
	i__3 = i__ + j * q_dim1;
	i__4 = i__ + jm1 * q_dim1;
	q__4.r = q[i__3].r - q[i__4].r, q__4.i = q[i__3].i - q[i__4].i;
	i__5 = i__ + jp1 * q_dim1;
	q__3.r = q__4.r - q[i__5].r, q__3.i = q__4.i - q[i__5].i;
	q__2.r = q__3.r * .5f, q__2.i = q__3.i * .5f;
	i__6 = i__;
	q__1.r = q__2.r + b[i__6].r, q__1.i = q__2.i + b[i__6].i;
	q[i__1].r = q__1.r, q[i__1].i = q__1.i;
/* L158: */
    }
    goto L164;
L159:
    switch (noddpr) {
	case 1:  goto L160;
	case 2:  goto L162;
    }
L160:
    i__2 = m;
    for (i__ = 1; i__ <= i__2; ++i__) {
	ip1 = ip + i__;
	i__1 = i__ + j * q_dim1;
	i__3 = ip1;
	i__4 = i__;
	q__1.r = p[i__3].r + b[i__4].r, q__1.i = p[i__3].i + b[i__4].i;
	q[i__1].r = q__1.r, q[i__1].i = q__1.i;
/* L161: */
    }
    ip -= m;
    goto L164;
L162:
    i__2 = m;
    for (i__ = 1; i__ <= i__2; ++i__) {
	i__1 = i__ + j * q_dim1;
	i__3 = i__ + j * q_dim1;
	i__4 = i__ + jm1 * q_dim1;
	q__2.r = q[i__3].r - q[i__4].r, q__2.i = q[i__3].i - q[i__4].i;
	i__5 = i__;
	q__1.r = q__2.r + b[i__5].r, q__1.i = q__2.i + b[i__5].i;
	q[i__1].r = q__1.r, q[i__1].i = q__1.i;
/* L163: */
    }
L164:

/*     START BACK SUBSTITUTION. */

    jst /= 2;
    jsh = jst / 2;
    nun <<= 1;
    if (nun > n) {
	goto L183;
    }
    i__2 = n;
    i__1 = l;
    for (j = jst; i__1 < 0 ? j >= i__2 : j <= i__2; j += i__1) {
	jm1 = j - jsh;
	jp1 = j + jsh;
	jm2 = j - jst;
	jp2 = j + jst;
	if (j > jst) {
	    goto L166;
	}
	i__3 = m;
	for (i__ = 1; i__ <= i__3; ++i__) {
	    i__4 = i__;
	    i__5 = i__ + j * q_dim1;
	    i__6 = i__ + jp2 * q_dim1;
	    q__1.r = q[i__5].r + q[i__6].r, q__1.i = q[i__5].i + q[i__6].i;
	    b[i__4].r = q__1.r, b[i__4].i = q__1.i;
/* L165: */
	}
	goto L170;
L166:
	if (jp2 <= n) {
	    goto L168;
	}
	i__3 = m;
	for (i__ = 1; i__ <= i__3; ++i__) {
	    i__4 = i__;
	    i__5 = i__ + j * q_dim1;
	    i__6 = i__ + jm2 * q_dim1;
	    q__1.r = q[i__5].r + q[i__6].r, q__1.i = q[i__5].i + q[i__6].i;
	    b[i__4].r = q__1.r, b[i__4].i = q__1.i;
/* L167: */
	}
	if (jst < jstsav) {
	    irreg = 1;
	}
	switch (irreg) {
	    case 1:  goto L170;
	    case 2:  goto L171;
	}
L168:
	i__3 = m;
	for (i__ = 1; i__ <= i__3; ++i__) {
	    i__4 = i__;
	    i__5 = i__ + j * q_dim1;
	    i__6 = i__ + jm2 * q_dim1;
	    q__2.r = q[i__5].r + q[i__6].r, q__2.i = q[i__5].i + q[i__6].i;
	    i__7 = i__ + jp2 * q_dim1;
	    q__1.r = q__2.r + q[i__7].r, q__1.i = q__2.i + q[i__7].i;
	    b[i__4].r = q__1.r, b[i__4].i = q__1.i;
/* L169: */
	}
L170:
	cmpcsg_(&jst, &c__1, &c_b20, &c_b21, &tcos[1]);
	ideg = jst;
	jdeg = 0;
	goto L172;
L171:
	if (j + l > n) {
	    lr -= jst;
	}
	kr = jst + lr;
	cmpcsg_(&kr, &jstsav, &c_b21, &fi, &tcos[1]);
	cmpcsg_(&lr, &jstsav, &c_b21, &fi, &tcos[kr + 1]);
	ideg = kr;
	jdeg = lr;
L172:
	cmptrx_(&ideg, &jdeg, &m, &ba[1], &bb[1], &bc[1], &b[1], &tcos[1], &
		d__[1], &w[1]);
	if (jst > 1) {
	    goto L174;
	}
	i__3 = m;
	for (i__ = 1; i__ <= i__3; ++i__) {
	    i__4 = i__ + j * q_dim1;
	    i__5 = i__;
	    q[i__4].r = b[i__5].r, q[i__4].i = b[i__5].i;
/* L173: */
	}
	goto L182;
L174:
	if (jp2 > n) {
	    goto L177;
	}
L175:
	i__3 = m;
	for (i__ = 1; i__ <= i__3; ++i__) {
	    i__4 = i__ + j * q_dim1;
	    i__5 = i__ + j * q_dim1;
	    i__6 = i__ + jm1 * q_dim1;
	    q__4.r = q[i__5].r - q[i__6].r, q__4.i = q[i__5].i - q[i__6].i;
	    i__7 = i__ + jp1 * q_dim1;
	    q__3.r = q__4.r - q[i__7].r, q__3.i = q__4.i - q[i__7].i;
	    q__2.r = q__3.r * .5f, q__2.i = q__3.i * .5f;
	    i__8 = i__;
	    q__1.r = q__2.r + b[i__8].r, q__1.i = q__2.i + b[i__8].i;
	    q[i__4].r = q__1.r, q[i__4].i = q__1.i;
/* L176: */
	}
	goto L182;
L177:
	switch (irreg) {
	    case 1:  goto L175;
	    case 2:  goto L178;
	}
L178:
	if (j + jsh > n) {
	    goto L180;
	}
	i__3 = m;
	for (i__ = 1; i__ <= i__3; ++i__) {
	    ip1 = ip + i__;
	    i__4 = i__ + j * q_dim1;
	    i__5 = i__;
	    i__6 = ip1;
	    q__1.r = b[i__5].r + p[i__6].r, q__1.i = b[i__5].i + p[i__6].i;
	    q[i__4].r = q__1.r, q[i__4].i = q__1.i;
/* L179: */
	}
	ip -= m;
	goto L182;
L180:
	i__3 = m;
	for (i__ = 1; i__ <= i__3; ++i__) {
	    i__4 = i__ + j * q_dim1;
	    i__5 = i__;
	    i__6 = i__ + j * q_dim1;
	    q__2.r = b[i__5].r + q[i__6].r, q__2.i = b[i__5].i + q[i__6].i;
	    i__7 = i__ + jm1 * q_dim1;
	    q__1.r = q__2.r - q[i__7].r, q__1.i = q__2.i - q[i__7].i;
	    q[i__4].r = q__1.r, q[i__4].i = q__1.i;
/* L181: */
	}
L182:
	;
    }
    l /= 2;
    goto L164;
L183:

/*     RETURN STORAGE REQUIREMENTS FOR P VECTORS. */

    r__1 = (real) ipstor;
    q__1.r = r__1, q__1.i = 0.f;
    w[1].r = q__1.r, w[1].i = q__1.i;
    return 0;
} /* cmposd_ */

