/* postg2.f -- translated by f2c (version 12.02.01).
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
static real c_b9 = .5f;
static real c_b18 = 0.f;
static integer c__0 = 0;
static integer c__3 = 3;
static real c_b98 = 1.f;

/* DECK POSTG2 */
/* Subroutine */ int postg2_(integer *nperod, integer *n, integer *m, real *a,
	 real *bb, real *c__, integer *idimq, real *q, real *b, real *b2, 
	real *b3, real *w, real *w2, real *w3, real *d__, real *tcos, real *p)
{
    /* System generated locals */
    integer q_dim1, q_offset, i__1, i__2, i__3;
    real r__1, r__2;
    static integer equiv_3[4];

    /* Local variables */
    static integer i__, j;
#define k (equiv_3)
    static real t;
#define k1 (equiv_3)
#define k2 (equiv_3 + 1)
#define k3 (equiv_3 + 2)
#define k4 (equiv_3 + 3)
    static real fi;
    static integer ii, ip, jr, kr, np, mr, nr, lr, jm1, jm2, jm3, jp1, jp2, 
	    i2r, jp3;
    extern /* Subroutine */ int tri3_(integer *, real *, real *, real *, 
	    integer *, real *, real *, real *, real *, real *, real *, real *,
	     real *);
    static integer nrod;
    static real fnum;
    extern /* Subroutine */ int trix_(integer *, integer *, integer *, real *,
	     real *, real *, real *, real *, real *, real *);
    static real fnum2;
    static integer i2rby2, nlast, ijump, jstep, jstop;
    extern /* Subroutine */ int s1merg_(real *, integer *, integer *, integer 
	    *, integer *, integer *), cosgen_(integer *, integer *, real *, 
	    real *, real *);
    static integer nlastp, nrodpr, jstart, ipstor;

/* ***BEGIN PROLOGUE  POSTG2 */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to POISTG */
/* ***LIBRARY   SLATEC */
/* ***TYPE      SINGLE PRECISION (POSTG2-S) */
/* ***AUTHOR  (UNKNOWN) */
/* ***DESCRIPTION */

/*     Subroutine to solve Poisson's equation on a staggered grid. */

/* ***SEE ALSO  POISTG */
/* ***ROUTINES CALLED  COSGEN, S1MERG, TRI3, TRIX */
/* ***REVISION HISTORY  (YYMMDD) */
/*   801001  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900402  Added TYPE section.  (WRB) */
/*   920130  Modified to use merge routine S1MERG rather than deleted */
/*           routine MERGE.  (WRB) */
/* ***END PROLOGUE  POSTG2 */

/* ***FIRST EXECUTABLE STATEMENT  POSTG2 */
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
    np = *nperod;
    fnum = np / 3 * .5f;
    fnum2 = np / 2 * .5f;
    mr = *m;
    ip = -mr;
    ipstor = 0;
    i2r = 1;
    jr = 2;
    nr = *n;
    nlast = *n;
    kr = 1;
    lr = 0;
    if (nr <= 3) {
	goto L142;
    }
L101:
    jr = i2r << 1;
    nrod = 1;
    if (nr / 2 << 1 == nr) {
	nrod = 0;
    }
    jstart = 1;
    jstop = nlast - jr;
    if (nrod == 0) {
	jstop -= i2r;
    }
    i2rby2 = i2r / 2;
    if (jstop >= jstart) {
	goto L102;
    }
    j = jr;
    goto L115;
L102:

/*     REGULAR REDUCTION. */

    ijump = 1;
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
	    goto L106;
	}
	cosgen_(&i2r, &c__1, &fnum, &c_b9, &tcos[1]);
	if (i2r != 1) {
	    goto L104;
	}
	i__3 = mr;
	for (i__ = 1; i__ <= i__3; ++i__) {
	    b[i__] = q[i__ + q_dim1];
	    q[i__ + q_dim1] = q[i__ + (q_dim1 << 1)];
/* L103: */
	}
	goto L112;
L104:
	i__3 = mr;
	for (i__ = 1; i__ <= i__3; ++i__) {
	    b[i__] = q[i__ + q_dim1] + (q[i__ + jp2 * q_dim1] - q[i__ + jp1 * 
		    q_dim1] - q[i__ + jp3 * q_dim1]) * .5f;
	    q[i__ + q_dim1] = q[i__ + jp2 * q_dim1] + q[i__ + q_dim1] - q[i__ 
		    + jp1 * q_dim1];
/* L105: */
	}
	goto L112;
L106:
	switch (ijump) {
	    case 1:  goto L107;
	    case 2:  goto L108;
	}
L107:
	ijump = 2;
	cosgen_(&i2r, &c__1, &c_b9, &c_b18, &tcos[1]);
L108:
	if (i2r != 1) {
	    goto L110;
	}
	i__3 = mr;
	for (i__ = 1; i__ <= i__3; ++i__) {
	    b[i__] = q[i__ + j * q_dim1] * 2.f;
	    q[i__ + j * q_dim1] = q[i__ + jm2 * q_dim1] + q[i__ + jp2 * 
		    q_dim1];
/* L109: */
	}
	goto L112;
L110:
	i__3 = mr;
	for (i__ = 1; i__ <= i__3; ++i__) {
	    fi = q[i__ + j * q_dim1];
	    q[i__ + j * q_dim1] = q[i__ + j * q_dim1] - q[i__ + jm1 * q_dim1] 
		    - q[i__ + jp1 * q_dim1] + q[i__ + jm2 * q_dim1] + q[i__ + 
		    jp2 * q_dim1];
	    b[i__] = fi + q[i__ + j * q_dim1] - q[i__ + jm3 * q_dim1] - q[i__ 
		    + jp3 * q_dim1];
/* L111: */
	}
L112:
	trix_(&i2r, &c__0, &mr, &a[1], &bb[1], &c__[1], &b[1], &tcos[1], &d__[
		1], &w[1]);
	i__3 = mr;
	for (i__ = 1; i__ <= i__3; ++i__) {
	    q[i__ + j * q_dim1] += b[i__];
/* L113: */
	}

/*     END OF REDUCTION FOR REGULAR UNKNOWNS. */

/* L114: */
    }

/*     BEGIN SPECIAL REDUCTION FOR LAST UNKNOWN. */

    j = jstop + jr;
L115:
    nlast = j;
    jm1 = j - i2rby2;
    jm2 = j - i2r;
    jm3 = jm2 - i2rby2;
    if (nrod == 0) {
	goto L125;
    }

/*     ODD NUMBER OF UNKNOWNS */

    if (i2r != 1) {
	goto L117;
    }
    i__2 = mr;
    for (i__ = 1; i__ <= i__2; ++i__) {
	b[i__] = q[i__ + j * q_dim1];
	q[i__ + j * q_dim1] = q[i__ + jm2 * q_dim1];
/* L116: */
    }
    goto L123;
L117:
    i__2 = mr;
    for (i__ = 1; i__ <= i__2; ++i__) {
	b[i__] = q[i__ + j * q_dim1] + (q[i__ + jm2 * q_dim1] - q[i__ + jm1 * 
		q_dim1] - q[i__ + jm3 * q_dim1]) * .5f;
/* L118: */
    }
    if (nrodpr != 0) {
	goto L120;
    }
    i__2 = mr;
    for (i__ = 1; i__ <= i__2; ++i__) {
	ii = ip + i__;
	q[i__ + j * q_dim1] = q[i__ + jm2 * q_dim1] + p[ii];
/* L119: */
    }
    ip -= mr;
    goto L122;
L120:
    i__2 = mr;
    for (i__ = 1; i__ <= i__2; ++i__) {
	q[i__ + j * q_dim1] = q[i__ + j * q_dim1] - q[i__ + jm1 * q_dim1] + q[
		i__ + jm2 * q_dim1];
/* L121: */
    }
L122:
    if (lr == 0) {
	goto L123;
    }
    cosgen_(&lr, &c__1, &fnum2, &c_b9, &tcos[kr + 1]);
L123:
    cosgen_(&kr, &c__1, &fnum2, &c_b9, &tcos[1]);
    trix_(&kr, &lr, &mr, &a[1], &bb[1], &c__[1], &b[1], &tcos[1], &d__[1], &w[
	    1]);
    i__2 = mr;
    for (i__ = 1; i__ <= i__2; ++i__) {
	q[i__ + j * q_dim1] += b[i__];
/* L124: */
    }
    kr += i2r;
    goto L141;
L125:

/*     EVEN NUMBER OF UNKNOWNS */

    jp1 = j + i2rby2;
    jp2 = j + i2r;
    if (i2r != 1) {
	goto L129;
    }
    i__2 = mr;
    for (i__ = 1; i__ <= i__2; ++i__) {
	b[i__] = q[i__ + j * q_dim1];
/* L126: */
    }
    tcos[1] = 0.f;
    trix_(&c__1, &c__0, &mr, &a[1], &bb[1], &c__[1], &b[1], &tcos[1], &d__[1],
	     &w[1]);
    ip = 0;
    ipstor = mr;
    i__2 = mr;
    for (i__ = 1; i__ <= i__2; ++i__) {
	p[i__] = b[i__];
	b[i__] += q[i__ + *n * q_dim1];
/* L127: */
    }
    tcos[1] = (np / 2 << 1) - 1.f;
    tcos[2] = 0.f;
    trix_(&c__1, &c__1, &mr, &a[1], &bb[1], &c__[1], &b[1], &tcos[1], &d__[1],
	     &w[1]);
    i__2 = mr;
    for (i__ = 1; i__ <= i__2; ++i__) {
	q[i__ + j * q_dim1] = q[i__ + jm2 * q_dim1] + p[i__] + b[i__];
/* L128: */
    }
    goto L140;
L129:
    i__2 = mr;
    for (i__ = 1; i__ <= i__2; ++i__) {
	b[i__] = q[i__ + j * q_dim1] + (q[i__ + jm2 * q_dim1] - q[i__ + jm1 * 
		q_dim1] - q[i__ + jm3 * q_dim1]) * .5f;
/* L130: */
    }
    if (nrodpr != 0) {
	goto L132;
    }
    i__2 = mr;
    for (i__ = 1; i__ <= i__2; ++i__) {
	ii = ip + i__;
	b[i__] += p[ii];
/* L131: */
    }
    goto L134;
L132:
    i__2 = mr;
    for (i__ = 1; i__ <= i__2; ++i__) {
	b[i__] = b[i__] + q[i__ + jp2 * q_dim1] - q[i__ + jp1 * q_dim1];
/* L133: */
    }
L134:
    cosgen_(&i2r, &c__1, &c_b9, &c_b18, &tcos[1]);
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
/* L135: */
    }
    if (lr == 0) {
	goto L136;
    }
    cosgen_(&lr, &c__1, &fnum2, &c_b9, &tcos[i2r + 1]);
    s1merg_(&tcos[1], &c__0, &i2r, &i2r, &lr, &kr);
    goto L138;
L136:
    i__2 = i2r;
    for (i__ = 1; i__ <= i__2; ++i__) {
	ii = kr + i__;
	tcos[ii] = tcos[i__];
/* L137: */
    }
L138:
    cosgen_(&kr, &c__1, &fnum2, &c_b9, &tcos[1]);
    trix_(&kr, &kr, &mr, &a[1], &bb[1], &c__[1], &b[1], &tcos[1], &d__[1], &w[
	    1]);
    i__2 = mr;
    for (i__ = 1; i__ <= i__2; ++i__) {
	ii = ip + i__;
	q[i__ + j * q_dim1] = q[i__ + jm2 * q_dim1] + p[ii] + b[i__];
/* L139: */
    }
L140:
    lr = kr;
    kr += jr;
L141:
    nr = (nlast - 1) / jr + 1;
    if (nr <= 3) {
	goto L142;
    }
    i2r = jr;
    nrodpr = nrod;
    goto L101;
L142:

/*      BEGIN SOLUTION */

    j = jr + 1;
    jm1 = j - i2r;
    jp1 = j + i2r;
    jm2 = nlast - i2r;
    if (nr == 2) {
	goto L180;
    }
    if (lr != 0) {
	goto L167;
    }
    if (*n != 3) {
	goto L156;
    }

/*     CASE N = 3. */

    switch (np) {
	case 1:  goto L143;
	case 2:  goto L148;
	case 3:  goto L143;
    }
L143:
    i__2 = mr;
    for (i__ = 1; i__ <= i__2; ++i__) {
	b[i__] = q[i__ + (q_dim1 << 1)];
	b2[i__] = q[i__ + q_dim1] + q[i__ + q_dim1 * 3];
	b3[i__] = 0.f;
/* L144: */
    }
    switch (np) {
	case 1:  goto L146;
	case 2:  goto L146;
	case 3:  goto L145;
    }
L145:
    tcos[1] = -1.f;
    tcos[2] = 1.f;
    *k1 = 1;
    goto L147;
L146:
    tcos[1] = -2.f;
    tcos[2] = 1.f;
    tcos[3] = -1.f;
    *k1 = 2;
L147:
    *k2 = 1;
    *k3 = 0;
    *k4 = 0;
    goto L150;
L148:
    i__2 = mr;
    for (i__ = 1; i__ <= i__2; ++i__) {
	b[i__] = q[i__ + (q_dim1 << 1)];
	b2[i__] = q[i__ + q_dim1 * 3];
	b3[i__] = q[i__ + q_dim1];
/* L149: */
    }
    cosgen_(&c__3, &c__1, &c_b9, &c_b18, &tcos[1]);
    tcos[4] = -1.f;
    tcos[5] = 1.f;
    tcos[6] = -1.f;
    tcos[7] = 1.f;
    *k1 = 3;
    *k2 = 2;
    *k3 = 1;
    *k4 = 1;
L150:
    tri3_(&mr, &a[1], &bb[1], &c__[1], k, &b[1], &b2[1], &b3[1], &tcos[1], &
	    d__[1], &w[1], &w2[1], &w3[1]);
    i__2 = mr;
    for (i__ = 1; i__ <= i__2; ++i__) {
	b[i__] = b[i__] + b2[i__] + b3[i__];
/* L151: */
    }
    switch (np) {
	case 1:  goto L153;
	case 2:  goto L153;
	case 3:  goto L152;
    }
L152:
    tcos[1] = 2.f;
    trix_(&c__1, &c__0, &mr, &a[1], &bb[1], &c__[1], &b[1], &tcos[1], &d__[1],
	     &w[1]);
L153:
    i__2 = mr;
    for (i__ = 1; i__ <= i__2; ++i__) {
	q[i__ + (q_dim1 << 1)] = b[i__];
	b[i__] = q[i__ + q_dim1] + b[i__];
/* L154: */
    }
    tcos[1] = fnum * 4.f - 1.f;
    trix_(&c__1, &c__0, &mr, &a[1], &bb[1], &c__[1], &b[1], &tcos[1], &d__[1],
	     &w[1]);
    i__2 = mr;
    for (i__ = 1; i__ <= i__2; ++i__) {
	q[i__ + q_dim1] = b[i__];
/* L155: */
    }
    jr = 1;
    i2r = 0;
    goto L188;

/*     CASE N = 2**P+1 */

L156:
    i__2 = mr;
    for (i__ = 1; i__ <= i__2; ++i__) {
	b[i__] = q[i__ + j * q_dim1] + q[i__ + q_dim1] - q[i__ + jm1 * q_dim1]
		 + q[i__ + nlast * q_dim1] - q[i__ + jm2 * q_dim1];
/* L157: */
    }
    switch (np) {
	case 1:  goto L158;
	case 2:  goto L160;
	case 3:  goto L158;
    }
L158:
    i__2 = mr;
    for (i__ = 1; i__ <= i__2; ++i__) {
	b2[i__] = q[i__ + q_dim1] + q[i__ + nlast * q_dim1] + q[i__ + j * 
		q_dim1] - q[i__ + jm1 * q_dim1] - q[i__ + jp1 * q_dim1];
	b3[i__] = 0.f;
/* L159: */
    }
    *k1 = nlast - 1;
    *k2 = nlast + jr - 1;
    i__2 = jr - 1;
    cosgen_(&i__2, &c__1, &c_b18, &c_b98, &tcos[nlast]);
    tcos[*k2] = (real) ((np << 1) - 4);
    r__1 = .5f - fnum;
    cosgen_(&jr, &c__1, &r__1, &c_b9, &tcos[*k2 + 1]);
    *k3 = (3 - np) / 2;
    i__2 = jr - *k3;
    i__1 = *k2 - *k3;
    i__3 = jr + *k3;
    s1merg_(&tcos[1], k1, &i__2, &i__1, &i__3, &c__0);
    *k1 = *k1 - 1 + *k3;
    cosgen_(&jr, &c__1, &fnum, &c_b9, &tcos[*k1 + 1]);
    *k2 = jr;
    *k3 = 0;
    *k4 = 0;
    goto L162;
L160:
    i__2 = mr;
    for (i__ = 1; i__ <= i__2; ++i__) {
	fi = (q[i__ + j * q_dim1] - q[i__ + jm1 * q_dim1] - q[i__ + jp1 * 
		q_dim1]) / 2.f;
	b2[i__] = q[i__ + q_dim1] + fi;
	b3[i__] = q[i__ + nlast * q_dim1] + fi;
/* L161: */
    }
    *k1 = nlast + jr - 1;
    *k2 = *k1 + jr - 1;
    i__2 = jr - 1;
    cosgen_(&i__2, &c__1, &c_b18, &c_b98, &tcos[*k1 + 1]);
    cosgen_(&nlast, &c__1, &c_b9, &c_b18, &tcos[*k2 + 1]);
    i__2 = jr - 1;
    s1merg_(&tcos[1], k1, &i__2, k2, &nlast, &c__0);
    *k3 = *k1 + nlast - 1;
    *k4 = *k3 + jr;
    cosgen_(&jr, &c__1, &c_b9, &c_b9, &tcos[*k3 + 1]);
    cosgen_(&jr, &c__1, &c_b18, &c_b9, &tcos[*k4 + 1]);
    s1merg_(&tcos[1], k3, &jr, k4, &jr, k1);
    *k2 = nlast - 1;
    *k3 = jr;
    *k4 = jr;
L162:
    tri3_(&mr, &a[1], &bb[1], &c__[1], k, &b[1], &b2[1], &b3[1], &tcos[1], &
	    d__[1], &w[1], &w2[1], &w3[1]);
    i__2 = mr;
    for (i__ = 1; i__ <= i__2; ++i__) {
	b[i__] = b[i__] + b2[i__] + b3[i__];
/* L163: */
    }
    if (np != 3) {
	goto L164;
    }
    tcos[1] = 2.f;
    trix_(&c__1, &c__0, &mr, &a[1], &bb[1], &c__[1], &b[1], &tcos[1], &d__[1],
	     &w[1]);
L164:
    i__2 = mr;
    for (i__ = 1; i__ <= i__2; ++i__) {
	q[i__ + j * q_dim1] = b[i__] + (q[i__ + j * q_dim1] - q[i__ + jm1 * 
		q_dim1] - q[i__ + jp1 * q_dim1]) * .5f;
	b[i__] = q[i__ + j * q_dim1] + q[i__ + q_dim1];
/* L165: */
    }
    cosgen_(&jr, &c__1, &fnum, &c_b9, &tcos[1]);
    trix_(&jr, &c__0, &mr, &a[1], &bb[1], &c__[1], &b[1], &tcos[1], &d__[1], &
	    w[1]);
    i__2 = mr;
    for (i__ = 1; i__ <= i__2; ++i__) {
	q[i__ + q_dim1] = q[i__ + q_dim1] - q[i__ + jm1 * q_dim1] + b[i__];
/* L166: */
    }
    goto L188;

/*     CASE OF GENERAL N WITH NR = 3 . */

L167:
    i__2 = mr;
    for (i__ = 1; i__ <= i__2; ++i__) {
	b[i__] = q[i__ + q_dim1] - q[i__ + jm1 * q_dim1] + q[i__ + j * q_dim1]
		;
/* L168: */
    }
    if (nrod != 0) {
	goto L170;
    }
    i__2 = mr;
    for (i__ = 1; i__ <= i__2; ++i__) {
	ii = ip + i__;
	b[i__] += p[ii];
/* L169: */
    }
    goto L172;
L170:
    i__2 = mr;
    for (i__ = 1; i__ <= i__2; ++i__) {
	b[i__] = b[i__] + q[i__ + nlast * q_dim1] - q[i__ + jm2 * q_dim1];
/* L171: */
    }
L172:
    i__2 = mr;
    for (i__ = 1; i__ <= i__2; ++i__) {
	t = (q[i__ + j * q_dim1] - q[i__ + jm1 * q_dim1] - q[i__ + jp1 * 
		q_dim1]) * .5f;
	q[i__ + j * q_dim1] = t;
	b2[i__] = q[i__ + nlast * q_dim1] + t;
	b3[i__] = q[i__ + q_dim1] + t;
/* L173: */
    }
    *k1 = kr + (jr << 1);
    i__2 = jr - 1;
    cosgen_(&i__2, &c__1, &c_b18, &c_b98, &tcos[*k1 + 1]);
    *k2 = *k1 + jr;
    tcos[*k2] = (real) ((np << 1) - 4);
    *k4 = (np - 1) * (3 - np);
    *k3 = *k2 + 1 - *k4;
    i__2 = kr + jr + *k4;
    r__1 = *k4 / 2.f;
    r__2 = 1.f - *k4;
    cosgen_(&i__2, &c__1, &r__1, &r__2, &tcos[*k3]);
    *k4 = 1 - np / 3;
    i__2 = jr - *k4;
    i__1 = *k2 - *k4;
    i__3 = kr + jr + *k4;
    s1merg_(&tcos[1], k1, &i__2, &i__1, &i__3, &c__0);
    if (np == 3) {
	--(*k1);
    }
    *k2 = kr + jr;
    *k4 = *k1 + *k2;
    cosgen_(&kr, &c__1, &fnum2, &c_b9, &tcos[*k4 + 1]);
    *k3 = *k4 + kr;
    cosgen_(&jr, &c__1, &fnum, &c_b9, &tcos[*k3 + 1]);
    s1merg_(&tcos[1], k4, &kr, k3, &jr, k1);
    *k4 = *k3 + jr;
    cosgen_(&lr, &c__1, &fnum2, &c_b9, &tcos[*k4 + 1]);
    i__2 = *k1 + *k2;
    s1merg_(&tcos[1], k3, &jr, k4, &lr, &i__2);
    cosgen_(&kr, &c__1, &fnum2, &c_b9, &tcos[*k3 + 1]);
    *k3 = kr;
    *k4 = kr;
    tri3_(&mr, &a[1], &bb[1], &c__[1], k, &b[1], &b2[1], &b3[1], &tcos[1], &
	    d__[1], &w[1], &w2[1], &w3[1]);
    i__2 = mr;
    for (i__ = 1; i__ <= i__2; ++i__) {
	b[i__] = b[i__] + b2[i__] + b3[i__];
/* L174: */
    }
    if (np != 3) {
	goto L175;
    }
    tcos[1] = 2.f;
    trix_(&c__1, &c__0, &mr, &a[1], &bb[1], &c__[1], &b[1], &tcos[1], &d__[1],
	     &w[1]);
L175:
    i__2 = mr;
    for (i__ = 1; i__ <= i__2; ++i__) {
	q[i__ + j * q_dim1] += b[i__];
	b[i__] = q[i__ + q_dim1] + q[i__ + j * q_dim1];
/* L176: */
    }
    cosgen_(&jr, &c__1, &fnum, &c_b9, &tcos[1]);
    trix_(&jr, &c__0, &mr, &a[1], &bb[1], &c__[1], &b[1], &tcos[1], &d__[1], &
	    w[1]);
    if (jr != 1) {
	goto L178;
    }
    i__2 = mr;
    for (i__ = 1; i__ <= i__2; ++i__) {
	q[i__ + q_dim1] = b[i__];
/* L177: */
    }
    goto L188;
L178:
    i__2 = mr;
    for (i__ = 1; i__ <= i__2; ++i__) {
	q[i__ + q_dim1] = q[i__ + q_dim1] - q[i__ + jm1 * q_dim1] + b[i__];
/* L179: */
    }
    goto L188;
L180:

/*     CASE OF GENERAL N AND NR = 2 . */

    i__2 = mr;
    for (i__ = 1; i__ <= i__2; ++i__) {
	ii = ip + i__;
	b3[i__] = 0.f;
	b[i__] = q[i__ + q_dim1] + p[ii];
	q[i__ + q_dim1] -= q[i__ + jm1 * q_dim1];
	b2[i__] = q[i__ + q_dim1] + q[i__ + nlast * q_dim1];
/* L181: */
    }
    *k1 = kr + jr;
    *k2 = *k1 + jr;
    i__2 = jr - 1;
    cosgen_(&i__2, &c__1, &c_b18, &c_b98, &tcos[*k1 + 1]);
    switch (np) {
	case 1:  goto L182;
	case 2:  goto L183;
	case 3:  goto L182;
    }
L182:
    tcos[*k2] = (real) ((np << 1) - 4);
    cosgen_(&kr, &c__1, &c_b18, &c_b98, &tcos[*k2 + 1]);
    goto L184;
L183:
    i__2 = kr + 1;
    cosgen_(&i__2, &c__1, &c_b9, &c_b18, &tcos[*k2]);
L184:
    *k4 = 1 - np / 3;
    i__2 = jr - *k4;
    i__1 = *k2 - *k4;
    i__3 = kr + *k4;
    s1merg_(&tcos[1], k1, &i__2, &i__1, &i__3, &c__0);
    if (np == 3) {
	--(*k1);
    }
    *k2 = kr;
    cosgen_(&kr, &c__1, &fnum2, &c_b9, &tcos[*k1 + 1]);
    *k4 = *k1 + kr;
    cosgen_(&lr, &c__1, &fnum2, &c_b9, &tcos[*k4 + 1]);
    *k3 = lr;
    *k4 = 0;
    tri3_(&mr, &a[1], &bb[1], &c__[1], k, &b[1], &b2[1], &b3[1], &tcos[1], &
	    d__[1], &w[1], &w2[1], &w3[1]);
    i__2 = mr;
    for (i__ = 1; i__ <= i__2; ++i__) {
	b[i__] += b2[i__];
/* L185: */
    }
    if (np != 3) {
	goto L186;
    }
    tcos[1] = 2.f;
    trix_(&c__1, &c__0, &mr, &a[1], &bb[1], &c__[1], &b[1], &tcos[1], &d__[1],
	     &w[1]);
L186:
    i__2 = mr;
    for (i__ = 1; i__ <= i__2; ++i__) {
	q[i__ + q_dim1] += b[i__];
/* L187: */
    }
L188:

/*     START BACK SUBSTITUTION. */

    j = nlast - jr;
    i__2 = mr;
    for (i__ = 1; i__ <= i__2; ++i__) {
	b[i__] = q[i__ + nlast * q_dim1] + q[i__ + j * q_dim1];
/* L189: */
    }
    jm2 = nlast - i2r;
    if (jr != 1) {
	goto L191;
    }
    i__2 = mr;
    for (i__ = 1; i__ <= i__2; ++i__) {
	q[i__ + nlast * q_dim1] = 0.f;
/* L190: */
    }
    goto L195;
L191:
    if (nrod != 0) {
	goto L193;
    }
    i__2 = mr;
    for (i__ = 1; i__ <= i__2; ++i__) {
	ii = ip + i__;
	q[i__ + nlast * q_dim1] = p[ii];
/* L192: */
    }
    ip -= mr;
    goto L195;
L193:
    i__2 = mr;
    for (i__ = 1; i__ <= i__2; ++i__) {
	q[i__ + nlast * q_dim1] -= q[i__ + jm2 * q_dim1];
/* L194: */
    }
L195:
    cosgen_(&kr, &c__1, &fnum2, &c_b9, &tcos[1]);
    cosgen_(&lr, &c__1, &fnum2, &c_b9, &tcos[kr + 1]);
    trix_(&kr, &lr, &mr, &a[1], &bb[1], &c__[1], &b[1], &tcos[1], &d__[1], &w[
	    1]);
    i__2 = mr;
    for (i__ = 1; i__ <= i__2; ++i__) {
	q[i__ + nlast * q_dim1] += b[i__];
/* L196: */
    }
    nlastp = nlast;
L197:
    jstep = jr;
    jr = i2r;
    i2r /= 2;
    if (jr == 0) {
	goto L210;
    }
    jstart = jr + 1;
    kr -= jr;
    if (nlast + jr > *n) {
	goto L198;
    }
    kr -= jr;
    nlast += jr;
    jstop = nlast - jstep;
    goto L199;
L198:
    jstop = nlast - jr;
L199:
    lr = kr - jr;
    cosgen_(&jr, &c__1, &c_b9, &c_b18, &tcos[1]);
    i__2 = jstop;
    i__1 = jstep;
    for (j = jstart; i__1 < 0 ? j >= i__2 : j <= i__2; j += i__1) {
	jm2 = j - jr;
	jp2 = j + jr;
	if (j != jr) {
	    goto L201;
	}
	i__3 = mr;
	for (i__ = 1; i__ <= i__3; ++i__) {
	    b[i__] = q[i__ + j * q_dim1] + q[i__ + jp2 * q_dim1];
/* L200: */
	}
	goto L203;
L201:
	i__3 = mr;
	for (i__ = 1; i__ <= i__3; ++i__) {
	    b[i__] = q[i__ + j * q_dim1] + q[i__ + jm2 * q_dim1] + q[i__ + 
		    jp2 * q_dim1];
/* L202: */
	}
L203:
	if (jr != 1) {
	    goto L205;
	}
	i__3 = mr;
	for (i__ = 1; i__ <= i__3; ++i__) {
	    q[i__ + j * q_dim1] = 0.f;
/* L204: */
	}
	goto L207;
L205:
	jm1 = j - i2r;
	jp1 = j + i2r;
	i__3 = mr;
	for (i__ = 1; i__ <= i__3; ++i__) {
	    q[i__ + j * q_dim1] = (q[i__ + j * q_dim1] - q[i__ + jm1 * q_dim1]
		     - q[i__ + jp1 * q_dim1]) * .5f;
/* L206: */
	}
L207:
	trix_(&jr, &c__0, &mr, &a[1], &bb[1], &c__[1], &b[1], &tcos[1], &d__[
		1], &w[1]);
	i__3 = mr;
	for (i__ = 1; i__ <= i__3; ++i__) {
	    q[i__ + j * q_dim1] += b[i__];
/* L208: */
	}
/* L209: */
    }
    nrod = 1;
    if (nlast + i2r <= *n) {
	nrod = 0;
    }
    if (nlastp != nlast) {
	goto L188;
    }
    goto L197;
L210:

/*     RETURN STORAGE REQUIREMENTS FOR P VECTORS. */

    w[1] = (real) ipstor;
    return 0;
} /* postg2_ */

#undef k4
#undef k3
#undef k2
#undef k1
#undef k


