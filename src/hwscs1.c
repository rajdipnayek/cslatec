/* hwscs1.f -- translated by f2c (version 12.02.01).
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

/* DECK HWSCS1 */
/* Subroutine */ int hwscs1_(integer *intl, real *ts, real *tf, integer *m, 
	integer *mbdcnd, real *bdts, real *bdtf, real *rs, real *rf, integer *
	n, integer *nbdcnd, real *bdrs, real *bdrf, real *elmbda, real *f, 
	integer *idimf, real *pertrb, real *w, real *s, real *an, real *bn, 
	real *cn, real *r__, real *am, real *bm, real *cm, real *sint, real *
	bmh)
{
    /* System generated locals */
    integer f_dim1, f_offset, i__1, i__2;
    real r__1;

    /* Local variables */
    static integer i__, j, l;
    static real r2, t1, ar, at, dr, ct, cr;
    static integer mp, np;
    static real xp, dr2, rf2;
    static integer mp1, np1;
    static real rs2, hne, hdr, dth;
    static integer jrf, itf;
    static real tdr, tdt, czr;
    static integer its;
    static real wtf;
    static integer jrs;
    static real wrf, yph, sum, rsq, xps, yps, wrs, wts, wrz;
    static integer iflg;
    static real hdth;
    static integer jrfm, itfm;
    static real yhld;
    static integer ictr, munk, nunk;
    static real sdts;
    static integer itsp, jrsp;
    static real wtnm, theta;
    static integer ising;
    extern /* Subroutine */ int blktri_(integer *, integer *, integer *, real 
	    *, real *, real *, integer *, integer *, real *, real *, real *, 
	    integer *, real *, integer *, real *);
    static integer ierror;

/* ***BEGIN PROLOGUE  HWSCS1 */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to HWSCSP */
/* ***LIBRARY   SLATEC */
/* ***TYPE      SINGLE PRECISION (HWSCS1-S) */
/* ***AUTHOR  (UNKNOWN) */
/* ***SEE ALSO  HWSCSP */
/* ***ROUTINES CALLED  BLKTRI */
/* ***REVISION HISTORY  (YYMMDD) */
/*   801001  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   891009  Removed unreferenced variables.  (WRB) */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900402  Added TYPE section.  (WRB) */
/* ***END PROLOGUE  HWSCS1 */
/* ***FIRST EXECUTABLE STATEMENT  HWSCS1 */
    /* Parameter adjustments */
    --bdts;
    --bdtf;
    --bdrs;
    --bdrf;
    f_dim1 = *idimf;
    f_offset = 1 + f_dim1;
    f -= f_offset;
    --w;
    --s;
    --an;
    --bn;
    --cn;
    --r__;
    --am;
    --bm;
    --cm;
    --sint;
    --bmh;

    /* Function Body */
    mp1 = *m + 1;
    dth = (*tf - *ts) / *m;
    tdt = dth + dth;
    hdth = dth / 2.f;
    sdts = 1.f / (dth * dth);
    i__1 = mp1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	theta = *ts + (i__ - 1) * dth;
	sint[i__] = sin(theta);
	if (sint[i__] != 0.f) {
	    goto L101;
	} else {
	    goto L102;
	}
L101:
	t1 = sdts / sint[i__];
	am[i__] = t1 * sin(theta - hdth);
	cm[i__] = t1 * sin(theta + hdth);
	bm[i__] = -(am[i__] + cm[i__]);
L102:
	;
    }
    np1 = *n + 1;
    dr = (*rf - *rs) / *n;
    hdr = dr / 2.f;
    tdr = dr + dr;
    dr2 = dr * dr;
    czr = dth * 6.f / (dr2 * (cos(*ts) - cos(*tf)));
    i__1 = np1;
    for (j = 1; j <= i__1; ++j) {
	r__[j] = *rs + (j - 1) * dr;
/* Computing 2nd power */
	r__1 = r__[j] - hdr;
	an[j] = r__1 * r__1 / dr2;
/* Computing 2nd power */
	r__1 = r__[j] + hdr;
	cn[j] = r__1 * r__1 / dr2;
	bn[j] = -(an[j] + cn[j]);
/* L103: */
    }
    mp = 1;
    np = 1;

/* BOUNDARY CONDITION AT PHI=PS */

    switch (*mbdcnd) {
	case 1:  goto L104;
	case 2:  goto L104;
	case 3:  goto L105;
	case 4:  goto L105;
	case 5:  goto L106;
	case 6:  goto L106;
	case 7:  goto L104;
	case 8:  goto L105;
	case 9:  goto L106;
    }
L104:
    at = am[2];
    its = 2;
    goto L107;
L105:
    at = am[1];
    its = 1;
    cm[1] += am[1];
    goto L107;
L106:
    its = 1;
    bm[1] = sdts * -4.f;
    cm[1] = -bm[1];

/* BOUNDARY CONDITION AT PHI=PF */

L107:
    switch (*mbdcnd) {
	case 1:  goto L108;
	case 2:  goto L109;
	case 3:  goto L109;
	case 4:  goto L108;
	case 5:  goto L108;
	case 6:  goto L109;
	case 7:  goto L110;
	case 8:  goto L110;
	case 9:  goto L110;
    }
L108:
    ct = cm[*m];
    itf = *m;
    goto L111;
L109:
    ct = cm[*m + 1];
    am[*m + 1] += cm[*m + 1];
    itf = *m + 1;
    goto L111;
L110:
    itf = *m + 1;
    am[*m + 1] = sdts * 4.f;
    bm[*m + 1] = -am[*m + 1];
L111:
    wts = sint[its + 1] * am[its + 1] / cm[its];
    wtf = sint[itf - 1] * cm[itf - 1] / am[itf];
    itsp = its + 1;
    itfm = itf - 1;

/* BOUNDARY CONDITION AT R=RS */

    ictr = 0;
    switch (*nbdcnd) {
	case 1:  goto L112;
	case 2:  goto L112;
	case 3:  goto L113;
	case 4:  goto L113;
	case 5:  goto L114;
	case 6:  goto L114;
    }
L112:
    ar = an[2];
    jrs = 2;
    goto L118;
L113:
    ar = an[1];
    jrs = 1;
    cn[1] += an[1];
    goto L118;
L114:
    jrs = 2;
    ictr = 1;
    s[*n] = an[*n] / bn[*n];
    i__1 = *n;
    for (j = 3; j <= i__1; ++j) {
	l = *n - j + 2;
	s[l] = an[l] / (bn[l] - cn[l] * s[l + 1]);
/* L115: */
    }
    s[2] = -s[2];
    i__1 = *n;
    for (j = 3; j <= i__1; ++j) {
	s[j] = -s[j] * s[j - 1];
/* L116: */
    }
    wtnm = wts + wtf;
    i__1 = itfm;
    for (i__ = itsp; i__ <= i__1; ++i__) {
	wtnm += sint[i__];
/* L117: */
    }
    yps = czr * wtnm * (s[2] - 1.f);

/* BOUNDARY CONDITION AT R=RF */

L118:
    switch (*nbdcnd) {
	case 1:  goto L119;
	case 2:  goto L120;
	case 3:  goto L120;
	case 4:  goto L119;
	case 5:  goto L119;
	case 6:  goto L120;
    }
L119:
    cr = cn[*n];
    jrf = *n;
    goto L121;
L120:
    cr = cn[*n + 1];
    an[*n + 1] += cn[*n + 1];
    jrf = *n + 1;
L121:
/* Computing 2nd power */
    r__1 = r__[jrs];
    wrs = an[jrs + 1] * (r__1 * r__1) / cn[jrs];
/* Computing 2nd power */
    r__1 = r__[jrf];
    wrf = cn[jrf - 1] * (r__1 * r__1) / an[jrf];
    wrz = an[jrs] / czr;
    jrsp = jrs + 1;
    jrfm = jrf - 1;
    munk = itf - its + 1;
    nunk = jrf - jrs + 1;
    i__1 = itf;
    for (i__ = its; i__ <= i__1; ++i__) {
	bmh[i__] = bm[i__];
/* L122: */
    }
    ising = 0;
    switch (*nbdcnd) {
	case 1:  goto L132;
	case 2:  goto L132;
	case 3:  goto L123;
	case 4:  goto L132;
	case 5:  goto L132;
	case 6:  goto L123;
    }
L123:
    switch (*mbdcnd) {
	case 1:  goto L132;
	case 2:  goto L132;
	case 3:  goto L124;
	case 4:  goto L132;
	case 5:  goto L132;
	case 6:  goto L124;
	case 7:  goto L132;
	case 8:  goto L124;
	case 9:  goto L124;
    }
L124:
    if (*elmbda >= 0.f) {
	goto L125;
    } else {
	goto L132;
    }
L125:
    ising = 1;
    sum = wts * wrs + wts * wrf + wtf * wrs + wtf * wrf;
    if (ictr != 0) {
	goto L126;
    } else {
	goto L127;
    }
L126:
    sum += wrz;
L127:
    i__1 = jrfm;
    for (j = jrsp; j <= i__1; ++j) {
/* Computing 2nd power */
	r__1 = r__[j];
	r2 = r__1 * r__1;
	i__2 = itfm;
	for (i__ = itsp; i__ <= i__2; ++i__) {
	    sum += r2 * sint[i__];
/* L128: */
	}
/* L129: */
    }
    i__1 = jrfm;
    for (j = jrsp; j <= i__1; ++j) {
/* Computing 2nd power */
	r__1 = r__[j];
	sum += (wts + wtf) * (r__1 * r__1);
/* L130: */
    }
    i__1 = itfm;
    for (i__ = itsp; i__ <= i__1; ++i__) {
	sum += (wrs + wrf) * sint[i__];
/* L131: */
    }
    hne = sum;
L132:
    switch (*mbdcnd) {
	case 1:  goto L133;
	case 2:  goto L133;
	case 3:  goto L133;
	case 4:  goto L133;
	case 5:  goto L134;
	case 6:  goto L134;
	case 7:  goto L133;
	case 8:  goto L133;
	case 9:  goto L134;
    }
L133:
/* Computing 2nd power */
    r__1 = sint[its];
    bm[its] = bmh[its] + *elmbda / (r__1 * r__1);
L134:
    switch (*mbdcnd) {
	case 1:  goto L135;
	case 2:  goto L135;
	case 3:  goto L135;
	case 4:  goto L135;
	case 5:  goto L135;
	case 6:  goto L135;
	case 7:  goto L136;
	case 8:  goto L136;
	case 9:  goto L136;
    }
L135:
/* Computing 2nd power */
    r__1 = sint[itf];
    bm[itf] = bmh[itf] + *elmbda / (r__1 * r__1);
L136:
    i__1 = itfm;
    for (i__ = itsp; i__ <= i__1; ++i__) {
/* Computing 2nd power */
	r__1 = sint[i__];
	bm[i__] = bmh[i__] + *elmbda / (r__1 * r__1);
/* L137: */
    }
    switch (*mbdcnd) {
	case 1:  goto L138;
	case 2:  goto L138;
	case 3:  goto L140;
	case 4:  goto L140;
	case 5:  goto L142;
	case 6:  goto L142;
	case 7:  goto L138;
	case 8:  goto L140;
	case 9:  goto L142;
    }
L138:
    i__1 = jrf;
    for (j = jrs; j <= i__1; ++j) {
/* Computing 2nd power */
	r__1 = r__[j];
	f[j * f_dim1 + 2] -= at * f[j * f_dim1 + 1] / (r__1 * r__1);
/* L139: */
    }
    goto L142;
L140:
    i__1 = jrf;
    for (j = jrs; j <= i__1; ++j) {
/* Computing 2nd power */
	r__1 = r__[j];
	f[j * f_dim1 + 1] += tdt * bdts[j] * at / (r__1 * r__1);
/* L141: */
    }
L142:
    switch (*mbdcnd) {
	case 1:  goto L143;
	case 2:  goto L145;
	case 3:  goto L145;
	case 4:  goto L143;
	case 5:  goto L143;
	case 6:  goto L145;
	case 7:  goto L147;
	case 8:  goto L147;
	case 9:  goto L147;
    }
L143:
    i__1 = jrf;
    for (j = jrs; j <= i__1; ++j) {
/* Computing 2nd power */
	r__1 = r__[j];
	f[*m + j * f_dim1] -= ct * f[*m + 1 + j * f_dim1] / (r__1 * r__1);
/* L144: */
    }
    goto L147;
L145:
    i__1 = jrf;
    for (j = jrs; j <= i__1; ++j) {
/* Computing 2nd power */
	r__1 = r__[j];
	f[*m + 1 + j * f_dim1] -= tdt * bdtf[j] * ct / (r__1 * r__1);
/* L146: */
    }
L147:
    switch (*nbdcnd) {
	case 1:  goto L151;
	case 2:  goto L151;
	case 3:  goto L153;
	case 4:  goto L153;
	case 5:  goto L148;
	case 6:  goto L148;
    }
L148:
    if (*mbdcnd - 3 != 0) {
	goto L155;
    } else {
	goto L149;
    }
L149:
    yhld = f[its + f_dim1] - czr / tdt * (sin(*tf) * bdtf[2] - sin(*ts) * 
	    bdts[2]);
    i__1 = mp1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	f[i__ + f_dim1] = yhld;
/* L150: */
    }
    goto L155;
L151:
/* Computing 2nd power */
    r__1 = *rs + dr;
    rs2 = r__1 * r__1;
    i__1 = itf;
    for (i__ = its; i__ <= i__1; ++i__) {
	f[i__ + (f_dim1 << 1)] -= ar * f[i__ + f_dim1] / rs2;
/* L152: */
    }
    goto L155;
L153:
    i__1 = itf;
    for (i__ = its; i__ <= i__1; ++i__) {
/* Computing 2nd power */
	r__1 = *rs;
	f[i__ + f_dim1] += tdr * bdrs[i__] * ar / (r__1 * r__1);
/* L154: */
    }
L155:
    switch (*nbdcnd) {
	case 1:  goto L156;
	case 2:  goto L158;
	case 3:  goto L158;
	case 4:  goto L156;
	case 5:  goto L156;
	case 6:  goto L158;
    }
L156:
/* Computing 2nd power */
    r__1 = *rf - dr;
    rf2 = r__1 * r__1;
    i__1 = itf;
    for (i__ = its; i__ <= i__1; ++i__) {
	f[i__ + *n * f_dim1] -= cr * f[i__ + (*n + 1) * f_dim1] / rf2;
/* L157: */
    }
    goto L160;
L158:
    i__1 = itf;
    for (i__ = its; i__ <= i__1; ++i__) {
/* Computing 2nd power */
	r__1 = *rf;
	f[i__ + (*n + 1) * f_dim1] -= tdr * bdrf[i__] * cr / (r__1 * r__1);
/* L159: */
    }
L160:
    *pertrb = 0.f;
    if (ising != 0) {
	goto L161;
    } else {
	goto L170;
    }
L161:
    sum = wts * wrs * f[its + jrs * f_dim1] + wts * wrf * f[its + jrf * 
	    f_dim1] + wtf * wrs * f[itf + jrs * f_dim1] + wtf * wrf * f[itf + 
	    jrf * f_dim1];
    if (ictr != 0) {
	goto L162;
    } else {
	goto L163;
    }
L162:
    sum += wrz * f[its + f_dim1];
L163:
    i__1 = jrfm;
    for (j = jrsp; j <= i__1; ++j) {
/* Computing 2nd power */
	r__1 = r__[j];
	r2 = r__1 * r__1;
	i__2 = itfm;
	for (i__ = itsp; i__ <= i__2; ++i__) {
	    sum += r2 * sint[i__] * f[i__ + j * f_dim1];
/* L164: */
	}
/* L165: */
    }
    i__1 = jrfm;
    for (j = jrsp; j <= i__1; ++j) {
/* Computing 2nd power */
	r__1 = r__[j];
	sum += r__1 * r__1 * (wts * f[its + j * f_dim1] + wtf * f[itf + j * 
		f_dim1]);
/* L166: */
    }
    i__1 = itfm;
    for (i__ = itsp; i__ <= i__1; ++i__) {
	sum += sint[i__] * (wrs * f[i__ + jrs * f_dim1] + wrf * f[i__ + jrf * 
		f_dim1]);
/* L167: */
    }
    *pertrb = sum / hne;
    i__1 = np1;
    for (j = 1; j <= i__1; ++j) {
	i__2 = mp1;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    f[i__ + j * f_dim1] -= *pertrb;
/* L168: */
	}
/* L169: */
    }
L170:
    i__1 = jrf;
    for (j = jrs; j <= i__1; ++j) {
/* Computing 2nd power */
	r__1 = r__[j];
	rsq = r__1 * r__1;
	i__2 = itf;
	for (i__ = its; i__ <= i__2; ++i__) {
	    f[i__ + j * f_dim1] = rsq * f[i__ + j * f_dim1];
/* L171: */
	}
/* L172: */
    }
    iflg = *intl;
L173:
    blktri_(&iflg, &np, &nunk, &an[jrs], &bn[jrs], &cn[jrs], &mp, &munk, &am[
	    its], &bm[its], &cm[its], idimf, &f[its + jrs * f_dim1], &ierror, 
	    &w[1]);
    ++iflg;
    if (iflg - 1 != 0) {
	goto L174;
    } else {
	goto L173;
    }
L174:
    if (*nbdcnd != 0) {
	goto L177;
    } else {
	goto L175;
    }
L175:
    i__1 = mp1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	f[i__ + (jrf + 1) * f_dim1] = f[i__ + jrs * f_dim1];
/* L176: */
    }
L177:
    if (*mbdcnd != 0) {
	goto L180;
    } else {
	goto L178;
    }
L178:
    i__1 = np1;
    for (j = 1; j <= i__1; ++j) {
	f[itf + 1 + j * f_dim1] = f[its + j * f_dim1];
/* L179: */
    }
L180:
    xp = 0.f;
    if (ictr != 0) {
	goto L181;
    } else {
	goto L188;
    }
L181:
    if (ising != 0) {
	goto L186;
    } else {
	goto L182;
    }
L182:
    sum = wts * f[its + (f_dim1 << 1)] + wtf * f[itf + (f_dim1 << 1)];
    i__1 = itfm;
    for (i__ = itsp; i__ <= i__1; ++i__) {
	sum += sint[i__] * f[i__ + (f_dim1 << 1)];
/* L183: */
    }
    yph = czr * sum;
    xp = (f[its + f_dim1] - yph) / yps;
    i__1 = jrf;
    for (j = jrs; j <= i__1; ++j) {
	xps = xp * s[j];
	i__2 = itf;
	for (i__ = its; i__ <= i__2; ++i__) {
	    f[i__ + j * f_dim1] += xps;
/* L184: */
	}
/* L185: */
    }
L186:
    i__1 = mp1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	f[i__ + f_dim1] = xp;
/* L187: */
    }
L188:
    return 0;
} /* hwscs1_ */

