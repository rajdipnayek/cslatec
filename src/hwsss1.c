/* hwsss1.f -- translated by f2c (version 12.02.01).
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

/* DECK HWSSS1 */
/* Subroutine */ int hwsss1_(real *ts, real *tf, integer *m, integer *mbdcnd, 
	real *bdts, real *bdtf, real *ps, real *pf, integer *n, integer *
	nbdcnd, real *bdps, real *bdpf, real *elmbda, real *f, integer *idimf,
	 real *pertrb, real *am, real *bm, real *cm, real *sn, real *ss, real 
	*sint, real *d__)
{
    /* System generated locals */
    integer f_dim1, f_offset, i__1, i__2;

    /* Local variables */
    static integer i__, j;
    static real t1, cf;
    static integer ii;
    static real fm, fn, cp, at, ct, wp;
    static integer mp1, np1, iid;
    static real dfn, hld, fjj, hne, den, dfs, dth;
    static integer mbr, nbr, itf, jpf;
    static real dnn, dsn;
    static integer inp;
    static real tdp, cnp, dss, tdt;
    static integer isp;
    static real wpf;
    static integer jps, its;
    static real wtf, dns, csp, rtn, sum, rts, wps, wts, fim1, dth2, sum1, 
	    sum2, dphi, hdth;
    static integer jpfm, itfm;
    static real yhld;
    static integer munk, nunk, jpsp, itsp;
    static real dphi2, theta;
    static integer ising;
    extern /* Subroutine */ int genbun_(integer *, integer *, integer *, 
	    integer *, real *, real *, real *, integer *, real *, integer *, 
	    real *);
    static integer ierror;

/* ***BEGIN PROLOGUE  HWSSS1 */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to HWSSSP */
/* ***LIBRARY   SLATEC */
/* ***TYPE      SINGLE PRECISION (HWSSS1-S) */
/* ***AUTHOR  (UNKNOWN) */
/* ***SEE ALSO  HWSSSP */
/* ***ROUTINES CALLED  GENBUN */
/* ***REVISION HISTORY  (YYMMDD) */
/*   801001  DATE WRITTEN */
/*   891009  Removed unreferenced variables.  (WRB) */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900402  Added TYPE section.  (WRB) */
/* ***END PROLOGUE  HWSSS1 */

/* ***FIRST EXECUTABLE STATEMENT  HWSSS1 */
    /* Parameter adjustments */
    --bdts;
    --bdtf;
    --bdps;
    --bdpf;
    f_dim1 = *idimf;
    f_offset = 1 + f_dim1;
    f -= f_offset;
    --am;
    --bm;
    --cm;
    --sn;
    --ss;
    --sint;
    --d__;

    /* Function Body */
    mp1 = *m + 1;
    np1 = *n + 1;
    fn = (real) (*n);
    fm = (real) (*m);
    dth = (*tf - *ts) / fm;
    hdth = dth / 2.f;
    tdt = dth + dth;
    dphi = (*pf - *ps) / fn;
    tdp = dphi + dphi;
    dphi2 = dphi * dphi;
    dth2 = dth * dth;
    cp = 4.f / (fn * dth2);
    wp = fn * sin(hdth) / 4.f;
    i__1 = mp1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	fim1 = (real) (i__ - 1);
	theta = fim1 * dth + *ts;
	sint[i__] = sin(theta);
	if (sint[i__] != 0.f) {
	    goto L101;
	} else {
	    goto L102;
	}
L101:
	t1 = 1.f / (dth2 * sint[i__]);
	am[i__] = t1 * sin(theta - hdth);
	cm[i__] = t1 * sin(theta + hdth);
	bm[i__] = -am[i__] - cm[i__] + *elmbda;
L102:
	;
    }
    inp = 0;
    isp = 0;

/* BOUNDARY CONDITION AT THETA=TS */

    mbr = *mbdcnd + 1;
    switch (mbr) {
	case 1:  goto L103;
	case 2:  goto L104;
	case 3:  goto L104;
	case 4:  goto L105;
	case 5:  goto L105;
	case 6:  goto L106;
	case 7:  goto L106;
	case 8:  goto L104;
	case 9:  goto L105;
	case 10:  goto L106;
    }
L103:
    its = 1;
    goto L107;
L104:
    at = am[2];
    its = 2;
    goto L107;
L105:
    at = am[1];
    its = 1;
    cm[1] = am[1] + cm[1];
    goto L107;
L106:
    at = am[2];
    inp = 1;
    its = 2;

/* BOUNDARY CONDITION THETA=TF */

L107:
    switch (mbr) {
	case 1:  goto L108;
	case 2:  goto L109;
	case 3:  goto L110;
	case 4:  goto L110;
	case 5:  goto L109;
	case 6:  goto L109;
	case 7:  goto L110;
	case 8:  goto L111;
	case 9:  goto L111;
	case 10:  goto L111;
    }
L108:
    itf = *m;
    goto L112;
L109:
    ct = cm[*m];
    itf = *m;
    goto L112;
L110:
    ct = cm[*m + 1];
    am[*m + 1] += cm[*m + 1];
    itf = *m + 1;
    goto L112;
L111:
    itf = *m;
    isp = 1;
    ct = cm[*m];

/* COMPUTE HOMOGENEOUS SOLUTION WITH SOLUTION AT POLE EQUAL TO ONE */

L112:
    itsp = its + 1;
    itfm = itf - 1;
    wts = sint[its + 1] * am[its + 1] / cm[its];
    wtf = sint[itf - 1] * cm[itf - 1] / am[itf];
    munk = itf - its + 1;
    if (isp <= 0) {
	goto L116;
    } else {
	goto L113;
    }
L113:
    d__[its] = cm[its] / bm[its];
    i__1 = *m;
    for (i__ = itsp; i__ <= i__1; ++i__) {
	d__[i__] = cm[i__] / (bm[i__] - am[i__] * d__[i__ - 1]);
/* L114: */
    }
    ss[*m] = -d__[*m];
    iid = *m - its;
    i__1 = iid;
    for (ii = 1; ii <= i__1; ++ii) {
	i__ = *m - ii;
	ss[i__] = -d__[i__] * ss[i__ + 1];
/* L115: */
    }
    ss[*m + 1] = 1.f;
L116:
    if (inp <= 0) {
	goto L120;
    } else {
	goto L117;
    }
L117:
    sn[1] = 1.f;
    d__[itf] = am[itf] / bm[itf];
    iid = itf - 2;
    i__1 = iid;
    for (ii = 1; ii <= i__1; ++ii) {
	i__ = itf - ii;
	d__[i__] = am[i__] / (bm[i__] - cm[i__] * d__[i__ + 1]);
/* L118: */
    }
    sn[2] = -d__[2];
    i__1 = itf;
    for (i__ = 3; i__ <= i__1; ++i__) {
	sn[i__] = -d__[i__] * sn[i__ - 1];
/* L119: */
    }

/* BOUNDARY CONDITIONS AT PHI=PS */

L120:
    nbr = *nbdcnd + 1;
    wps = 1.f;
    wpf = 1.f;
    switch (nbr) {
	case 1:  goto L121;
	case 2:  goto L122;
	case 3:  goto L122;
	case 4:  goto L123;
	case 5:  goto L123;
    }
L121:
    jps = 1;
    goto L124;
L122:
    jps = 2;
    goto L124;
L123:
    jps = 1;
    wps = .5f;

/* BOUNDARY CONDITION AT PHI=PF */

L124:
    switch (nbr) {
	case 1:  goto L125;
	case 2:  goto L126;
	case 3:  goto L127;
	case 4:  goto L127;
	case 5:  goto L126;
    }
L125:
    jpf = *n;
    goto L128;
L126:
    jpf = *n;
    goto L128;
L127:
    wpf = .5f;
    jpf = *n + 1;
L128:
    jpsp = jps + 1;
    jpfm = jpf - 1;
    nunk = jpf - jps + 1;
    fjj = (real) (jpfm - jpsp + 1);

/* SCALE COEFFICIENTS FOR SUBROUTINE GENBUN */

    i__1 = itf;
    for (i__ = its; i__ <= i__1; ++i__) {
	cf = dphi2 * sint[i__] * sint[i__];
	am[i__] = cf * am[i__];
	bm[i__] = cf * bm[i__];
	cm[i__] = cf * cm[i__];
/* L129: */
    }
    am[its] = 0.f;
    cm[itf] = 0.f;
    ising = 0;
    switch (mbr) {
	case 1:  goto L130;
	case 2:  goto L138;
	case 3:  goto L138;
	case 4:  goto L130;
	case 5:  goto L138;
	case 6:  goto L138;
	case 7:  goto L130;
	case 8:  goto L138;
	case 9:  goto L130;
	case 10:  goto L130;
    }
L130:
    switch (nbr) {
	case 1:  goto L131;
	case 2:  goto L138;
	case 3:  goto L138;
	case 4:  goto L131;
	case 5:  goto L138;
    }
L131:
    if (*elmbda >= 0.f) {
	goto L132;
    } else {
	goto L138;
    }
L132:
    ising = 1;
    sum = wts * wps + wts * wpf + wtf * wps + wtf * wpf;
    if (inp <= 0) {
	goto L134;
    } else {
	goto L133;
    }
L133:
    sum += wp;
L134:
    if (isp <= 0) {
	goto L136;
    } else {
	goto L135;
    }
L135:
    sum += wp;
L136:
    sum1 = 0.f;
    i__1 = itfm;
    for (i__ = itsp; i__ <= i__1; ++i__) {
	sum1 += sint[i__];
/* L137: */
    }
    sum += fjj * (sum1 + wts + wtf);
    sum += (wps + wpf) * sum1;
    hne = sum;
L138:
    switch (mbr) {
	case 1:  goto L146;
	case 2:  goto L142;
	case 3:  goto L142;
	case 4:  goto L144;
	case 5:  goto L144;
	case 6:  goto L139;
	case 7:  goto L139;
	case 8:  goto L142;
	case 9:  goto L144;
	case 10:  goto L139;
    }
L139:
    if (*nbdcnd - 3 != 0) {
	goto L146;
    } else {
	goto L140;
    }
L140:
    yhld = f[jps * f_dim1 + 1] - 4.f / (fn * dphi * dth2) * (bdpf[2] - bdps[2]
	    );
    i__1 = np1;
    for (j = 1; j <= i__1; ++j) {
	f[j * f_dim1 + 1] = yhld;
/* L141: */
    }
    goto L146;
L142:
    i__1 = jpf;
    for (j = jps; j <= i__1; ++j) {
	f[j * f_dim1 + 2] -= at * f[j * f_dim1 + 1];
/* L143: */
    }
    goto L146;
L144:
    i__1 = jpf;
    for (j = jps; j <= i__1; ++j) {
	f[j * f_dim1 + 1] += tdt * bdts[j] * at;
/* L145: */
    }
L146:
    switch (mbr) {
	case 1:  goto L154;
	case 2:  goto L150;
	case 3:  goto L152;
	case 4:  goto L152;
	case 5:  goto L150;
	case 6:  goto L150;
	case 7:  goto L152;
	case 8:  goto L147;
	case 9:  goto L147;
	case 10:  goto L147;
    }
L147:
    if (*nbdcnd - 3 != 0) {
	goto L154;
    } else {
	goto L148;
    }
L148:
    yhld = f[*m + 1 + jps * f_dim1] - 4.f / (fn * dphi * dth2) * (bdpf[*m] - 
	    bdps[*m]);
    i__1 = np1;
    for (j = 1; j <= i__1; ++j) {
	f[*m + 1 + j * f_dim1] = yhld;
/* L149: */
    }
    goto L154;
L150:
    i__1 = jpf;
    for (j = jps; j <= i__1; ++j) {
	f[*m + j * f_dim1] -= ct * f[*m + 1 + j * f_dim1];
/* L151: */
    }
    goto L154;
L152:
    i__1 = jpf;
    for (j = jps; j <= i__1; ++j) {
	f[*m + 1 + j * f_dim1] -= tdt * bdtf[j] * ct;
/* L153: */
    }
L154:
    switch (nbr) {
	case 1:  goto L159;
	case 2:  goto L155;
	case 3:  goto L155;
	case 4:  goto L157;
	case 5:  goto L157;
    }
L155:
    i__1 = itf;
    for (i__ = its; i__ <= i__1; ++i__) {
	f[i__ + (f_dim1 << 1)] -= f[i__ + f_dim1] / (dphi2 * sint[i__] * sint[
		i__]);
/* L156: */
    }
    goto L159;
L157:
    i__1 = itf;
    for (i__ = its; i__ <= i__1; ++i__) {
	f[i__ + f_dim1] += tdp * bdps[i__] / (dphi2 * sint[i__] * sint[i__]);
/* L158: */
    }
L159:
    switch (nbr) {
	case 1:  goto L164;
	case 2:  goto L160;
	case 3:  goto L162;
	case 4:  goto L162;
	case 5:  goto L160;
    }
L160:
    i__1 = itf;
    for (i__ = its; i__ <= i__1; ++i__) {
	f[i__ + *n * f_dim1] -= f[i__ + (*n + 1) * f_dim1] / (dphi2 * sint[
		i__] * sint[i__]);
/* L161: */
    }
    goto L164;
L162:
    i__1 = itf;
    for (i__ = its; i__ <= i__1; ++i__) {
	f[i__ + (*n + 1) * f_dim1] -= tdp * bdpf[i__] / (dphi2 * sint[i__] * 
		sint[i__]);
/* L163: */
    }
L164:
    *pertrb = 0.f;
    if (ising != 0) {
	goto L165;
    } else {
	goto L176;
    }
L165:
    sum = wts * wps * f[its + jps * f_dim1] + wts * wpf * f[its + jpf * 
	    f_dim1] + wtf * wps * f[itf + jps * f_dim1] + wtf * wpf * f[itf + 
	    jpf * f_dim1];
    if (inp <= 0) {
	goto L167;
    } else {
	goto L166;
    }
L166:
    sum += wp * f[jps * f_dim1 + 1];
L167:
    if (isp <= 0) {
	goto L169;
    } else {
	goto L168;
    }
L168:
    sum += wp * f[*m + 1 + jps * f_dim1];
L169:
    i__1 = itfm;
    for (i__ = itsp; i__ <= i__1; ++i__) {
	sum1 = 0.f;
	i__2 = jpfm;
	for (j = jpsp; j <= i__2; ++j) {
	    sum1 += f[i__ + j * f_dim1];
/* L170: */
	}
	sum += sint[i__] * sum1;
/* L171: */
    }
    sum1 = 0.f;
    sum2 = 0.f;
    i__1 = jpfm;
    for (j = jpsp; j <= i__1; ++j) {
	sum1 += f[its + j * f_dim1];
	sum2 += f[itf + j * f_dim1];
/* L172: */
    }
    sum = sum + wts * sum1 + wtf * sum2;
    sum1 = 0.f;
    sum2 = 0.f;
    i__1 = itfm;
    for (i__ = itsp; i__ <= i__1; ++i__) {
	sum1 += sint[i__] * f[i__ + jps * f_dim1];
	sum2 += sint[i__] * f[i__ + jpf * f_dim1];
/* L173: */
    }
    sum = sum + wps * sum1 + wpf * sum2;
    *pertrb = sum / hne;
    i__1 = np1;
    for (j = 1; j <= i__1; ++j) {
	i__2 = mp1;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    f[i__ + j * f_dim1] -= *pertrb;
/* L174: */
	}
/* L175: */
    }

/* SCALE RIGHT SIDE FOR SUBROUTINE GENBUN */

L176:
    i__1 = itf;
    for (i__ = its; i__ <= i__1; ++i__) {
	cf = dphi2 * sint[i__] * sint[i__];
	i__2 = jpf;
	for (j = jps; j <= i__2; ++j) {
	    f[i__ + j * f_dim1] = cf * f[i__ + j * f_dim1];
/* L177: */
	}
/* L178: */
    }
    genbun_(nbdcnd, &nunk, &c__1, &munk, &am[its], &bm[its], &cm[its], idimf, 
	    &f[its + jps * f_dim1], &ierror, &d__[1]);
    if (ising <= 0) {
	goto L186;
    } else {
	goto L179;
    }
L179:
    if (inp <= 0) {
	goto L183;
    } else {
	goto L180;
    }
L180:
    if (isp <= 0) {
	goto L181;
    } else {
	goto L186;
    }
L181:
    i__1 = np1;
    for (j = 1; j <= i__1; ++j) {
	f[j * f_dim1 + 1] = 0.f;
/* L182: */
    }
    goto L209;
L183:
    if (isp <= 0) {
	goto L186;
    } else {
	goto L184;
    }
L184:
    i__1 = np1;
    for (j = 1; j <= i__1; ++j) {
	f[*m + 1 + j * f_dim1] = 0.f;
/* L185: */
    }
    goto L209;
L186:
    if (inp <= 0) {
	goto L193;
    } else {
	goto L187;
    }
L187:
    sum = wps * f[its + jps * f_dim1] + wpf * f[its + jpf * f_dim1];
    i__1 = jpfm;
    for (j = jpsp; j <= i__1; ++j) {
	sum += f[its + j * f_dim1];
/* L188: */
    }
    dfn = cp * sum;
    dnn = cp * ((wps + wpf + fjj) * (sn[2] - 1.f)) + *elmbda;
    dsn = cp * (wps + wpf + fjj) * sn[*m];
    if (isp <= 0) {
	goto L189;
    } else {
	goto L194;
    }
L189:
    cnp = (f[f_dim1 + 1] - dfn) / dnn;
    i__1 = itf;
    for (i__ = its; i__ <= i__1; ++i__) {
	hld = cnp * sn[i__];
	i__2 = jpf;
	for (j = jps; j <= i__2; ++j) {
	    f[i__ + j * f_dim1] += hld;
/* L190: */
	}
/* L191: */
    }
    i__1 = np1;
    for (j = 1; j <= i__1; ++j) {
	f[j * f_dim1 + 1] = cnp;
/* L192: */
    }
    goto L209;
L193:
    if (isp <= 0) {
	goto L209;
    } else {
	goto L194;
    }
L194:
    sum = wps * f[itf + jps * f_dim1] + wpf * f[itf + jpf * f_dim1];
    i__1 = jpfm;
    for (j = jpsp; j <= i__1; ++j) {
	sum += f[itf + j * f_dim1];
/* L195: */
    }
    dfs = cp * sum;
    dss = cp * ((wps + wpf + fjj) * (ss[*m] - 1.f)) + *elmbda;
    dns = cp * (wps + wpf + fjj) * ss[2];
    if (inp <= 0) {
	goto L196;
    } else {
	goto L200;
    }
L196:
    csp = (f[*m + 1 + f_dim1] - dfs) / dss;
    i__1 = itf;
    for (i__ = its; i__ <= i__1; ++i__) {
	hld = csp * ss[i__];
	i__2 = jpf;
	for (j = jps; j <= i__2; ++j) {
	    f[i__ + j * f_dim1] += hld;
/* L197: */
	}
/* L198: */
    }
    i__1 = np1;
    for (j = 1; j <= i__1; ++j) {
	f[*m + 1 + j * f_dim1] = csp;
/* L199: */
    }
    goto L209;
L200:
    rtn = f[f_dim1 + 1] - dfn;
    rts = f[*m + 1 + f_dim1] - dfs;
    if (ising <= 0) {
	goto L202;
    } else {
	goto L201;
    }
L201:
    csp = 0.f;
    cnp = rtn / dnn;
    goto L205;
L202:
    if (dabs(dnn) - dabs(dsn) <= 0.f) {
	goto L204;
    } else {
	goto L203;
    }
L203:
    den = dss - dns * dsn / dnn;
    rts -= rtn * dsn / dnn;
    csp = rts / den;
    cnp = (rtn - csp * dns) / dnn;
    goto L205;
L204:
    den = dns - dss * dnn / dsn;
    rtn -= rts * dnn / dsn;
    csp = rtn / den;
    cnp = (rts - dss * csp) / dsn;
L205:
    i__1 = itf;
    for (i__ = its; i__ <= i__1; ++i__) {
	hld = cnp * sn[i__] + csp * ss[i__];
	i__2 = jpf;
	for (j = jps; j <= i__2; ++j) {
	    f[i__ + j * f_dim1] += hld;
/* L206: */
	}
/* L207: */
    }
    i__1 = np1;
    for (j = 1; j <= i__1; ++j) {
	f[j * f_dim1 + 1] = cnp;
	f[*m + 1 + j * f_dim1] = csp;
/* L208: */
    }
L209:
    if (*nbdcnd != 0) {
	goto L212;
    } else {
	goto L210;
    }
L210:
    i__1 = mp1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	f[i__ + (jpf + 1) * f_dim1] = f[i__ + jps * f_dim1];
/* L211: */
    }
L212:
    return 0;
} /* hwsss1_ */

