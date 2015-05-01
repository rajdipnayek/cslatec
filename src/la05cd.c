/* la05cd.f -- translated by f2c (version 12.02.01).
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

/* Common Block Declarations */

struct {
    doublereal small;
    integer lp, lenl, lenu, ncp, lrow, lcol;
} la05dd_;

#define la05dd_1 la05dd_

/* Table of constant values */

static logical c_false = FALSE_;
static logical c_true = TRUE_;
static integer c__1 = 1;
static integer c_n6 = -6;
static integer c__2 = 2;
static integer c_n7 = -7;
static integer c_n8 = -8;

/* DECK LA05CD */
/* Subroutine */ int la05cd_(doublereal *a, integer *ind, integer *ia, 
	integer *n, integer *ip, integer *iw, doublereal *w, doublereal *g, 
	doublereal *u, integer *mm)
{
    /* System generated locals */
    address a__1[2];
    integer ind_dim1, ind_offset, iw_dim1, iw_offset, ip_dim1, ip_offset, 
	    i__1, i__2, i__3, i__4[2];
    doublereal d__1, d__2;
    char ch__1[62];

    /* Local variables */
    static integer i__, j, k, l, m, m1;
    static doublereal am;
    static integer ii, ij, kj;
    static doublereal au;
    static integer jm, im, kp, kl, kr, km, is, in, ir, jp, kq, ks, kk, nz, 
	    mcp, kpl, krl, ins, jns, knp, ipp, last, last1, last2;
    static char xern1[8];
    extern /* Subroutine */ int la05ed_(doublereal *, integer *, integer *, 
	    integer *, integer *, integer *, logical *), xermsg_(char *, char 
	    *, char *, integer *, integer *, ftnlen, ftnlen, ftnlen), xsetun_(
	    integer *);

    /* Fortran I/O blocks */
    static icilist io___37 = { 0, xern1, 0, "(I8)", 8, 1 };


/* ***BEGIN PROLOGUE  LA05CD */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to DSPLP */
/* ***LIBRARY   SLATEC */
/* ***TYPE      DOUBLE PRECISION (LA05CS-D, LA05CD-D) */
/* ***AUTHOR  (UNKNOWN) */
/* ***DESCRIPTION */

/*     THIS SUBPROGRAM IS A SLIGHT MODIFICATION OF A SUBPROGRAM */
/*     FROM THE C. 1979 AERE HARWELL LIBRARY.  THE NAME OF THE */
/*     CORRESPONDING HARWELL CODE CAN BE OBTAINED BY DELETING */
/*     THE FINAL LETTER =D= IN THE NAMES USED HERE. */
/*     REVISED SEP. 13, 1979. */

/*     ROYALTIES HAVE BEEN PAID TO AERE-UK FOR USE OF THEIR CODES */
/*     IN THE PACKAGE GIVEN HERE.  ANY PRIMARY USAGE OF THE HARWELL */
/*     SUBROUTINES REQUIRES A ROYALTY AGREEMENT AND PAYMENT BETWEEN */
/*     THE USER AND AERE-UK.  ANY USAGE OF THE SANDIA WRITTEN CODES */
/*     DSPLP( ) (WHICH USES THE HARWELL SUBROUTINES) IS PERMITTED. */

/* ***SEE ALSO  DSPLP */
/* ***ROUTINES CALLED  LA05ED, XERMSG, XSETUN */
/* ***COMMON BLOCKS    LA05DD */
/* ***REVISION HISTORY  (YYMMDD) */
/*   811215  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890831  Modified array declarations.  (WRB) */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ) */
/*   900402  Added TYPE section.  (WRB) */
/*   900510  Convert XERRWV calls to XERMSG calls.  (RWC) */
/*   920410  Corrected second dimension on IW declaration.  (WRB) */
/*   920422  Changed upper limit on DO from LAST to LAST-1.  (WRB) */
/* ***END PROLOGUE  LA05CD */

/* ***FIRST EXECUTABLE STATEMENT  LA05CD */
    /* Parameter adjustments */
    --a;
    ind_dim1 = *ia;
    ind_offset = 1 + ind_dim1;
    ind -= ind_offset;
    iw_dim1 = *n;
    iw_offset = 1 + iw_dim1;
    iw -= iw_offset;
    ip_dim1 = *n;
    ip_offset = 1 + ip_dim1;
    ip -= ip_offset;
    --w;

    /* Function Body */
    xsetun_(&la05dd_1.lp);
    if (*g < 0.) {
	goto L620;
    }
    jm = *mm;
/* MCP LIMITS THE VALUE OF NCP PERMITTED BEFORE AN ERROR RETURN RESULTS. */
    mcp = la05dd_1.ncp + 20;
/* REMOVE OLD COLUMN */
    la05dd_1.lenu -= iw[jm + (iw_dim1 << 1)];
    kp = ip[jm + (ip_dim1 << 1)];
    im = ind[kp + ind_dim1];
    kl = kp + iw[jm + (iw_dim1 << 1)] - 1;
    iw[jm + (iw_dim1 << 1)] = 0;
    i__1 = kl;
    for (k = kp; k <= i__1; ++k) {
	i__ = ind[k + ind_dim1];
	ind[k + ind_dim1] = 0;
	kr = ip[i__ + ip_dim1];
	nz = iw[i__ + iw_dim1] - 1;
	iw[i__ + iw_dim1] = nz;
	krl = kr + nz;
	i__2 = krl;
	for (km = kr; km <= i__2; ++km) {
	    if (ind[km + (ind_dim1 << 1)] == jm) {
		goto L20;
	    }
/* L10: */
	}
L20:
	a[km] = a[krl];
	ind[km + (ind_dim1 << 1)] = ind[krl + (ind_dim1 << 1)];
	ind[krl + (ind_dim1 << 1)] = 0;
/* L30: */
    }

/* INSERT NEW COLUMN */
    i__1 = *n;
    for (ii = 1; ii <= i__1; ++ii) {
	i__ = iw[ii + iw_dim1 * 3];
	if (i__ == im) {
	    m = ii;
	}
	if ((d__1 = w[i__], abs(d__1)) <= la05dd_1.small) {
	    goto L100;
	}
	++la05dd_1.lenu;
	last = ii;
	if (la05dd_1.lcol + la05dd_1.lenl < *ia) {
	    goto L40;
	}
/* COMPRESS COLUMN FILE IF NECESSARY. */
	if (la05dd_1.ncp >= mcp || la05dd_1.lenl + la05dd_1.lenu >= *ia) {
	    goto L610;
	}
	la05ed_(&a[1], &ind[ind_offset], &ip[(ip_dim1 << 1) + 1], n, &iw[(
		iw_dim1 << 1) + 1], ia, &c_false);
L40:
	++la05dd_1.lcol;
	nz = iw[jm + (iw_dim1 << 1)];
	if (nz == 0) {
	    ip[jm + (ip_dim1 << 1)] = la05dd_1.lcol;
	}
	iw[jm + (iw_dim1 << 1)] = nz + 1;
	ind[la05dd_1.lcol + ind_dim1] = i__;
	nz = iw[i__ + iw_dim1];
	kpl = ip[i__ + ip_dim1] + nz;
	if (kpl > la05dd_1.lrow) {
	    goto L50;
	}
	if (ind[kpl + (ind_dim1 << 1)] == 0) {
	    goto L90;
	}
/* NEW ENTRY HAS TO BE CREATED. */
L50:
	if (la05dd_1.lenl + la05dd_1.lrow + nz < *ia) {
	    goto L60;
	}
	if (la05dd_1.ncp >= mcp || la05dd_1.lenl + la05dd_1.lenu + nz >= *ia) 
		{
	    goto L610;
	}
/* COMPRESS ROW FILE IF NECESSARY. */
	la05ed_(&a[1], &ind[(ind_dim1 << 1) + 1], &ip[ip_offset], n, &iw[
		iw_offset], ia, &c_true);
L60:
	kp = ip[i__ + ip_dim1];
	ip[i__ + ip_dim1] = la05dd_1.lrow + 1;
	if (nz == 0) {
	    goto L80;
	}
	kpl = kp + nz - 1;
	i__2 = kpl;
	for (k = kp; k <= i__2; ++k) {
	    ++la05dd_1.lrow;
	    a[la05dd_1.lrow] = a[k];
	    ind[la05dd_1.lrow + (ind_dim1 << 1)] = ind[k + (ind_dim1 << 1)];
	    ind[k + (ind_dim1 << 1)] = 0;
/* L70: */
	}
L80:
	++la05dd_1.lrow;
	kpl = la05dd_1.lrow;
/* PLACE NEW ELEMENT AT END OF ROW. */
L90:
	iw[i__ + iw_dim1] = nz + 1;
	a[kpl] = w[i__];
	ind[kpl + (ind_dim1 << 1)] = jm;
L100:
	w[i__] = 0.;
/* L110: */
    }
    if (iw[im + iw_dim1] == 0 || iw[jm + (iw_dim1 << 1)] == 0 || m > last) {
	goto L590;
    }

/* FIND COLUMN SINGLETONS, OTHER THAN THE SPIKE. NON-SINGLETONS ARE */
/*     MARKED WITH W(J)=1. ONLY IW(.,3) IS REVISED AND IW(.,4) IS USED */
/*     FOR WORKSPACE. */
    ins = m;
    m1 = m;
    w[jm] = 1.;
    i__1 = last;
    for (ii = m; ii <= i__1; ++ii) {
	i__ = iw[ii + iw_dim1 * 3];
	j = iw[ii + (iw_dim1 << 2)];
	if (w[j] == 0.f) {
	    goto L130;
	}
	kp = ip[i__ + ip_dim1];
	kl = kp + iw[i__ + iw_dim1] - 1;
	i__2 = kl;
	for (k = kp; k <= i__2; ++k) {
	    j = ind[k + (ind_dim1 << 1)];
	    w[j] = 1.;
/* L120: */
	}
	iw[ins + (iw_dim1 << 2)] = i__;
	++ins;
	goto L140;
/* PLACE SINGLETONS IN NEW POSITION. */
L130:
	iw[m1 + iw_dim1 * 3] = i__;
	++m1;
L140:
	;
    }
/* PLACE NON-SINGLETONS IN NEW POSITION. */
    ij = m + 1;
    i__1 = last - 1;
    for (ii = m1; ii <= i__1; ++ii) {
	iw[ii + iw_dim1 * 3] = iw[ij + (iw_dim1 << 2)];
	++ij;
/* L150: */
    }
/* PLACE SPIKE AT END. */
    iw[last + iw_dim1 * 3] = im;

/* FIND ROW SINGLETONS, APART FROM SPIKE ROW. NON-SINGLETONS ARE MARKED */
/*     WITH W(I)=2. AGAIN ONLY IW(.,3) IS REVISED AND IW(.,4) IS USED */
/*     FOR WORKSPACE. */
    last1 = last;
    jns = last;
    w[im] = 2.;
    j = jm;
    i__1 = last;
    for (ij = m1; ij <= i__1; ++ij) {
	ii = last + m1 - ij;
	i__ = iw[ii + iw_dim1 * 3];
	if (w[i__] != 2.) {
	    goto L170;
	}
	k = ip[i__ + ip_dim1];
	if (ii != last) {
	    j = ind[k + (ind_dim1 << 1)];
	}
	kp = ip[j + (ip_dim1 << 1)];
	kl = kp + iw[j + (iw_dim1 << 1)] - 1;
	iw[jns + (iw_dim1 << 2)] = i__;
	--jns;
	i__2 = kl;
	for (k = kp; k <= i__2; ++k) {
	    i__ = ind[k + ind_dim1];
	    w[i__] = 2.;
/* L160: */
	}
	goto L180;
L170:
	iw[last1 + iw_dim1 * 3] = i__;
	--last1;
L180:
	;
    }
    i__1 = last1;
    for (ii = m1; ii <= i__1; ++ii) {
	++jns;
	i__ = iw[jns + (iw_dim1 << 2)];
	w[i__] = 3.;
	iw[ii + iw_dim1 * 3] = i__;
/* L190: */
    }

/* DEAL WITH SINGLETON SPIKE COLUMN. NOTE THAT BUMP ROWS ARE MARKED BY */
/*    W(I)=3. */
    i__1 = last1;
    for (ii = m1; ii <= i__1; ++ii) {
	kp = ip[jm + (ip_dim1 << 1)];
	kl = kp + iw[jm + (iw_dim1 << 1)] - 1;
	is = 0;
	i__2 = kl;
	for (k = kp; k <= i__2; ++k) {
	    l = ind[k + ind_dim1];
	    if (w[l] != 3.) {
		goto L200;
	    }
	    if (is != 0) {
		goto L240;
	    }
	    i__ = l;
	    knp = k;
	    is = 1;
L200:
	    ;
	}
	if (is == 0) {
	    goto L590;
	}
/* MAKE A(I,JM) A PIVOT. */
	ind[knp + ind_dim1] = ind[kp + ind_dim1];
	ind[kp + ind_dim1] = i__;
	kp = ip[i__ + ip_dim1];
	i__2 = *ia;
	for (k = kp; k <= i__2; ++k) {
	    if (ind[k + (ind_dim1 << 1)] == jm) {
		goto L220;
	    }
/* L210: */
	}
L220:
	am = a[kp];
	a[kp] = a[k];
	a[k] = am;
	ind[k + (ind_dim1 << 1)] = ind[kp + (ind_dim1 << 1)];
	ind[kp + (ind_dim1 << 1)] = jm;
	jm = ind[k + (ind_dim1 << 1)];
	iw[ii + (iw_dim1 << 2)] = i__;
	w[i__] = 2.;
/* L230: */
    }
    ii = last1;
    goto L260;
L240:
    in = m1;
    i__1 = last1;
    for (ij = ii; ij <= i__1; ++ij) {
	iw[ij + (iw_dim1 << 2)] = iw[in + iw_dim1 * 3];
	++in;
/* L250: */
    }
L260:
    last2 = last1 - 1;
    if (m1 == last1) {
	goto L570;
    }
    i__1 = last2;
    for (i__ = m1; i__ <= i__1; ++i__) {
	iw[i__ + iw_dim1 * 3] = iw[i__ + (iw_dim1 << 2)];
/* L270: */
    }
    m1 = ii;
    if (m1 == last1) {
	goto L570;
    }

/* CLEAR W */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	w[i__] = 0.;
/* L280: */
    }

/* PERFORM ELIMINATION */
    ir = iw[last1 + iw_dim1 * 3];
    i__1 = last1;
    for (ii = m1; ii <= i__1; ++ii) {
	ipp = iw[ii + iw_dim1 * 3];
	kp = ip[ipp + ip_dim1];
	kr = ip[ir + ip_dim1];
	jp = ind[kp + (ind_dim1 << 1)];
	if (ii == last1) {
	    jp = jm;
	}
/* SEARCH NON-PIVOT ROW FOR ELEMENT TO BE ELIMINATED. */
/*  AND BRING IT TO FRONT OF ITS ROW */
	krl = kr + iw[ir + iw_dim1] - 1;
	i__2 = krl;
	for (knp = kr; knp <= i__2; ++knp) {
	    if (jp == ind[knp + (ind_dim1 << 1)]) {
		goto L300;
	    }
/* L290: */
	}
	if (ii - last1 != 0) {
	    goto L560;
	} else {
	    goto L590;
	}
/* BRING ELEMENT TO BE ELIMINATED TO FRONT OF ITS ROW. */
L300:
	am = a[knp];
	a[knp] = a[kr];
	a[kr] = am;
	ind[knp + (ind_dim1 << 1)] = ind[kr + (ind_dim1 << 1)];
	ind[kr + (ind_dim1 << 1)] = jp;
	if (ii == last1) {
	    goto L310;
	}
	if ((d__1 = a[kp], abs(d__1)) < *u * abs(am)) {
	    goto L310;
	}
	if (abs(am) < *u * (d__1 = a[kp], abs(d__1))) {
	    goto L340;
	}
	if (iw[ipp + iw_dim1] <= iw[ir + iw_dim1]) {
	    goto L340;
	}
/* PERFORM INTERCHANGE */
L310:
	iw[last1 + iw_dim1 * 3] = ipp;
	iw[ii + iw_dim1 * 3] = ir;
	ir = ipp;
	ipp = iw[ii + iw_dim1 * 3];
	k = kr;
	kr = kp;
	kp = k;
	kj = ip[jp + (ip_dim1 << 1)];
	i__2 = *ia;
	for (k = kj; k <= i__2; ++k) {
	    if (ind[k + ind_dim1] == ipp) {
		goto L330;
	    }
/* L320: */
	}
L330:
	ind[k + ind_dim1] = ind[kj + ind_dim1];
	ind[kj + ind_dim1] = ipp;
L340:
	if (a[kp] == 0.) {
	    goto L590;
	}
	if (ii == last1) {
	    goto L560;
	}
	am = -a[kr] / a[kp];
/* COMPRESS ROW FILE UNLESS IT IS CERTAIN THAT THERE IS ROOM FOR NEW ROW. */
	if (la05dd_1.lrow + iw[ir + iw_dim1] + iw[ipp + iw_dim1] + 
		la05dd_1.lenl <= *ia) {
	    goto L350;
	}
	if (la05dd_1.ncp >= mcp || la05dd_1.lenu + iw[ir + iw_dim1] + iw[ipp 
		+ iw_dim1] + la05dd_1.lenl > *ia) {
	    goto L610;
	}
	la05ed_(&a[1], &ind[(ind_dim1 << 1) + 1], &ip[ip_offset], n, &iw[
		iw_offset], ia, &c_true);
	kp = ip[ipp + ip_dim1];
	kr = ip[ir + ip_dim1];
L350:
	krl = kr + iw[ir + iw_dim1] - 1;
	kq = kp + 1;
	kpl = kp + iw[ipp + iw_dim1] - 1;
/* PLACE PIVOT ROW (EXCLUDING PIVOT ITSELF) IN W. */
	if (kq > kpl) {
	    goto L370;
	}
	i__2 = kpl;
	for (k = kq; k <= i__2; ++k) {
	    j = ind[k + (ind_dim1 << 1)];
	    w[j] = a[k];
/* L360: */
	}
L370:
	ip[ir + ip_dim1] = la05dd_1.lrow + 1;

/* TRANSFER MODIFIED ELEMENTS. */
	ind[kr + (ind_dim1 << 1)] = 0;
	++kr;
	if (kr > krl) {
	    goto L430;
	}
	i__2 = krl;
	for (ks = kr; ks <= i__2; ++ks) {
	    j = ind[ks + (ind_dim1 << 1)];
	    au = a[ks] + am * w[j];
	    ind[ks + (ind_dim1 << 1)] = 0;
/* IF ELEMENT IS VERY SMALL REMOVE IT FROM U. */
	    if (abs(au) <= la05dd_1.small) {
		goto L380;
	    }
/* Computing MAX */
	    d__1 = *g, d__2 = abs(au);
	    *g = max(d__1,d__2);
	    ++la05dd_1.lrow;
	    a[la05dd_1.lrow] = au;
	    ind[la05dd_1.lrow + (ind_dim1 << 1)] = j;
	    goto L410;
L380:
	    --la05dd_1.lenu;
/* REMOVE ELEMENT FROM COL FILE. */
	    k = ip[j + (ip_dim1 << 1)];
	    kl = k + iw[j + (iw_dim1 << 1)] - 1;
	    iw[j + (iw_dim1 << 1)] = kl - k;
	    i__3 = kl;
	    for (kk = k; kk <= i__3; ++kk) {
		if (ind[kk + ind_dim1] == ir) {
		    goto L400;
		}
/* L390: */
	    }
L400:
	    ind[kk + ind_dim1] = ind[kl + ind_dim1];
	    ind[kl + ind_dim1] = 0;
L410:
	    w[j] = 0.;
/* L420: */
	}

/* SCAN PIVOT ROW FOR FILLS. */
L430:
	if (kq > kpl) {
	    goto L520;
	}
	i__2 = kpl;
	for (ks = kq; ks <= i__2; ++ks) {
	    j = ind[ks + (ind_dim1 << 1)];
	    au = am * w[j];
	    if (abs(au) <= la05dd_1.small) {
		goto L500;
	    }
	    ++la05dd_1.lrow;
	    a[la05dd_1.lrow] = au;
	    ind[la05dd_1.lrow + (ind_dim1 << 1)] = j;
	    ++la05dd_1.lenu;

/* CREATE FILL IN COLUMN FILE. */
	    nz = iw[j + (iw_dim1 << 1)];
	    k = ip[j + (ip_dim1 << 1)];
	    kl = k + nz - 1;
/* IF POSSIBLE PLACE NEW ELEMENT AT END OF PRESENT ENTRY. */
	    if (kl != la05dd_1.lcol) {
		goto L440;
	    }
	    if (la05dd_1.lcol + la05dd_1.lenl >= *ia) {
		goto L460;
	    }
	    ++la05dd_1.lcol;
	    goto L450;
L440:
	    if (ind[kl + 1 + ind_dim1] != 0) {
		goto L460;
	    }
L450:
	    ind[kl + 1 + ind_dim1] = ir;
	    goto L490;
/* NEW ENTRY HAS TO BE CREATED. */
L460:
	    if (la05dd_1.lcol + la05dd_1.lenl + nz + 1 < *ia) {
		goto L470;
	    }
/* COMPRESS COLUMN FILE IF THERE IS NOT ROOM FOR NEW ENTRY. */
	    if (la05dd_1.ncp >= mcp || la05dd_1.lenu + la05dd_1.lenl + nz + 1 
		    >= *ia) {
		goto L610;
	    }
	    la05ed_(&a[1], &ind[ind_offset], &ip[(ip_dim1 << 1) + 1], n, &iw[(
		    iw_dim1 << 1) + 1], ia, &c_false);
	    k = ip[j + (ip_dim1 << 1)];
	    kl = k + nz - 1;
/* TRANSFER OLD ENTRY INTO NEW. */
L470:
	    ip[j + (ip_dim1 << 1)] = la05dd_1.lcol + 1;
	    i__3 = kl;
	    for (kk = k; kk <= i__3; ++kk) {
		++la05dd_1.lcol;
		ind[la05dd_1.lcol + ind_dim1] = ind[kk + ind_dim1];
		ind[kk + ind_dim1] = 0;
/* L480: */
	    }
/* ADD NEW ELEMENT. */
	    ++la05dd_1.lcol;
	    ind[la05dd_1.lcol + ind_dim1] = ir;
L490:
/* Computing MAX */
	    d__1 = *g, d__2 = abs(au);
	    *g = max(d__1,d__2);
	    iw[j + (iw_dim1 << 1)] = nz + 1;
L500:
	    w[j] = 0.;
/* L510: */
	}
L520:
	iw[ir + iw_dim1] = la05dd_1.lrow + 1 - ip[ir + ip_dim1];

/* STORE MULTIPLIER */
	if (la05dd_1.lenl + la05dd_1.lcol + 1 <= *ia) {
	    goto L530;
	}
/* COMPRESS COL FILE IF NECESSARY. */
	if (la05dd_1.ncp >= mcp) {
	    goto L610;
	}
	la05ed_(&a[1], &ind[ind_offset], &ip[(ip_dim1 << 1) + 1], n, &iw[(
		iw_dim1 << 1) + 1], ia, &c_false);
L530:
	k = *ia - la05dd_1.lenl;
	++la05dd_1.lenl;
	a[k] = am;
	ind[k + ind_dim1] = ipp;
	ind[k + (ind_dim1 << 1)] = ir;
/* CREATE BLANK IN PIVOTAL COLUMN. */
	kp = ip[jp + (ip_dim1 << 1)];
	nz = iw[jp + (iw_dim1 << 1)] - 1;
	kl = kp + nz;
	i__2 = kl;
	for (k = kp; k <= i__2; ++k) {
	    if (ind[k + ind_dim1] == ir) {
		goto L550;
	    }
/* L540: */
	}
L550:
	ind[k + ind_dim1] = ind[kl + ind_dim1];
	iw[jp + (iw_dim1 << 1)] = nz;
	ind[kl + ind_dim1] = 0;
	--la05dd_1.lenu;
L560:
	;
    }

/* CONSTRUCT COLUMN PERMUTATION AND STORE IT IN IW(.,4) */
L570:
    i__1 = last;
    for (ii = m; ii <= i__1; ++ii) {
	i__ = iw[ii + iw_dim1 * 3];
	k = ip[i__ + ip_dim1];
	j = ind[k + (ind_dim1 << 1)];
	iw[ii + (iw_dim1 << 2)] = j;
/* L580: */
    }
    return 0;

/*     THE FOLLOWING INSTRUCTIONS IMPLEMENT THE FAILURE EXITS. */

L590:
    if (la05dd_1.lp > 0) {
	s_wsfi(&io___37);
	do_fio(&c__1, (char *)&(*mm), (ftnlen)sizeof(integer));
	e_wsfi();
/* Writing concatenation */
	i__4[0] = 54, a__1[0] = "SINGULAR MATRIX AFTER REPLACEMENT OF COLUMN"
		".  INDEX = ";
	i__4[1] = 8, a__1[1] = xern1;
	s_cat(ch__1, a__1, i__4, &c__2, (ftnlen)62);
	xermsg_("SLATEC", "LA05CD", ch__1, &c_n6, &c__1, (ftnlen)6, (ftnlen)6,
		 (ftnlen)62);
    }
    *g = -6.;
    return 0;

L610:
    if (la05dd_1.lp > 0) {
	xermsg_("SLATEC", "LA05CD", "LENGTHS OF ARRAYS A(*) AND IND(*,2) ARE"
		" TOO SMALL.", &c_n7, &c__1, (ftnlen)6, (ftnlen)6, (ftnlen)50);
    }
    *g = -7.;
    return 0;

L620:
    if (la05dd_1.lp > 0) {
	xermsg_("SLATEC", "LA05CD", "EARLIER ENTRY GAVE ERROR RETURN.", &c_n8,
		 &c__2, (ftnlen)6, (ftnlen)6, (ftnlen)32);
    }
    *g = -8.;
    return 0;
} /* la05cd_ */

