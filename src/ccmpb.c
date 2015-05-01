/* ccmpb.f -- translated by f2c (version 12.02.01).
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
    integer npp, k;
    real eps, cnv;
    integer nm, ncmplx, ik;
} ccblk_;

#define ccblk_1 ccblk_

/* Table of constant values */

static integer c__4 = 4;
static integer c__2 = 2;

/* DECK CCMPB */
/* Subroutine */ int ccmpb_(integer *n, integer *ierror, real *an, real *bn, 
	real *cn, real *b, real *ah, real *bh)
{
    /* System generated locals */
    integer i__1, i__2, i__3, i__4;
    real r__1, r__2, r__3;

    /* Local variables */
    static integer i__, j, l;
    static real d1, d2, d3;
    static integer i2, j2, i4, l1, l2, j1, ib, if__, nb, jf, lh, ir, js, ls, 
	    ifd;
    static real arg;
    static integer kdo, n2m2, ipl, nmp;
    extern /* Subroutine */ int cpadd_(integer *, integer *, real *, real *, 
	    real *, real *, real *), inxcb_(integer *, integer *, integer *, 
	    integer *);
    static real bnorm;
    extern /* Subroutine */ int tevlc_(integer *, real *, real *, integer *);
    extern doublereal r1mach_(integer *);

/* ***BEGIN PROLOGUE  CCMPB */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to CBLKTR */
/* ***LIBRARY   SLATEC */
/* ***TYPE      COMPLEX (COMPB-S, CCMPB-C) */
/* ***AUTHOR  (UNKNOWN) */
/* ***DESCRIPTION */

/*     CCMPB computes the roots of the B polynomials using subroutine */
/*     TEVLC which is a modification the EISPACK program TQLRAT. */
/*     IERROR is set to 4 if either TEVLC fails or if A(J+1)*C(J) is */
/*     less than zero for some J.  AH,BH are temporary work arrays. */

/* ***SEE ALSO  CBLKTR */
/* ***ROUTINES CALLED  CPADD, INXCB, R1MACH, TEVLC */
/* ***COMMON BLOCKS    CCBLK */
/* ***REVISION HISTORY  (YYMMDD) */
/*   801001  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900402  Added TYPE section.  (WRB) */
/* ***END PROLOGUE  CCMPB */

/* ***FIRST EXECUTABLE STATEMENT  CCMPB */
    /* Parameter adjustments */
    --bh;
    --ah;
    --b;
    --cn;
    --bn;
    --an;

    /* Function Body */
    ccblk_1.eps = r1mach_(&c__4);
    bnorm = dabs(bn[1]);
    i__1 = ccblk_1.nm;
    for (j = 2; j <= i__1; ++j) {
/* Computing MAX */
	r__2 = bnorm, r__3 = (r__1 = bn[j], dabs(r__1));
	bnorm = dmax(r__2,r__3);
	arg = an[j] * cn[j - 1];
	if (arg >= 0.f) {
	    goto L101;
	} else {
	    goto L119;
	}
L101:
	r__1 = sqrt(arg);
	b[j] = r_sign(&r__1, &an[j]);
/* L102: */
    }
    ccblk_1.cnv = ccblk_1.eps * bnorm;
    if__ = pow_ii(&c__2, &ccblk_1.k);
    kdo = ccblk_1.k - 1;
    i__1 = kdo;
    for (l = 1; l <= i__1; ++l) {
	ir = l - 1;
	i2 = pow_ii(&c__2, &ir);
	i4 = i2 + i2;
	ipl = i4 - 1;
	ifd = if__ - i4;
	i__2 = ifd;
	i__3 = i4;
	for (i__ = i4; i__3 < 0 ? i__ >= i__2 : i__ <= i__2; i__ += i__3) {
	    inxcb_(&i__, &l, &ib, &nb);
	    if (nb <= 0) {
		goto L108;
	    } else {
		goto L103;
	    }
L103:
	    js = i__ - ipl;
	    jf = js + nb - 1;
	    ls = 0;
	    i__4 = jf;
	    for (j = js; j <= i__4; ++j) {
		++ls;
		bh[ls] = bn[j];
		ah[ls] = b[j];
/* L104: */
	    }
	    tevlc_(&nb, &bh[1], &ah[1], ierror);
	    if (*ierror != 0) {
		goto L118;
	    } else {
		goto L105;
	    }
L105:
	    lh = ib - 1;
	    i__4 = nb;
	    for (j = 1; j <= i__4; ++j) {
		++lh;
		b[lh] = -bh[j];
/* L106: */
	    }
/* L107: */
	}
L108:
	;
    }
    i__1 = ccblk_1.nm;
    for (j = 1; j <= i__1; ++j) {
	b[j] = -bn[j];
/* L109: */
    }
    if (ccblk_1.npp != 0) {
	goto L117;
    } else {
	goto L110;
    }
L110:
    nmp = ccblk_1.nm + 1;
    nb = ccblk_1.nm + nmp;
    i__1 = nb;
    for (j = 1; j <= i__1; ++j) {
	l1 = (j - 1) % nmp + 1;
	l2 = (j + ccblk_1.nm - 1) % nmp + 1;
	arg = an[l1] * cn[l2];
	if (arg >= 0.f) {
	    goto L111;
	} else {
	    goto L119;
	}
L111:
	r__1 = sqrt(arg);
	r__2 = -an[l1];
	bh[j] = r_sign(&r__1, &r__2);
	ah[j] = -bn[l1];
/* L112: */
    }
    tevlc_(&nb, &ah[1], &bh[1], ierror);
    if (*ierror != 0) {
	goto L118;
    } else {
	goto L113;
    }
L113:
    i__1 = ccblk_1.k - 1;
    inxcb_(&if__, &i__1, &j2, &lh);
    i__1 = if__ / 2;
    i__3 = ccblk_1.k - 1;
    inxcb_(&i__1, &i__3, &j1, &lh);
    ++j2;
    lh = j2;
    n2m2 = j2 + ccblk_1.nm + ccblk_1.nm - 2;
L114:
    d1 = (r__1 = b[j1] - b[j2 - 1], dabs(r__1));
    d2 = (r__1 = b[j1] - b[j2], dabs(r__1));
    d3 = (r__1 = b[j1] - b[j2 + 1], dabs(r__1));
    if (d2 < d1 && d2 < d3) {
	goto L115;
    }
    b[lh] = b[j2];
    ++j2;
    ++lh;
    if (j2 - n2m2 <= 0) {
	goto L114;
    } else {
	goto L116;
    }
L115:
    ++j2;
    ++j1;
    if (j2 - n2m2 <= 0) {
	goto L114;
    } else {
	goto L116;
    }
L116:
    b[lh] = b[n2m2 + 1];
    i__1 = ccblk_1.k - 1;
    inxcb_(&if__, &i__1, &j1, &j2);
    j2 = j1 + nmp + nmp;
    i__1 = ccblk_1.nm + 1;
    cpadd_(&i__1, ierror, &an[1], &cn[1], &b[j1], &b[j1], &b[j2]);
L117:
    return 0;
L118:
    *ierror = 4;
    return 0;
L119:
    *ierror = 5;
    return 0;
} /* ccmpb_ */

