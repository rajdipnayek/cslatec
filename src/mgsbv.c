/* mgsbv.f -- translated by f2c (version 12.02.01).
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
    real ae, re, tol;
    integer nxpts, nic, nopg, mxnon, ndisk, ntape, neq, indpvt, integ, nps, 
	    ntp, neqivp, numort, nfcc, icoco;
} ml18jr_;

#define ml18jr_1 ml18jr_

struct {
    real uro, sru, eps, sqovfl, twou, fouru;
    integer lpar;
} ml5mco_;

#define ml5mco_1 ml5mco_

/* Table of constant values */

static integer c__1 = 1;

/* DECK MGSBV */
/* Subroutine */ int mgsbv_(integer *m, integer *n, real *a, integer *ia, 
	integer *niv, integer *iflag, real *s, real *p, integer *ip, integer *
	inhomo, real *v, real *w, real *wcnd)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3;
    real r__1, r__2;

    /* Local variables */
    static integer i__, j, k, l;
    static real t, y;
    static integer m2, kd, kj, jk, jp, kp, nn, jq, lr, nr, ix;
    static real vl;
    static integer jy, jz, iz;
    static real sv, ry;
    static integer ip1, np1;
    static real dot, pjp;
    static integer lix, nrm1;
    extern doublereal sdot_(integer *, real *, integer *, real *, integer *);
    static integer nivn, nmnr;
    static real psave;
    extern doublereal prvec_(integer *, real *, real *);
    static real vnorm;

/* ***BEGIN PROLOGUE  MGSBV */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to BVSUP */
/* ***LIBRARY   SLATEC */
/* ***TYPE      SINGLE PRECISION (MGSBV-S, DMGSBV-D) */
/* ***AUTHOR  Watts, H. A., (SNLA) */
/* ***DESCRIPTION */

/* ********************************************************************** */
/* Orthogonalize a set of N real vectors and determine their rank */

/* ********************************************************************** */
/* INPUT */
/* ********************************************************************** */
/*   M = Dimension of vectors */
/*   N = No. of vectors */
/*   A = Array whose first N cols contain the vectors */
/*   IA = First dimension of array A (col length) */
/*   NIV = Number of independent vectors needed */
/*   INHOMO = 1 Corresponds to having a non-zero particular solution */
/*   V = Particular solution vector (not included in the pivoting) */
/*   INDPVT = 1 Means pivoting will not be used */

/* ********************************************************************** */
/* OUTPUT */
/* ********************************************************************** */
/*   NIV = No. of linear independent vectors in input set */
/*     A = Matrix whose first NIV cols. contain NIV orthogonal vectors */
/*         which span the vector space determined by the input vectors */
/*   IFLAG */
/*          = 0 success */
/*          = 1 incorrect input */
/*          = 2 rank of new vectors less than N */
/*   P = Decomposition matrix.  P is upper triangular and */
/*             (old vectors) = (new vectors) * P. */
/*         The old vectors will be reordered due to pivoting */
/*         The dimension of p must be .GE. N*(N+1)/2. */
/*             (  N*(2*N+1) when N .NE. NFCC ) */
/*   IP = Pivoting vector. The dimension of IP must be .GE. N. */
/*             (  2*N when N .NE. NFCC ) */
/*   S = Square of norms of incoming vectors */
/*   V = Vector which is orthogonal to the vectors of A */
/*   W = Orthogonalization information for the vector V */
/*   WCND = Worst case (smallest) norm decrement value of the */
/*          vectors being orthogonalized  (represents a test */
/*          for linear dependence of the vectors) */
/* ********************************************************************** */

/* ***SEE ALSO  BVSUP */
/* ***ROUTINES CALLED  PRVEC, SDOT */
/* ***COMMON BLOCKS    ML18JR, ML5MCO */
/* ***REVISION HISTORY  (YYMMDD) */
/*   750601  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890831  Modified array declarations.  (WRB) */
/*   890921  Realigned order of variables in certain COMMON blocks. */
/*           (WRB) */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900328  Added TYPE section.  (WRB) */
/*   910722  Updated AUTHOR section.  (ALS) */
/* ***END PROLOGUE  MGSBV */





/* ***FIRST EXECUTABLE STATEMENT  MGSBV */
    /* Parameter adjustments */
    a_dim1 = *ia;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --s;
    --p;
    --ip;
    --v;
    --w;

    /* Function Body */
    if (*m > 0 && *n > 0 && *ia >= *m) {
	goto L10;
    }
    *iflag = 1;
    return 0;

L10:
    jp = 0;
    *iflag = 0;
    np1 = *n + 1;
    y = 0.f;
    m2 = *m / 2;

/*     CALCULATE SQUARE OF NORMS OF INCOMING VECTORS AND SEARCH FOR */
/*     VECTOR WITH LARGEST MAGNITUDE */

    j = 0;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	vl = sdot_(m, &a[i__ * a_dim1 + 1], &c__1, &a[i__ * a_dim1 + 1], &
		c__1);
	s[i__] = vl;
	if (*n == ml18jr_1.nfcc) {
	    goto L25;
	}
	j = (i__ << 1) - 1;
	p[j] = vl;
	ip[j] = j;
L25:
	++j;
	p[j] = vl;
	ip[j] = j;
	if (vl <= y) {
	    goto L30;
	}
	y = vl;
	ix = i__;
L30:
	;
    }
    if (ml18jr_1.indpvt != 1) {
	goto L33;
    }
    ix = 1;
    y = p[1];
L33:
    lix = ix;
    if (*n != ml18jr_1.nfcc) {
	lix = (ix << 1) - 1;
    }
    p[lix] = p[1];
    s[np1] = 0.f;
    if (*inhomo == 1) {
	s[np1] = sdot_(m, &v[1], &c__1, &v[1], &c__1);
    }
    *wcnd = 1.f;
    nivn = *niv;
    *niv = 0;

    if (y == 0.f) {
	goto L170;
    }
/* ********************************************************************** */
    i__1 = *n;
    for (nr = 1; nr <= i__1; ++nr) {
	if (nivn == *niv) {
	    goto L150;
	}
	*niv = nr;
	if (ix == nr) {
	    goto L80;
	}

/*     PIVOTING OF COLUMNS OF P MATRIX */

	nn = *n;
	lix = ix;
	lr = nr;
	if (*n == ml18jr_1.nfcc) {
	    goto L40;
	}
	nn = ml18jr_1.nfcc;
	lix = (ix << 1) - 1;
	lr = (nr << 1) - 1;
L40:
	if (nr == 1) {
	    goto L60;
	}
	kd = lix - lr;
	kj = lr;
	nrm1 = lr - 1;
	i__2 = nrm1;
	for (j = 1; j <= i__2; ++j) {
	    psave = p[kj];
	    jk = kj + kd;
	    p[kj] = p[jk];
	    p[jk] = psave;
/* L50: */
	    kj = kj + nn - j;
	}
	jy = jk + nmnr;
	jz = jy - kd;
	p[jy] = p[jz];
L60:
	iz = ip[lix];
	ip[lix] = ip[lr];
	ip[lr] = iz;
	sv = s[ix];
	s[ix] = s[nr];
	s[nr] = sv;
	if (*n == ml18jr_1.nfcc) {
	    goto L69;
	}
	if (nr == 1) {
	    goto L67;
	}
	kj = lr + 1;
	i__2 = nrm1;
	for (k = 1; k <= i__2; ++k) {
	    psave = p[kj];
	    jk = kj + kd;
	    p[kj] = p[jk];
	    p[jk] = psave;
/* L65: */
	    kj = kj + ml18jr_1.nfcc - k;
	}
L67:
	iz = ip[lix + 1];
	ip[lix + 1] = ip[lr + 1];
	ip[lr + 1] = iz;

/*     PIVOTING OF COLUMNS OF VECTORS */

L69:
	i__2 = *m;
	for (l = 1; l <= i__2; ++l) {
	    t = a[l + ix * a_dim1];
	    a[l + ix * a_dim1] = a[l + nr * a_dim1];
/* L70: */
	    a[l + nr * a_dim1] = t;
	}

/*     CALCULATE P(NR,NR) AS NORM SQUARED OF PIVOTAL VECTOR */

L80:
	++jp;
	p[jp] = y;
	ry = 1.f / y;
	nmnr = *n - nr;
	if (*n == ml18jr_1.nfcc) {
	    goto L85;
	}
	nmnr = ml18jr_1.nfcc - ((nr << 1) - 1);
	++jp;
	p[jp] = 0.f;
	kp = jp + nmnr;
	p[kp] = y;
L85:
	if (nr == *n || nivn == *niv) {
	    goto L125;
	}

/*    CALCULATE ORTHOGONAL PROJECTION VECTORS AND SEARCH FOR LARGEST NORM */

	y = 0.f;
	ip1 = nr + 1;
	ix = ip1;
/*     **************************************** */
	i__2 = *n;
	for (j = ip1; j <= i__2; ++j) {
	    dot = sdot_(m, &a[nr * a_dim1 + 1], &c__1, &a[j * a_dim1 + 1], &
		    c__1);
	    ++jp;
	    jq = jp + nmnr;
	    if (*n != ml18jr_1.nfcc) {
		jq = jq + nmnr - 1;
	    }
	    p[jq] = p[jp] - dot * (dot * ry);
	    p[jp] = dot * ry;
	    i__3 = *m;
	    for (i__ = 1; i__ <= i__3; ++i__) {
/* L90: */
		a[i__ + j * a_dim1] -= p[jp] * a[i__ + nr * a_dim1];
	    }
	    if (*n == ml18jr_1.nfcc) {
		goto L99;
	    }
	    kp = jp + nmnr;
	    ++jp;
	    pjp = ry * prvec_(m, &a[nr * a_dim1 + 1], &a[j * a_dim1 + 1]);
	    p[jp] = pjp;
	    p[kp] = -pjp;
	    ++kp;
	    p[kp] = ry * dot;
	    i__3 = m2;
	    for (k = 1; k <= i__3; ++k) {
		l = m2 + k;
		a[k + j * a_dim1] -= pjp * a[l + nr * a_dim1];
/* L95: */
		a[l + j * a_dim1] += pjp * a[k + nr * a_dim1];
	    }
	    p[jq] -= pjp * (pjp / ry);

/*     TEST FOR CANCELLATION IN RECURRENCE RELATION */

L99:
	    if (p[jq] > s[j] * ml5mco_1.sru) {
		goto L100;
	    }
	    p[jq] = sdot_(m, &a[j * a_dim1 + 1], &c__1, &a[j * a_dim1 + 1], &
		    c__1);
L100:
	    if (p[jq] <= y) {
		goto L120;
	    }
	    y = p[jq];
	    ix = j;
L120:
	    ;
	}
	if (*n != ml18jr_1.nfcc) {
	    jp = kp;
	}
/*     **************************************** */
	if (ml18jr_1.indpvt == 1) {
	    ix = ip1;
	}

/*     RECOMPUTE NORM SQUARED OF PIVOTAL VECTOR WITH SCALAR PRODUCT */

	y = sdot_(m, &a[ix * a_dim1 + 1], &c__1, &a[ix * a_dim1 + 1], &c__1);
	if (y <= ml5mco_1.eps * s[ix]) {
	    goto L170;
	}
/* Computing MIN */
	r__1 = *wcnd, r__2 = y / s[ix];
	*wcnd = dmin(r__1,r__2);

/*     COMPUTE ORTHOGONAL PROJECTION OF PARTICULAR SOLUTION */

L125:
	if (*inhomo != 1) {
	    goto L140;
	}
	lr = nr;
	if (*n != ml18jr_1.nfcc) {
	    lr = (nr << 1) - 1;
	}
	w[lr] = sdot_(m, &a[nr * a_dim1 + 1], &c__1, &v[1], &c__1) * ry;
	i__2 = *m;
	for (i__ = 1; i__ <= i__2; ++i__) {
/* L130: */
	    v[i__] -= w[lr] * a[i__ + nr * a_dim1];
	}
	if (*n == ml18jr_1.nfcc) {
	    goto L140;
	}
	lr = nr << 1;
	w[lr] = ry * prvec_(m, &v[1], &a[nr * a_dim1 + 1]);
	i__2 = m2;
	for (k = 1; k <= i__2; ++k) {
	    l = m2 + k;
	    v[k] += w[lr] * a[l + nr * a_dim1];
/* L135: */
	    v[l] -= w[lr] * a[k + nr * a_dim1];
	}
L140:
	;
    }
/* ********************************************************************** */

/*     TEST FOR LINEAR DEPENDENCE OF PARTICULAR SOLUTION */

L150:
    if (*inhomo != 1) {
	return 0;
    }
    if (*n > 1 && s[np1] < 1.f) {
	return 0;
    }
    vnorm = sdot_(m, &v[1], &c__1, &v[1], &c__1);
    if (s[np1] != 0.f) {
/* Computing MIN */
	r__1 = *wcnd, r__2 = vnorm / s[np1];
	*wcnd = dmin(r__1,r__2);
    }
    if (vnorm >= ml5mco_1.eps * s[np1]) {
	return 0;
    }
L170:
    *iflag = 2;
    *wcnd = ml5mco_1.eps;
    return 0;
} /* mgsbv_ */

