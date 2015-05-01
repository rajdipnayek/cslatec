/* dmgsbv.f -- translated by f2c (version 12.02.01).
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
    doublereal ae, re, tol;
    integer nxpts, nic, nopg, mxnon, ndisk, ntape, neq, indpvt, integ, nps, 
	    ntp, neqivp, numort, nfcc, icoco;
} dml18j_;

#define dml18j_1 dml18j_

struct {
    doublereal uro, sru, eps, sqovfl, twou, fouru;
    integer lpar;
} dml5mc_;

#define dml5mc_1 dml5mc_

/* Table of constant values */

static integer c__1 = 1;

/* DECK DMGSBV */
/* Subroutine */ int dmgsbv_(integer *m, integer *n, doublereal *a, integer *
	ia, integer *niv, integer *iflag, doublereal *s, doublereal *p, 
	integer *ip, integer *inhomo, doublereal *v, doublereal *w, 
	doublereal *wcnd)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3;
    doublereal d__1, d__2;

    /* Local variables */
    static integer i__, j, k, l;
    static doublereal t, y;
    static integer m2, kd, jk, kj, jp, jq, kp, nn, lr, nr, ix;
    static doublereal vl;
    static integer iz, jy, jz;
    static doublereal sv, ry;
    static integer ip1, np1;
    static doublereal dot, pjp;
    static integer lix, nrm1;
    extern doublereal ddot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    static integer nivn, nmnr;
    static doublereal psave, vnorm;
    extern doublereal dprvec_(integer *, doublereal *, doublereal *);

/* ***BEGIN PROLOGUE  DMGSBV */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to DBVSUP */
/* ***LIBRARY   SLATEC */
/* ***TYPE      DOUBLE PRECISION (MGSBV-S, DMGSBV-D) */
/* ***AUTHOR  Watts, H. A., (SNLA) */
/* ***DESCRIPTION */

/* ********************************************************************** */
/* Orthogonalize a set of N double precision vectors and determine their */
/* rank. */

/* ********************************************************************** */
/* INPUT */
/* ********************************************************************** */
/*   M = dimension of vectors. */
/*   N = no. of vectors. */
/*   A = array whose first N cols contain the vectors. */
/*   IA = first dimension of array A (col length). */
/*   NIV = number of independent vectors needed. */
/*   INHOMO = 1 corresponds to having a non-zero particular solution. */
/*   V = particular solution vector (not included in the pivoting). */
/*   INDPVT = 1 means pivoting will not be used. */

/* ********************************************************************** */
/* OUTPUT */
/* ********************************************************************** */
/*   NIV = no. of linear independent vectors in input set. */
/*     A = matrix whose first NIV cols. contain NIV orthogonal vectors */
/*         which span the vector space determined by the input vectors. */
/*   IFLAG */
/*          = 0 success */
/*          = 1 incorrect input */
/*          = 2 rank of new vectors less than N */
/*   P = decomposition matrix.  P is upper triangular and */
/*             (old vectors) = (new vectors) * P. */
/*         The old vectors will be reordered due to pivoting. */
/*         The dimension of P must be .GE. N*(N+1)/2. */
/*             (  N*(2*N+1) when N .NE. NFCC ) */
/*   IP = pivoting vector. The dimension of IP must be .GE. N. */
/*             (  2*N when N .NE. NFCC ) */
/*   S = square of norms of incoming vectors. */
/*   V = vector which is orthogonal to the vectors of A. */
/*   W = orthogonalization information for the vector V. */
/*   WCND = worst case (smallest) norm decrement value of the */
/*          vectors being orthogonalized  (represents a test */
/*          for linear dependence of the vectors). */
/* ********************************************************************** */

/* ***SEE ALSO  DBVSUP */
/* ***ROUTINES CALLED  DDOT, DPRVEC */
/* ***COMMON BLOCKS    DML18J, DML5MC */
/* ***REVISION HISTORY  (YYMMDD) */
/*   750601  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890831  Modified array declarations.  (WRB) */
/*   890921  Realigned order of variables in certain COMMON blocks. */
/*           (WRB) */
/*   890921  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900328  Added TYPE section.  (WRB) */
/*   910722  Updated AUTHOR section.  (ALS) */
/* ***END PROLOGUE  DMGSBV */





/* ***FIRST EXECUTABLE STATEMENT  DMGSBV */
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
    goto L280;
L10:
/*        BEGIN BLOCK PERMITTING ...EXITS TO 270 */
/*           BEGIN BLOCK PERMITTING ...EXITS TO 260 */

    jp = 0;
    *iflag = 0;
    np1 = *n + 1;
    y = 0.;
    m2 = *m / 2;

/*              CALCULATE SQUARE OF NORMS OF INCOMING VECTORS AND SEARCH */
/*              FOR VECTOR WITH LARGEST MAGNITUDE */

    j = 0;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	vl = ddot_(m, &a[i__ * a_dim1 + 1], &c__1, &a[i__ * a_dim1 + 1], &
		c__1);
	s[i__] = vl;
	if (*n == dml18j_1.nfcc) {
	    goto L20;
	}
	j = (i__ << 1) - 1;
	p[j] = vl;
	ip[j] = j;
L20:
	++j;
	p[j] = vl;
	ip[j] = j;
	if (vl <= y) {
	    goto L30;
	}
	y = vl;
	ix = i__;
L30:
/* L40: */
	;
    }
    if (dml18j_1.indpvt != 1) {
	goto L50;
    }
    ix = 1;
    y = p[1];
L50:
    lix = ix;
    if (*n != dml18j_1.nfcc) {
	lix = (ix << 1) - 1;
    }
    p[lix] = p[1];
    s[np1] = 0.;
    if (*inhomo == 1) {
	s[np1] = ddot_(m, &v[1], &c__1, &v[1], &c__1);
    }
    *wcnd = 1.;
    nivn = *niv;
    *niv = 0;

/*           ...EXIT */
    if (y == 0.) {
	goto L260;
    }
/*              ********************************************************* */
    i__1 = *n;
    for (nr = 1; nr <= i__1; ++nr) {
/*                 BEGIN BLOCK PERMITTING ...EXITS TO 230 */
/*              ......EXIT */
	if (nivn == *niv) {
	    goto L250;
	}
	*niv = nr;
	if (ix == nr) {
	    goto L130;
	}

/*                       PIVOTING OF COLUMNS OF P MATRIX */

	nn = *n;
	lix = ix;
	lr = nr;
	if (*n == dml18j_1.nfcc) {
	    goto L60;
	}
	nn = dml18j_1.nfcc;
	lix = (ix << 1) - 1;
	lr = (nr << 1) - 1;
L60:
	if (nr == 1) {
	    goto L80;
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
	    kj = kj + nn - j;
/* L70: */
	}
	jy = jk + nmnr;
	jz = jy - kd;
	p[jy] = p[jz];
L80:
	iz = ip[lix];
	ip[lix] = ip[lr];
	ip[lr] = iz;
	sv = s[ix];
	s[ix] = s[nr];
	s[nr] = sv;
	if (*n == dml18j_1.nfcc) {
	    goto L110;
	}
	if (nr == 1) {
	    goto L100;
	}
	kj = lr + 1;
	i__2 = nrm1;
	for (k = 1; k <= i__2; ++k) {
	    psave = p[kj];
	    jk = kj + kd;
	    p[kj] = p[jk];
	    p[jk] = psave;
	    kj = kj + dml18j_1.nfcc - k;
/* L90: */
	}
L100:
	iz = ip[lix + 1];
	ip[lix + 1] = ip[lr + 1];
	ip[lr + 1] = iz;
L110:

/*                       PIVOTING OF COLUMNS OF VECTORS */

	i__2 = *m;
	for (l = 1; l <= i__2; ++l) {
	    t = a[l + ix * a_dim1];
	    a[l + ix * a_dim1] = a[l + nr * a_dim1];
	    a[l + nr * a_dim1] = t;
/* L120: */
	}
L130:

/*                    CALCULATE P(NR,NR) AS NORM SQUARED OF PIVOTAL */
/*                    VECTOR */

	++jp;
	p[jp] = y;
	ry = 1. / y;
	nmnr = *n - nr;
	if (*n == dml18j_1.nfcc) {
	    goto L140;
	}
	nmnr = dml18j_1.nfcc - ((nr << 1) - 1);
	++jp;
	p[jp] = 0.;
	kp = jp + nmnr;
	p[kp] = y;
L140:
	if (nr == *n || nivn == *niv) {
	    goto L200;
	}

/*                       CALCULATE ORTHOGONAL PROJECTION VECTORS AND */
/*                       SEARCH FOR LARGEST NORM */

	y = 0.;
	ip1 = nr + 1;
	ix = ip1;
/*                       ************************************************ */
	i__2 = *n;
	for (j = ip1; j <= i__2; ++j) {
	    dot = ddot_(m, &a[nr * a_dim1 + 1], &c__1, &a[j * a_dim1 + 1], &
		    c__1);
	    ++jp;
	    jq = jp + nmnr;
	    if (*n != dml18j_1.nfcc) {
		jq = jq + nmnr - 1;
	    }
	    p[jq] = p[jp] - dot * (dot * ry);
	    p[jp] = dot * ry;
	    i__3 = *m;
	    for (i__ = 1; i__ <= i__3; ++i__) {
		a[i__ + j * a_dim1] -= p[jp] * a[i__ + nr * a_dim1];
/* L150: */
	    }
	    if (*n == dml18j_1.nfcc) {
		goto L170;
	    }
	    kp = jp + nmnr;
	    ++jp;
	    pjp = ry * dprvec_(m, &a[nr * a_dim1 + 1], &a[j * a_dim1 + 1]);
	    p[jp] = pjp;
	    p[kp] = -pjp;
	    ++kp;
	    p[kp] = ry * dot;
	    i__3 = m2;
	    for (k = 1; k <= i__3; ++k) {
		l = m2 + k;
		a[k + j * a_dim1] -= pjp * a[l + nr * a_dim1];
		a[l + j * a_dim1] += pjp * a[k + nr * a_dim1];
/* L160: */
	    }
	    p[jq] -= pjp * (pjp / ry);
L170:

/*                          TEST FOR CANCELLATION IN RECURRENCE RELATION */

	    if (p[jq] <= s[j] * dml5mc_1.sru) {
		p[jq] = ddot_(m, &a[j * a_dim1 + 1], &c__1, &a[j * a_dim1 + 1]
			, &c__1);
	    }
	    if (p[jq] <= y) {
		goto L180;
	    }
	    y = p[jq];
	    ix = j;
L180:
/* L190: */
	    ;
	}
	if (*n != dml18j_1.nfcc) {
	    jp = kp;
	}
/*                       ************************************************ */
	if (dml18j_1.indpvt == 1) {
	    ix = ip1;
	}

/*                       RECOMPUTE NORM SQUARED OF PIVOTAL VECTOR WITH */
/*                       SCALAR PRODUCT */

	y = ddot_(m, &a[ix * a_dim1 + 1], &c__1, &a[ix * a_dim1 + 1], &c__1);
/*           ............EXIT */
	if (y <= dml5mc_1.eps * s[ix]) {
	    goto L260;
	}
/* Computing MIN */
	d__1 = *wcnd, d__2 = y / s[ix];
	*wcnd = min(d__1,d__2);
L200:

/*                    COMPUTE ORTHOGONAL PROJECTION OF PARTICULAR */
/*                    SOLUTION */

/*                 ...EXIT */
	if (*inhomo != 1) {
	    goto L230;
	}
	lr = nr;
	if (*n != dml18j_1.nfcc) {
	    lr = (nr << 1) - 1;
	}
	w[lr] = ddot_(m, &a[nr * a_dim1 + 1], &c__1, &v[1], &c__1) * ry;
	i__2 = *m;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    v[i__] -= w[lr] * a[i__ + nr * a_dim1];
/* L210: */
	}
/*                 ...EXIT */
	if (*n == dml18j_1.nfcc) {
	    goto L230;
	}
	lr = nr << 1;
	w[lr] = ry * dprvec_(m, &v[1], &a[nr * a_dim1 + 1]);
	i__2 = m2;
	for (k = 1; k <= i__2; ++k) {
	    l = m2 + k;
	    v[k] += w[lr] * a[l + nr * a_dim1];
	    v[l] -= w[lr] * a[k + nr * a_dim1];
/* L220: */
	}
L230:
/* L240: */
	;
    }
L250:
/*              ********************************************************* */

/*                  TEST FOR LINEAR DEPENDENCE OF PARTICULAR SOLUTION */

/*        ......EXIT */
    if (*inhomo != 1) {
	goto L270;
    }
    if (*n > 1 && s[np1] < 1.f) {
	goto L270;
    }
    vnorm = ddot_(m, &v[1], &c__1, &v[1], &c__1);
    if (s[np1] != 0.) {
/* Computing MIN */
	d__1 = *wcnd, d__2 = vnorm / s[np1];
	*wcnd = min(d__1,d__2);
    }
/*        ......EXIT */
    if (vnorm >= dml5mc_1.eps * s[np1]) {
	goto L270;
    }
L260:
    *iflag = 2;
    *wcnd = dml5mc_1.eps;
L270:
L280:
    return 0;
} /* dmgsbv_ */

