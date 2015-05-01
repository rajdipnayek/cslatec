/* poisp2.f -- translated by f2c (version 12.02.01).
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
static integer c__2 = 2;

/* DECK POISP2 */
/* Subroutine */ int poisp2_(integer *m, integer *n, real *a, real *bb, real *
	c__, real *q, integer *idimq, real *b, real *b2, real *b3, real *w, 
	real *w2, real *w3, real *d__, real *tcos, real *p)
{
    /* System generated locals */
    integer q_dim1, q_offset, i__1, i__2;

    /* Local variables */
    static integer i__, j;
    static real s, t;
    static integer lh, mr, nr, nrm1, nrmj, nrpj;
    extern /* Subroutine */ int poisd2_(integer *, integer *, integer *, real 
	    *, real *, real *, real *, integer *, real *, real *, real *, 
	    real *, real *), poisn2_(integer *, integer *, integer *, integer 
	    *, real *, real *, real *, real *, integer *, real *, real *, 
	    real *, real *, real *, real *, real *, real *, real *);
    static integer ipstor;

/* ***BEGIN PROLOGUE  POISP2 */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to GENBUN */
/* ***LIBRARY   SLATEC */
/* ***TYPE      SINGLE PRECISION (POISP2-S, CMPOSP-C) */
/* ***AUTHOR  (UNKNOWN) */
/* ***DESCRIPTION */

/*     Subroutine to solve Poisson equation with periodic boundary */
/*     conditions. */

/* ***SEE ALSO  GENBUN */
/* ***ROUTINES CALLED  POISD2, POISN2 */
/* ***REVISION HISTORY  (YYMMDD) */
/*   801001  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900402  Added TYPE section.  (WRB) */
/* ***END PROLOGUE  POISP2 */

/* ***FIRST EXECUTABLE STATEMENT  POISP2 */
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
    mr = *m;
    nr = (*n + 1) / 2;
    nrm1 = nr - 1;
    if (nr << 1 != *n) {
	goto L107;
    }

/*     EVEN NUMBER OF UNKNOWNS */

    i__1 = nrm1;
    for (j = 1; j <= i__1; ++j) {
	nrmj = nr - j;
	nrpj = nr + j;
	i__2 = mr;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    s = q[i__ + nrmj * q_dim1] - q[i__ + nrpj * q_dim1];
	    t = q[i__ + nrmj * q_dim1] + q[i__ + nrpj * q_dim1];
	    q[i__ + nrmj * q_dim1] = s;
	    q[i__ + nrpj * q_dim1] = t;
/* L101: */
	}
/* L102: */
    }
    i__1 = mr;
    for (i__ = 1; i__ <= i__1; ++i__) {
	q[i__ + nr * q_dim1] *= 2.f;
	q[i__ + *n * q_dim1] *= 2.f;
/* L103: */
    }
    poisd2_(&mr, &nrm1, &c__1, &a[1], &bb[1], &c__[1], &q[q_offset], idimq, &
	    b[1], &w[1], &d__[1], &tcos[1], &p[1]);
    ipstor = w[1];
    i__1 = nr + 1;
    poisn2_(&mr, &i__1, &c__1, &c__1, &a[1], &bb[1], &c__[1], &q[nr * q_dim1 
	    + 1], idimq, &b[1], &b2[1], &b3[1], &w[1], &w2[1], &w3[1], &d__[1]
	    , &tcos[1], &p[1]);
/* Computing MAX */
    i__1 = ipstor, i__2 = (integer) w[1];
    ipstor = max(i__1,i__2);
    i__1 = nrm1;
    for (j = 1; j <= i__1; ++j) {
	nrmj = nr - j;
	nrpj = nr + j;
	i__2 = mr;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    s = (q[i__ + nrpj * q_dim1] + q[i__ + nrmj * q_dim1]) * .5f;
	    t = (q[i__ + nrpj * q_dim1] - q[i__ + nrmj * q_dim1]) * .5f;
	    q[i__ + nrmj * q_dim1] = s;
	    q[i__ + nrpj * q_dim1] = t;
/* L104: */
	}
/* L105: */
    }
    i__1 = mr;
    for (i__ = 1; i__ <= i__1; ++i__) {
	q[i__ + nr * q_dim1] *= .5f;
	q[i__ + *n * q_dim1] *= .5f;
/* L106: */
    }
    goto L118;
L107:

/*     ODD  NUMBER OF UNKNOWNS */

    i__1 = nrm1;
    for (j = 1; j <= i__1; ++j) {
	nrpj = *n + 1 - j;
	i__2 = mr;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    s = q[i__ + j * q_dim1] - q[i__ + nrpj * q_dim1];
	    t = q[i__ + j * q_dim1] + q[i__ + nrpj * q_dim1];
	    q[i__ + j * q_dim1] = s;
	    q[i__ + nrpj * q_dim1] = t;
/* L108: */
	}
/* L109: */
    }
    i__1 = mr;
    for (i__ = 1; i__ <= i__1; ++i__) {
	q[i__ + nr * q_dim1] *= 2.f;
/* L110: */
    }
    lh = nrm1 / 2;
    i__1 = lh;
    for (j = 1; j <= i__1; ++j) {
	nrmj = nr - j;
	i__2 = mr;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    s = q[i__ + j * q_dim1];
	    q[i__ + j * q_dim1] = q[i__ + nrmj * q_dim1];
	    q[i__ + nrmj * q_dim1] = s;
/* L111: */
	}
/* L112: */
    }
    poisd2_(&mr, &nrm1, &c__2, &a[1], &bb[1], &c__[1], &q[q_offset], idimq, &
	    b[1], &w[1], &d__[1], &tcos[1], &p[1]);
    ipstor = w[1];
    poisn2_(&mr, &nr, &c__2, &c__1, &a[1], &bb[1], &c__[1], &q[nr * q_dim1 + 
	    1], idimq, &b[1], &b2[1], &b3[1], &w[1], &w2[1], &w3[1], &d__[1], 
	    &tcos[1], &p[1]);
/* Computing MAX */
    i__1 = ipstor, i__2 = (integer) w[1];
    ipstor = max(i__1,i__2);
    i__1 = nrm1;
    for (j = 1; j <= i__1; ++j) {
	nrpj = nr + j;
	i__2 = mr;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    s = (q[i__ + nrpj * q_dim1] + q[i__ + j * q_dim1]) * .5f;
	    t = (q[i__ + nrpj * q_dim1] - q[i__ + j * q_dim1]) * .5f;
	    q[i__ + nrpj * q_dim1] = t;
	    q[i__ + j * q_dim1] = s;
/* L113: */
	}
/* L114: */
    }
    i__1 = mr;
    for (i__ = 1; i__ <= i__1; ++i__) {
	q[i__ + nr * q_dim1] *= .5f;
/* L115: */
    }
    i__1 = lh;
    for (j = 1; j <= i__1; ++j) {
	nrmj = nr - j;
	i__2 = mr;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    s = q[i__ + j * q_dim1];
	    q[i__ + j * q_dim1] = q[i__ + nrmj * q_dim1];
	    q[i__ + nrmj * q_dim1] = s;
/* L116: */
	}
/* L117: */
    }
L118:

/*     RETURN STORAGE REQUIREMENTS FOR P VECTORS. */

    w[1] = (real) ipstor;
    return 0;
} /* poisp2_ */

