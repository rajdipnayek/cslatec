/* cmposp.f -- translated by f2c (version 12.02.01).
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

/* DECK CMPOSP */
/* Subroutine */ int cmposp_(integer *m, integer *n, complex *a, complex *bb, 
	complex *c__, complex *q, integer *idimq, complex *b, complex *b2, 
	complex *b3, complex *w, complex *w2, complex *w3, complex *d__, 
	complex *tcos, complex *p)
{
    /* System generated locals */
    integer q_dim1, q_offset, i__1, i__2, i__3, i__4;
    real r__1;
    complex q__1, q__2;

    /* Local variables */
    static integer i__, j;
    static complex s, t;
    static integer lh, mr, nr, nrm1, nrmj, nrpj;
    extern /* Subroutine */ int cmposd_(integer *, integer *, integer *, 
	    complex *, complex *, complex *, complex *, integer *, complex *, 
	    complex *, complex *, complex *, complex *), cmposn_(integer *, 
	    integer *, integer *, integer *, complex *, complex *, complex *, 
	    complex *, integer *, complex *, complex *, complex *, complex *, 
	    complex *, complex *, complex *, complex *, complex *);
    static integer ipstor;

/* ***BEGIN PROLOGUE  CMPOSP */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to CMGNBN */
/* ***LIBRARY   SLATEC */
/* ***TYPE      COMPLEX (POISP2-S, CMPOSP-C) */
/* ***AUTHOR  (UNKNOWN) */
/* ***DESCRIPTION */

/*     Subroutine to solve Poisson's equation with periodic boundary */
/*     conditions. */

/* ***SEE ALSO  CMGNBN */
/* ***ROUTINES CALLED  CMPOSD, CMPOSN */
/* ***REVISION HISTORY  (YYMMDD) */
/*   801001  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900402  Added TYPE section.  (WRB) */
/* ***END PROLOGUE  CMPOSP */

/* ***FIRST EXECUTABLE STATEMENT  CMPOSP */
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
	    i__3 = i__ + nrmj * q_dim1;
	    i__4 = i__ + nrpj * q_dim1;
	    q__1.r = q[i__3].r - q[i__4].r, q__1.i = q[i__3].i - q[i__4].i;
	    s.r = q__1.r, s.i = q__1.i;
	    i__3 = i__ + nrmj * q_dim1;
	    i__4 = i__ + nrpj * q_dim1;
	    q__1.r = q[i__3].r + q[i__4].r, q__1.i = q[i__3].i + q[i__4].i;
	    t.r = q__1.r, t.i = q__1.i;
	    i__3 = i__ + nrmj * q_dim1;
	    q[i__3].r = s.r, q[i__3].i = s.i;
	    i__3 = i__ + nrpj * q_dim1;
	    q[i__3].r = t.r, q[i__3].i = t.i;
/* L101: */
	}
/* L102: */
    }
    i__1 = mr;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = i__ + nr * q_dim1;
	i__3 = i__ + nr * q_dim1;
	q__1.r = q[i__3].r * 2.f, q__1.i = q[i__3].i * 2.f;
	q[i__2].r = q__1.r, q[i__2].i = q__1.i;
	i__2 = i__ + *n * q_dim1;
	i__3 = i__ + *n * q_dim1;
	q__1.r = q[i__3].r * 2.f, q__1.i = q[i__3].i * 2.f;
	q[i__2].r = q__1.r, q[i__2].i = q__1.i;
/* L103: */
    }
    cmposd_(&mr, &nrm1, &c__1, &a[1], &bb[1], &c__[1], &q[q_offset], idimq, &
	    b[1], &w[1], &d__[1], &tcos[1], &p[1]);
    ipstor = w[1].r;
    i__1 = nr + 1;
    cmposn_(&mr, &i__1, &c__1, &c__1, &a[1], &bb[1], &c__[1], &q[nr * q_dim1 
	    + 1], idimq, &b[1], &b2[1], &b3[1], &w[1], &w2[1], &w3[1], &d__[1]
	    , &tcos[1], &p[1]);
/* Computing MAX */
    i__1 = ipstor, i__2 = (integer) w[1].r;
    ipstor = max(i__1,i__2);
    i__1 = nrm1;
    for (j = 1; j <= i__1; ++j) {
	nrmj = nr - j;
	nrpj = nr + j;
	i__2 = mr;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    i__3 = i__ + nrpj * q_dim1;
	    i__4 = i__ + nrmj * q_dim1;
	    q__2.r = q[i__3].r + q[i__4].r, q__2.i = q[i__3].i + q[i__4].i;
	    q__1.r = q__2.r * .5f, q__1.i = q__2.i * .5f;
	    s.r = q__1.r, s.i = q__1.i;
	    i__3 = i__ + nrpj * q_dim1;
	    i__4 = i__ + nrmj * q_dim1;
	    q__2.r = q[i__3].r - q[i__4].r, q__2.i = q[i__3].i - q[i__4].i;
	    q__1.r = q__2.r * .5f, q__1.i = q__2.i * .5f;
	    t.r = q__1.r, t.i = q__1.i;
	    i__3 = i__ + nrmj * q_dim1;
	    q[i__3].r = s.r, q[i__3].i = s.i;
	    i__3 = i__ + nrpj * q_dim1;
	    q[i__3].r = t.r, q[i__3].i = t.i;
/* L104: */
	}
/* L105: */
    }
    i__1 = mr;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = i__ + nr * q_dim1;
	i__3 = i__ + nr * q_dim1;
	q__1.r = q[i__3].r * .5f, q__1.i = q[i__3].i * .5f;
	q[i__2].r = q__1.r, q[i__2].i = q__1.i;
	i__2 = i__ + *n * q_dim1;
	i__3 = i__ + *n * q_dim1;
	q__1.r = q[i__3].r * .5f, q__1.i = q[i__3].i * .5f;
	q[i__2].r = q__1.r, q[i__2].i = q__1.i;
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
	    i__3 = i__ + j * q_dim1;
	    i__4 = i__ + nrpj * q_dim1;
	    q__1.r = q[i__3].r - q[i__4].r, q__1.i = q[i__3].i - q[i__4].i;
	    s.r = q__1.r, s.i = q__1.i;
	    i__3 = i__ + j * q_dim1;
	    i__4 = i__ + nrpj * q_dim1;
	    q__1.r = q[i__3].r + q[i__4].r, q__1.i = q[i__3].i + q[i__4].i;
	    t.r = q__1.r, t.i = q__1.i;
	    i__3 = i__ + j * q_dim1;
	    q[i__3].r = s.r, q[i__3].i = s.i;
	    i__3 = i__ + nrpj * q_dim1;
	    q[i__3].r = t.r, q[i__3].i = t.i;
/* L108: */
	}
/* L109: */
    }
    i__1 = mr;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = i__ + nr * q_dim1;
	i__3 = i__ + nr * q_dim1;
	q__1.r = q[i__3].r * 2.f, q__1.i = q[i__3].i * 2.f;
	q[i__2].r = q__1.r, q[i__2].i = q__1.i;
/* L110: */
    }
    lh = nrm1 / 2;
    i__1 = lh;
    for (j = 1; j <= i__1; ++j) {
	nrmj = nr - j;
	i__2 = mr;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    i__3 = i__ + j * q_dim1;
	    s.r = q[i__3].r, s.i = q[i__3].i;
	    i__3 = i__ + j * q_dim1;
	    i__4 = i__ + nrmj * q_dim1;
	    q[i__3].r = q[i__4].r, q[i__3].i = q[i__4].i;
	    i__3 = i__ + nrmj * q_dim1;
	    q[i__3].r = s.r, q[i__3].i = s.i;
/* L111: */
	}
/* L112: */
    }
    cmposd_(&mr, &nrm1, &c__2, &a[1], &bb[1], &c__[1], &q[q_offset], idimq, &
	    b[1], &w[1], &d__[1], &tcos[1], &p[1]);
    ipstor = w[1].r;
    cmposn_(&mr, &nr, &c__2, &c__1, &a[1], &bb[1], &c__[1], &q[nr * q_dim1 + 
	    1], idimq, &b[1], &b2[1], &b3[1], &w[1], &w2[1], &w3[1], &d__[1], 
	    &tcos[1], &p[1]);
/* Computing MAX */
    i__1 = ipstor, i__2 = (integer) w[1].r;
    ipstor = max(i__1,i__2);
    i__1 = nrm1;
    for (j = 1; j <= i__1; ++j) {
	nrpj = nr + j;
	i__2 = mr;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    i__3 = i__ + nrpj * q_dim1;
	    i__4 = i__ + j * q_dim1;
	    q__2.r = q[i__3].r + q[i__4].r, q__2.i = q[i__3].i + q[i__4].i;
	    q__1.r = q__2.r * .5f, q__1.i = q__2.i * .5f;
	    s.r = q__1.r, s.i = q__1.i;
	    i__3 = i__ + nrpj * q_dim1;
	    i__4 = i__ + j * q_dim1;
	    q__2.r = q[i__3].r - q[i__4].r, q__2.i = q[i__3].i - q[i__4].i;
	    q__1.r = q__2.r * .5f, q__1.i = q__2.i * .5f;
	    t.r = q__1.r, t.i = q__1.i;
	    i__3 = i__ + nrpj * q_dim1;
	    q[i__3].r = t.r, q[i__3].i = t.i;
	    i__3 = i__ + j * q_dim1;
	    q[i__3].r = s.r, q[i__3].i = s.i;
/* L113: */
	}
/* L114: */
    }
    i__1 = mr;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = i__ + nr * q_dim1;
	i__3 = i__ + nr * q_dim1;
	q__1.r = q[i__3].r * .5f, q__1.i = q[i__3].i * .5f;
	q[i__2].r = q__1.r, q[i__2].i = q__1.i;
/* L115: */
    }
    i__1 = lh;
    for (j = 1; j <= i__1; ++j) {
	nrmj = nr - j;
	i__2 = mr;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    i__3 = i__ + j * q_dim1;
	    s.r = q[i__3].r, s.i = q[i__3].i;
	    i__3 = i__ + j * q_dim1;
	    i__4 = i__ + nrmj * q_dim1;
	    q[i__3].r = q[i__4].r, q[i__3].i = q[i__4].i;
	    i__3 = i__ + nrmj * q_dim1;
	    q[i__3].r = s.r, q[i__3].i = s.i;
/* L116: */
	}
/* L117: */
    }
L118:

/*     RETURN STORAGE REQUIREMENTS FOR P VECTORS. */

    r__1 = (real) ipstor;
    q__1.r = r__1, q__1.i = 0.f;
    w[1].r = q__1.r, w[1].i = q__1.i;
    return 0;
} /* cmposp_ */

