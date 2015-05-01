/* dh12.f -- translated by f2c (version 12.02.01).
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

/* DECK DH12 */
/* Subroutine */ int dh12_(integer *mode, integer *lpivot, integer *l1, 
	integer *m, doublereal *u, integer *iue, doublereal *up, doublereal *
	c__, integer *ice, integer *icv, integer *ncv)
{
    /* System generated locals */
    integer u_dim1, u_offset, i__1, i__2;
    doublereal d__1, d__2;

    /* Local variables */
    static doublereal b;
    static integer i__, j, i2, i3, i4;
    static doublereal cl, sm;
    static integer kl1, kl2, l1m1;
    static doublereal one;
    static integer klp;
    extern doublereal ddot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    static integer incr;
    static doublereal ul1m1;
    static integer mml1p2;
    static doublereal clinv;
    extern /* Subroutine */ int dswap_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), daxpy_(integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *);

/* ***BEGIN PROLOGUE  DH12 */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to DHFTI, DLSEI and DWNNLS */
/* ***LIBRARY   SLATEC */
/* ***TYPE      DOUBLE PRECISION (H12-S, DH12-D) */
/* ***AUTHOR  (UNKNOWN) */
/* ***DESCRIPTION */

/*      *** DOUBLE PRECISION VERSION OF H12 ****** */

/*     C.L.Lawson and R.J.Hanson, Jet Propulsion Laboratory, 1973 Jun 12 */
/*     to appear in 'Solving Least Squares Problems', Prentice-Hall, 1974 */

/*     Construction and/or application of a single */
/*     Householder transformation..     Q = I + U*(U**T)/B */

/*     MODE    = 1 or 2   to select algorithm  H1  or  H2 . */
/*     LPIVOT is the index of the pivot element. */
/*     L1,M   If L1 .LE. M   the transformation will be constructed to */
/*            zero elements indexed from L1 through M.   If L1 GT. M */
/*            THE SUBROUTINE DOES AN IDENTITY TRANSFORMATION. */
/*     U(),IUE,UP    On entry to H1 U() contains the pivot vector. */
/*                   IUE is the storage increment between elements. */
/*                                       On exit from H1 U() and UP */
/*                   contain quantities defining the vector U of the */
/*                   Householder transformation.   On entry to H2 U() */
/*                   and UP should contain quantities previously computed */
/*                   by H1.  These will not be modified by H2. */
/*     C()    On entry to H1 or H2 C() contains a matrix which will be */
/*            regarded as a set of vectors to which the Householder */
/*            transformation is to be applied.  On exit C() contains the */
/*            set of transformed vectors. */
/*     ICE    Storage increment between elements of vectors in C(). */
/*     ICV    Storage increment between vectors in C(). */
/*     NCV    Number of vectors in C() to be transformed. If NCV .LE. 0 */
/*            no operations will be done on C(). */

/* ***SEE ALSO  DHFTI, DLSEI, DWNNLS */
/* ***ROUTINES CALLED  DAXPY, DDOT, DSWAP */
/* ***REVISION HISTORY  (YYMMDD) */
/*   790101  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890831  Modified array declarations.  (WRB) */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900328  Added TYPE section.  (WRB) */
/*   900911  Added DDOT to DOUBLE PRECISION statement.  (WRB) */
/* ***END PROLOGUE  DH12 */
/*     BEGIN BLOCK PERMITTING ...EXITS TO 140 */
/* ***FIRST EXECUTABLE STATEMENT  DH12 */
    /* Parameter adjustments */
    u_dim1 = *iue;
    u_offset = 1 + u_dim1;
    u -= u_offset;
    --c__;

    /* Function Body */
    one = 1.;

/*     ...EXIT */
    if (0 >= *lpivot || *lpivot >= *l1 || *l1 > *m) {
	goto L140;
    }
    cl = (d__1 = u[*lpivot * u_dim1 + 1], abs(d__1));
    if (*mode == 2) {
	goto L40;
    }
/*           ****** CONSTRUCT THE TRANSFORMATION. ****** */
    i__1 = *m;
    for (j = *l1; j <= i__1; ++j) {
/* Computing MAX */
	d__2 = (d__1 = u[j * u_dim1 + 1], abs(d__1));
	cl = max(d__2,cl);
/* L10: */
    }
    if (cl > 0.) {
	goto L20;
    }
/*     .........EXIT */
    goto L140;
L20:
    clinv = one / cl;
/* Computing 2nd power */
    d__1 = u[*lpivot * u_dim1 + 1] * clinv;
    sm = d__1 * d__1;
    i__1 = *m;
    for (j = *l1; j <= i__1; ++j) {
/* Computing 2nd power */
	d__1 = u[j * u_dim1 + 1] * clinv;
	sm += d__1 * d__1;
/* L30: */
    }
    cl *= sqrt(sm);
    if (u[*lpivot * u_dim1 + 1] > 0.) {
	cl = -cl;
    }
    *up = u[*lpivot * u_dim1 + 1] - cl;
    u[*lpivot * u_dim1 + 1] = cl;
    goto L50;
L40:
/*        ****** APPLY THE TRANSFORMATION  I+U*(U**T)/B  TO C. ****** */

    if (cl > 0.) {
	goto L50;
    }
/*     ......EXIT */
    goto L140;
L50:
/*     ...EXIT */
    if (*ncv <= 0) {
	goto L140;
    }
    b = *up * u[*lpivot * u_dim1 + 1];
/*        B  MUST BE NONPOSITIVE HERE.  IF B = 0., RETURN. */

    if (b < 0.) {
	goto L60;
    }
/*     ......EXIT */
    goto L140;
L60:
    b = one / b;
    mml1p2 = *m - *l1 + 2;
    if (mml1p2 <= 20) {
	goto L80;
    }
    l1m1 = *l1 - 1;
    kl1 = (l1m1 - 1) * *ice + 1;
    kl2 = kl1;
    klp = (*lpivot - 1) * *ice + 1;
    ul1m1 = u[l1m1 * u_dim1 + 1];
    u[l1m1 * u_dim1 + 1] = *up;
    if (*lpivot != l1m1) {
	dswap_(ncv, &c__[kl1], icv, &c__[klp], icv);
    }
    i__1 = *ncv;
    for (j = 1; j <= i__1; ++j) {
	sm = ddot_(&mml1p2, &u[l1m1 * u_dim1 + 1], iue, &c__[kl1], ice);
	sm *= b;
	daxpy_(&mml1p2, &sm, &u[l1m1 * u_dim1 + 1], iue, &c__[kl1], ice);
	kl1 += *icv;
/* L70: */
    }
    u[l1m1 * u_dim1 + 1] = ul1m1;
/*     ......EXIT */
    if (*lpivot == l1m1) {
	goto L140;
    }
    kl1 = kl2;
    dswap_(ncv, &c__[kl1], icv, &c__[klp], icv);
    goto L130;
L80:
    i2 = 1 - *icv + *ice * (*lpivot - 1);
    incr = *ice * (*l1 - *lpivot);
    i__1 = *ncv;
    for (j = 1; j <= i__1; ++j) {
	i2 += *icv;
	i3 = i2 + incr;
	i4 = i3;
	sm = c__[i2] * *up;
	i__2 = *m;
	for (i__ = *l1; i__ <= i__2; ++i__) {
	    sm += c__[i3] * u[i__ * u_dim1 + 1];
	    i3 += *ice;
/* L90: */
	}
	if (sm == 0.) {
	    goto L110;
	}
	sm *= b;
	c__[i2] += sm * *up;
	i__2 = *m;
	for (i__ = *l1; i__ <= i__2; ++i__) {
	    c__[i4] += sm * u[i__ * u_dim1 + 1];
	    i4 += *ice;
/* L100: */
	}
L110:
/* L120: */
	;
    }
L130:
L140:
    return 0;
} /* dh12_ */

