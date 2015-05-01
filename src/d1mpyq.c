/* d1mpyq.f -- translated by f2c (version 12.02.01).
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

/* DECK D1MPYQ */
/* Subroutine */ int d1mpyq_(integer *m, integer *n, doublereal *a, integer *
	lda, doublereal *v, doublereal *w)
{
    /* Initialized data */

    static doublereal one = 1.;

    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;
    doublereal d__1, d__2;

    /* Local variables */
    static integer i__, j, nm1, nmj;
    static doublereal cos__, sin__, temp;

/* ***BEGIN PROLOGUE  D1MPYQ */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to DNSQ and DNSQE */
/* ***LIBRARY   SLATEC */
/* ***TYPE      DOUBLE PRECISION (R1MPYQ-S, D1MPYQ-D) */
/* ***AUTHOR  (UNKNOWN) */
/* ***DESCRIPTION */

/*     Given an M by N matrix A, this subroutine computes A*Q where */
/*     Q is the product of 2*(N - 1) transformations */

/*           GV(N-1)*...*GV(1)*GW(1)*...*GW(N-1) */

/*     and GV(I), GW(I) are Givens rotations in the (I,N) plane which */
/*     eliminate elements in the I-th and N-th planes, respectively. */
/*     Q itself is not given, rather the information to recover the */
/*     GV, GW rotations is supplied. */

/*     The SUBROUTINE statement is */

/*       SUBROUTINE D1MPYQ(M,N,A,LDA,V,W) */

/*     where */

/*       M is a positive integer input variable set to the number */
/*         of rows of A. */

/*       N IS a positive integer input variable set to the number */
/*         of columns of A. */

/*       A is an M by N array. On input A must contain the matrix */
/*         to be postmultiplied by the orthogonal matrix Q */
/*         described above. On output A*Q has replaced A. */

/*       LDA is a positive integer input variable not less than M */
/*         which specifies the leading dimension of the array A. */

/*       V is an input array of length N. V(I) must contain the */
/*         information necessary to recover the Givens rotation GV(I) */
/*         described above. */

/*       W is an input array of length N. W(I) must contain the */
/*         information necessary to recover the Givens rotation GW(I) */
/*         described above. */

/* ***SEE ALSO  DNSQ, DNSQE */
/* ***ROUTINES CALLED  (NONE) */
/* ***REVISION HISTORY  (YYMMDD) */
/*   800301  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890831  Modified array declarations.  (WRB) */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900326  Removed duplicate information from DESCRIPTION section. */
/*           (WRB) */
/*   900328  Added TYPE section.  (WRB) */
/* ***END PROLOGUE  D1MPYQ */
    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --v;
    --w;

    /* Function Body */

/*     APPLY THE FIRST SET OF GIVENS ROTATIONS TO A. */

/* ***FIRST EXECUTABLE STATEMENT  D1MPYQ */
    nm1 = *n - 1;
    if (nm1 < 1) {
	goto L50;
    }
    i__1 = nm1;
    for (nmj = 1; nmj <= i__1; ++nmj) {
	j = *n - nmj;
	if ((d__1 = v[j], abs(d__1)) > one) {
	    cos__ = one / v[j];
	}
	if ((d__1 = v[j], abs(d__1)) > one) {
/* Computing 2nd power */
	    d__2 = cos__;
	    sin__ = sqrt(one - d__2 * d__2);
	}
	if ((d__1 = v[j], abs(d__1)) <= one) {
	    sin__ = v[j];
	}
	if ((d__1 = v[j], abs(d__1)) <= one) {
/* Computing 2nd power */
	    d__2 = sin__;
	    cos__ = sqrt(one - d__2 * d__2);
	}
	i__2 = *m;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    temp = cos__ * a[i__ + j * a_dim1] - sin__ * a[i__ + *n * a_dim1];
	    a[i__ + *n * a_dim1] = sin__ * a[i__ + j * a_dim1] + cos__ * a[
		    i__ + *n * a_dim1];
	    a[i__ + j * a_dim1] = temp;
/* L10: */
	}
/* L20: */
    }

/*     APPLY THE SECOND SET OF GIVENS ROTATIONS TO A. */

    i__1 = nm1;
    for (j = 1; j <= i__1; ++j) {
	if ((d__1 = w[j], abs(d__1)) > one) {
	    cos__ = one / w[j];
	}
	if ((d__1 = w[j], abs(d__1)) > one) {
/* Computing 2nd power */
	    d__2 = cos__;
	    sin__ = sqrt(one - d__2 * d__2);
	}
	if ((d__1 = w[j], abs(d__1)) <= one) {
	    sin__ = w[j];
	}
	if ((d__1 = w[j], abs(d__1)) <= one) {
/* Computing 2nd power */
	    d__2 = sin__;
	    cos__ = sqrt(one - d__2 * d__2);
	}
	i__2 = *m;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    temp = cos__ * a[i__ + j * a_dim1] + sin__ * a[i__ + *n * a_dim1];
	    a[i__ + *n * a_dim1] = -sin__ * a[i__ + j * a_dim1] + cos__ * a[
		    i__ + *n * a_dim1];
	    a[i__ + j * a_dim1] = temp;
/* L30: */
	}
/* L40: */
    }
L50:
    return 0;

/*     LAST CARD OF SUBROUTINE D1MPYQ. */

} /* d1mpyq_ */

