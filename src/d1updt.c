/* d1updt.f -- translated by f2c (version 12.02.01).
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

static integer c__2 = 2;

/* DECK D1UPDT */
/* Subroutine */ int d1updt_(integer *m, integer *n, doublereal *s, integer *
	ls, doublereal *u, doublereal *v, doublereal *w, logical *sing)
{
    /* Initialized data */

    static doublereal one = 1.;
    static doublereal p5 = .5;
    static doublereal p25 = .25;
    static doublereal zero = 0.;

    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1, d__2;

    /* Local variables */
    static integer i__, j, l, jj, nm1;
    static doublereal tan__;
    static integer nmj;
    static doublereal cos__, sin__, tau, temp, giant, cotan;
    extern doublereal d1mach_(integer *);

/* ***BEGIN PROLOGUE  D1UPDT */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to DNSQ and DNSQE */
/* ***LIBRARY   SLATEC */
/* ***TYPE      DOUBLE PRECISION (R1UPDT-S, D1UPDT-D) */
/* ***AUTHOR  (UNKNOWN) */
/* ***DESCRIPTION */

/*     Given an M by N lower trapezoidal matrix S, an M-vector U, */
/*     and an N-vector V, the problem is to determine an */
/*     orthogonal matrix Q such that */

/*                   t */
/*           (S + U*V )*Q */

/*     is again lower trapezoidal. */

/*     This subroutine determines Q as the product of 2*(N - 1) */
/*     transformations */

/*           GV(N-1)*...*GV(1)*GW(1)*...*GW(N-1) */

/*     where GV(I), GW(I) are Givens rotations in the (I,N) plane */
/*     which eliminate elements in the I-th and N-th planes, */
/*     respectively. Q itself is not accumulated, rather the */
/*     information to recover the GV, GW rotations is returned. */

/*     The SUBROUTINE statement is */

/*       SUBROUTINE D1UPDT(M,N,S,LS,U,V,W,SING) */

/*     where */

/*       M is a positive integer input variable set to the number */
/*         of rows of S. */

/*       N is a positive integer input variable set to the number */
/*         of columns of S. N must not exceed M. */

/*       S is an array of length LS. On input S must contain the lower */
/*         trapezoidal matrix S stored by columns. On output S contains */
/*         the lower trapezoidal matrix produced as described above. */

/*       LS is a positive integer input variable not less than */
/*         (N*(2*M-N+1))/2. */

/*       U is an input array of length M which must contain the */
/*         vector U. */

/*       V is an array of length N. On input V must contain the vector */
/*         V. On output V(I) contains the information necessary to */
/*         recover the Givens rotation GV(I) described above. */

/*       W is an output array of length M. W(I) contains information */
/*         necessary to recover the Givens rotation GW(I) described */
/*         above. */

/*       SING is a LOGICAL output variable. SING is set TRUE if any */
/*         of the diagonal elements of the output S are zero. Otherwise */
/*         SING is set FALSE. */

/* ***SEE ALSO  DNSQ, DNSQE */
/* ***ROUTINES CALLED  D1MACH */
/* ***REVISION HISTORY  (YYMMDD) */
/*   800301  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890831  Modified array declarations.  (WRB) */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900326  Removed duplicate information from DESCRIPTION section. */
/*           (WRB) */
/*   900328  Added TYPE section.  (WRB) */
/* ***END PROLOGUE  D1UPDT */
    /* Parameter adjustments */
    --w;
    --v;
    --u;
    --s;

    /* Function Body */

/*     GIANT IS THE LARGEST MAGNITUDE. */

/* ***FIRST EXECUTABLE STATEMENT  D1UPDT */
    giant = d1mach_(&c__2);

/*     INITIALIZE THE DIAGONAL ELEMENT POINTER. */

    jj = *n * ((*m << 1) - *n + 1) / 2 - (*m - *n);

/*     MOVE THE NONTRIVIAL PART OF THE LAST COLUMN OF S INTO W. */

    l = jj;
    i__1 = *m;
    for (i__ = *n; i__ <= i__1; ++i__) {
	w[i__] = s[l];
	++l;
/* L10: */
    }

/*     ROTATE THE VECTOR V INTO A MULTIPLE OF THE N-TH UNIT VECTOR */
/*     IN SUCH A WAY THAT A SPIKE IS INTRODUCED INTO W. */

    nm1 = *n - 1;
    if (nm1 < 1) {
	goto L70;
    }
    i__1 = nm1;
    for (nmj = 1; nmj <= i__1; ++nmj) {
	j = *n - nmj;
	jj -= *m - j + 1;
	w[j] = zero;
	if (v[j] == zero) {
	    goto L50;
	}

/*        DETERMINE A GIVENS ROTATION WHICH ELIMINATES THE */
/*        J-TH ELEMENT OF V. */

	if ((d__1 = v[*n], abs(d__1)) >= (d__2 = v[j], abs(d__2))) {
	    goto L20;
	}
	cotan = v[*n] / v[j];
/* Computing 2nd power */
	d__1 = cotan;
	sin__ = p5 / sqrt(p25 + p25 * (d__1 * d__1));
	cos__ = sin__ * cotan;
	tau = one;
	if (abs(cos__) * giant > one) {
	    tau = one / cos__;
	}
	goto L30;
L20:
	tan__ = v[j] / v[*n];
/* Computing 2nd power */
	d__1 = tan__;
	cos__ = p5 / sqrt(p25 + p25 * (d__1 * d__1));
	sin__ = cos__ * tan__;
	tau = sin__;
L30:

/*        APPLY THE TRANSFORMATION TO V AND STORE THE INFORMATION */
/*        NECESSARY TO RECOVER THE GIVENS ROTATION. */

	v[*n] = sin__ * v[j] + cos__ * v[*n];
	v[j] = tau;

/*        APPLY THE TRANSFORMATION TO S AND EXTEND THE SPIKE IN W. */

	l = jj;
	i__2 = *m;
	for (i__ = j; i__ <= i__2; ++i__) {
	    temp = cos__ * s[l] - sin__ * w[i__];
	    w[i__] = sin__ * s[l] + cos__ * w[i__];
	    s[l] = temp;
	    ++l;
/* L40: */
	}
L50:
/* L60: */
	;
    }
L70:

/*     ADD THE SPIKE FROM THE RANK 1 UPDATE TO W. */

    i__1 = *m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	w[i__] += v[*n] * u[i__];
/* L80: */
    }

/*     ELIMINATE THE SPIKE. */

    *sing = FALSE_;
    if (nm1 < 1) {
	goto L140;
    }
    i__1 = nm1;
    for (j = 1; j <= i__1; ++j) {
	if (w[j] == zero) {
	    goto L120;
	}

/*        DETERMINE A GIVENS ROTATION WHICH ELIMINATES THE */
/*        J-TH ELEMENT OF THE SPIKE. */

	if ((d__1 = s[jj], abs(d__1)) >= (d__2 = w[j], abs(d__2))) {
	    goto L90;
	}
	cotan = s[jj] / w[j];
/* Computing 2nd power */
	d__1 = cotan;
	sin__ = p5 / sqrt(p25 + p25 * (d__1 * d__1));
	cos__ = sin__ * cotan;
	tau = one;
	if (abs(cos__) * giant > one) {
	    tau = one / cos__;
	}
	goto L100;
L90:
	tan__ = w[j] / s[jj];
/* Computing 2nd power */
	d__1 = tan__;
	cos__ = p5 / sqrt(p25 + p25 * (d__1 * d__1));
	sin__ = cos__ * tan__;
	tau = sin__;
L100:

/*        APPLY THE TRANSFORMATION TO S AND REDUCE THE SPIKE IN W. */

	l = jj;
	i__2 = *m;
	for (i__ = j; i__ <= i__2; ++i__) {
	    temp = cos__ * s[l] + sin__ * w[i__];
	    w[i__] = -sin__ * s[l] + cos__ * w[i__];
	    s[l] = temp;
	    ++l;
/* L110: */
	}

/*        STORE THE INFORMATION NECESSARY TO RECOVER THE */
/*        GIVENS ROTATION. */

	w[j] = tau;
L120:

/*        TEST FOR ZERO DIAGONAL ELEMENTS IN THE OUTPUT S. */

	if (s[jj] == zero) {
	    *sing = TRUE_;
	}
	jj += *m - j + 1;
/* L130: */
    }
L140:

/*     MOVE W BACK INTO THE LAST COLUMN OF THE OUTPUT S. */

    l = jj;
    i__1 = *m;
    for (i__ = *n; i__ <= i__1; ++i__) {
	s[l] = w[i__];
	++l;
/* L150: */
    }
    if (s[jj] == zero) {
	*sing = TRUE_;
    }
    return 0;

/*     LAST CARD OF SUBROUTINE D1UPDT. */

} /* d1updt_ */

