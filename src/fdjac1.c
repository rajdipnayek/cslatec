/* fdjac1.f -- translated by f2c (version 12.02.01).
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

static integer c__4 = 4;

/* DECK FDJAC1 */
/* Subroutine */ int fdjac1_(S_fp fcn, integer *n, real *x, real *fvec, real *
	fjac, integer *ldfjac, integer *iflag, integer *ml, integer *mu, real 
	*epsfcn, real *wa1, real *wa2)
{
    /* Initialized data */

    static real zero = 0.f;

    /* System generated locals */
    integer fjac_dim1, fjac_offset, i__1, i__2, i__3, i__4;
    real r__1;

    /* Local variables */
    static real h__;
    static integer i__, j, k;
    static real eps, temp;
    static integer msum;
    extern doublereal r1mach_(integer *);
    static real epsmch;

/* ***BEGIN PROLOGUE  FDJAC1 */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to SNSQ and SNSQE */
/* ***LIBRARY   SLATEC */
/* ***TYPE      SINGLE PRECISION (FDJAC1-S, DFDJC1-D) */
/* ***AUTHOR  (UNKNOWN) */
/* ***DESCRIPTION */

/*     This subroutine computes a forward-difference approximation */
/*     to the N by N Jacobian matrix associated with a specified */
/*     problem of N functions in N VARIABLES. If the Jacobian has */
/*     a banded form, then function evaluations are saved by only */
/*     approximating the nonzero terms. */

/*     The subroutine statement is */

/*       SUBROUTINE FDJAC1(FCN,N,X,FVEC,FJAC,LDFJAC,IFLAG,ML,MU,EPSFCN, */
/*                         WA1,WA2) */

/*     where */

/*       FCN is the name of the user-supplied subroutine which */
/*         calculates the functions. FCN must be declared */
/*         in an external statement in the user calling */
/*         program, and should be written as follows. */

/*         SUBROUTINE FCN(N,X,FVEC,IFLAG) */
/*         INTEGER N,IFLAG */
/*         REAL X(N),FVEC(N) */
/*         ---------- */
/*         Calculate the functions at X and */
/*         return this vector in FVEC. */
/*         ---------- */
/*         RETURN */
/*         END */

/*         The value of IFLAG should not be changed by FCN unless */
/*         the user wants to terminate execution of FDJAC1. */
/*         In this case set IFLAG to a negative integer. */

/*       N Is a positive integer input variable set to the number */
/*         of functions and variables. */

/*       X is an input array of length N. */

/*       FVEC is an input array of length N which must contain the */
/*         functions evaluated at X. */

/*       FJAC is an output N by N array which contains the */
/*         approximation to the Jacobian matrix evaluated at X. */

/*       LDFJAC is a positive integer input variable not less than N */
/*         which specifies the leading dimension of the array FJAC. */

/*       IFLAG is an integer variable which can be used to terminate */
/*         the execution of FDJAC1. See description of FCN. */

/*       ML is a nonnegative integer input variable which specifies */
/*         the number of subdiagonals within the band of the */
/*         Jacobian matrix. If the Jacobian is not banded, set */
/*         ML to at least N - 1. */

/*       EPSFCN is an input variable used in determining a suitable */
/*         step length for the forward-difference approximation. This */
/*         approximation assumes that the relative errors in the */
/*         functions are of the order of EPSFCN. If EPSFCN is less */
/*         than the machine precision, it is assumed that the relative */
/*         errors in the functions are of the order of the machine */
/*         precision. */

/*       MU is a nonnegative integer input variable which specifies */
/*         the number of superdiagonals within the band of the */
/*         Jacobian matrix. If the Jacobian is not banded, set */
/*         MU to at least N - 1. */

/*       WA1 and WA2 are work arrays of length N. If ML + MU + 1 is at */
/*         least N, then the Jacobian is considered dense, and WA2 is */
/*         not referenced. */

/* ***SEE ALSO  SNSQ, SNSQE */
/* ***ROUTINES CALLED  R1MACH */
/* ***REVISION HISTORY  (YYMMDD) */
/*   800301  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890831  Modified array declarations.  (WRB) */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900326  Removed duplicate information from DESCRIPTION section. */
/*           (WRB) */
/*   900328  Added TYPE section.  (WRB) */
/* ***END PROLOGUE  FDJAC1 */
    /* Parameter adjustments */
    --x;
    --fvec;
    fjac_dim1 = *ldfjac;
    fjac_offset = 1 + fjac_dim1;
    fjac -= fjac_offset;
    --wa1;
    --wa2;

    /* Function Body */
/* ***FIRST EXECUTABLE STATEMENT  FDJAC1 */
    epsmch = r1mach_(&c__4);

    eps = sqrt((dmax(*epsfcn,epsmch)));
    msum = *ml + *mu + 1;
    if (msum < *n) {
	goto L40;
    }

/*        COMPUTATION OF DENSE APPROXIMATE JACOBIAN. */

    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	temp = x[j];
	h__ = eps * dabs(temp);
	if (h__ == zero) {
	    h__ = eps;
	}
	x[j] = temp + h__;
	(*fcn)(n, &x[1], &wa1[1], iflag);
	if (*iflag < 0) {
	    goto L30;
	}
	x[j] = temp;
	i__2 = *n;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    fjac[i__ + j * fjac_dim1] = (wa1[i__] - fvec[i__]) / h__;
/* L10: */
	}
/* L20: */
    }
L30:
    goto L110;
L40:

/*        COMPUTATION OF BANDED APPROXIMATE JACOBIAN. */

    i__1 = msum;
    for (k = 1; k <= i__1; ++k) {
	i__2 = *n;
	i__3 = msum;
	for (j = k; i__3 < 0 ? j >= i__2 : j <= i__2; j += i__3) {
	    wa2[j] = x[j];
	    h__ = eps * (r__1 = wa2[j], dabs(r__1));
	    if (h__ == zero) {
		h__ = eps;
	    }
	    x[j] = wa2[j] + h__;
/* L60: */
	}
	(*fcn)(n, &x[1], &wa1[1], iflag);
	if (*iflag < 0) {
	    goto L100;
	}
	i__3 = *n;
	i__2 = msum;
	for (j = k; i__2 < 0 ? j >= i__3 : j <= i__3; j += i__2) {
	    x[j] = wa2[j];
	    h__ = eps * (r__1 = wa2[j], dabs(r__1));
	    if (h__ == zero) {
		h__ = eps;
	    }
	    i__4 = *n;
	    for (i__ = 1; i__ <= i__4; ++i__) {
		fjac[i__ + j * fjac_dim1] = zero;
		if (i__ >= j - *mu && i__ <= j + *ml) {
		    fjac[i__ + j * fjac_dim1] = (wa1[i__] - fvec[i__]) / h__;
		}
/* L70: */
	    }
/* L80: */
	}
/* L90: */
    }
L100:
L110:
    return 0;

/*     LAST CARD OF SUBROUTINE FDJAC1. */

} /* fdjac1_ */

