/* fdjac3.f -- translated by f2c (version 12.02.01).
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

/* DECK FDJAC3 */
/* Subroutine */ int fdjac3_(S_fp fcn, integer *m, integer *n, real *x, real *
	fvec, real *fjac, integer *ldfjac, integer *iflag, real *epsfcn, real 
	*wa)
{
    /* Initialized data */

    static real zero = 0.f;

    /* System generated locals */
    integer fjac_dim1, fjac_offset, i__1, i__2;

    /* Local variables */
    static real h__;
    static integer i__, j;
    static real eps, temp;
    extern doublereal r1mach_(integer *);
    static real epsmch;

/* ***BEGIN PROLOGUE  FDJAC3 */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to SNLS1 and SNLS1E */
/* ***LIBRARY   SLATEC */
/* ***TYPE      SINGLE PRECISION (FDJAC3-S, DFDJC3-D) */
/* ***AUTHOR  (UNKNOWN) */
/* ***DESCRIPTION */

/*     This subroutine computes a forward-difference approximation */
/*     to the M by N Jacobian matrix associated with a specified */
/*     problem of M functions in N variables. */

/*     The subroutine statement is */

/*       SUBROUTINE FDJAC3(FCN,M,N,X,FVEC,FJAC,LDFJAC,IFLAG,EPSFCN,WA) */

/*     where */

/*       FCN is the name of the user-supplied subroutine which */
/*         calculates the functions. FCN must be declared */
/*         in an external statement in the user calling */
/*         program, and should be written as follows. */

/*         SUBROUTINE FCN(IFLAG,M,N,X,FVEC,FJAC,LDFJAC) */
/*         INTEGER LDFJAC,M,N,IFLAG */
/*         REAL X(N),FVEC(M),FJAC(LDFJAC,N) */
/*         ---------- */
/*         When IFLAG.EQ.1 calculate the functions at X and */
/*         return this vector in FVEC. */
/*         ---------- */
/*         RETURN */
/*         END */

/*         The value of IFLAG should not be changed by FCN unless */
/*         the user wants to terminate execution of FDJAC3. */
/*         In this case set IFLAG to a negative integer. */

/*       M is a positive integer input variable set to the number */
/*         of functions. */

/*       N is a positive integer input variable set to the number */
/*         of variables. N must not exceed M. */

/*       X is an input array of length N. */

/*       FVEC is an input array of length M which must contain the */
/*         functions evaluated at X. */

/*       FJAC is an output M by N array which contains the */
/*         approximation to the Jacobian matrix evaluated at X. */

/*       LDFJAC is a positive integer input variable not less than M */
/*         which specifies the leading dimension of the array FJAC. */

/*       IFLAG is an integer variable which can be used to terminate */
/*         THE EXECUTION OF FDJAC3. See description of FCN. */

/*       EPSFCN is an input variable used in determining a suitable */
/*         step length for the forward-difference approximation. This */
/*         approximation assumes that the relative errors in the */
/*         functions are of the order of EPSFCN. If EPSFCN is less */
/*         than the machine precision, it is assumed that the relative */
/*         errors in the functions are of the order of the machine */
/*         precision. */

/*       WA is a work array of length M. */

/* ***SEE ALSO  SNLS1, SNLS1E */
/* ***ROUTINES CALLED  R1MACH */
/* ***REVISION HISTORY  (YYMMDD) */
/*   800301  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890831  Modified array declarations.  (WRB) */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900326  Removed duplicate information from DESCRIPTION section. */
/*           (WRB) */
/*   900328  Added TYPE section.  (WRB) */
/* ***END PROLOGUE  FDJAC3 */
    /* Parameter adjustments */
    --x;
    --fvec;
    fjac_dim1 = *ldfjac;
    fjac_offset = 1 + fjac_dim1;
    fjac -= fjac_offset;
    --wa;

    /* Function Body */
/* ***FIRST EXECUTABLE STATEMENT  FDJAC3 */
    epsmch = r1mach_(&c__4);

    eps = sqrt((dmax(*epsfcn,epsmch)));
/*      SET IFLAG=1 TO INDICATE THAT FUNCTION VALUES */
/*           ARE TO BE RETURNED BY FCN. */
    *iflag = 1;
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	temp = x[j];
	h__ = eps * dabs(temp);
	if (h__ == zero) {
	    h__ = eps;
	}
	x[j] = temp + h__;
	(*fcn)(iflag, m, n, &x[1], &wa[1], &fjac[fjac_offset], ldfjac);
	if (*iflag < 0) {
	    goto L30;
	}
	x[j] = temp;
	i__2 = *m;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    fjac[i__ + j * fjac_dim1] = (wa[i__] - fvec[i__]) / h__;
/* L10: */
	}
/* L20: */
    }
L30:
    return 0;

/*     LAST CARD OF SUBROUTINE FDJAC3. */

} /* fdjac3_ */

