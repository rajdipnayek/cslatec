/* chkder.f -- translated by f2c (version 12.02.01).
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

/* DECK CHKDER */
/* Subroutine */ int chkder_(integer *m, integer *n, real *x, real *fvec, 
	real *fjac, integer *ldfjac, real *xp, real *fvecp, integer *mode, 
	real *err)
{
    /* Initialized data */

    static real factor = 100.f;
    static real one = 1.f;
    static real zero = 0.f;

    /* System generated locals */
    integer fjac_dim1, fjac_offset, i__1, i__2;
    real r__1, r__2, r__3, r__4, r__5;

    /* Local variables */
    static integer i__, j;
    static real eps, epsf, temp;
    extern doublereal r1mach_(integer *);
    static real epsmch, epslog;

/* ***BEGIN PROLOGUE  CHKDER */
/* ***PURPOSE  Check the gradients of M nonlinear functions in N */
/*            variables, evaluated at a point X, for consistency */
/*            with the functions themselves. */
/* ***LIBRARY   SLATEC */
/* ***CATEGORY  F3, G4C */
/* ***TYPE      SINGLE PRECISION (CHKDER-S, DCKDER-D) */
/* ***KEYWORDS  GRADIENTS, JACOBIAN, MINPACK, NONLINEAR */
/* ***AUTHOR  Hiebert, K. L. (SNLA) */
/* ***DESCRIPTION */

/*   This subroutine is a companion routine to SNLS1,SNLS1E,SNSQ,and */
/*   SNSQE which may be used to check the calculation of the Jacobian. */

/*     SUBROUTINE CHKDER */

/*     This subroutine checks the gradients of M nonlinear functions */
/*     in N variables, evaluated at a point X, for consistency with */
/*     the functions themselves. The user must call CKDER twice, */
/*     first with MODE = 1 and then with MODE = 2. */

/*     MODE = 1. On input, X must contain the point of evaluation. */
/*               On output, XP is set to a neighboring point. */

/*     MODE = 2. On input, FVEC must contain the functions and the */
/*                         rows of FJAC must contain the gradients */
/*                         of the respective functions each evaluated */
/*                         at X, and FVECP must contain the functions */
/*                         evaluated at XP. */
/*               On output, ERR contains measures of correctness of */
/*                          the respective gradients. */

/*     The subroutine does not perform reliably if cancellation or */
/*     rounding errors cause a severe loss of significance in the */
/*     evaluation of a function. Therefore, none of the components */
/*     of X should be unusually small (in particular, zero) or any */
/*     other value which may cause loss of significance. */

/*     The SUBROUTINE statement is */

/*       SUBROUTINE CHKDER(M,N,X,FVEC,FJAC,LDFJAC,XP,FVECP,MODE,ERR) */

/*     where */

/*       M is a positive integer input variable set to the number */
/*         of functions. */

/*       N is a positive integer input variable set to the number */
/*         of variables. */

/*       X is an input array of length N. */

/*       FVEC is an array of length M. On input when MODE = 2, */
/*         FVEC must contain the functions evaluated at X. */

/*       FJAC is an M by N array. On input when MODE = 2, */
/*         the rows of FJAC must contain the gradients of */
/*         the respective functions evaluated at X. */

/*       LDFJAC is a positive integer input parameter not less than M */
/*         which specifies the leading dimension of the array FJAC. */

/*       XP is an array of length N. On output when MODE = 1, */
/*         XP is set to a neighboring point of X. */

/*       FVECP is an array of length M. On input when MODE = 2, */
/*         FVECP must contain the functions evaluated at XP. */

/*       MODE is an integer input variable set to 1 on the first call */
/*         and 2 on the second. Other values of MODE are equivalent */
/*         to MODE = 1. */

/*       ERR is an array of length M. On output when MODE = 2, */
/*         ERR contains measures of correctness of the respective */
/*         gradients. If there is no severe loss of significance, */
/*         then if ERR(I) is 1.0 the I-th gradient is correct, */
/*         while if ERR(I) is 0.0 the I-th gradient is incorrect. */
/*         For values of ERR between 0.0 and 1.0, the categorization */
/*         is less certain. In general, a value of ERR(I) greater */
/*         than 0.5 indicates that the I-th gradient is probably */
/*         correct, while a value of ERR(I) less than 0.5 indicates */
/*         that the I-th gradient is probably incorrect. */

/* ***REFERENCES  M. J. D. Powell, A hybrid method for nonlinear equa- */
/*                 tions. In Numerical Methods for Nonlinear Algebraic */
/*                 Equations, P. Rabinowitz, Editor.  Gordon and Breach, */
/*                 1988. */
/* ***ROUTINES CALLED  R1MACH */
/* ***REVISION HISTORY  (YYMMDD) */
/*   800301  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890831  Modified array declarations.  (WRB) */
/*   890831  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900326  Removed duplicate information from DESCRIPTION section. */
/*           (WRB) */
/*   920501  Reformatted the REFERENCES section.  (WRB) */
/* ***END PROLOGUE  CHKDER */

    /* Parameter adjustments */
    --x;
    --fvec;
    fjac_dim1 = *ldfjac;
    fjac_offset = 1 + fjac_dim1;
    fjac -= fjac_offset;
    --xp;
    --fvecp;
    --err;

    /* Function Body */
/* ***FIRST EXECUTABLE STATEMENT  CHKDER */
    epsmch = r1mach_(&c__4);

    eps = sqrt(epsmch);

    if (*mode == 2) {
	goto L20;
    }

/*        MODE = 1. */

    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	temp = eps * (r__1 = x[j], dabs(r__1));
	if (temp == zero) {
	    temp = eps;
	}
	xp[j] = x[j] + temp;
/* L10: */
    }
    goto L70;
L20:

/*        MODE = 2. */

    epsf = factor * epsmch;
    epslog = r_lg10(&eps);
    i__1 = *m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	err[i__] = zero;
/* L30: */
    }
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	temp = (r__1 = x[j], dabs(r__1));
	if (temp == zero) {
	    temp = one;
	}
	i__2 = *m;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    err[i__] += temp * fjac[i__ + j * fjac_dim1];
/* L40: */
	}
/* L50: */
    }
    i__1 = *m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	temp = one;
	if (fvec[i__] != zero && fvecp[i__] != zero && (r__2 = fvecp[i__] - 
		fvec[i__], dabs(r__2)) >= epsf * (r__1 = fvec[i__], dabs(r__1)
		)) {
	    temp = eps * (r__3 = (fvecp[i__] - fvec[i__]) / eps - err[i__], 
		    dabs(r__3)) / ((r__4 = fvec[i__], dabs(r__4)) + (r__5 = 
		    fvecp[i__], dabs(r__5)));
	}
	err[i__] = one;
	if (temp > epsmch && temp < eps) {
	    err[i__] = (r_lg10(&temp) - epslog) / epslog;
	}
	if (temp >= eps) {
	    err[i__] = zero;
	}
/* L60: */
    }
L70:

    return 0;

/*     LAST CARD OF SUBROUTINE CHKDER. */

} /* chkder_ */

