/* qelg.f -- translated by f2c (version 12.02.01).
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
static integer c__2 = 2;

/* DECK QELG */
/* Subroutine */ int qelg_(integer *n, real *epstab, real *result, real *
	abserr, real *res3la, integer *nres)
{
    /* System generated locals */
    integer i__1;
    real r__1, r__2, r__3;

    /* Local variables */
    static integer i__;
    static real e0, e1, e2, e3;
    static integer k1, k2, k3, ib, ie;
    static real ss;
    static integer ib2;
    static real res;
    static integer num;
    static real err1, err2, err3, tol1, tol2, tol3;
    static integer indx;
    static real e1abs, oflow, error, delta1, delta2, delta3;
    extern doublereal r1mach_(integer *);
    static real epmach, epsinf;
    static integer newelm, limexp;

/* ***BEGIN PROLOGUE  QELG */
/* ***SUBSIDIARY */
/* ***PURPOSE  The routine determines the limit of a given sequence of */
/*            approximations, by means of the Epsilon algorithm of */
/*            P. Wynn. An estimate of the absolute error is also given. */
/*            The condensed Epsilon table is computed. Only those */
/*            elements needed for the computation of the next diagonal */
/*            are preserved. */
/* ***LIBRARY   SLATEC */
/* ***TYPE      SINGLE PRECISION (QELG-S, DQELG-D) */
/* ***KEYWORDS  CONVERGENCE ACCELERATION, EPSILON ALGORITHM, EXTRAPOLATION */
/* ***AUTHOR  Piessens, Robert */
/*             Applied Mathematics and Programming Division */
/*             K. U. Leuven */
/*           de Doncker, Elise */
/*             Applied Mathematics and Programming Division */
/*             K. U. Leuven */
/* ***DESCRIPTION */

/*           Epsilon algorithm */
/*           Standard fortran subroutine */
/*           Real version */

/*           PARAMETERS */
/*              N      - Integer */
/*                       EPSTAB(N) contains the new element in the */
/*                       first column of the epsilon table. */

/*              EPSTAB - Real */
/*                       Vector of dimension 52 containing the elements */
/*                       of the two lower diagonals of the triangular */
/*                       epsilon table. The elements are numbered */
/*                       starting at the right-hand corner of the */
/*                       triangle. */

/*              RESULT - Real */
/*                       Resulting approximation to the integral */

/*              ABSERR - Real */
/*                       Estimate of the absolute error computed from */
/*                       RESULT and the 3 previous results */

/*              RES3LA - Real */
/*                       Vector of dimension 3 containing the last 3 */
/*                       results */

/*              NRES   - Integer */
/*                       Number of calls to the routine */
/*                       (should be zero at first call) */

/* ***SEE ALSO  QAGIE, QAGOE, QAGPE, QAGSE */
/* ***ROUTINES CALLED  R1MACH */
/* ***REVISION HISTORY  (YYMMDD) */
/*   800101  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890531  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900328  Added TYPE section.  (WRB) */
/* ***END PROLOGUE  QELG */


/*           LIST OF MAJOR VARIABLES */
/*           ----------------------- */

/*           E0     - THE 4 ELEMENTS ON WHICH THE */
/*           E1       COMPUTATION OF A NEW ELEMENT IN */
/*           E2       THE EPSILON TABLE IS BASED */
/*           E3                 E0 */
/*                        E3    E1    NEW */
/*                              E2 */
/*           NEWELM - NUMBER OF ELEMENTS TO BE COMPUTED IN THE NEW */
/*                    DIAGONAL */
/*           ERROR  - ERROR = ABS(E1-E0)+ABS(E2-E1)+ABS(NEW-E2) */
/*           RESULT - THE ELEMENT IN THE NEW DIAGONAL WITH LEAST VALUE */
/*                    OF ERROR */

/*           MACHINE DEPENDENT CONSTANTS */
/*           --------------------------- */

/*           EPMACH IS THE LARGEST RELATIVE SPACING. */
/*           OFLOW IS THE LARGEST POSITIVE MAGNITUDE. */
/*           LIMEXP IS THE MAXIMUM NUMBER OF ELEMENTS THE EPSILON */
/*           TABLE CAN CONTAIN. IF THIS NUMBER IS REACHED, THE UPPER */
/*           DIAGONAL OF THE EPSILON TABLE IS DELETED. */

/* ***FIRST EXECUTABLE STATEMENT  QELG */
    /* Parameter adjustments */
    --res3la;
    --epstab;

    /* Function Body */
    epmach = r1mach_(&c__4);
    oflow = r1mach_(&c__2);
    ++(*nres);
    *abserr = oflow;
    *result = epstab[*n];
    if (*n < 3) {
	goto L100;
    }
    limexp = 50;
    epstab[*n + 2] = epstab[*n];
    newelm = (*n - 1) / 2;
    epstab[*n] = oflow;
    num = *n;
    k1 = *n;
    i__1 = newelm;
    for (i__ = 1; i__ <= i__1; ++i__) {
	k2 = k1 - 1;
	k3 = k1 - 2;
	res = epstab[k1 + 2];
	e0 = epstab[k3];
	e1 = epstab[k2];
	e2 = res;
	e1abs = dabs(e1);
	delta2 = e2 - e1;
	err2 = dabs(delta2);
/* Computing MAX */
	r__1 = dabs(e2);
	tol2 = dmax(r__1,e1abs) * epmach;
	delta3 = e1 - e0;
	err3 = dabs(delta3);
/* Computing MAX */
	r__1 = e1abs, r__2 = dabs(e0);
	tol3 = dmax(r__1,r__2) * epmach;
	if (err2 > tol2 || err3 > tol3) {
	    goto L10;
	}

/*           IF E0, E1 AND E2 ARE EQUAL TO WITHIN MACHINE */
/*           ACCURACY, CONVERGENCE IS ASSUMED. */
/*           RESULT = E2 */
/*           ABSERR = ABS(E1-E0)+ABS(E2-E1) */

	*result = res;
	*abserr = err2 + err3;
/* ***JUMP OUT OF DO-LOOP */
	goto L100;
L10:
	e3 = epstab[k1];
	epstab[k1] = e1;
	delta1 = e1 - e3;
	err1 = dabs(delta1);
/* Computing MAX */
	r__1 = e1abs, r__2 = dabs(e3);
	tol1 = dmax(r__1,r__2) * epmach;

/*           IF TWO ELEMENTS ARE VERY CLOSE TO EACH OTHER, OMIT */
/*           A PART OF THE TABLE BY ADJUSTING THE VALUE OF N */

	if (err1 <= tol1 || err2 <= tol2 || err3 <= tol3) {
	    goto L20;
	}
	ss = 1.f / delta1 + 1.f / delta2 - 1.f / delta3;
	epsinf = (r__1 = ss * e1, dabs(r__1));

/*           TEST TO DETECT IRREGULAR BEHAVIOUR IN THE TABLE, AND */
/*           EVENTUALLY OMIT A PART OF THE TABLE ADJUSTING THE VALUE */
/*           OF N. */

	if (epsinf > 1e-4f) {
	    goto L30;
	}
L20:
	*n = i__ + i__ - 1;
/* ***JUMP OUT OF DO-LOOP */
	goto L50;

/*           COMPUTE A NEW ELEMENT AND EVENTUALLY ADJUST */
/*           THE VALUE OF RESULT. */

L30:
	res = e1 + 1.f / ss;
	epstab[k1] = res;
	k1 += -2;
	error = err2 + (r__1 = res - e2, dabs(r__1)) + err3;
	if (error > *abserr) {
	    goto L40;
	}
	*abserr = error;
	*result = res;
L40:
	;
    }

/*           SHIFT THE TABLE. */

L50:
    if (*n == limexp) {
	*n = (limexp / 2 << 1) - 1;
    }
    ib = 1;
    if (num / 2 << 1 == num) {
	ib = 2;
    }
    ie = newelm + 1;
    i__1 = ie;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ib2 = ib + 2;
	epstab[ib] = epstab[ib2];
	ib = ib2;
/* L60: */
    }
    if (num == *n) {
	goto L80;
    }
    indx = num - *n + 1;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	epstab[i__] = epstab[indx];
	++indx;
/* L70: */
    }
L80:
    if (*nres >= 4) {
	goto L90;
    }
    res3la[*nres] = *result;
    *abserr = oflow;
    goto L100;

/*           COMPUTE ERROR ESTIMATE */

L90:
    *abserr = (r__1 = *result - res3la[3], dabs(r__1)) + (r__2 = *result - 
	    res3la[2], dabs(r__2)) + (r__3 = *result - res3la[1], dabs(r__3));
    res3la[1] = res3la[2];
    res3la[2] = res3la[3];
    res3la[3] = *result;
L100:
/* Computing MAX */
    r__1 = *abserr, r__2 = epmach * 5.f * dabs(*result);
    *abserr = dmax(r__1,r__2);
    return 0;
} /* qelg_ */

