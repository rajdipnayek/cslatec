/* cmpcsg.f -- translated by f2c (version 12.02.01).
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

/* DECK CMPCSG */
/* Subroutine */ int cmpcsg_(integer *n, integer *ijump, real *fnum, real *
	fden, complex *a)
{
    /* System generated locals */
    integer i__1, i__2, i__3;
    real r__1;
    complex q__1;

    /* Local variables */
    static integer i__, k;
    static real x, y;
    static integer k1, k2, k3, k4, k5;
    static real pi;
    static integer np1;
    static real dum, pibyn;
    extern doublereal pimach_(real *);

/* ***BEGIN PROLOGUE  CMPCSG */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to CMGNBN */
/* ***LIBRARY   SLATEC */
/* ***TYPE      COMPLEX (COSGEN-S, CMPCSG-C) */
/* ***AUTHOR  (UNKNOWN) */
/* ***DESCRIPTION */

/*     This subroutine computes required cosine values in ascending */
/*     order.  When IJUMP .GT. 1 the routine computes values */

/*        2*COS(J*PI/L) , J=1,2,...,L and J .NE. 0(MOD N/IJUMP+1) */

/*     where L = IJUMP*(N/IJUMP+1). */


/*     when IJUMP = 1 it computes */

/*            2*COS((J-FNUM)*PI/(N+FDEN)) ,  J=1, 2, ... ,N */

/*     where */
/*        FNUM = 0.5, FDEN = 0.0,  for regular reduction values. */
/*        FNUM = 0.0, FDEN = 1.0, for B-R and C-R when ISTAG = 1 */
/*        FNUM = 0.0, FDEN = 0.5, for B-R and C-R when ISTAG = 2 */
/*        FNUM = 0.5, FDEN = 0.5, for B-R and C-R when ISTAG = 2 */
/*                                in CMPOSN only. */

/* ***SEE ALSO  CMGNBN */
/* ***ROUTINES CALLED  PIMACH */
/* ***REVISION HISTORY  (YYMMDD) */
/*   801001  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900402  Added TYPE section.  (WRB) */
/* ***END PROLOGUE  CMPCSG */


/* ***FIRST EXECUTABLE STATEMENT  CMPCSG */
    /* Parameter adjustments */
    --a;

    /* Function Body */
    pi = pimach_(&dum);
    if (*n == 0) {
	goto L105;
    }
    if (*ijump == 1) {
	goto L103;
    }
    k3 = *n / *ijump + 1;
    k4 = k3 - 1;
    pibyn = pi / (*n + *ijump);
    i__1 = *ijump;
    for (k = 1; k <= i__1; ++k) {
	k1 = (k - 1) * k3;
	k5 = (k - 1) * k4;
	i__2 = k4;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    x = (real) (k1 + i__);
	    k2 = k5 + i__;
	    i__3 = k2;
	    r__1 = cos(x * pibyn) * -2.f;
	    q__1.r = r__1, q__1.i = 0.f;
	    a[i__3].r = q__1.r, a[i__3].i = q__1.i;
/* L101: */
	}
/* L102: */
    }
    goto L105;
L103:
    np1 = *n + 1;
    y = pi / (*n + *fden);
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	x = np1 - i__ - *fnum;
	i__2 = i__;
	r__1 = cos(x * y) * 2.f;
	q__1.r = r__1, q__1.i = 0.f;
	a[i__2].r = q__1.r, a[i__2].i = q__1.i;
/* L104: */
    }
L105:
    return 0;
} /* cmpcsg_ */

