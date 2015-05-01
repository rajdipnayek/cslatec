/* ssort.f -- translated by f2c (version 12.02.01).
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

/* DECK SSORT */
/* Subroutine */ int ssort_(real *x, real *y, integer *n, integer *kflag)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__, j, k, l, m;
    static real r__, t;
    static integer ij, il[21], kk, nn, iu[21];
    static real tt, ty, tty;
    extern /* Subroutine */ int xermsg_(char *, char *, char *, integer *, 
	    integer *, ftnlen, ftnlen, ftnlen);

/* ***BEGIN PROLOGUE  SSORT */
/* ***PURPOSE  Sort an array and optionally make the same interchanges in */
/*            an auxiliary array.  The array may be sorted in increasing */
/*            or decreasing order.  A slightly modified QUICKSORT */
/*            algorithm is used. */
/* ***LIBRARY   SLATEC */
/* ***CATEGORY  N6A2B */
/* ***TYPE      SINGLE PRECISION (SSORT-S, DSORT-D, ISORT-I) */
/* ***KEYWORDS  SINGLETON QUICKSORT, SORT, SORTING */
/* ***AUTHOR  Jones, R. E., (SNLA) */
/*           Wisniewski, J. A., (SNLA) */
/* ***DESCRIPTION */

/*   SSORT sorts array X and optionally makes the same interchanges in */
/*   array Y.  The array X may be sorted in increasing order or */
/*   decreasing order.  A slightly modified quicksort algorithm is used. */

/*   Description of Parameters */
/*      X - array of values to be sorted   (usually abscissas) */
/*      Y - array to be (optionally) carried along */
/*      N - number of values in array X to be sorted */
/*      KFLAG - control parameter */
/*            =  2  means sort X in increasing order and carry Y along. */
/*            =  1  means sort X in increasing order (ignoring Y) */
/*            = -1  means sort X in decreasing order (ignoring Y) */
/*            = -2  means sort X in decreasing order and carry Y along. */

/* ***REFERENCES  R. C. Singleton, Algorithm 347, An efficient algorithm */
/*                 for sorting with minimal storage, Communications of */
/*                 the ACM, 12, 3 (1969), pp. 185-187. */
/* ***ROUTINES CALLED  XERMSG */
/* ***REVISION HISTORY  (YYMMDD) */
/*   761101  DATE WRITTEN */
/*   761118  Modified to use the Singleton quicksort algorithm.  (JAW) */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890831  Modified array declarations.  (WRB) */
/*   891009  Removed unreferenced statement labels.  (WRB) */
/*   891024  Changed category.  (WRB) */
/*   891024  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ) */
/*   901012  Declared all variables; changed X,Y to SX,SY. (M. McClain) */
/*   920501  Reformatted the REFERENCES section.  (DWL, WRB) */
/*   920519  Clarified error messages.  (DWL) */
/*   920801  Declarations section rebuilt and code restructured to use */
/*           IF-THEN-ELSE-ENDIF.  (RWC, WRB) */
/* ***END PROLOGUE  SSORT */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. Local Arrays .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */
/* ***FIRST EXECUTABLE STATEMENT  SSORT */
    /* Parameter adjustments */
    --y;
    --x;

    /* Function Body */
    nn = *n;
    if (nn < 1) {
	xermsg_("SLATEC", "SSORT", "The number of values to be sorted is not"
		" positive.", &c__1, &c__1, (ftnlen)6, (ftnlen)5, (ftnlen)50);
	return 0;
    }

    kk = abs(*kflag);
    if (kk != 1 && kk != 2) {
	xermsg_("SLATEC", "SSORT", "The sort control parameter, K, is not 2,"
		" 1, -1, or -2.", &c__2, &c__1, (ftnlen)6, (ftnlen)5, (ftnlen)
		54);
	return 0;
    }

/*     Alter array X to get decreasing order if needed */

    if (*kflag <= -1) {
	i__1 = nn;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    x[i__] = -x[i__];
/* L10: */
	}
    }

    if (kk == 2) {
	goto L100;
    }

/*     Sort X only */

    m = 1;
    i__ = 1;
    j = nn;
    r__ = .375f;

L20:
    if (i__ == j) {
	goto L60;
    }
    if (r__ <= .5898437f) {
	r__ += .0390625f;
    } else {
	r__ += -.21875f;
    }

L30:
    k = i__;

/*     Select a central element of the array and save it in location T */

    ij = i__ + (integer) ((j - i__) * r__);
    t = x[ij];

/*     If first element of array is greater than T, interchange with T */

    if (x[i__] > t) {
	x[ij] = x[i__];
	x[i__] = t;
	t = x[ij];
    }
    l = j;

/*     If last element of array is less than than T, interchange with T */

    if (x[j] < t) {
	x[ij] = x[j];
	x[j] = t;
	t = x[ij];

/*        If first element of array is greater than T, interchange with T */

	if (x[i__] > t) {
	    x[ij] = x[i__];
	    x[i__] = t;
	    t = x[ij];
	}
    }

/*     Find an element in the second half of the array which is smaller */
/*     than T */

L40:
    --l;
    if (x[l] > t) {
	goto L40;
    }

/*     Find an element in the first half of the array which is greater */
/*     than T */

L50:
    ++k;
    if (x[k] < t) {
	goto L50;
    }

/*     Interchange these elements */

    if (k <= l) {
	tt = x[l];
	x[l] = x[k];
	x[k] = tt;
	goto L40;
    }

/*     Save upper and lower subscripts of the array yet to be sorted */

    if (l - i__ > j - k) {
	il[m - 1] = i__;
	iu[m - 1] = l;
	i__ = k;
	++m;
    } else {
	il[m - 1] = k;
	iu[m - 1] = j;
	j = l;
	++m;
    }
    goto L70;

/*     Begin again on another portion of the unsorted array */

L60:
    --m;
    if (m == 0) {
	goto L190;
    }
    i__ = il[m - 1];
    j = iu[m - 1];

L70:
    if (j - i__ >= 1) {
	goto L30;
    }
    if (i__ == 1) {
	goto L20;
    }
    --i__;

L80:
    ++i__;
    if (i__ == j) {
	goto L60;
    }
    t = x[i__ + 1];
    if (x[i__] <= t) {
	goto L80;
    }
    k = i__;

L90:
    x[k + 1] = x[k];
    --k;
    if (t < x[k]) {
	goto L90;
    }
    x[k + 1] = t;
    goto L80;

/*     Sort X and carry Y along */

L100:
    m = 1;
    i__ = 1;
    j = nn;
    r__ = .375f;

L110:
    if (i__ == j) {
	goto L150;
    }
    if (r__ <= .5898437f) {
	r__ += .0390625f;
    } else {
	r__ += -.21875f;
    }

L120:
    k = i__;

/*     Select a central element of the array and save it in location T */

    ij = i__ + (integer) ((j - i__) * r__);
    t = x[ij];
    ty = y[ij];

/*     If first element of array is greater than T, interchange with T */

    if (x[i__] > t) {
	x[ij] = x[i__];
	x[i__] = t;
	t = x[ij];
	y[ij] = y[i__];
	y[i__] = ty;
	ty = y[ij];
    }
    l = j;

/*     If last element of array is less than T, interchange with T */

    if (x[j] < t) {
	x[ij] = x[j];
	x[j] = t;
	t = x[ij];
	y[ij] = y[j];
	y[j] = ty;
	ty = y[ij];

/*        If first element of array is greater than T, interchange with T */

	if (x[i__] > t) {
	    x[ij] = x[i__];
	    x[i__] = t;
	    t = x[ij];
	    y[ij] = y[i__];
	    y[i__] = ty;
	    ty = y[ij];
	}
    }

/*     Find an element in the second half of the array which is smaller */
/*     than T */

L130:
    --l;
    if (x[l] > t) {
	goto L130;
    }

/*     Find an element in the first half of the array which is greater */
/*     than T */

L140:
    ++k;
    if (x[k] < t) {
	goto L140;
    }

/*     Interchange these elements */

    if (k <= l) {
	tt = x[l];
	x[l] = x[k];
	x[k] = tt;
	tty = y[l];
	y[l] = y[k];
	y[k] = tty;
	goto L130;
    }

/*     Save upper and lower subscripts of the array yet to be sorted */

    if (l - i__ > j - k) {
	il[m - 1] = i__;
	iu[m - 1] = l;
	i__ = k;
	++m;
    } else {
	il[m - 1] = k;
	iu[m - 1] = j;
	j = l;
	++m;
    }
    goto L160;

/*     Begin again on another portion of the unsorted array */

L150:
    --m;
    if (m == 0) {
	goto L190;
    }
    i__ = il[m - 1];
    j = iu[m - 1];

L160:
    if (j - i__ >= 1) {
	goto L120;
    }
    if (i__ == 1) {
	goto L110;
    }
    --i__;

L170:
    ++i__;
    if (i__ == j) {
	goto L150;
    }
    t = x[i__ + 1];
    ty = y[i__ + 1];
    if (x[i__] <= t) {
	goto L170;
    }
    k = i__;

L180:
    x[k + 1] = x[k];
    y[k + 1] = y[k];
    --k;
    if (t < x[k]) {
	goto L180;
    }
    x[k + 1] = t;
    y[k + 1] = ty;
    goto L170;

/*     Clean up */

L190:
    if (*kflag <= -1) {
	i__1 = nn;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    x[i__] = -x[i__];
/* L200: */
	}
    }
    return 0;
} /* ssort_ */

