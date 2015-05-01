/* qs2i1r.f -- translated by f2c (version 12.02.01).
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

/* DECK QS2I1R */
/* Subroutine */ int qs2i1r_(integer *ia, integer *ja, real *a, integer *n, 
	integer *kflag)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__, j, k, l, m;
    static real r__;
    static integer ij;
    static real ta;
    static integer kk, il[21], nn, it, jt, iu[21], iit, jjt;
    static real tta;
    extern /* Subroutine */ int xermsg_(char *, char *, char *, integer *, 
	    integer *, ftnlen, ftnlen, ftnlen);

/* ***BEGIN PROLOGUE  QS2I1R */
/* ***SUBSIDIARY */
/* ***PURPOSE  Sort an integer array, moving an integer and real array. */
/*            This routine sorts the integer array IA and makes the same */
/*            interchanges in the integer array JA and the real array A. */
/*            The array IA may be sorted in increasing order or decreas- */
/*            ing order.  A slightly modified QUICKSORT algorithm is */
/*            used. */
/* ***LIBRARY   SLATEC (SLAP) */
/* ***CATEGORY  N6A2A */
/* ***TYPE      SINGLE PRECISION (QS2I1R-S, QS2I1D-D) */
/* ***KEYWORDS  SINGLETON QUICKSORT, SLAP, SORT, SORTING */
/* ***AUTHOR  Jones, R. E., (SNLA) */
/*           Kahaner, D. K., (NBS) */
/*           Seager, M. K., (LLNL) seager@llnl.gov */
/*           Wisniewski, J. A., (SNLA) */
/* ***DESCRIPTION */
/*     Written by Rondall E Jones */
/*     Modified by John A. Wisniewski to use the Singleton QUICKSORT */
/*     algorithm. date 18 November 1976. */

/*     Further modified by David K. Kahaner */
/*     National Bureau of Standards */
/*     August, 1981 */

/*     Even further modification made to bring the code up to the */
/*     Fortran 77 level and make it more readable and to carry */
/*     along one integer array and one real array during the sort by */
/*     Mark K. Seager */
/*     Lawrence Livermore National Laboratory */
/*     November, 1987 */
/*     This routine was adapted from the ISORT routine. */

/*     ABSTRACT */
/*         This routine sorts an integer array IA and makes the same */
/*         interchanges in the integer array JA and the real array A. */
/*         The array IA may be sorted in increasing order or decreasing */
/*         order.  A slightly modified quicksort algorithm is used. */

/*     DESCRIPTION OF PARAMETERS */
/*        IA - Integer array of values to be sorted. */
/*        JA - Integer array to be carried along. */
/*         A - Real array to be carried along. */
/*         N - Number of values in integer array IA to be sorted. */
/*     KFLAG - Control parameter */
/*           = 1 means sort IA in INCREASING order. */
/*           =-1 means sort IA in DECREASING order. */

/* ***SEE ALSO  SS2Y */
/* ***REFERENCES  R. C. Singleton, Algorithm 347, An Efficient Algorithm */
/*                 for Sorting With Minimal Storage, Communications ACM */
/*                 12:3 (1969), pp.185-7. */
/* ***ROUTINES CALLED  XERMSG */
/* ***REVISION HISTORY  (YYMMDD) */
/*   761118  DATE WRITTEN */
/*   890125  Previous REVISION DATE */
/*   890915  Made changes requested at July 1989 CML Meeting.  (MKS) */
/*   890922  Numerous changes to prologue to make closer to SLATEC */
/*           standard.  (FNF) */
/*   890929  Numerous changes to reduce SP/DP differences.  (FNF) */
/*   900805  Changed XERROR calls to calls to XERMSG.  (RWC) */
/*   910411  Prologue converted to Version 4.0 format.  (BAB) */
/*   910506  Made subsidiary to SS2Y and corrected reference.  (FNF) */
/*   920511  Added complete declaration section.  (WRB) */
/*   920929  Corrected format of reference.  (FNF) */
/*   921012  Added E0's to f.p. constants.  (FNF) */
/* ***END PROLOGUE  QS2I1R */
/* VD$R NOVECTOR */
/* VD$R NOCONCUR */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. Local Arrays .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */
/* ***FIRST EXECUTABLE STATEMENT  QS2I1R */
    /* Parameter adjustments */
    --a;
    --ja;
    --ia;

    /* Function Body */
    nn = *n;
    if (nn < 1) {
	xermsg_("SLATEC", "QS2I1R", "The number of values to be sorted was n"
		"ot positive.", &c__1, &c__1, (ftnlen)6, (ftnlen)6, (ftnlen)51)
		;
	return 0;
    }
    if (*n == 1) {
	return 0;
    }
    kk = abs(*kflag);
    if (kk != 1) {
	xermsg_("SLATEC", "QS2I1R", "The sort control parameter, K, was not "
		"1 or -1.", &c__2, &c__1, (ftnlen)6, (ftnlen)6, (ftnlen)47);
	return 0;
    }

/*     Alter array IA to get decreasing order if needed. */

    if (*kflag < 1) {
	i__1 = nn;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    ia[i__] = -ia[i__];
/* L20: */
	}
    }

/*     Sort IA and carry JA and A along. */
/*     And now...Just a little black magic... */
    m = 1;
    i__ = 1;
    j = nn;
    r__ = .375f;
L210:
    if (r__ <= .5898437f) {
	r__ += .0390625f;
    } else {
	r__ += -.21875f;
    }
L225:
    k = i__;

/*     Select a central element of the array and save it in location */
/*     it, jt, at. */

    ij = i__ + (integer) ((j - i__) * r__);
    it = ia[ij];
    jt = ja[ij];
    ta = a[ij];

/*     If first element of array is greater than it, interchange with it. */

    if (ia[i__] > it) {
	ia[ij] = ia[i__];
	ia[i__] = it;
	it = ia[ij];
	ja[ij] = ja[i__];
	ja[i__] = jt;
	jt = ja[ij];
	a[ij] = a[i__];
	a[i__] = ta;
	ta = a[ij];
    }
    l = j;

/*     If last element of array is less than it, swap with it. */

    if (ia[j] < it) {
	ia[ij] = ia[j];
	ia[j] = it;
	it = ia[ij];
	ja[ij] = ja[j];
	ja[j] = jt;
	jt = ja[ij];
	a[ij] = a[j];
	a[j] = ta;
	ta = a[ij];

/*     If first element of array is greater than it, swap with it. */

	if (ia[i__] > it) {
	    ia[ij] = ia[i__];
	    ia[i__] = it;
	    it = ia[ij];
	    ja[ij] = ja[i__];
	    ja[i__] = jt;
	    jt = ja[ij];
	    a[ij] = a[i__];
	    a[i__] = ta;
	    ta = a[ij];
	}
    }

/*     Find an element in the second half of the array which is */
/*     smaller than it. */

L240:
    --l;
    if (ia[l] > it) {
	goto L240;
    }

/*     Find an element in the first half of the array which is */
/*     greater than it. */

L245:
    ++k;
    if (ia[k] < it) {
	goto L245;
    }

/*     Interchange these elements. */

    if (k <= l) {
	iit = ia[l];
	ia[l] = ia[k];
	ia[k] = iit;
	jjt = ja[l];
	ja[l] = ja[k];
	ja[k] = jjt;
	tta = a[l];
	a[l] = a[k];
	a[k] = tta;
	goto L240;
    }

/*     Save upper and lower subscripts of the array yet to be sorted. */

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
    goto L260;

/*     Begin again on another portion of the unsorted array. */

L255:
    --m;
    if (m == 0) {
	goto L300;
    }
    i__ = il[m - 1];
    j = iu[m - 1];
L260:
    if (j - i__ >= 1) {
	goto L225;
    }
    if (i__ == j) {
	goto L255;
    }
    if (i__ == 1) {
	goto L210;
    }
    --i__;
L265:
    ++i__;
    if (i__ == j) {
	goto L255;
    }
    it = ia[i__ + 1];
    jt = ja[i__ + 1];
    ta = a[i__ + 1];
    if (ia[i__] <= it) {
	goto L265;
    }
    k = i__;
L270:
    ia[k + 1] = ia[k];
    ja[k + 1] = ja[k];
    a[k + 1] = a[k];
    --k;
    if (it < ia[k]) {
	goto L270;
    }
    ia[k + 1] = it;
    ja[k + 1] = jt;
    a[k + 1] = ta;
    goto L265;

/*     Clean up, if necessary. */

L300:
    if (*kflag < 1) {
	i__1 = nn;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    ia[i__] = -ia[i__];
/* L310: */
	}
    }
    return 0;
/* ------------- LAST LINE OF QS2I1R FOLLOWS ---------------------------- */
} /* qs2i1r_ */

