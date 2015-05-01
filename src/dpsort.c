/* dpsort.f -- translated by f2c (version 12.02.01).
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

/* DECK DPSORT */
/* Subroutine */ int dpsort_(doublereal *dx, integer *n, integer *iperm, 
	integer *kflag, integer *ier)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer i__, j, k, l, m;
    static doublereal r__;
    static integer ij, il[21], kk, lm, nn, iu[21], lmt, indx;
    static doublereal temp;
    static integer indx0, istrt;
    extern /* Subroutine */ int xermsg_(char *, char *, char *, integer *, 
	    integer *, ftnlen, ftnlen, ftnlen);

/* ***BEGIN PROLOGUE  DPSORT */
/* ***PURPOSE  Return the permutation vector generated by sorting a given */
/*            array and, optionally, rearrange the elements of the array. */
/*            The array may be sorted in increasing or decreasing order. */
/*            A slightly modified quicksort algorithm is used. */
/* ***LIBRARY   SLATEC */
/* ***CATEGORY  N6A1B, N6A2B */
/* ***TYPE      DOUBLE PRECISION (SPSORT-S, DPSORT-D, IPSORT-I, HPSORT-H) */
/* ***KEYWORDS  NUMBER SORTING, PASSIVE SORTING, SINGLETON QUICKSORT, SORT */
/* ***AUTHOR  Jones, R. E., (SNLA) */
/*           Rhoads, G. S., (NBS) */
/*           Wisniewski, J. A., (SNLA) */
/* ***DESCRIPTION */

/*   DPSORT returns the permutation vector IPERM generated by sorting */
/*   the array DX and, optionally, rearranges the values in DX.  DX may */
/*   be sorted in increasing or decreasing order.  A slightly modified */
/*   quicksort algorithm is used. */

/*   IPERM is such that DX(IPERM(I)) is the Ith value in the */
/*   rearrangement of DX.  IPERM may be applied to another array by */
/*   calling IPPERM, SPPERM, DPPERM or HPPERM. */

/*   The main difference between DPSORT and its active sorting equivalent */
/*   DSORT is that the data are referenced indirectly rather than */
/*   directly.  Therefore, DPSORT should require approximately twice as */
/*   long to execute as DSORT.  However, DPSORT is more general. */

/*   Description of Parameters */
/*      DX - input/output -- double precision array of values to be */
/*           sorted.  If ABS(KFLAG) = 2, then the values in DX will be */
/*           rearranged on output; otherwise, they are unchanged. */
/*      N  - input -- number of values in array DX to be sorted. */
/*      IPERM - output -- permutation array such that IPERM(I) is the */
/*              index of the value in the original order of the */
/*              DX array that is in the Ith location in the sorted */
/*              order. */
/*      KFLAG - input -- control parameter: */
/*            =  2  means return the permutation vector resulting from */
/*                  sorting DX in increasing order and sort DX also. */
/*            =  1  means return the permutation vector resulting from */
/*                  sorting DX in increasing order and do not sort DX. */
/*            = -1  means return the permutation vector resulting from */
/*                  sorting DX in decreasing order and do not sort DX. */
/*            = -2  means return the permutation vector resulting from */
/*                  sorting DX in decreasing order and sort DX also. */
/*      IER - output -- error indicator: */
/*          =  0  if no error, */
/*          =  1  if N is zero or negative, */
/*          =  2  if KFLAG is not 2, 1, -1, or -2. */
/* ***REFERENCES  R. C. Singleton, Algorithm 347, An efficient algorithm */
/*                 for sorting with minimal storage, Communications of */
/*                 the ACM, 12, 3 (1969), pp. 185-187. */
/* ***ROUTINES CALLED  XERMSG */
/* ***REVISION HISTORY  (YYMMDD) */
/*   761101  DATE WRITTEN */
/*   761118  Modified by John A. Wisniewski to use the Singleton */
/*           quicksort algorithm. */
/*   870423  Modified by Gregory S. Rhoads for passive sorting with the */
/*           option for the rearrangement of the original data. */
/*   890619  Double precision version of SPSORT created by D. W. Lozier. */
/*   890620  Algorithm for rearranging the data vector corrected by R. */
/*           Boisvert. */
/*   890622  Prologue upgraded to Version 4.0 style by D. Lozier. */
/*   891128  Error when KFLAG.LT.0 and N=1 corrected by R. Boisvert. */
/*   920507  Modified by M. McClain to revise prologue text. */
/*   920818  Declarations section rebuilt and code restructured to use */
/*           IF-THEN-ELSE-ENDIF.  (SMR, WRB) */
/* ***END PROLOGUE  DPSORT */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. Local Arrays .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */
/* ***FIRST EXECUTABLE STATEMENT  DPSORT */
    /* Parameter adjustments */
    --iperm;
    --dx;

    /* Function Body */
    *ier = 0;
    nn = *n;
    if (nn < 1) {
	*ier = 1;
	xermsg_("SLATEC", "DPSORT", "The number of values to be sorted, N, i"
		"s not positive.", ier, &c__1, (ftnlen)6, (ftnlen)6, (ftnlen)
		54);
	return 0;
    }

    kk = abs(*kflag);
    if (kk != 1 && kk != 2) {
	*ier = 2;
	xermsg_("SLATEC", "DPSORT", "The sort control parameter, KFLAG, is n"
		"ot 2, 1, -1, or -2.", ier, &c__1, (ftnlen)6, (ftnlen)6, (
		ftnlen)58);
	return 0;
    }

/*     Initialize permutation vector */

    i__1 = nn;
    for (i__ = 1; i__ <= i__1; ++i__) {
	iperm[i__] = i__;
/* L10: */
    }

/*     Return if only one value is to be sorted */

    if (nn == 1) {
	return 0;
    }

/*     Alter array DX to get decreasing order if needed */

    if (*kflag <= -1) {
	i__1 = nn;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    dx[i__] = -dx[i__];
/* L20: */
	}
    }

/*     Sort DX only */

    m = 1;
    i__ = 1;
    j = nn;
    r__ = .375;

L30:
    if (i__ == j) {
	goto L80;
    }
    if (r__ <= .5898437) {
	r__ += .0390625;
    } else {
	r__ += -.21875;
    }

L40:
    k = i__;

/*     Select a central element of the array and save it in location L */

    ij = i__ + (integer) ((j - i__) * r__);
    lm = iperm[ij];

/*     If first element of array is greater than LM, interchange with LM */

    if (dx[iperm[i__]] > dx[lm]) {
	iperm[ij] = iperm[i__];
	iperm[i__] = lm;
	lm = iperm[ij];
    }
    l = j;

/*     If last element of array is less than LM, interchange with LM */

    if (dx[iperm[j]] < dx[lm]) {
	iperm[ij] = iperm[j];
	iperm[j] = lm;
	lm = iperm[ij];

/*        If first element of array is greater than LM, interchange */
/*        with LM */

	if (dx[iperm[i__]] > dx[lm]) {
	    iperm[ij] = iperm[i__];
	    iperm[i__] = lm;
	    lm = iperm[ij];
	}
    }
    goto L60;
L50:
    lmt = iperm[l];
    iperm[l] = iperm[k];
    iperm[k] = lmt;

/*     Find an element in the second half of the array which is smaller */
/*     than LM */

L60:
    --l;
    if (dx[iperm[l]] > dx[lm]) {
	goto L60;
    }

/*     Find an element in the first half of the array which is greater */
/*     than LM */

L70:
    ++k;
    if (dx[iperm[k]] < dx[lm]) {
	goto L70;
    }

/*     Interchange these elements */

    if (k <= l) {
	goto L50;
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
    goto L90;

/*     Begin again on another portion of the unsorted array */

L80:
    --m;
    if (m == 0) {
	goto L120;
    }
    i__ = il[m - 1];
    j = iu[m - 1];

L90:
    if (j - i__ >= 1) {
	goto L40;
    }
    if (i__ == 1) {
	goto L30;
    }
    --i__;

L100:
    ++i__;
    if (i__ == j) {
	goto L80;
    }
    lm = iperm[i__ + 1];
    if (dx[iperm[i__]] <= dx[lm]) {
	goto L100;
    }
    k = i__;

L110:
    iperm[k + 1] = iperm[k];
    --k;
    if (dx[lm] < dx[iperm[k]]) {
	goto L110;
    }
    iperm[k + 1] = lm;
    goto L100;

/*     Clean up */

L120:
    if (*kflag <= -1) {
	i__1 = nn;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    dx[i__] = -dx[i__];
/* L130: */
	}
    }

/*     Rearrange the values of DX if desired */

    if (kk == 2) {

/*        Use the IPERM vector as a flag. */
/*        If IPERM(I) < 0, then the I-th value is in correct location */

	i__1 = nn;
	for (istrt = 1; istrt <= i__1; ++istrt) {
	    if (iperm[istrt] >= 0) {
		indx = istrt;
		indx0 = indx;
		temp = dx[istrt];
L140:
		if (iperm[indx] > 0) {
		    dx[indx] = dx[iperm[indx]];
		    indx0 = indx;
		    iperm[indx] = -iperm[indx];
		    indx = (i__2 = iperm[indx], abs(i__2));
		    goto L140;
		}
		dx[indx0] = temp;
	    }
/* L150: */
	}

/*        Revert the signs of the IPERM values */

	i__1 = nn;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    iperm[i__] = -iperm[i__];
/* L160: */
	}

    }

    return 0;
} /* dpsort_ */

