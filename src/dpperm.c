/* dpperm.f -- translated by f2c (version 12.02.01).
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

/* DECK DPPERM */
/* Subroutine */ int dpperm_(doublereal *dx, integer *n, integer *iperm, 
	integer *ier)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer i__, indx, indx0;
    static doublereal dtemp;
    static integer istrt;
    extern /* Subroutine */ int xermsg_(char *, char *, char *, integer *, 
	    integer *, ftnlen, ftnlen, ftnlen);

/* ***BEGIN PROLOGUE  DPPERM */
/* ***PURPOSE  Rearrange a given array according to a prescribed */
/*            permutation vector. */
/* ***LIBRARY   SLATEC */
/* ***CATEGORY  N8 */
/* ***TYPE      DOUBLE PRECISION (SPPERM-S, DPPERM-D, IPPERM-I, HPPERM-H) */
/* ***KEYWORDS  PERMUTATION, REARRANGEMENT */
/* ***AUTHOR  McClain, M. A., (NIST) */
/*           Rhoads, G. S., (NBS) */
/* ***DESCRIPTION */

/*         DPPERM rearranges the data vector DX according to the */
/*         permutation IPERM: DX(I) <--- DX(IPERM(I)).  IPERM could come */
/*         from one of the sorting routines IPSORT, SPSORT, DPSORT or */
/*         HPSORT. */

/*     Description of Parameters */
/*         DX - input/output -- double precision array of values to be */
/*                   rearranged. */
/*         N - input -- number of values in double precision array DX. */
/*         IPERM - input -- permutation vector. */
/*         IER - output -- error indicator: */
/*             =  0  if no error, */
/*             =  1  if N is zero or negative, */
/*             =  2  if IPERM is not a valid permutation. */

/* ***REFERENCES  (NONE) */
/* ***ROUTINES CALLED  XERMSG */
/* ***REVISION HISTORY  (YYMMDD) */
/*   901004  DATE WRITTEN */
/*   920507  Modified by M. McClain to revise prologue text. */
/* ***END PROLOGUE  DPPERM */
/* ***FIRST EXECUTABLE STATEMENT  DPPERM */
    /* Parameter adjustments */
    --iperm;
    --dx;

    /* Function Body */
    *ier = 0;
    if (*n < 1) {
	*ier = 1;
	xermsg_("SLATEC", "DPPERM", "The number of values to be rearranged, "
		"N, is not positive.", ier, &c__1, (ftnlen)6, (ftnlen)6, (
		ftnlen)58);
	return 0;
    }

/*     CHECK WHETHER IPERM IS A VALID PERMUTATION */

    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	indx = (i__2 = iperm[i__], abs(i__2));
	if (indx >= 1 && indx <= *n) {
	    if (iperm[indx] > 0) {
		iperm[indx] = -iperm[indx];
		goto L100;
	    }
	}
	*ier = 2;
	xermsg_("SLATEC", "DPPERM", "The permutation vector, IPERM, is not v"
		"alid.", ier, &c__1, (ftnlen)6, (ftnlen)6, (ftnlen)44);
	return 0;
L100:
	;
    }

/*     REARRANGE THE VALUES OF DX */

/*     USE THE IPERM VECTOR AS A FLAG. */
/*     IF IPERM(I) > 0, THEN THE I-TH VALUE IS IN CORRECT LOCATION */

    i__1 = *n;
    for (istrt = 1; istrt <= i__1; ++istrt) {
	if (iperm[istrt] > 0) {
	    goto L330;
	}
	indx = istrt;
	indx0 = indx;
	dtemp = dx[istrt];
L320:
	if (iperm[indx] >= 0) {
	    goto L325;
	}
	dx[indx] = dx[-iperm[indx]];
	indx0 = indx;
	iperm[indx] = -iperm[indx];
	indx = iperm[indx];
	goto L320;
L325:
	dx[indx0] = dtemp;
L330:
	;
    }

    return 0;
} /* dpperm_ */

