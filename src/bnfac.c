/* bnfac.f -- translated by f2c (version 12.02.01).
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

/* DECK BNFAC */
/* Subroutine */ int bnfac_(real *w, integer *nroww, integer *nrow, integer *
	nbandl, integer *nbandu, integer *iflag)
{
    /* System generated locals */
    integer w_dim1, w_offset, i__1, i__2, i__3;

    /* Local variables */
    static integer i__, j, k, ipk, jmax, kmax, midmk;
    static real pivot;
    static integer nrowm1, middle;
    static real factor;

/* ***BEGIN PROLOGUE  BNFAC */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to BINT4 and BINTK */
/* ***LIBRARY   SLATEC */
/* ***TYPE      SINGLE PRECISION (BNFAC-S, DBNFAC-D) */
/* ***AUTHOR  (UNKNOWN) */
/* ***DESCRIPTION */

/*  BNFAC is the BANFAC routine from */
/*        * A Practical Guide to Splines *  by C. de Boor */

/*  Returns in  W  the lu-factorization (without pivoting) of the banded */
/*  matrix  A  of order  NROW  with  (NBANDL + 1 + NBANDU) bands or diag- */
/*  onals in the work array  W . */

/* *****  I N P U T  ****** */
/*  W.....Work array of size  (NROWW,NROW)  containing the interesting */
/*        part of a banded matrix  A , with the diagonals or bands of  A */
/*        stored in the rows of  W , while columns of  A  correspond to */
/*        columns of  W . This is the storage mode used in  LINPACK  and */
/*        results in efficient innermost loops. */
/*           Explicitly,  A  has  NBANDL  bands below the diagonal */
/*                            +     1     (main) diagonal */
/*                            +   NBANDU  bands above the diagonal */
/*        and thus, with    MIDDLE = NBANDU + 1, */
/*          A(I+J,J)  is in  W(I+MIDDLE,J)  for I=-NBANDU,...,NBANDL */
/*                                              J=1,...,NROW . */
/*        For example, the interesting entries of A (1,2)-banded matrix */
/*        of order  9  would appear in the first  1+1+2 = 4  rows of  W */
/*        as follows. */
/*                          13 24 35 46 57 68 79 */
/*                       12 23 34 45 56 67 78 89 */
/*                    11 22 33 44 55 66 77 88 99 */
/*                    21 32 43 54 65 76 87 98 */

/*        All other entries of  W  not identified in this way with an en- */
/*        try of  A  are never referenced . */
/*  NROWW.....Row dimension of the work array  W . */
/*        must be  .GE.  NBANDL + 1 + NBANDU  . */
/*  NBANDL.....Number of bands of  A  below the main diagonal */
/*  NBANDU.....Number of bands of  A  above the main diagonal . */

/* *****  O U T P U T  ****** */
/*  IFLAG.....Integer indicating success( = 1) or failure ( = 2) . */
/*     If  IFLAG = 1, then */
/*  W.....contains the LU-factorization of  A  into a unit lower triangu- */
/*        lar matrix  L  and an upper triangular matrix  U (both banded) */
/*        and stored in customary fashion over the corresponding entries */
/*        of  A . This makes it possible to solve any particular linear */
/*        system  A*X = B  for  X  by A */
/*              CALL BNSLV ( W, NROWW, NROW, NBANDL, NBANDU, B ) */
/*        with the solution X  contained in  B  on return . */
/*     If  IFLAG = 2, then */
/*        one of  NROW-1, NBANDL,NBANDU failed to be nonnegative, or else */
/*        one of the potential pivots was found to be zero indicating */
/*        that  A  does not have an LU-factorization. This implies that */
/*        A  is singular in case it is totally positive . */

/* *****  M E T H O D  ****** */
/*     Gauss elimination  W I T H O U T  pivoting is used. The routine is */
/*  intended for use with matrices  A  which do not require row inter- */
/*  changes during factorization, especially for the  T O T A L L Y */
/*  P O S I T I V E  matrices which occur in spline calculations. */
/*     The routine should not be used for an arbitrary banded matrix. */

/* ***SEE ALSO  BINT4, BINTK */
/* ***ROUTINES CALLED  (NONE) */
/* ***REVISION HISTORY  (YYMMDD) */
/*   800901  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890831  Modified array declarations.  (WRB) */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900328  Added TYPE section.  (WRB) */
/* ***END PROLOGUE  BNFAC */


/* ***FIRST EXECUTABLE STATEMENT  BNFAC */
    /* Parameter adjustments */
    w_dim1 = *nroww;
    w_offset = 1 + w_dim1;
    w -= w_offset;

    /* Function Body */
    *iflag = 1;
    middle = *nbandu + 1;
/*                         W(MIDDLE,.) CONTAINS THE MAIN DIAGONAL OF  A . */
    nrowm1 = *nrow - 1;
    if (nrowm1 < 0) {
	goto L120;
    } else if (nrowm1 == 0) {
	goto L110;
    } else {
	goto L10;
    }
L10:
    if (*nbandl > 0) {
	goto L30;
    }
/*                A IS UPPER TRIANGULAR. CHECK THAT DIAGONAL IS NONZERO . */
    i__1 = nrowm1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (w[middle + i__ * w_dim1] == 0.f) {
	    goto L120;
	}
/* L20: */
    }
    goto L110;
L30:
    if (*nbandu > 0) {
	goto L60;
    }
/*              A IS LOWER TRIANGULAR. CHECK THAT DIAGONAL IS NONZERO AND */
/*                 DIVIDE EACH COLUMN BY ITS DIAGONAL . */
    i__1 = nrowm1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	pivot = w[middle + i__ * w_dim1];
	if (pivot == 0.f) {
	    goto L120;
	}
/* Computing MIN */
	i__2 = *nbandl, i__3 = *nrow - i__;
	jmax = min(i__2,i__3);
	i__2 = jmax;
	for (j = 1; j <= i__2; ++j) {
	    w[middle + j + i__ * w_dim1] /= pivot;
/* L40: */
	}
/* L50: */
    }
    return 0;

/*        A  IS NOT JUST A TRIANGULAR MATRIX. CONSTRUCT LU FACTORIZATION */
L60:
    i__1 = nrowm1;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*                                  W(MIDDLE,I)  IS PIVOT FOR I-TH STEP . */
	pivot = w[middle + i__ * w_dim1];
	if (pivot == 0.f) {
	    goto L120;
	}
/*                 JMAX  IS THE NUMBER OF (NONZERO) ENTRIES IN COLUMN  I */
/*                     BELOW THE DIAGONAL . */
/* Computing MIN */
	i__2 = *nbandl, i__3 = *nrow - i__;
	jmax = min(i__2,i__3);
/*              DIVIDE EACH ENTRY IN COLUMN  I  BELOW DIAGONAL BY PIVOT . */
	i__2 = jmax;
	for (j = 1; j <= i__2; ++j) {
	    w[middle + j + i__ * w_dim1] /= pivot;
/* L70: */
	}
/*                 KMAX  IS THE NUMBER OF (NONZERO) ENTRIES IN ROW  I  TO */
/*                     THE RIGHT OF THE DIAGONAL . */
/* Computing MIN */
	i__2 = *nbandu, i__3 = *nrow - i__;
	kmax = min(i__2,i__3);
/*                  SUBTRACT  A(I,I+K)*(I-TH COLUMN) FROM (I+K)-TH COLUMN */
/*                  (BELOW ROW  I ) . */
	i__2 = kmax;
	for (k = 1; k <= i__2; ++k) {
	    ipk = i__ + k;
	    midmk = middle - k;
	    factor = w[midmk + ipk * w_dim1];
	    i__3 = jmax;
	    for (j = 1; j <= i__3; ++j) {
		w[midmk + j + ipk * w_dim1] -= w[middle + j + i__ * w_dim1] * 
			factor;
/* L80: */
	    }
/* L90: */
	}
/* L100: */
    }
/*                                       CHECK THE LAST DIAGONAL ENTRY . */
L110:
    if (w[middle + *nrow * w_dim1] != 0.f) {
	return 0;
    }
L120:
    *iflag = 2;
    return 0;
} /* bnfac_ */

