/* bnslv.f -- translated by f2c (version 12.02.01).
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

/* DECK BNSLV */
/* Subroutine */ int bnslv_(real *w, integer *nroww, integer *nrow, integer *
	nbandl, integer *nbandu, real *b)
{
    /* System generated locals */
    integer w_dim1, w_offset, i__1, i__2, i__3;

    /* Local variables */
    static integer i__, j, jmax, nrowm1, middle;

/* ***BEGIN PROLOGUE  BNSLV */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to BINT4 and BINTK */
/* ***LIBRARY   SLATEC */
/* ***TYPE      SINGLE PRECISION (BNSLV-S, DBNSLV-D) */
/* ***AUTHOR  (UNKNOWN) */
/* ***DESCRIPTION */

/*  BNSLV is the BANSLV routine from */
/*        * A Practical Guide to Splines *  by C. de Boor */

/*  Companion routine to  BNFAC . It returns the solution  X  of the */
/*  linear system  A*X = B  in place of  B , given the LU-factorization */
/*  for  A  in the work array  W from BNFAC. */

/* *****  I N P U T  ****** */
/*  W, NROWW,NROW,NBANDL,NBANDU.....Describe the LU-factorization of a */
/*        banded matrix  A  of order  NROW  as constructed in  BNFAC . */
/*        For details, see  BNFAC . */
/*  B.....Right side of the system to be solved . */

/* *****  O U T P U T  ****** */
/*  B.....Contains the solution  X , of order  NROW . */

/* *****  M E T H O D  ****** */
/*     (With  A = L*U, as stored in  W,) the unit lower triangular system */
/*  L(U*X) = B  is solved for  Y = U*X, and  Y  stored in  B . Then the */
/*  upper triangular system  U*X = Y  is solved for  X  . The calcul- */
/*  ations are so arranged that the innermost loops stay within columns. */

/* ***SEE ALSO  BINT4, BINTK */
/* ***ROUTINES CALLED  (NONE) */
/* ***REVISION HISTORY  (YYMMDD) */
/*   800901  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890831  Modified array declarations.  (WRB) */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900328  Added TYPE section.  (WRB) */
/* ***END PROLOGUE  BNSLV */

/* ***FIRST EXECUTABLE STATEMENT  BNSLV */
    /* Parameter adjustments */
    w_dim1 = *nroww;
    w_offset = 1 + w_dim1;
    w -= w_offset;
    --b;

    /* Function Body */
    middle = *nbandu + 1;
    if (*nrow == 1) {
	goto L80;
    }
    nrowm1 = *nrow - 1;
    if (*nbandl == 0) {
	goto L30;
    }
/*                                 FORWARD PASS */
/*            FOR I=1,2,...,NROW-1, SUBTRACT  RIGHT SIDE(I)*(I-TH COLUMN */
/*            OF  L )  FROM RIGHT SIDE  (BELOW I-TH ROW) . */
    i__1 = nrowm1;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing MIN */
	i__2 = *nbandl, i__3 = *nrow - i__;
	jmax = min(i__2,i__3);
	i__2 = jmax;
	for (j = 1; j <= i__2; ++j) {
	    b[i__ + j] -= b[i__] * w[middle + j + i__ * w_dim1];
/* L10: */
	}
/* L20: */
    }
/*                                 BACKWARD PASS */
/*            FOR I=NROW,NROW-1,...,1, DIVIDE RIGHT SIDE(I) BY I-TH DIAG- */
/*            ONAL ENTRY OF  U, THEN SUBTRACT  RIGHT SIDE(I)*(I-TH COLUMN */
/*            OF  U)  FROM RIGHT SIDE  (ABOVE I-TH ROW). */
L30:
    if (*nbandu > 0) {
	goto L50;
    }
/*                                A  IS LOWER TRIANGULAR . */
    i__1 = *nrow;
    for (i__ = 1; i__ <= i__1; ++i__) {
	b[i__] /= w[i__ * w_dim1 + 1];
/* L40: */
    }
    return 0;
L50:
    i__ = *nrow;
L60:
    b[i__] /= w[middle + i__ * w_dim1];
/* Computing MIN */
    i__1 = *nbandu, i__2 = i__ - 1;
    jmax = min(i__1,i__2);
    i__1 = jmax;
    for (j = 1; j <= i__1; ++j) {
	b[i__ - j] -= b[i__] * w[middle - j + i__ * w_dim1];
/* L70: */
    }
    --i__;
    if (i__ > 1) {
	goto L60;
    }
L80:
    b[1] /= w[middle + w_dim1];
    return 0;
} /* bnslv_ */

