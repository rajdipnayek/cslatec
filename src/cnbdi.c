/* cnbdi.f -- translated by f2c (version 12.02.01).
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

/* DECK CNBDI */
/* Subroutine */ int cnbdi_(complex *abe, integer *lda, integer *n, integer *
	ml, integer *mu, integer *ipvt, complex *det)
{
    /* System generated locals */
    integer abe_dim1, abe_offset, i__1, i__2;
    real r__1, r__2;
    complex q__1, q__2;

    /* Local variables */
    static integer i__;
    static real ten;

/* ***BEGIN PROLOGUE  CNBDI */
/* ***PURPOSE  Compute the determinant of a band matrix using the factors */
/*            computed by CNBCO or CNBFA. */
/* ***LIBRARY   SLATEC */
/* ***CATEGORY  D3C2 */
/* ***TYPE      COMPLEX (SNBDI-S, DNBDI-D, CNBDI-C) */
/* ***KEYWORDS  BANDED, DETERMINANT, LINEAR EQUATIONS, NONSYMMETRIC */
/* ***AUTHOR  Voorhees, E. A., (LANL) */
/* ***DESCRIPTION */

/*     CNBDI computes the determinant of a band matrix */
/*     using the factors computed by CNBCO or CNBFA. */
/*     If the inverse is needed, use CNBSL  N  times. */

/*     On Entry */

/*        ABE     COMPLEX(LDA, NC) */
/*                the output from CNBCO or CNBFA. */
/*                NC must be .GE. 2*ML+MU+1 . */

/*        LDA     INTEGER */
/*                the leading dimension of the array  ABE . */

/*        N       INTEGER */
/*                the order of the original matrix. */

/*        ML      INTEGER */
/*                number of diagonals below the main diagonal. */

/*        MU      INTEGER */
/*                number of diagonals above the main diagonal. */

/*        IPVT    INTEGER(N) */
/*                the pivot vector from CNBCO or CNBFA. */

/*     On Return */

/*        DET     COMPLEX(2) */
/*                determinant of original matrix. */
/*                Determinant = DET(1) * 10.0**DET(2) */
/*                with  1.0 .LE. CABS1(DET(1)) .LT. 10.0 */
/*                or  DET(1) = 0.0 . */

/* ***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W. */
/*                 Stewart, LINPACK Users' Guide, SIAM, 1979. */
/* ***ROUTINES CALLED  (NONE) */
/* ***REVISION HISTORY  (YYMMDD) */
/*   800730  DATE WRITTEN */
/*   890831  Modified array declarations.  (WRB) */
/*   890831  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   920501  Reformatted the REFERENCES section.  (WRB) */
/* ***END PROLOGUE  CNBDI */


/* ***FIRST EXECUTABLE STATEMENT  CNBDI */
    /* Parameter adjustments */
    abe_dim1 = *lda;
    abe_offset = 1 + abe_dim1;
    abe -= abe_offset;
    --ipvt;
    --det;

    /* Function Body */
    det[1].r = 1.f, det[1].i = 0.f;
    det[2].r = 0.f, det[2].i = 0.f;
    ten = 10.f;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (ipvt[i__] != i__) {
	    q__1.r = -det[1].r, q__1.i = -det[1].i;
	    det[1].r = q__1.r, det[1].i = q__1.i;
	}
	i__2 = i__ + (*ml + 1) * abe_dim1;
	q__1.r = abe[i__2].r * det[1].r - abe[i__2].i * det[1].i, q__1.i = 
		abe[i__2].r * det[1].i + abe[i__2].i * det[1].r;
	det[1].r = q__1.r, det[1].i = q__1.i;
	if ((r__1 = det[1].r, dabs(r__1)) + (r__2 = r_imag(&det[1]), dabs(
		r__2)) == 0.f) {
	    goto L60;
	}
L10:
	if ((r__1 = det[1].r, dabs(r__1)) + (r__2 = r_imag(&det[1]), dabs(
		r__2)) >= 1.f) {
	    goto L20;
	}
	q__2.r = ten, q__2.i = 0.f;
	q__1.r = q__2.r * det[1].r - q__2.i * det[1].i, q__1.i = q__2.r * det[
		1].i + q__2.i * det[1].r;
	det[1].r = q__1.r, det[1].i = q__1.i;
	q__1.r = det[2].r - 1.f, q__1.i = det[2].i - 0.f;
	det[2].r = q__1.r, det[2].i = q__1.i;
	goto L10;
L20:
L30:
	if ((r__1 = det[1].r, dabs(r__1)) + (r__2 = r_imag(&det[1]), dabs(
		r__2)) < ten) {
	    goto L40;
	}
	q__2.r = ten, q__2.i = 0.f;
	c_div(&q__1, &det[1], &q__2);
	det[1].r = q__1.r, det[1].i = q__1.i;
	q__1.r = det[2].r + 1.f, q__1.i = det[2].i + 0.f;
	det[2].r = q__1.r, det[2].i = q__1.i;
	goto L30;
L40:
/* L50: */
	;
    }
L60:
    return 0;
} /* cnbdi_ */

