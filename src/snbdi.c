/* snbdi.f -- translated by f2c (version 12.02.01).
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

/* DECK SNBDI */
/* Subroutine */ int snbdi_(real *abe, integer *lda, integer *n, integer *ml, 
	integer *mu, integer *ipvt, real *det)
{
    /* System generated locals */
    integer abe_dim1, abe_offset, i__1;

    /* Local variables */
    static integer i__;
    static real ten;

/* ***BEGIN PROLOGUE  SNBDI */
/* ***PURPOSE  Compute the determinant of a band matrix using the factors */
/*            computed by SNBCO or SNBFA. */
/* ***LIBRARY   SLATEC */
/* ***CATEGORY  D3A2 */
/* ***TYPE      SINGLE PRECISION (SNBDI-S, DNBDI-D, CNBDI-C) */
/* ***KEYWORDS  BANDED, DETERMINANT, LINEAR EQUATIONS, NONSYMMETRIC */
/* ***AUTHOR  Voorhees, E. A., (LANL) */
/* ***DESCRIPTION */

/*     SNBDI computes the determinant of a band matrix */
/*     using the factors computed by SNBCO or SNBFA. */
/*     If the inverse is needed, use SNBSL  N  times. */

/*     On Entry */

/*        ABE     REAL(LDA, NC) */
/*                the output from SNBCO or SNBFA. */
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
/*                the pivot vector from SNBCO or SNBFA. */

/*     On Return */

/*        DET     REAL(2) */
/*                determinant of original matrix. */
/*                Determinant = DET(1) * 10.0**DET(2) */
/*                with  1.0 .LE. ABS(DET(1)) .LT. 10.0 */
/*                or  DET(1) = 0.0 . */

/* ***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W. */
/*                 Stewart, LINPACK Users' Guide, SIAM, 1979. */
/* ***ROUTINES CALLED  (NONE) */
/* ***REVISION HISTORY  (YYMMDD) */
/*   800725  DATE WRITTEN */
/*   890831  Modified array declarations.  (WRB) */
/*   890831  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   920501  Reformatted the REFERENCES section.  (WRB) */
/* ***END PROLOGUE  SNBDI */

/* ***FIRST EXECUTABLE STATEMENT  SNBDI */
    /* Parameter adjustments */
    abe_dim1 = *lda;
    abe_offset = 1 + abe_dim1;
    abe -= abe_offset;
    --ipvt;
    --det;

    /* Function Body */
    det[1] = 1.f;
    det[2] = 0.f;
    ten = 10.f;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (ipvt[i__] != i__) {
	    det[1] = -det[1];
	}
	det[1] = abe[i__ + (*ml + 1) * abe_dim1] * det[1];
	if (det[1] == 0.f) {
	    goto L60;
	}
L10:
	if (dabs(det[1]) >= 1.f) {
	    goto L20;
	}
	det[1] = ten * det[1];
	det[2] += -1.f;
	goto L10;
L20:
L30:
	if (dabs(det[1]) < ten) {
	    goto L40;
	}
	det[1] /= ten;
	det[2] += 1.f;
	goto L30;
L40:
/* L50: */
	;
    }
L60:
    return 0;
} /* snbdi_ */

