/* dpbdi.f -- translated by f2c (version 12.02.01).
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

/* DECK DPBDI */
/* Subroutine */ int dpbdi_(doublereal *abd, integer *lda, integer *n, 
	integer *m, doublereal *det)
{
    /* System generated locals */
    integer abd_dim1, abd_offset, i__1;
    doublereal d__1;

    /* Local variables */
    static integer i__;
    static doublereal s;

/* ***BEGIN PROLOGUE  DPBDI */
/* ***PURPOSE  Compute the determinant of a symmetric positive definite */
/*            band matrix using the factors computed by DPBCO or DPBFA. */
/* ***LIBRARY   SLATEC (LINPACK) */
/* ***CATEGORY  D3B2 */
/* ***TYPE      DOUBLE PRECISION (SPBDI-S, DPBDI-D, CPBDI-C) */
/* ***KEYWORDS  BANDED, DETERMINANT, INVERSE, LINEAR ALGEBRA, LINPACK, */
/*             MATRIX, POSITIVE DEFINITE */
/* ***AUTHOR  Moler, C. B., (U. of New Mexico) */
/* ***DESCRIPTION */

/*     DPBDI computes the determinant */
/*     of a double precision symmetric positive definite band matrix */
/*     using the factors computed by DPBCO or DPBFA. */
/*     If the inverse is needed, use DPBSL  N  times. */

/*     On Entry */

/*        ABD     DOUBLE PRECISION(LDA, N) */
/*                the output from DPBCO or DPBFA. */

/*        LDA     INTEGER */
/*                the leading dimension of the array  ABD . */

/*        N       INTEGER */
/*                the order of the matrix  A . */

/*        M       INTEGER */
/*                the number of diagonals above the main diagonal. */

/*     On Return */

/*        DET     DOUBLE PRECISION(2) */
/*                determinant of original matrix in the form */
/*                DETERMINANT = DET(1) * 10.0**DET(2) */
/*                with  1.0 .LE. DET(1) .LT. 10.0 */
/*                or  DET(1) .EQ. 0.0 . */

/* ***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W. */
/*                 Stewart, LINPACK Users' Guide, SIAM, 1979. */
/* ***ROUTINES CALLED  (NONE) */
/* ***REVISION HISTORY  (YYMMDD) */
/*   780814  DATE WRITTEN */
/*   890831  Modified array declarations.  (WRB) */
/*   890831  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900326  Removed duplicate information from DESCRIPTION section. */
/*           (WRB) */
/*   920501  Reformatted the REFERENCES section.  (WRB) */
/* ***END PROLOGUE  DPBDI */

/* ***FIRST EXECUTABLE STATEMENT  DPBDI */

/*     COMPUTE DETERMINANT */

    /* Parameter adjustments */
    abd_dim1 = *lda;
    abd_offset = 1 + abd_dim1;
    abd -= abd_offset;
    --det;

    /* Function Body */
    det[1] = 1.;
    det[2] = 0.;
    s = 10.;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing 2nd power */
	d__1 = abd[*m + 1 + i__ * abd_dim1];
	det[1] = d__1 * d__1 * det[1];
	if (det[1] == 0.) {
	    goto L60;
	}
L10:
	if (det[1] >= 1.) {
	    goto L20;
	}
	det[1] = s * det[1];
	det[2] += -1.;
	goto L10;
L20:
L30:
	if (det[1] < s) {
	    goto L40;
	}
	det[1] /= s;
	det[2] += 1.;
	goto L30;
L40:
/* L50: */
	;
    }
L60:
    return 0;
} /* dpbdi_ */

