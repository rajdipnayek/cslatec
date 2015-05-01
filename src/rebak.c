/* rebak.f -- translated by f2c (version 12.02.01).
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

/* DECK REBAK */
/* Subroutine */ int rebak_(integer *nm, integer *n, real *b, real *dl, 
	integer *m, real *z__)
{
    /* System generated locals */
    integer b_dim1, b_offset, z_dim1, z_offset, i__1, i__2, i__3;

    /* Local variables */
    static integer i__, j, k;
    static real x;
    static integer i1, ii;

/* ***BEGIN PROLOGUE  REBAK */
/* ***PURPOSE  Form the eigenvectors of a generalized symmetric */
/*            eigensystem from the eigenvectors of derived matrix output */
/*            from REDUC or REDUC2. */
/* ***LIBRARY   SLATEC (EISPACK) */
/* ***CATEGORY  D4C4 */
/* ***TYPE      SINGLE PRECISION (REBAK-S) */
/* ***KEYWORDS  EIGENVALUES, EIGENVECTORS, EISPACK */
/* ***AUTHOR  Smith, B. T., et al. */
/* ***DESCRIPTION */

/*     This subroutine is a translation of the ALGOL procedure REBAKA, */
/*     NUM. MATH. 11, 99-110(1968) by Martin and Wilkinson. */
/*     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 303-314(1971). */

/*     This subroutine forms the eigenvectors of a generalized */
/*     SYMMETRIC eigensystem by back transforming those of the */
/*     derived symmetric matrix determined by  REDUC  or  REDUC2. */

/*     On Input */

/*        NM must be set to the row dimension of the two-dimensional */
/*          array parameters, B and Z, as declared in the calling */
/*          program dimension statement.  NM is an INTEGER variable. */

/*        N is the order of the matrix system.  N is an INTEGER */
/*          variable.  N must be less than or equal to NM. */

/*        B contains information about the similarity transformation */
/*          (Cholesky decomposition) used in the reduction by  REDUC */
/*          or  REDUC2  in its strict lower triangle.  B is a two- */
/*          dimensional REAL array, dimensioned B(NM,N). */

/*        DL contains further information about the transformation. */
/*          DL is a one-dimensional REAL array, dimensioned DL(N). */

/*        M is the number of eigenvectors to be back transformed. */
/*          M is an INTEGER variable. */

/*        Z contains the eigenvectors to be back transformed in its */
/*          first M columns.  Z is a two-dimensional REAL array */
/*          dimensioned Z(NM,M). */

/*     On Output */

/*        Z contains the transformed eigenvectors in its first */
/*          M columns. */

/*     Questions and comments should be directed to B. S. Garbow, */
/*     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY */
/*     ------------------------------------------------------------------ */

/* ***REFERENCES  B. T. Smith, J. M. Boyle, J. J. Dongarra, B. S. Garbow, */
/*                 Y. Ikebe, V. C. Klema and C. B. Moler, Matrix Eigen- */
/*                 system Routines - EISPACK Guide, Springer-Verlag, */
/*                 1976. */
/* ***ROUTINES CALLED  (NONE) */
/* ***REVISION HISTORY  (YYMMDD) */
/*   760101  DATE WRITTEN */
/*   890831  Modified array declarations.  (WRB) */
/*   890831  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   920501  Reformatted the REFERENCES section.  (WRB) */
/* ***END PROLOGUE  REBAK */


/* ***FIRST EXECUTABLE STATEMENT  REBAK */
    /* Parameter adjustments */
    z_dim1 = *nm;
    z_offset = 1 + z_dim1;
    z__ -= z_offset;
    b_dim1 = *nm;
    b_offset = 1 + b_dim1;
    b -= b_offset;
    --dl;

    /* Function Body */
    if (*m == 0) {
	goto L200;
    }

    i__1 = *m;
    for (j = 1; j <= i__1; ++j) {
/*     .......... FOR I=N STEP -1 UNTIL 1 DO -- .......... */
	i__2 = *n;
	for (ii = 1; ii <= i__2; ++ii) {
	    i__ = *n + 1 - ii;
	    i1 = i__ + 1;
	    x = z__[i__ + j * z_dim1];
	    if (i__ == *n) {
		goto L80;
	    }

	    i__3 = *n;
	    for (k = i1; k <= i__3; ++k) {
/* L60: */
		x -= b[k + i__ * b_dim1] * z__[k + j * z_dim1];
	    }

L80:
	    z__[i__ + j * z_dim1] = x / dl[i__];
/* L100: */
	}
    }

L200:
    return 0;
} /* rebak_ */

