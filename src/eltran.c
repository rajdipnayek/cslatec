/* eltran.f -- translated by f2c (version 12.02.01).
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

/* DECK ELTRAN */
/* Subroutine */ int eltran_(integer *nm, integer *n, integer *low, integer *
	igh, real *a, integer *int__, real *z__)
{
    /* System generated locals */
    integer a_dim1, a_offset, z_dim1, z_offset, i__1, i__2;

    /* Local variables */
    static integer i__, j, kl, mm, mp, mp1;

/* ***BEGIN PROLOGUE  ELTRAN */
/* ***PURPOSE  Accumulates the stabilized elementary similarity */
/*            transformations used in the reduction of a real general */
/*            matrix to upper Hessenberg form by ELMHES. */
/* ***LIBRARY   SLATEC (EISPACK) */
/* ***CATEGORY  D4C4 */
/* ***TYPE      SINGLE PRECISION (ELTRAN-S) */
/* ***KEYWORDS  EIGENVALUES, EIGENVECTORS, EISPACK */
/* ***AUTHOR  Smith, B. T., et al. */
/* ***DESCRIPTION */

/*     This subroutine is a translation of the ALGOL procedure ELMTRANS, */
/*     NUM. MATH. 16, 181-204(1970) by Peters and Wilkinson. */
/*     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 372-395(1971). */

/*     This subroutine accumulates the stabilized elementary */
/*     similarity transformations used in the reduction of a */
/*     REAL GENERAL matrix to upper Hessenberg form by  ELMHES. */

/*     On INPUT */

/*        NM must be set to the row dimension of the two-dimensional */
/*          array parameters, A and Z, as declared in the calling */
/*          program dimension statement.  NM is an INTEGER variable. */

/*        N is the order of the matrix A.  N is an INTEGER variable. */
/*          N must be less than or equal to NM. */

/*        LOW and IGH are two INTEGER variables determined by the */
/*          balancing subroutine  BALANC.  If  BALANC  has not been */
/*          used, set LOW=1 and IGH equal to the order of the matrix, N. */

/*        A contains the multipliers which were used in the reduction */
/*          by  ELMHES  in its lower triangle below the subdiagonal. */
/*          A is a two-dimensional REAL array, dimensioned A(NM,IGH). */

/*        INT contains information on the rows and columns interchanged */
/*          in the reduction by  ELMHES.  Only elements LOW through IGH */
/*          are used.  INT is a one-dimensional INTEGER array, */
/*          dimensioned INT(IGH). */

/*     On OUTPUT */

/*        Z contains the transformation matrix produced in the reduction */
/*          by  ELMHES.  Z is a two-dimensional REAL array, dimensioned */
/*          Z(NM,N). */

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
/* ***END PROLOGUE  ELTRAN */


/* ***FIRST EXECUTABLE STATEMENT  ELTRAN */
    /* Parameter adjustments */
    z_dim1 = *nm;
    z_offset = 1 + z_dim1;
    z__ -= z_offset;
    a_dim1 = *nm;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --int__;

    /* Function Body */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {

	i__2 = *n;
	for (j = 1; j <= i__2; ++j) {
/* L60: */
	    z__[i__ + j * z_dim1] = 0.f;
	}

	z__[i__ + i__ * z_dim1] = 1.f;
/* L80: */
    }

    kl = *igh - *low - 1;
    if (kl < 1) {
	goto L200;
    }
/*     .......... FOR MP=IGH-1 STEP -1 UNTIL LOW+1 DO -- .......... */
    i__1 = kl;
    for (mm = 1; mm <= i__1; ++mm) {
	mp = *igh - mm;
	mp1 = mp + 1;

	i__2 = *igh;
	for (i__ = mp1; i__ <= i__2; ++i__) {
/* L100: */
	    z__[i__ + mp * z_dim1] = a[i__ + (mp - 1) * a_dim1];
	}

	i__ = int__[mp];
	if (i__ == mp) {
	    goto L140;
	}

	i__2 = *igh;
	for (j = mp; j <= i__2; ++j) {
	    z__[mp + j * z_dim1] = z__[i__ + j * z_dim1];
	    z__[i__ + j * z_dim1] = 0.f;
/* L130: */
	}

	z__[i__ + mp * z_dim1] = 1.f;
L140:
	;
    }

L200:
    return 0;
} /* eltran_ */

