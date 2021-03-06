/* elmbak.f -- translated by f2c (version 12.02.01).
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

/* DECK ELMBAK */
/* Subroutine */ int elmbak_(integer *nm, integer *low, integer *igh, real *a,
	 integer *int__, integer *m, real *z__)
{
    /* System generated locals */
    integer a_dim1, a_offset, z_dim1, z_offset, i__1, i__2, i__3;

    /* Local variables */
    static integer i__, j;
    static real x;
    static integer la, mm, mp, kp1, mp1;

/* ***BEGIN PROLOGUE  ELMBAK */
/* ***PURPOSE  Form the eigenvectors of a real general matrix from the */
/*            eigenvectors of the upper Hessenberg matrix output from */
/*            ELMHES. */
/* ***LIBRARY   SLATEC (EISPACK) */
/* ***CATEGORY  D4C4 */
/* ***TYPE      SINGLE PRECISION (ELMBAK-S, COMBAK-C) */
/* ***KEYWORDS  EIGENVALUES, EIGENVECTORS, EISPACK */
/* ***AUTHOR  Smith, B. T., et al. */
/* ***DESCRIPTION */

/*     This subroutine is a translation of the ALGOL procedure ELMBAK, */
/*     NUM. MATH. 12, 349-368(1968) by Martin and Wilkinson. */
/*     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 339-358(1971). */

/*     This subroutine forms the eigenvectors of a REAL GENERAL */
/*     matrix by back transforming those of the corresponding */
/*     upper Hessenberg matrix determined by  ELMHES. */

/*     On INPUT */

/*        NM must be set to the row dimension of the two-dimensional */
/*          array parameters, A and Z, as declared in the calling */
/*          program dimension statement.  NM is an INTEGER variable. */

/*        LOW and IGH are two INTEGER variables determined by the */
/*          balancing subroutine  BALANC.  If  BALANC  has not been */
/*          used, set LOW=1 and IGH equal to the order of the matrix. */

/*        A contains the multipliers which were used in the reduction */
/*          by  ELMHES  in its lower triangle below the subdiagonal. */
/*          A is a two-dimensional REAL array, dimensioned A(NM,IGH). */

/*        INT contains information on the rows and columns interchanged */
/*          in the reduction by  ELMHES.  Only elements LOW through IGH */
/*          are used.  INT is a one-dimensional INTEGER array, */
/*          dimensioned INT(IGH). */

/*        M is the number of columns of Z to be back transformed. */
/*          M is an INTEGER variable. */

/*        Z contains the real and imaginary parts of the eigenvectors */
/*          to be back transformed in its first M columns.  Z is a */
/*          two-dimensional REAL array, dimensioned Z(NM,M). */

/*     On OUTPUT */

/*        Z contains the real and imaginary parts of the transformed */
/*          eigenvectors in its first M columns. */

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
/* ***END PROLOGUE  ELMBAK */


/* ***FIRST EXECUTABLE STATEMENT  ELMBAK */
    /* Parameter adjustments */
    z_dim1 = *nm;
    z_offset = 1 + z_dim1;
    z__ -= z_offset;
    a_dim1 = *nm;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --int__;

    /* Function Body */
    if (*m == 0) {
	goto L200;
    }
    la = *igh - 1;
    kp1 = *low + 1;
    if (la < kp1) {
	goto L200;
    }
/*     .......... FOR MP=IGH-1 STEP -1 UNTIL LOW+1 DO -- .......... */
    i__1 = la;
    for (mm = kp1; mm <= i__1; ++mm) {
	mp = *low + *igh - mm;
	mp1 = mp + 1;

	i__2 = *igh;
	for (i__ = mp1; i__ <= i__2; ++i__) {
	    x = a[i__ + (mp - 1) * a_dim1];
	    if (x == 0.f) {
		goto L110;
	    }

	    i__3 = *m;
	    for (j = 1; j <= i__3; ++j) {
/* L100: */
		z__[i__ + j * z_dim1] += x * z__[mp + j * z_dim1];
	    }

L110:
	    ;
	}

	i__ = int__[mp];
	if (i__ == mp) {
	    goto L140;
	}

	i__2 = *m;
	for (j = 1; j <= i__2; ++j) {
	    x = z__[i__ + j * z_dim1];
	    z__[i__ + j * z_dim1] = z__[mp + j * z_dim1];
	    z__[mp + j * z_dim1] = x;
/* L130: */
	}

L140:
	;
    }

L200:
    return 0;
} /* elmbak_ */

