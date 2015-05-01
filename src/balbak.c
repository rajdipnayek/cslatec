/* balbak.f -- translated by f2c (version 12.02.01).
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

/* DECK BALBAK */
/* Subroutine */ int balbak_(integer *nm, integer *n, integer *low, integer *
	igh, real *scale, integer *m, real *z__)
{
    /* System generated locals */
    integer z_dim1, z_offset, i__1, i__2;

    /* Local variables */
    static integer i__, j, k;
    static real s;
    static integer ii;

/* ***BEGIN PROLOGUE  BALBAK */
/* ***PURPOSE  Form the eigenvectors of a real general matrix from the */
/*            eigenvectors of matrix output from BALANC. */
/* ***LIBRARY   SLATEC (EISPACK) */
/* ***CATEGORY  D4C4 */
/* ***TYPE      SINGLE PRECISION (BALBAK-S, CBABK2-C) */
/* ***KEYWORDS  EIGENVECTORS, EISPACK */
/* ***AUTHOR  Smith, B. T., et al. */
/* ***DESCRIPTION */

/*     This subroutine is a translation of the ALGOL procedure BALBAK, */
/*     NUM. MATH. 13, 293-304(1969) by Parlett and Reinsch. */
/*     HANDBOOK FOR AUTO. COMP., Vol.II-LINEAR ALGEBRA, 315-326(1971). */

/*     This subroutine forms the eigenvectors of a REAL GENERAL */
/*     matrix by back transforming those of the corresponding */
/*     balanced matrix determined by  BALANC. */

/*     On INPUT */

/*        NM must be set to the row dimension of the two-dimensional */
/*          array parameter, Z, as declared in the calling program */
/*          dimension statement.  NM is an INTEGER variable. */

/*        N is the number of components of the vectors in matrix Z. */
/*          N is an INTEGER variable.  N must be less than or equal */
/*          to NM. */

/*        LOW and IGH are INTEGER variables determined by  BALANC. */

/*        SCALE contains information determining the permutations and */
/*          scaling factors used by  BALANC.  SCALE is a one-dimensional */
/*          REAL array, dimensioned SCALE(N). */

/*        M is the number of columns of Z to be back transformed. */
/*          M is an INTEGER variable. */

/*        Z contains the real and imaginary parts of the eigen- */
/*          vectors to be back transformed in its first M columns. */
/*          Z is a two-dimensional REAL array, dimensioned Z(NM,M). */

/*     On OUTPUT */

/*        Z contains the real and imaginary parts of the */
/*          transformed eigenvectors in its first M columns. */

/*     Questions and comments should be directed to B. S. Garbow, */
/*     Applied Mathematics Division, ARGONNE NATIONAL LABORATORY */
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
/* ***END PROLOGUE  BALBAK */


/* ***FIRST EXECUTABLE STATEMENT  BALBAK */
    /* Parameter adjustments */
    z_dim1 = *nm;
    z_offset = 1 + z_dim1;
    z__ -= z_offset;
    --scale;

    /* Function Body */
    if (*m == 0) {
	goto L200;
    }
    if (*igh == *low) {
	goto L120;
    }

    i__1 = *igh;
    for (i__ = *low; i__ <= i__1; ++i__) {
	s = scale[i__];
/*     .......... LEFT HAND EIGENVECTORS ARE BACK TRANSFORMED */
/*                IF THE FOREGOING STATEMENT IS REPLACED BY */
/*                S=1.0E0/SCALE(I). .......... */
	i__2 = *m;
	for (j = 1; j <= i__2; ++j) {
/* L100: */
	    z__[i__ + j * z_dim1] *= s;
	}

/* L110: */
    }
/*     ......... FOR I=LOW-1 STEP -1 UNTIL 1, */
/*               IGH+1 STEP 1 UNTIL N DO -- .......... */
L120:
    i__1 = *n;
    for (ii = 1; ii <= i__1; ++ii) {
	i__ = ii;
	if (i__ >= *low && i__ <= *igh) {
	    goto L140;
	}
	if (i__ < *low) {
	    i__ = *low - ii;
	}
	k = scale[i__];
	if (k == i__) {
	    goto L140;
	}

	i__2 = *m;
	for (j = 1; j <= i__2; ++j) {
	    s = z__[i__ + j * z_dim1];
	    z__[i__ + j * z_dim1] = z__[k + j * z_dim1];
	    z__[k + j * z_dim1] = s;
/* L130: */
	}

L140:
	;
    }

L200:
    return 0;
} /* balbak_ */

