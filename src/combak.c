/* combak.f -- translated by f2c (version 12.02.01).
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

/* DECK COMBAK */
/* Subroutine */ int combak_(integer *nm, integer *low, integer *igh, real *
	ar, real *ai, integer *int__, integer *m, real *zr, real *zi)
{
    /* System generated locals */
    integer ar_dim1, ar_offset, ai_dim1, ai_offset, zr_dim1, zr_offset, 
	    zi_dim1, zi_offset, i__1, i__2, i__3;

    /* Local variables */
    static integer i__, j, la, mm, mp;
    static real xi, xr;
    static integer kp1, mp1;

/* ***BEGIN PROLOGUE  COMBAK */
/* ***PURPOSE  Form the eigenvectors of a complex general matrix from the */
/*            eigenvectors of a upper Hessenberg matrix output from */
/*            COMHES. */
/* ***LIBRARY   SLATEC (EISPACK) */
/* ***CATEGORY  D4C4 */
/* ***TYPE      COMPLEX (ELMBAK-S, COMBAK-C) */
/* ***KEYWORDS  EIGENVALUES, EIGENVECTORS, EISPACK */
/* ***AUTHOR  Smith, B. T., et al. */
/* ***DESCRIPTION */

/*     This subroutine is a translation of the ALGOL procedure COMBAK, */
/*     NUM. MATH. 12, 349-368(1968) by Martin and Wilkinson. */
/*     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 339-358(1971). */

/*     This subroutine forms the eigenvectors of a COMPLEX GENERAL */
/*     matrix by back transforming those of the corresponding */
/*     upper Hessenberg matrix determined by  COMHES. */

/*     On INPUT */

/*        NM must be set to the row dimension of the two-dimensional */
/*          array parameters, AR, AI, ZR and ZI, as declared in the */
/*          calling program dimension statement.  NM is an INTEGER */
/*          variable. */

/*        LOW and IGH are two INTEGER variables determined by the */
/*          balancing subroutine  CBAL.  If  CBAL  has not been used, */
/*          set LOW=1 and IGH equal to the order of the matrix. */

/*        AR and AI contain the multipliers which were used in the */
/*           reduction by  COMHES  in their lower triangles below */
/*           the subdiagonal.  AR and AI are two-dimensional REAL */
/*           arrays, dimensioned AR(NM,IGH) and AI(NM,IGH). */

/*        INT contains information on the rows and columns */
/*          interchanged in the reduction by  COMHES.  Only */
/*          elements LOW through IGH are used.  INT is a */
/*          one-dimensional INTEGER array, dimensioned INT(IGH). */

/*        M is the number of eigenvectors to be back transformed. */
/*          M is an INTEGER variable. */

/*        ZR and ZI contain the real and imaginary parts, respectively, */
/*          of the eigenvectors to be back transformed in their first M */
/*          columns.  ZR and ZI are two-dimensional REAL arrays, */
/*          dimensioned ZR(NM,M) and ZI(NM,M). */

/*     On OUTPUT */

/*        ZR and ZI contain the real and imaginary parts, respectively, */
/*          of the transformed eigenvectors in their first M columns. */

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
/* ***END PROLOGUE  COMBAK */


/* ***FIRST EXECUTABLE STATEMENT  COMBAK */
    /* Parameter adjustments */
    zi_dim1 = *nm;
    zi_offset = 1 + zi_dim1;
    zi -= zi_offset;
    zr_dim1 = *nm;
    zr_offset = 1 + zr_dim1;
    zr -= zr_offset;
    ai_dim1 = *nm;
    ai_offset = 1 + ai_dim1;
    ai -= ai_offset;
    ar_dim1 = *nm;
    ar_offset = 1 + ar_dim1;
    ar -= ar_offset;
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
	    xr = ar[i__ + (mp - 1) * ar_dim1];
	    xi = ai[i__ + (mp - 1) * ai_dim1];
	    if (xr == 0.f && xi == 0.f) {
		goto L110;
	    }

	    i__3 = *m;
	    for (j = 1; j <= i__3; ++j) {
		zr[i__ + j * zr_dim1] = zr[i__ + j * zr_dim1] + xr * zr[mp + 
			j * zr_dim1] - xi * zi[mp + j * zi_dim1];
		zi[i__ + j * zi_dim1] = zi[i__ + j * zi_dim1] + xr * zi[mp + 
			j * zi_dim1] + xi * zr[mp + j * zr_dim1];
/* L100: */
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
	    xr = zr[i__ + j * zr_dim1];
	    zr[i__ + j * zr_dim1] = zr[mp + j * zr_dim1];
	    zr[mp + j * zr_dim1] = xr;
	    xi = zi[i__ + j * zi_dim1];
	    zi[i__ + j * zi_dim1] = zi[mp + j * zi_dim1];
	    zi[mp + j * zi_dim1] = xi;
/* L130: */
	}

L140:
	;
    }

L200:
    return 0;
} /* combak_ */

