/* cortb.f -- translated by f2c (version 12.02.01).
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

/* DECK CORTB */
/* Subroutine */ int cortb_(integer *nm, integer *low, integer *igh, real *ar,
	 real *ai, real *ortr, real *orti, integer *m, real *zr, real *zi)
{
    /* System generated locals */
    integer ar_dim1, ar_offset, ai_dim1, ai_offset, zr_dim1, zr_offset, 
	    zi_dim1, zi_offset, i__1, i__2, i__3;

    /* Local variables */
    static real h__;
    static integer i__, j, la;
    static real gi, gr;
    static integer mm, mp, kp1, mp1;

/* ***BEGIN PROLOGUE  CORTB */
/* ***PURPOSE  Form the eigenvectors of a complex general matrix from */
/*            eigenvectors of upper Hessenberg matrix output from */
/*            CORTH. */
/* ***LIBRARY   SLATEC (EISPACK) */
/* ***CATEGORY  D4C4 */
/* ***TYPE      COMPLEX (ORTBAK-S, CORTB-C) */
/* ***KEYWORDS  EIGENVALUES, EIGENVECTORS, EISPACK */
/* ***AUTHOR  Smith, B. T., et al. */
/* ***DESCRIPTION */

/*     This subroutine is a translation of a complex analogue of */
/*     the ALGOL procedure ORTBAK, NUM. MATH. 12, 349-368(1968) */
/*     by Martin and Wilkinson. */
/*     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 339-358(1971). */

/*     This subroutine forms the eigenvectors of a COMPLEX GENERAL */
/*     matrix by back transforming those of the corresponding */
/*     upper Hessenberg matrix determined by  CORTH. */

/*     On INPUT */

/*        NM must be set to the row dimension of the two-dimensional */
/*          array parameters, AR, AI, ZR, and ZI, as declared in the */
/*          calling program dimension statement.  NM is an INTEGER */
/*          variable. */

/*        LOW and IGH are two INTEGER variables determined by the */
/*          balancing subroutine  CBAL.  If  CBAL  has not been used, */
/*          set LOW=1 and IGH equal to the order of the matrix. */

/*        AR and AI contain information about the unitary trans- */
/*          formations used in the reduction by  CORTH  in their */
/*          strict lower triangles.  AR and AI are two-dimensional */
/*          REAL arrays, dimensioned AR(NM,IGH) and AI(NM,IGH). */

/*        ORTR and ORTI contain further information about the unitary */
/*          transformations used in the reduction by  CORTH.  Only */
/*          elements LOW through IGH are used.  ORTR and ORTI are */
/*          one-dimensional REAL arrays, dimensioned ORTR(IGH) and */
/*          ORTI(IGH). */

/*        M is the number of columns of Z=(ZR,ZI) to be back transformed. */
/*          M is an INTEGER variable. */

/*        ZR and ZI contain the real and imaginary parts, respectively, */
/*          of the eigenvectors to be back transformed in their first */
/*          M columns.  ZR and ZI are two-dimensional REAL arrays, */
/*          dimensioned ZR(NM,M) and ZI(NM,M). */

/*     On OUTPUT */

/*        ZR and ZI contain the real and imaginary parts, respectively, */
/*          of the transformed eigenvectors in their first M columns. */

/*        ORTR and ORTI have been altered. */

/*     Note that CORTB preserves vector Euclidean norms. */

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
/* ***END PROLOGUE  CORTB */


/* ***FIRST EXECUTABLE STATEMENT  CORTB */
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
    --ortr;
    --orti;

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
	if (ar[mp + (mp - 1) * ar_dim1] == 0.f && ai[mp + (mp - 1) * ai_dim1] 
		== 0.f) {
	    goto L140;
	}
/*     .......... H BELOW IS NEGATIVE OF H FORMED IN CORTH .......... */
	h__ = ar[mp + (mp - 1) * ar_dim1] * ortr[mp] + ai[mp + (mp - 1) * 
		ai_dim1] * orti[mp];
	mp1 = mp + 1;

	i__2 = *igh;
	for (i__ = mp1; i__ <= i__2; ++i__) {
	    ortr[i__] = ar[i__ + (mp - 1) * ar_dim1];
	    orti[i__] = ai[i__ + (mp - 1) * ai_dim1];
/* L100: */
	}

	i__2 = *m;
	for (j = 1; j <= i__2; ++j) {
	    gr = 0.f;
	    gi = 0.f;

	    i__3 = *igh;
	    for (i__ = mp; i__ <= i__3; ++i__) {
		gr = gr + ortr[i__] * zr[i__ + j * zr_dim1] + orti[i__] * zi[
			i__ + j * zi_dim1];
		gi = gi + ortr[i__] * zi[i__ + j * zi_dim1] - orti[i__] * zr[
			i__ + j * zr_dim1];
/* L110: */
	    }

	    gr /= h__;
	    gi /= h__;

	    i__3 = *igh;
	    for (i__ = mp; i__ <= i__3; ++i__) {
		zr[i__ + j * zr_dim1] = zr[i__ + j * zr_dim1] + gr * ortr[i__]
			 - gi * orti[i__];
		zi[i__ + j * zi_dim1] = zi[i__ + j * zi_dim1] + gr * orti[i__]
			 + gi * ortr[i__];
/* L120: */
	    }

/* L130: */
	}

L140:
	;
    }

L200:
    return 0;
} /* cortb_ */

