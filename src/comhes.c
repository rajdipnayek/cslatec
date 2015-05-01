/* comhes.f -- translated by f2c (version 12.02.01).
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

/* DECK COMHES */
/* Subroutine */ int comhes_(integer *nm, integer *n, integer *low, integer *
	igh, real *ar, real *ai, integer *int__)
{
    /* System generated locals */
    integer ar_dim1, ar_offset, ai_dim1, ai_offset, i__1, i__2, i__3;
    real r__1, r__2;

    /* Local variables */
    static integer i__, j, m, la;
    static real xi, yi, xr, yr;
    static integer mm1, kp1, mp1;
    extern /* Subroutine */ int cdiv_(real *, real *, real *, real *, real *, 
	    real *);

/* ***BEGIN PROLOGUE  COMHES */
/* ***PURPOSE  Reduce a complex general matrix to complex upper Hessenberg */
/*            form using stabilized elementary similarity */
/*            transformations. */
/* ***LIBRARY   SLATEC (EISPACK) */
/* ***CATEGORY  D4C1B2 */
/* ***TYPE      COMPLEX (ELMHES-S, COMHES-C) */
/* ***KEYWORDS  EIGENVALUES, EIGENVECTORS, EISPACK */
/* ***AUTHOR  Smith, B. T., et al. */
/* ***DESCRIPTION */

/*     This subroutine is a translation of the ALGOL procedure COMHES, */
/*     NUM. MATH. 12, 349-368(1968) by Martin and Wilkinson. */
/*     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 339-358(1971). */

/*     Given a COMPLEX GENERAL matrix, this subroutine */
/*     reduces a submatrix situated in rows and columns */
/*     LOW through IGH to upper Hessenberg form by */
/*     stabilized elementary similarity transformations. */

/*     On INPUT */

/*        NM must be set to the row dimension of the two-dimensional */
/*          array parameters, AR and AI, as declared in the calling */
/*          program dimension statement.  NM is an INTEGER variable. */

/*        N is the order of the matrix A=(AR,AI).  N is an INTEGER */
/*          variable.  N must be less than or equal to NM. */

/*        LOW and IGH are two INTEGER variables determined by the */
/*          balancing subroutine  CBAL.  If  CBAL  has not been used, */
/*          set LOW=1 and IGH equal to the order of the matrix, N. */

/*        AR and AI contain the real and imaginary parts, respectively, */
/*          of the complex input matrix.  AR and AI are two-dimensional */
/*          REAL arrays, dimensioned AR(NM,N) and AI(NM,N). */

/*     On OUTPUT */

/*        AR and AI contain the real and imaginary parts, respectively, */
/*          of the upper Hessenberg matrix.  The multipliers which */
/*          were used in the reduction are stored in the remaining */
/*          triangles under the Hessenberg matrix. */

/*        INT contains information on the rows and columns */
/*          interchanged in the reduction.  Only elements LOW through */
/*          IGH are used.  INT is a one-dimensional INTEGER array, */
/*          dimensioned INT(IGH). */

/*     Calls CDIV for complex division. */

/*     Questions and comments should be directed to B. S. Garbow, */
/*     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY */
/*     ------------------------------------------------------------------ */

/* ***REFERENCES  B. T. Smith, J. M. Boyle, J. J. Dongarra, B. S. Garbow, */
/*                 Y. Ikebe, V. C. Klema and C. B. Moler, Matrix Eigen- */
/*                 system Routines - EISPACK Guide, Springer-Verlag, */
/*                 1976. */
/* ***ROUTINES CALLED  CDIV */
/* ***REVISION HISTORY  (YYMMDD) */
/*   760101  DATE WRITTEN */
/*   890831  Modified array declarations.  (WRB) */
/*   890831  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   920501  Reformatted the REFERENCES section.  (WRB) */
/* ***END PROLOGUE  COMHES */


/* ***FIRST EXECUTABLE STATEMENT  COMHES */
    /* Parameter adjustments */
    ai_dim1 = *nm;
    ai_offset = 1 + ai_dim1;
    ai -= ai_offset;
    ar_dim1 = *nm;
    ar_offset = 1 + ar_dim1;
    ar -= ar_offset;
    --int__;

    /* Function Body */
    la = *igh - 1;
    kp1 = *low + 1;
    if (la < kp1) {
	goto L200;
    }

    i__1 = la;
    for (m = kp1; m <= i__1; ++m) {
	mm1 = m - 1;
	xr = 0.f;
	xi = 0.f;
	i__ = m;

	i__2 = *igh;
	for (j = m; j <= i__2; ++j) {
	    if ((r__1 = ar[j + mm1 * ar_dim1], dabs(r__1)) + (r__2 = ai[j + 
		    mm1 * ai_dim1], dabs(r__2)) <= dabs(xr) + dabs(xi)) {
		goto L100;
	    }
	    xr = ar[j + mm1 * ar_dim1];
	    xi = ai[j + mm1 * ai_dim1];
	    i__ = j;
L100:
	    ;
	}

	int__[m] = i__;
	if (i__ == m) {
	    goto L130;
	}
/*     .......... INTERCHANGE ROWS AND COLUMNS OF AR AND AI .......... */
	i__2 = *n;
	for (j = mm1; j <= i__2; ++j) {
	    yr = ar[i__ + j * ar_dim1];
	    ar[i__ + j * ar_dim1] = ar[m + j * ar_dim1];
	    ar[m + j * ar_dim1] = yr;
	    yi = ai[i__ + j * ai_dim1];
	    ai[i__ + j * ai_dim1] = ai[m + j * ai_dim1];
	    ai[m + j * ai_dim1] = yi;
/* L110: */
	}

	i__2 = *igh;
	for (j = 1; j <= i__2; ++j) {
	    yr = ar[j + i__ * ar_dim1];
	    ar[j + i__ * ar_dim1] = ar[j + m * ar_dim1];
	    ar[j + m * ar_dim1] = yr;
	    yi = ai[j + i__ * ai_dim1];
	    ai[j + i__ * ai_dim1] = ai[j + m * ai_dim1];
	    ai[j + m * ai_dim1] = yi;
/* L120: */
	}
/*     .......... END INTERCHANGE .......... */
L130:
	if (xr == 0.f && xi == 0.f) {
	    goto L180;
	}
	mp1 = m + 1;

	i__2 = *igh;
	for (i__ = mp1; i__ <= i__2; ++i__) {
	    yr = ar[i__ + mm1 * ar_dim1];
	    yi = ai[i__ + mm1 * ai_dim1];
	    if (yr == 0.f && yi == 0.f) {
		goto L160;
	    }
	    cdiv_(&yr, &yi, &xr, &xi, &yr, &yi);
	    ar[i__ + mm1 * ar_dim1] = yr;
	    ai[i__ + mm1 * ai_dim1] = yi;

	    i__3 = *n;
	    for (j = m; j <= i__3; ++j) {
		ar[i__ + j * ar_dim1] = ar[i__ + j * ar_dim1] - yr * ar[m + j 
			* ar_dim1] + yi * ai[m + j * ai_dim1];
		ai[i__ + j * ai_dim1] = ai[i__ + j * ai_dim1] - yr * ai[m + j 
			* ai_dim1] - yi * ar[m + j * ar_dim1];
/* L140: */
	    }

	    i__3 = *igh;
	    for (j = 1; j <= i__3; ++j) {
		ar[j + m * ar_dim1] = ar[j + m * ar_dim1] + yr * ar[j + i__ * 
			ar_dim1] - yi * ai[j + i__ * ai_dim1];
		ai[j + m * ai_dim1] = ai[j + m * ai_dim1] + yr * ai[j + i__ * 
			ai_dim1] + yi * ar[j + i__ * ar_dim1];
/* L150: */
	    }

L160:
	    ;
	}

L180:
	;
    }

L200:
    return 0;
} /* comhes_ */

