/* cscale.f -- translated by f2c (version 12.02.01).
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

/* Table of constant values */

static integer c__1 = 1;
static real c_b12 = 2.f;

/* DECK CSCALE */
/* Subroutine */ int cscale_(real *a, integer *nrda, integer *nrow, integer *
	ncol, real *cols, real *colsav, real *rows, real *rowsav, real *anorm,
	 real *scales, integer *iscale, integer *ic)
{
    /* Initialized data */

    static real ten4 = 1e4f;
    static real ten20 = 1e20f;

    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;

    /* Local variables */
    static integer j, k;
    static real p, s, cs;
    static integer ip;
    extern doublereal sdot_(integer *, real *, integer *, real *, integer *);
    static real alog2, ascale;

/* ***BEGIN PROLOGUE  CSCALE */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to BVSUP */
/* ***LIBRARY   SLATEC */
/* ***TYPE      SINGLE PRECISION (CSCALE-S, DCSCAL-D) */
/* ***AUTHOR  Watts, H. A., (SNLA) */
/* ***DESCRIPTION */

/*     This routine scales the matrix A by columns when needed */

/* ***SEE ALSO  BVSUP */
/* ***ROUTINES CALLED  SDOT */
/* ***REVISION HISTORY  (YYMMDD) */
/*   750601  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890831  Modified array declarations.  (WRB) */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900328  Added TYPE section.  (WRB) */
/*   910722  Updated AUTHOR section.  (ALS) */
/* ***END PROLOGUE  CSCALE */

    /* Parameter adjustments */
    a_dim1 = *nrda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --cols;
    --colsav;
    --rows;
    --rowsav;
    --scales;

    /* Function Body */

/* ***FIRST EXECUTABLE STATEMENT  CSCALE */
    if (*iscale != -1) {
	goto L25;
    }

    if (*ic == 0) {
	goto L10;
    }
    i__1 = *ncol;
    for (k = 1; k <= i__1; ++k) {
/* L5: */
	cols[k] = sdot_(nrow, &a[k * a_dim1 + 1], &c__1, &a[k * a_dim1 + 1], &
		c__1);
    }

L10:
    ascale = *anorm / *ncol;
    i__1 = *ncol;
    for (k = 1; k <= i__1; ++k) {
	cs = cols[k];
	if (cs > ten4 * ascale || ten4 * cs < ascale) {
	    goto L50;
	}
	if (cs < 1.f / ten20 || cs > ten20) {
	    goto L50;
	}
/* L20: */
    }

L25:
    i__1 = *ncol;
    for (k = 1; k <= i__1; ++k) {
/* L30: */
	scales[k] = 1.f;
    }
    return 0;

L50:
    alog2 = log(2.f);
    *anorm = 0.f;
    i__1 = *ncol;
    for (k = 1; k <= i__1; ++k) {
	cs = cols[k];
	if (cs != 0.f) {
	    goto L60;
	}
	scales[k] = 1.f;
	goto L100;
L60:
	p = log(cs) / alog2;
	ip = p * -.5f;
	s = pow_ri(&c_b12, &ip);
	scales[k] = s;
	if (*ic == 1) {
	    goto L70;
	}
	cols[k] = s * s * cols[k];
	*anorm += cols[k];
	colsav[k] = cols[k];
L70:
	i__2 = *nrow;
	for (j = 1; j <= i__2; ++j) {
/* L80: */
	    a[j + k * a_dim1] = s * a[j + k * a_dim1];
	}
L100:
	;
    }

    if (*ic == 0) {
	return 0;
    }

    i__1 = *nrow;
    for (k = 1; k <= i__1; ++k) {
	rows[k] = sdot_(ncol, &a[k + a_dim1], nrda, &a[k + a_dim1], nrda);
	rowsav[k] = rows[k];
/* L200: */
	*anorm += rows[k];
    }
    return 0;
} /* cscale_ */

