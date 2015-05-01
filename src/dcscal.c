/* dcscal.f -- translated by f2c (version 12.02.01).
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
static doublereal c_b14 = 2.;

/* DECK DCSCAL */
/* Subroutine */ int dcscal_(doublereal *a, integer *nrda, integer *nrow, 
	integer *ncol, doublereal *cols, doublereal *colsav, doublereal *rows,
	 doublereal *rowsav, doublereal *anorm, doublereal *scales, integer *
	iscale, integer *ic)
{
    /* Initialized data */

    static doublereal ten4 = 1e4;
    static doublereal ten20 = 1e20;

    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;

    /* Local variables */
    static integer j, k;
    static doublereal p, s, cs;
    static integer ip;
    extern doublereal ddot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    static doublereal alog2, ascale;

/* ***BEGIN PROLOGUE  DCSCAL */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to DBVSUP and DSUDS */
/* ***LIBRARY   SLATEC */
/* ***TYPE      DOUBLE PRECISION (CSCALE-S, DCSCAL-D) */
/* ***AUTHOR  Watts, H. A., (SNLA) */
/* ***DESCRIPTION */

/*     This routine scales the matrix A by columns when needed. */

/* ***SEE ALSO  DBVSUP, DSUDS */
/* ***ROUTINES CALLED  DDOT */
/* ***REVISION HISTORY  (YYMMDD) */
/*   750601  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890831  Modified array declarations.  (WRB) */
/*   890911  Removed unnecessary intrinsics.  (WRB) */
/*   890911  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900328  Added TYPE section.  (WRB) */
/*   910722  Updated AUTHOR section.  (ALS) */
/* ***END PROLOGUE  DCSCAL */

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

/*     BEGIN BLOCK PERMITTING ...EXITS TO 130 */
/*        BEGIN BLOCK PERMITTING ...EXITS TO 60 */
/* ***FIRST EXECUTABLE STATEMENT  DCSCAL */
    if (*iscale != -1) {
	goto L40;
    }

    if (*ic == 0) {
	goto L20;
    }
    i__1 = *ncol;
    for (k = 1; k <= i__1; ++k) {
	cols[k] = ddot_(nrow, &a[k * a_dim1 + 1], &c__1, &a[k * a_dim1 + 1], &
		c__1);
/* L10: */
    }
L20:

    ascale = *anorm / *ncol;
    i__1 = *ncol;
    for (k = 1; k <= i__1; ++k) {
	cs = cols[k];
/*        .........EXIT */
	if (cs > ten4 * ascale || ten4 * cs < ascale) {
	    goto L60;
	}
/*        .........EXIT */
	if (cs < 1. / ten20 || cs > ten20) {
	    goto L60;
	}
/* L30: */
    }
L40:

    i__1 = *ncol;
    for (k = 1; k <= i__1; ++k) {
	scales[k] = 1.;
/* L50: */
    }
/*     ......EXIT */
    goto L130;
L60:

    alog2 = log(2.);
    *anorm = 0.;
    i__1 = *ncol;
    for (k = 1; k <= i__1; ++k) {
	cs = cols[k];
	if (cs != 0.) {
	    goto L70;
	}
	scales[k] = 1.;
	goto L100;
L70:
	p = log(cs) / alog2;
	ip = (integer) (p * -.5);
	s = pow_di(&c_b14, &ip);
	scales[k] = s;
	if (*ic == 1) {
	    goto L80;
	}
	cols[k] = s * s * cols[k];
	*anorm += cols[k];
	colsav[k] = cols[k];
L80:
	i__2 = *nrow;
	for (j = 1; j <= i__2; ++j) {
	    a[j + k * a_dim1] = s * a[j + k * a_dim1];
/* L90: */
	}
L100:
/* L110: */
	;
    }

/*     ...EXIT */
    if (*ic == 0) {
	goto L130;
    }

    i__1 = *nrow;
    for (k = 1; k <= i__1; ++k) {
	rows[k] = ddot_(ncol, &a[k + a_dim1], nrda, &a[k + a_dim1], nrda);
	rowsav[k] = rows[k];
	*anorm += rows[k];
/* L120: */
    }
L130:
    return 0;
} /* dcscal_ */

