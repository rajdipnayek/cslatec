/* chkpr4.f -- translated by f2c (version 12.02.01).
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

/* DECK CHKPR4 */
/* Subroutine */ int chkpr4_(integer *iorder, real *a, real *b, integer *m, 
	integer *mbdcnd, real *c__, real *d__, integer *n, integer *nbdcnd, 
	S_fp cofx, integer *idmn, integer *ierror)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__;
    static real ai, bi, ci, xi, dlx;

/* ***BEGIN PROLOGUE  CHKPR4 */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to SEPX4 */
/* ***LIBRARY   SLATEC */
/* ***TYPE      SINGLE PRECISION (CHKPR4-S) */
/* ***AUTHOR  (UNKNOWN) */
/* ***DESCRIPTION */

/*     This program checks the input parameters for errors. */

/* ***SEE ALSO  SEPX4 */
/* ***ROUTINES CALLED  (NONE) */
/* ***REVISION HISTORY  (YYMMDD) */
/*   801001  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900402  Added TYPE section.  (WRB) */
/* ***END PROLOGUE  CHKPR4 */
/* ***FIRST EXECUTABLE STATEMENT  CHKPR4 */
    *ierror = 1;
    if (*a >= *b || *c__ >= *d__) {
	return 0;
    }

/*     CHECK BOUNDARY SWITCHES */

    *ierror = 2;
    if (*mbdcnd < 0 || *mbdcnd > 4) {
	return 0;
    }
    *ierror = 3;
    if (*nbdcnd < 0 || *nbdcnd > 4) {
	return 0;
    }

/*     CHECK FIRST DIMENSION IN CALLING ROUTINE */

    *ierror = 5;
    if (*idmn < 7) {
	return 0;
    }

/*     CHECK M */

    *ierror = 6;
    if (*m > *idmn - 1 || *m < 6) {
	return 0;
    }

/*     CHECK N */

    *ierror = 7;
    if (*n < 5) {
	return 0;
    }

/*     CHECK IORDER */

    *ierror = 8;
    if (*iorder != 2 && *iorder != 4) {
	return 0;
    }

/*     CHECK THAT EQUATION IS ELLIPTIC */

    dlx = (*b - *a) / *m;
    i__1 = *m;
    for (i__ = 2; i__ <= i__1; ++i__) {
	xi = *a + (i__ - 1) * dlx;
	(*cofx)(&xi, &ai, &bi, &ci);
	if (ai > 0.f) {
	    goto L10;
	}
	*ierror = 10;
	return 0;
L10:
/* L30: */
	;
    }

/*     NO ERROR FOUND */

    *ierror = 0;
    return 0;
} /* chkpr4_ */

