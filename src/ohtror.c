/* ohtror.f -- translated by f2c (version 12.02.01).
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

/* DECK OHTROR */
/* Subroutine */ int ohtror_(real *q, integer *n, integer *nrda, real *diag, 
	integer *irank, real *div, real *td)
{
    /* System generated locals */
    integer q_dim1, q_offset, i__1, i__2, i__3;
    real r__1, r__2;

    /* Local variables */
    static integer j, k, l;
    static real dd, qs, sig;
    static integer kir;
    static real sqd;
    static integer irp;
    static real tdv;
    static integer kirm, nmir;
    extern doublereal sdot_(integer *, real *, integer *, real *, integer *);
    static real diagk;

/* ***BEGIN PROLOGUE  OHTROR */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to BVSUP */
/* ***LIBRARY   SLATEC */
/* ***TYPE      SINGLE PRECISION (OHTROR-S) */
/* ***AUTHOR  Watts, H. A., (SNLA) */
/* ***DESCRIPTION */

/*     For a rank deficient problem, additional orthogonal */
/*     HOUSEHOLDER transformations are applied to the right side */
/*     of Q to further reduce the triangular form. */
/*     Thus, after application of the routines ORTHOL and OHTROR */
/*     to the original matrix, the result is a nonsingular */
/*     triangular matrix while the remainder of the matrix */
/*     has been zeroed out. */

/* ***SEE ALSO  BVSUP */
/* ***ROUTINES CALLED  SDOT */
/* ***REVISION HISTORY  (YYMMDD) */
/*   750601  DATE WRITTEN */
/*   890831  Modified array declarations.  (WRB) */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900402  Added TYPE section.  (WRB) */
/*   910722  Updated AUTHOR section.  (ALS) */
/* ***END PROLOGUE  OHTROR */
/* ***FIRST EXECUTABLE STATEMENT  OHTROR */
    /* Parameter adjustments */
    q_dim1 = *nrda;
    q_offset = 1 + q_dim1;
    q -= q_offset;
    --diag;
    --div;
    --td;

    /* Function Body */
    nmir = *n - *irank;
    irp = *irank + 1;
    i__1 = *irank;
    for (k = 1; k <= i__1; ++k) {
	kir = irp - k;
	diagk = diag[kir];
	sig = diagk * diagk + sdot_(&nmir, &q[kir + irp * q_dim1], nrda, &q[
		kir + irp * q_dim1], nrda);
	r__1 = sqrt(sig);
	r__2 = -diagk;
	dd = r_sign(&r__1, &r__2);
	div[kir] = dd;
	tdv = diagk - dd;
	td[kir] = tdv;
	if (k == *irank) {
	    goto L30;
	}
	kirm = kir - 1;
	sqd = dd * diagk - sig;
	i__2 = kirm;
	for (j = 1; j <= i__2; ++j) {
	    qs = (tdv * q[j + kir * q_dim1] + sdot_(&nmir, &q[j + irp * 
		    q_dim1], nrda, &q[kir + irp * q_dim1], nrda)) / sqd;
	    q[j + kir * q_dim1] += qs * tdv;
	    i__3 = *n;
	    for (l = irp; l <= i__3; ++l) {
/* L10: */
		q[j + l * q_dim1] += qs * q[kir + l * q_dim1];
	    }
/* L20: */
	}
L30:
	;
    }
    return 0;
} /* ohtror_ */

