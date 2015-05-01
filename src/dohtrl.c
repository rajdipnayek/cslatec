/* dohtrl.f -- translated by f2c (version 12.02.01).
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

/* DECK DOHTRL */
/* Subroutine */ int dohtrl_(doublereal *q, integer *n, integer *nrda, 
	doublereal *diag, integer *irank, doublereal *div, doublereal *td)
{
    /* System generated locals */
    integer q_dim1, q_offset, i__1, i__2, i__3;
    doublereal d__1, d__2;

    /* Local variables */
    static integer j, k, l;
    static doublereal dd, qs, sig;
    static integer kir;
    static doublereal sqd;
    static integer irp;
    static doublereal tdv;
    extern doublereal ddot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    static integer kirm, nmir;
    static doublereal diagk;

/* ***BEGIN PROLOGUE  DOHTRL */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to DBVSUP and DSUDS */
/* ***LIBRARY   SLATEC */
/* ***TYPE      DOUBLE PRECISION (OHTROL-S, DOHTRL-D) */
/* ***AUTHOR  Watts, H. A., (SNLA) */
/* ***DESCRIPTION */

/*     For a rank deficient problem, additional orthogonal */
/*     HOUSEHOLDER transformations are applied to the left side */
/*     of Q to further reduce the triangular form. */
/*     Thus, after application of the routines DORTHR and DOHTRL */
/*     to the original matrix, the result is a nonsingular */
/*     triangular matrix while the remainder of the matrix */
/*     has been zeroed out. */

/* ***SEE ALSO  DBVSUP, DSUDS */
/* ***ROUTINES CALLED  DDOT */
/* ***REVISION HISTORY  (YYMMDD) */
/*   750601  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890831  Modified array declarations.  (WRB) */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900328  Added TYPE section.  (WRB) */
/*   910722  Updated AUTHOR section.  (ALS) */
/* ***END PROLOGUE  DOHTRL */
/* ***FIRST EXECUTABLE STATEMENT  DOHTRL */
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
	sig = diagk * diagk + ddot_(&nmir, &q[irp + kir * q_dim1], &c__1, &q[
		irp + kir * q_dim1], &c__1);
	d__1 = sqrt(sig);
	d__2 = -diagk;
	dd = d_sign(&d__1, &d__2);
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
	    qs = (tdv * q[kir + j * q_dim1] + ddot_(&nmir, &q[irp + j * 
		    q_dim1], &c__1, &q[irp + kir * q_dim1], &c__1)) / sqd;
	    q[kir + j * q_dim1] += qs * tdv;
	    i__3 = *n;
	    for (l = irp; l <= i__3; ++l) {
		q[l + j * q_dim1] += qs * q[l + kir * q_dim1];
/* L10: */
	    }
/* L20: */
	}
L30:
/* L40: */
	;
    }
    return 0;
} /* dohtrl_ */

