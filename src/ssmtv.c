/* ssmtv.f -- translated by f2c (version 12.02.01).
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

/* DECK SSMTV */
/* Subroutine */ int ssmtv_(integer *n, real *x, real *y, integer *nelt, 
	integer *ia, integer *ja, real *a, integer *isym)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer i__, j, ibgn, iend, jbgn, jend, icol, irow;

/* ***BEGIN PROLOGUE  SSMTV */
/* ***PURPOSE  SLAP Column Format Sparse Matrix Transpose Vector Product. */
/*            Routine to calculate the sparse matrix vector product: */
/*            Y = A'*X, where ' denotes transpose. */
/* ***LIBRARY   SLATEC (SLAP) */
/* ***CATEGORY  D1B4 */
/* ***TYPE      SINGLE PRECISION (SSMTV-S, DSMTV-D) */
/* ***KEYWORDS  MATRIX TRANSPOSE VECTOR MULTIPLY, SLAP, SPARSE */
/* ***AUTHOR  Greenbaum, Anne, (Courant Institute) */
/*           Seager, Mark K., (LLNL) */
/*             Lawrence Livermore National Laboratory */
/*             PO BOX 808, L-60 */
/*             Livermore, CA 94550 (510) 423-3141 */
/*             seager@llnl.gov */
/* ***DESCRIPTION */

/* *Usage: */
/*     INTEGER  N, NELT, IA(NELT), JA(NELT), ISYM */
/*     REAL     X(N), Y(N), A(NELT) */

/*     CALL SSMTV(N, X, Y, NELT, IA, JA, A, ISYM ) */

/* *Arguments: */
/* N      :IN       Integer. */
/*         Order of the Matrix. */
/* X      :IN       Real X(N). */
/*         The vector that should be multiplied by the transpose of */
/*         the matrix. */
/* Y      :OUT      Real Y(N). */
/*         The product of the transpose of the matrix and the vector. */
/* NELT   :IN       Integer. */
/*         Number of Non-Zeros stored in A. */
/* IA     :IN       Integer IA(NELT). */
/* JA     :IN       Integer JA(NELT). */
/* A      :IN       Real A(NELT). */
/*         These arrays should hold the matrix A in the SLAP Column */
/*         format.  See "Description", below. */
/* ISYM   :IN       Integer. */
/*         Flag to indicate symmetric storage format. */
/*         If ISYM=0, all non-zero entries of the matrix are stored. */
/*         If ISYM=1, the matrix is symmetric, and only the upper */
/*         or lower triangle of the matrix is stored. */

/* *Description */
/*       =================== S L A P Column format ================== */
/*       This routine  requires that  the matrix A  be stored in  the */
/*       SLAP Column format.  In this format the non-zeros are stored */
/*       counting down columns (except for  the diagonal entry, which */
/*       must appear first in each  "column")  and are stored  in the */
/*       real array A.  In other words, for each column in the matrix */
/*       put the diagonal entry in A.  Then put in the other non-zero */
/*       elements going down   the  column (except  the diagonal)  in */
/*       order.  The IA array holds the row  index for each non-zero. */
/*       The JA array holds the offsets into the IA, A arrays for the */
/*       beginning of   each    column.    That  is,    IA(JA(ICOL)), */
/*       A(JA(ICOL)) points to the beginning of the ICOL-th column in */
/*       IA and  A.  IA(JA(ICOL+1)-1),  A(JA(ICOL+1)-1) points to the */
/*       end  of   the ICOL-th  column.  Note   that  we  always have */
/*       JA(N+1) = NELT+1, where  N  is the number of columns in  the */
/*       matrix and  NELT   is the number of non-zeros in the matrix. */

/*       Here is an example of the  SLAP Column  storage format for a */
/*       5x5 Matrix (in the A and IA arrays '|'  denotes the end of a */
/*       column): */

/*           5x5 Matrix      SLAP Column format for 5x5 matrix on left. */
/*                              1  2  3    4  5    6  7    8    9 10 11 */
/*       |11 12  0  0 15|   A: 11 21 51 | 22 12 | 33 53 | 44 | 55 15 35 */
/*       |21 22  0  0  0|  IA:  1  2  5 |  2  1 |  3  5 |  4 |  5  1  3 */
/*       | 0  0 33  0 35|  JA:  1  4  6    8  9   12 */
/*       | 0  0  0 44  0| */
/*       |51  0 53  0 55| */

/*       With  the SLAP  format  the "inner  loops" of  this  routine */
/*       should vectorize   on machines with   hardware  support  for */
/*       vector gather/scatter operations.  Your compiler may require */
/*       a  compiler directive  to  convince   it that there  are  no */
/*       implicit vector  dependencies.  Compiler directives  for the */
/*       Alliant FX/Fortran and CRI CFT/CFT77 compilers  are supplied */
/*       with the standard SLAP distribution. */

/* *Cautions: */
/*     This   routine   assumes  that  the matrix A is stored in SLAP */
/*     Column format.  It does not check  for  this (for  speed)  and */
/*     evil, ugly, ornery and nasty things  will happen if the matrix */
/*     data  structure  is,  in fact, not SLAP Column.  Beware of the */
/*     wrong data structure!!! */

/* ***SEE ALSO  SSMV */
/* ***REFERENCES  (NONE) */
/* ***ROUTINES CALLED  (NONE) */
/* ***REVISION HISTORY  (YYMMDD) */
/*   871119  DATE WRITTEN */
/*   881213  Previous REVISION DATE */
/*   890915  Made changes requested at July 1989 CML Meeting.  (MKS) */
/*   890922  Numerous changes to prologue to make closer to SLATEC */
/*           standard.  (FNF) */
/*   890929  Numerous changes to reduce SP/DP differences.  (FNF) */
/*   910411  Prologue converted to Version 4.0 format.  (BAB) */
/*   920511  Added complete declaration section.  (WRB) */
/*   930701  Updated CATEGORY section.  (FNF, WRB) */
/* ***END PROLOGUE  SSMTV */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/* ***FIRST EXECUTABLE STATEMENT  SSMTV */

/*         Zero out the result vector. */

    /* Parameter adjustments */
    --y;
    --x;
    --a;
    --ja;
    --ia;

    /* Function Body */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	y[i__] = 0.f;
/* L10: */
    }

/*         Multiply by A-Transpose. */
/*         A-Transpose is stored by rows... */
/* VD$R NOCONCUR */
    i__1 = *n;
    for (irow = 1; irow <= i__1; ++irow) {
	ibgn = ja[irow];
	iend = ja[irow + 1] - 1;
/* VD$ ASSOC */
	i__2 = iend;
	for (i__ = ibgn; i__ <= i__2; ++i__) {
	    y[irow] += a[i__] * x[ia[i__]];
/* L20: */
	}
/* L30: */
    }

    if (*isym == 1) {

/*         The matrix is non-symmetric.  Need to get the other half in... */
/*         This loops assumes that the diagonal is the first entry in */
/*         each column. */

	i__1 = *n;
	for (icol = 1; icol <= i__1; ++icol) {
	    jbgn = ja[icol] + 1;
	    jend = ja[icol + 1] - 1;
	    if (jbgn > jend) {
		goto L50;
	    }
/* LLL. OPTION ASSERT (NOHAZARD) */
/* DIR$ IVDEP */
/* VD$ NODEPCHK */
	    i__2 = jend;
	    for (j = jbgn; j <= i__2; ++j) {
		y[ia[j]] += a[j] * x[icol];
/* L40: */
	    }
L50:
	    ;
	}
    }
    return 0;
/* ------------- LAST LINE OF SSMTV FOLLOWS ---------------------------- */
} /* ssmtv_ */

