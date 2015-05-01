/* ssds.f -- translated by f2c (version 12.02.01).
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

/* DECK SSDS */
/* Subroutine */ int ssds_(integer *n, integer *nelt, integer *ia, integer *
	ja, real *a, integer *isym, real *dinv)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer icol;

/* ***BEGIN PROLOGUE  SSDS */
/* ***PURPOSE  Diagonal Scaling Preconditioner SLAP Set Up. */
/*            Routine to compute the inverse of the diagonal of a matrix */
/*            stored in the SLAP Column format. */
/* ***LIBRARY   SLATEC (SLAP) */
/* ***CATEGORY  D2E */
/* ***TYPE      SINGLE PRECISION (SSDS-S, DSDS-D) */
/* ***KEYWORDS  DIAGONAL, SLAP SPARSE */
/* ***AUTHOR  Greenbaum, Anne, (Courant Institute) */
/*           Seager, Mark K., (LLNL) */
/*             Lawrence Livermore National Laboratory */
/*             PO BOX 808, L-60 */
/*             Livermore, CA 94550 (510) 423-3141 */
/*             seager@llnl.gov */
/* ***DESCRIPTION */

/* *Usage: */
/*     INTEGER N, NELT, IA(NELT), JA(NELT), ISYM */
/*     REAL    A(NELT), DINV(N) */

/*     CALL SSDS( N, NELT, IA, JA, A, ISYM, DINV ) */

/* *Arguments: */
/* N      :IN       Integer. */
/*         Order of the Matrix. */
/* NELT   :IN       Integer. */
/*         Number of elements in arrays IA, JA, and A. */
/* IA     :INOUT    Integer IA(NELT). */
/* JA     :INOUT    Integer JA(NELT). */
/* A      :INOUT    Real A(NELT). */
/*         These arrays should hold the matrix A in the SLAP Column */
/*         format.  See "Description", below. */
/* ISYM   :IN       Integer. */
/*         Flag to indicate symmetric storage format. */
/*         If ISYM=0, all non-zero entries of the matrix are stored. */
/*         If ISYM=1, the matrix is symmetric, and only the upper */
/*         or lower triangle of the matrix is stored. */
/* DINV   :OUT      Real DINV(N). */
/*         Upon return this array holds 1./DIAG(A). */

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

/*       With the SLAP  format  all  of  the   "inner  loops" of this */
/*       routine should vectorize  on  machines with hardware support */
/*       for vector   gather/scatter  operations.  Your compiler  may */
/*       require a compiler directive to  convince it that  there are */
/*       no  implicit  vector  dependencies.  Compiler directives for */
/*       the Alliant    FX/Fortran and CRI   CFT/CFT77 compilers  are */
/*       supplied with the standard SLAP distribution. */


/* *Cautions: */
/*       This routine assumes that the diagonal of A is all  non-zero */
/*       and that the operation DINV = 1.0/DIAG(A) will not underflow */
/*       or overflow.    This  is done so that the  loop  vectorizes. */
/*       Matrices  with zero or near zero or very  large entries will */
/*       have numerical difficulties  and  must  be fixed before this */
/*       routine is called. */
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
/* ***END PROLOGUE  SSDS */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/* ***FIRST EXECUTABLE STATEMENT  SSDS */

/*         Assume the Diagonal elements are the first in each column. */
/*         This loop should *VECTORIZE*.  If it does not you may have */
/*         to add a compiler directive.  We do not check for a zero */
/*         (or near zero) diagonal element since this would interfere */
/*         with vectorization.  If this makes you nervous put a check */
/*         in!  It will run much slower. */

    /* Parameter adjustments */
    --dinv;
    --a;
    --ja;
    --ia;

    /* Function Body */
    i__1 = *n;
    for (icol = 1; icol <= i__1; ++icol) {
	dinv[icol] = 1.f / a[ja[icol]];
/* L10: */
    }

    return 0;
/* ------------- LAST LINE OF SSDS FOLLOWS ---------------------------- */
} /* ssds_ */

