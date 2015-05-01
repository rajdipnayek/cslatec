/* dsdi.f -- translated by f2c (version 12.02.01).
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

/* DECK DSDI */
/* Subroutine */ int dsdi_(integer *n, doublereal *b, doublereal *x, integer *
	nelt, integer *ia, integer *ja, doublereal *a, integer *isym, 
	doublereal *rwork, integer *iwork)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__, locd;

/* ***BEGIN PROLOGUE  DSDI */
/* ***PURPOSE  Diagonal Matrix Vector Multiply. */
/*            Routine to calculate the product  X = DIAG*B, where DIAG */
/*            is a diagonal matrix. */
/* ***LIBRARY   SLATEC (SLAP) */
/* ***CATEGORY  D1B4 */
/* ***TYPE      DOUBLE PRECISION (SSDI-S, DSDI-D) */
/* ***KEYWORDS  ITERATIVE PRECONDITION, LINEAR SYSTEM SOLVE, SLAP, SPARSE */
/* ***AUTHOR  Greenbaum, Anne, (Courant Institute) */
/*           Seager, Mark K., (LLNL) */
/*             Lawrence Livermore National Laboratory */
/*             PO BOX 808, L-60 */
/*             Livermore, CA 94550 (510) 423-3141 */
/*             seager@llnl.gov */
/* ***DESCRIPTION */

/* *Usage: */
/*     INTEGER  N, NELT, IA(NELT), JA(NELT), ISYM, IWORK(10) */
/*     DOUBLE PRECISION B(N), X(N), A(NELT), RWORK(USER DEFINED) */

/*     CALL DSDI (N, B, X, NELT, IA, JA, A, ISYM, RWORK, IWORK) */

/* *Arguments: */
/* N      :IN       Integer */
/*         Order of the Matrix. */
/* B      :IN       Double Precision B(N). */
/*         Vector to multiply the diagonal by. */
/* X      :OUT      Double Precision X(N). */
/*         Result of DIAG*B. */
/* NELT   :DUMMY    Integer. */
/* IA     :DUMMY    Integer IA(NELT). */
/* JA     :DUMMY    Integer JA(NELT). */
/* A      :DUMMY    Double Precision A(NELT). */
/* ISYM   :DUMMY    Integer. */
/*         These are for compatibility with SLAP MSOLVE calling sequence. */
/* RWORK  :IN       Double Precision RWORK(USER DEFINED). */
/*         Work array holding the diagonal of some matrix to scale */
/*         B by.  This array must be set by the user or by a call */
/*         to the SLAP routine DSDS or DSD2S.  The length of RWORK */
/*         must be >= IWORK(4)+N. */
/* IWORK  :IN       Integer IWORK(10). */
/*         IWORK(4) holds the offset into RWORK for the diagonal matrix */
/*         to scale B by.  This is usually set up by the SLAP pre- */
/*         conditioner setup routines DSDS or DSD2S. */

/* *Description: */
/*         This routine is supplied with the SLAP package to perform */
/*         the  MSOLVE  operation for iterative drivers that require */
/*         diagonal  Scaling  (e.g., DSDCG, DSDBCG).   It  conforms */
/*         to the SLAP MSOLVE CALLING CONVENTION  and hence does not */
/*         require an interface routine as do some of the other pre- */
/*         conditioners supplied with SLAP. */

/* ***SEE ALSO  DSDS, DSD2S */
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
/* ***END PROLOGUE  DSDI */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/* ***FIRST EXECUTABLE STATEMENT  DSDI */

/*         Determine where the inverse of the diagonal */
/*         is in the work array and then scale by it. */

    /* Parameter adjustments */
    --x;
    --b;
    --a;
    --ja;
    --ia;
    --rwork;
    --iwork;

    /* Function Body */
    locd = iwork[4] - 1;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	x[i__] = rwork[locd + i__] * b[i__];
/* L10: */
    }
    return 0;
/* ------------- LAST LINE OF DSDI FOLLOWS ---------------------------- */
} /* dsdi_ */

