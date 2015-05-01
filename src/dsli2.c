/* dsli2.f -- translated by f2c (version 12.02.01).
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

/* DECK DSLI2 */
/* Subroutine */ int dsli2_(integer *n, doublereal *b, doublereal *x, integer 
	*nel, integer *iel, integer *jel, doublereal *el)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer i__, j, jbgn, jend, icol;

/* ***BEGIN PROLOGUE  DSLI2 */
/* ***PURPOSE  SLAP Lower Triangle Matrix Backsolve. */
/*            Routine to solve a system of the form  Lx = b , where L */
/*            is a lower triangular matrix. */
/* ***LIBRARY   SLATEC (SLAP) */
/* ***CATEGORY  D2A3 */
/* ***TYPE      DOUBLE PRECISION (SSLI2-S, DSLI2-D) */
/* ***KEYWORDS  ITERATIVE PRECONDITION, LINEAR SYSTEM SOLVE, SLAP, SPARSE */
/* ***AUTHOR  Greenbaum, Anne, (Courant Institute) */
/*           Seager, Mark K., (LLNL) */
/*             Lawrence Livermore National Laboratory */
/*             PO BOX 808, L-60 */
/*             Livermore, CA 94550 (510) 423-3141 */
/*             seager@llnl.gov */
/* ***DESCRIPTION */

/* *Usage: */
/*     INTEGER N, NEL, IEL(NEL), JEL(NEL) */
/*     DOUBLE PRECISION B(N), X(N), EL(NEL) */

/*     CALL DSLI2( N, B, X, NEL, IEL, JEL, EL ) */

/* *Arguments: */
/* N      :IN       Integer */
/*         Order of the Matrix. */
/* B      :IN       Double Precision B(N). */
/*         Right hand side vector. */
/* X      :OUT      Double Precision X(N). */
/*         Solution to Lx = b. */
/* NEL    :IN       Integer. */
/*         Number of non-zeros in the EL array. */
/* IEL    :IN       Integer IEL(NEL). */
/* JEL    :IN       Integer JEL(NEL). */
/* EL     :IN       Double Precision EL(NEL). */
/*         IEL, JEL, EL contain the unit lower triangular factor   of */
/*         the incomplete decomposition   of the A  matrix  stored in */
/*         SLAP Row format.  The diagonal of  ones *IS* stored.  This */
/*         structure can be  set up by the  DS2LT  routine.  See  the */
/*         "Description", below, for more details about the  SLAP Row */
/*         format. */

/* *Description: */
/*       This routine is supplied with the SLAP package  as a routine */
/*       to perform the MSOLVE operation in the DIR iteration routine */
/*       for the driver routine DSGS.  It must be called via the SLAP */
/*       MSOLVE calling sequence convention interface routine DSLI. */
/*         **** THIS ROUTINE ITSELF DOES NOT CONFORM TO THE **** */
/*               **** SLAP MSOLVE CALLING CONVENTION **** */

/*       ==================== S L A P Row format ==================== */

/*       This routine requires  that the matrix A  be  stored  in the */
/*       SLAP  Row format.   In this format  the non-zeros are stored */
/*       counting across  rows (except for the diagonal  entry, which */
/*       must  appear first  in each  "row")  and  are stored  in the */
/*       double precision  array A.  In other words, for each row  in */
/*       the matrix  put the diagonal  entry in A.   Then put in  the */
/*       other  non-zero elements  going across  the row  (except the */
/*       diagonal) in order.  The JA array holds the column index for */
/*       each non-zero.  The IA array holds the offsets  into the JA, */
/*       A  arrays  for  the   beginning  of  each  row.    That  is, */
/*       JA(IA(IROW)),A(IA(IROW)) are the first elements of the IROW- */
/*       th row in  JA and A,  and  JA(IA(IROW+1)-1), A(IA(IROW+1)-1) */
/*       are  the last elements  of the  IROW-th row.   Note  that we */
/*       always have  IA(N+1) = NELT+1, where N is the number of rows */
/*       in the matrix  and  NELT is the  number of non-zeros  in the */
/*       matrix. */

/*       Here is an example of the SLAP Row storage format for a  5x5 */
/*       Matrix (in the A and JA arrays '|' denotes the end of a row): */

/*           5x5 Matrix         SLAP Row format for 5x5 matrix on left. */
/*                              1  2  3    4  5    6  7    8    9 10 11 */
/*       |11 12  0  0 15|   A: 11 12 15 | 22 21 | 33 35 | 44 | 55 51 53 */
/*       |21 22  0  0  0|  JA:  1  2  5 |  2  1 |  3  5 |  4 |  5  1  3 */
/*       | 0  0 33  0 35|  IA:  1  4  6    8  9   12 */
/*       | 0  0  0 44  0| */
/*       |51  0 53  0 55| */

/*       With  the SLAP  Row format  the "inner loop" of this routine */
/*       should vectorize   on machines with   hardware  support  for */
/*       vector gather/scatter operations.  Your compiler may require */
/*       a  compiler directive  to  convince   it that there  are  no */
/*       implicit vector  dependencies.  Compiler directives  for the */
/*       Alliant FX/Fortran and CRI CFT/CFT77 compilers  are supplied */
/*       with the standard SLAP distribution. */

/* ***SEE ALSO  DSLI */
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
/*   921113  Corrected C***CATEGORY line.  (FNF) */
/*   930701  Updated CATEGORY section.  (FNF, WRB) */
/* ***END PROLOGUE  DSLI2 */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/* ***FIRST EXECUTABLE STATEMENT  DSLI2 */

/*         Initialize the solution by copying the right hands side */
/*         into it. */

    /* Parameter adjustments */
    --x;
    --b;
    --el;
    --jel;
    --iel;

    /* Function Body */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	x[i__] = b[i__];
/* L10: */
    }

/* VD$ NOCONCUR */
    i__1 = *n;
    for (icol = 1; icol <= i__1; ++icol) {
	x[icol] /= el[jel[icol]];
	jbgn = jel[icol] + 1;
	jend = jel[icol + 1] - 1;
	if (jbgn <= jend) {
/* LLL. OPTION ASSERT (NOHAZARD) */
/* DIR$ IVDEP */
/* VD$ NOCONCUR */
/* VD$ NODEPCHK */
	    i__2 = jend;
	    for (j = jbgn; j <= i__2; ++j) {
		x[iel[j]] -= el[j] * x[icol];
/* L20: */
	    }
	}
/* L30: */
    }

    return 0;
/* ------------- LAST LINE OF DSLI2 FOLLOWS ---------------------------- */
} /* dsli2_ */

