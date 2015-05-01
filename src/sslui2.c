/* sslui2.f -- translated by f2c (version 12.02.01).
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

/* DECK SSLUI2 */
/* Subroutine */ int sslui2_(integer *n, real *b, real *x, integer *il, 
	integer *jl, real *l, real *dinv, integer *iu, integer *ju, real *u)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer i__, j, jbgn, jend, icol, irow;

/* ***BEGIN PROLOGUE  SSLUI2 */
/* ***PURPOSE  SLAP Backsolve for LDU Factorization. */
/*            Routine to solve a system of the form  L*D*U X = B, */
/*            where L is a unit lower triangular matrix, D is a diagonal */
/*            matrix, and U is a unit upper triangular matrix. */
/* ***LIBRARY   SLATEC (SLAP) */
/* ***CATEGORY  D2E */
/* ***TYPE      SINGLE PRECISION (SSLUI2-S, DSLUI2-D) */
/* ***KEYWORDS  ITERATIVE PRECONDITION, NON-SYMMETRIC LINEAR SYSTEM SOLVE, */
/*             SLAP, SPARSE */
/* ***AUTHOR  Greenbaum, Anne, (Courant Institute) */
/*           Seager, Mark K., (LLNL) */
/*             Lawrence Livermore National Laboratory */
/*             PO BOX 808, L-60 */
/*             Livermore, CA 94550 (510) 423-3141 */
/*             seager@llnl.gov */
/* ***DESCRIPTION */

/* *Usage: */
/*     INTEGER N, IL(NL), JL(NL), IU(NU), JU(NU) */
/*     REAL    B(N), X(N), L(NL), DINV(N), U(NU) */

/*     CALL SSLUI2( N, B, X, IL, JL, L, DINV, IU, JU, U ) */

/* *Arguments: */
/* N      :IN       Integer */
/*         Order of the Matrix. */
/* B      :IN       Real B(N). */
/*         Right hand side. */
/* X      :OUT      Real X(N). */
/*         Solution of L*D*U x = b. */
/* IL     :IN       Integer IL(NL). */
/* JL     :IN       Integer JL(NL). */
/* L      :IN       Real     L(NL). */
/*         IL, JL, L contain the unit  lower triangular factor of the */
/*         incomplete decomposition of some matrix stored in SLAP Row */
/*         format.  The diagonal of ones *IS* stored.  This structure */
/*         can   be   set  up  by   the  SSILUS  routine.   See   the */
/*         "Description", below  for more   details about   the  SLAP */
/*         format.  (NL is the number of non-zeros in the L array.) */
/* DINV   :IN       Real DINV(N). */
/*         Inverse of the diagonal matrix D. */
/* IU     :IN       Integer IU(NU). */
/* JU     :IN       Integer JU(NU). */
/* U      :IN       Real     U(NU). */
/*         IU, JU, U contain the unit upper triangular factor  of the */
/*         incomplete decomposition  of  some  matrix stored in  SLAP */
/*         Column format.   The diagonal of ones  *IS* stored.   This */
/*         structure can be set up  by the SSILUS routine.  See   the */
/*         "Description", below   for  more   details about  the SLAP */
/*         format.  (NU is the number of non-zeros in the U array.) */

/* *Description: */
/*       This routine is supplied with  the SLAP package as a routine */
/*       to  perform  the  MSOLVE operation  in   the  SIR and   SBCG */
/*       iteration routines for  the  drivers SSILUR and SSLUBC.   It */
/*       must  be called  via   the  SLAP  MSOLVE  calling   sequence */
/*       convention interface routine SSLUI. */
/*         **** THIS ROUTINE ITSELF DOES NOT CONFORM TO THE **** */
/*               **** SLAP MSOLVE CALLING CONVENTION **** */

/*       IL, JL, L should contain the unit lower triangular factor of */
/*       the incomplete decomposition of the A matrix  stored in SLAP */
/*       Row format.  IU, JU, U should contain  the unit upper factor */
/*       of the  incomplete decomposition of  the A matrix  stored in */
/*       SLAP Column format This ILU factorization can be computed by */
/*       the SSILUS routine. The diagonals (which are all one's) are */
/*       stored. */

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

/*       ==================== S L A P Row format ==================== */

/*       This routine requires  that the matrix A  be  stored  in the */
/*       SLAP  Row format.   In this format  the non-zeros are stored */
/*       counting across  rows (except for the diagonal  entry, which */
/*       must appear first in each "row") and  are stored in the real */
/*       array A.  In other words, for each row in the matrix put the */
/*       diagonal entry in  A.   Then   put  in the   other  non-zero */
/*       elements   going  across the  row (except   the diagonal) in */
/*       order.   The  JA array  holds   the column   index for  each */
/*       non-zero.   The IA  array holds the  offsets into  the JA, A */
/*       arrays  for   the   beginning  of   each  row.   That    is, */
/*       JA(IA(IROW)),  A(IA(IROW)) points  to  the beginning  of the */
/*       IROW-th row in JA and A.   JA(IA(IROW+1)-1), A(IA(IROW+1)-1) */
/*       points to the  end of the  IROW-th row.  Note that we always */
/*       have IA(N+1) =  NELT+1, where  N  is  the number of rows  in */
/*       the matrix  and NELT  is the  number   of  non-zeros in  the */
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

/*       With  the SLAP  format  the "inner  loops" of  this  routine */
/*       should vectorize   on machines with   hardware  support  for */
/*       vector gather/scatter operations.  Your compiler may require */
/*       a  compiler directive  to  convince   it that there  are  no */
/*       implicit vector  dependencies.  Compiler directives  for the */
/*       Alliant FX/Fortran and CRI CFT/CFT77 compilers  are supplied */
/*       with the standard SLAP distribution. */

/* ***SEE ALSO  SSILUS */
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
/* ***END PROLOGUE  SSLUI2 */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/* ***FIRST EXECUTABLE STATEMENT  SSLUI2 */

/*         Solve  L*Y = B,  storing result in X, L stored by rows. */

    /* Parameter adjustments */
    --dinv;
    --x;
    --b;
    --il;
    --jl;
    --l;
    --iu;
    --ju;
    --u;

    /* Function Body */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	x[i__] = b[i__];
/* L10: */
    }
    i__1 = *n;
    for (irow = 2; irow <= i__1; ++irow) {
	jbgn = il[irow];
	jend = il[irow + 1] - 1;
	if (jbgn <= jend) {
/* LLL. OPTION ASSERT (NOHAZARD) */
/* DIR$ IVDEP */
/* VD$ ASSOC */
/* VD$ NODEPCHK */
	    i__2 = jend;
	    for (j = jbgn; j <= i__2; ++j) {
		x[irow] -= l[j] * x[jl[j]];
/* L20: */
	    }
	}
/* L30: */
    }

/*         Solve  D*Z = Y,  storing result in X. */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	x[i__] *= dinv[i__];
/* L40: */
    }

/*         Solve  U*X = Z, U stored by columns. */
    for (icol = *n; icol >= 2; --icol) {
	jbgn = ju[icol];
	jend = ju[icol + 1] - 1;
	if (jbgn <= jend) {
/* LLL. OPTION ASSERT (NOHAZARD) */
/* DIR$ IVDEP */
/* VD$ NODEPCHK */
	    i__1 = jend;
	    for (j = jbgn; j <= i__1; ++j) {
		x[iu[j]] -= u[j] * x[icol];
/* L50: */
	    }
	}
/* L60: */
    }

    return 0;
/* ------------- LAST LINE OF SSLUI2 FOLLOWS ---------------------------- */
} /* sslui2_ */

