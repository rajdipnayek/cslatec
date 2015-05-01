/* ssics.f -- translated by f2c (version 12.02.01).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/

#include <stdio.h>
#include <stdlib.h> /* For exit() */
#include <f2c.h>

/* Table of constant values */

static integer c__1 = 1;
static integer c__2 = 2;

/* DECK SSICS */
/* Subroutine */ int ssics_(integer *n, integer *nelt, integer *ia, integer *
	ja, real *a, integer *isym, integer *nel, integer *iel, integer *jel, 
	real *el, real *d__, real *r__, integer *iwarn)
{
    /* System generated locals */
    integer i__1, i__2, i__3, i__4;

    /* Local variables */
    static integer i__, j, ic, ir, irr, ibgn, iend, jbgn, jend, icol, irow;
    static char xern1[8];
    static integer icbgn, icend, irbgn, irend;
    static real eltmp;
    static integer jeltmp;
    extern /* Subroutine */ int xermsg_(char *, char *, char *, integer *, 
	    integer *, ftnlen, ftnlen, ftnlen);

/* ***BEGIN PROLOGUE  SSICS */
/* ***PURPOSE  Incompl. Cholesky Decomposition Preconditioner SLAP Set Up. */
/*            Routine to generate the Incomplete Cholesky decomposition, */
/*            L*D*L-trans, of a symmetric positive definite matrix, A, */
/*            which is stored in SLAP Column format.  The unit lower */
/*            triangular matrix L is stored by rows, and the inverse of */
/*            the diagonal matrix D is stored. */
/* ***LIBRARY   SLATEC (SLAP) */
/* ***CATEGORY  D2E */
/* ***TYPE      SINGLE PRECISION (SSICS-S, DSICS-D) */
/* ***KEYWORDS  INCOMPLETE CHOLESKY FACTORIZATION, */
/*             ITERATIVE PRECONDITION, LINEAR SYSTEM, SLAP SPARSE */
/* ***AUTHOR  Greenbaum, Anne, (Courant Institute) */
/*           Seager, Mark K., (LLNL) */
/*             Lawrence Livermore National Laboratory */
/*             PO BOX 808, L-60 */
/*             Livermore, CA 94550 (510) 423-3141 */
/*             seager@llnl.gov */
/* ***DESCRIPTION */

/* *Usage: */
/*     INTEGER N, NELT, IA(NELT), JA(NELT), ISYM */
/*     INTEGER NEL, IEL(NEL), JEL(NEL), IWARN */
/*     REAL    A(NELT), EL(NEL), D(N), R(N) */

/*     CALL SSICS( N, NELT, IA, JA, A, ISYM, NEL, IEL, JEL, EL, D, R, */
/*    $    IWARN ) */

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
/*         If ISYM=1, the matrix is symmetric, and only the lower */
/*         triangle of the matrix is stored. */
/* NEL    :OUT      Integer. */
/*         Number of non-zeros in the lower triangle of A.   Also */
/*         corresponds to the length of the IEL, JEL, EL arrays. */
/* IEL    :OUT      Integer IEL(NEL). */
/* JEL    :OUT      Integer JEL(NEL). */
/* EL     :OUT      Real     EL(NEL). */
/*         IEL, JEL, EL contain the unit lower triangular factor  of the */
/*         incomplete decomposition   of the A  matrix  stored  in  SLAP */
/*         Row format.   The Diagonal of   ones   *IS*   stored.     See */
/*         "Description", below for more details about the SLAP Row fmt. */
/* D      :OUT      Real D(N) */
/*         Upon return this array holds D(I) = 1./DIAG(A). */
/* R      :WORK     Real R(N). */
/*         Temporary real workspace needed for the factorization. */
/* IWARN  :OUT      Integer. */
/*         This is a warning variable and is zero if the IC factoriza- */
/*         tion goes well.  It is set to the row index corresponding to */
/*         the last zero pivot found.  See "Description", below. */

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

/*       With the SLAP  format some  of  the   "inner  loops" of this */
/*       routine should vectorize  on  machines with hardware support */
/*       for vector   gather/scatter  operations.  Your compiler  may */
/*       require a compiler directive to  convince it that  there are */
/*       no  implicit  vector  dependencies.  Compiler directives for */
/*       the Alliant    FX/Fortran and CRI   CFT/CFT77 compilers  are */
/*       supplied with the standard SLAP distribution. */

/*       The IC factorization does not always exist for SPD matrices. */
/*       In the event that a zero pivot is found it is set  to be 1.0 */
/*       and the factorization proceeds.   The integer variable IWARN */
/*       is set to the last row where the Diagonal was fudged.  This */
/*       eventuality hardly ever occurs in practice. */

/* ***SEE ALSO  SCG, SSICCG */
/* ***REFERENCES  1. Gene Golub and Charles Van Loan, Matrix Computations, */
/*                  Johns Hopkins University Press, Baltimore, Maryland, */
/*                  1983. */
/* ***ROUTINES CALLED  XERMSG */
/* ***REVISION HISTORY  (YYMMDD) */
/*   871119  DATE WRITTEN */
/*   881213  Previous REVISION DATE */
/*   890915  Made changes requested at July 1989 CML Meeting.  (MKS) */
/*   890922  Numerous changes to prologue to make closer to SLATEC */
/*           standard.  (FNF) */
/*   890929  Numerous changes to reduce SP/DP differences.  (FNF) */
/*   900805  Changed XERRWV calls to calls to XERMSG.  (RWC) */
/*   910411  Prologue converted to Version 4.0 format.  (BAB) */
/*   920511  Added complete declaration section.  (WRB) */
/*   920929  Corrected format of reference.  (FNF) */
/*   930701  Updated CATEGORY section.  (FNF, WRB) */
/* ***END PROLOGUE  SSICS */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. External Subroutines .. */
/* ***FIRST EXECUTABLE STATEMENT  SSICS */

/*         Set the lower triangle in IEL, JEL, EL */

    /* Parameter adjustments */
    --r__;
    --d__;
    --a;
    --ja;
    --ia;
    --el;
    --jel;
    --iel;

    /* Function Body */
    *iwarn = 0;

/*         All matrix elements stored in IA, JA, A.  Pick out the lower */
/*         triangle (making sure that the Diagonal of EL is one) and */
/*         store by rows. */

    *nel = 1;
    iel[1] = 1;
    jel[1] = 1;
    el[1] = 1.f;
    d__[1] = a[1];
/* VD$R NOCONCUR */
    i__1 = *n;
    for (irow = 2; irow <= i__1; ++irow) {
/*         Put in the Diagonal. */
	++(*nel);
	iel[irow] = *nel;
	jel[*nel] = irow;
	el[*nel] = 1.f;
	d__[irow] = a[ja[irow]];

/*         Look in all the lower triangle columns for a matching row. */
/*         Since the matrix is symmetric, we can look across the */
/*         ITOW-th row by looking down the IROW-th column (if it is */
/*         stored ISYM=0)... */
	if (*isym == 0) {
	    icbgn = ja[irow];
	    icend = ja[irow + 1] - 1;
	} else {
	    icbgn = 1;
	    icend = irow - 1;
	}
	i__2 = icend;
	for (ic = icbgn; ic <= i__2; ++ic) {
	    if (*isym == 0) {
		icol = ia[ic];
		if (icol >= irow) {
		    goto L20;
		}
	    } else {
		icol = ic;
	    }
	    jbgn = ja[icol] + 1;
	    jend = ja[icol + 1] - 1;
	    if (jbgn <= jend && ia[jend] >= irow) {
/* VD$ NOVECTOR */
		i__3 = jend;
		for (j = jbgn; j <= i__3; ++j) {
		    if (ia[j] == irow) {
			++(*nel);
			jel[*nel] = icol;
			el[*nel] = a[j];
			goto L20;
		    }
/* L10: */
		}
	    }
L20:
	    ;
	}
/* L30: */
    }
    iel[*n + 1] = *nel + 1;

/*         Sort ROWS of lower triangle into descending order (count out */
/*         along rows out from Diagonal). */

    i__1 = *n;
    for (irow = 2; irow <= i__1; ++irow) {
	ibgn = iel[irow] + 1;
	iend = iel[irow + 1] - 1;
	if (ibgn < iend) {
	    i__2 = iend - 1;
	    for (i__ = ibgn; i__ <= i__2; ++i__) {
/* VD$ NOVECTOR */
		i__3 = iend;
		for (j = i__ + 1; j <= i__3; ++j) {
		    if (jel[i__] > jel[j]) {
			jeltmp = jel[j];
			jel[j] = jel[i__];
			jel[i__] = jeltmp;
			eltmp = el[j];
			el[j] = el[i__];
			el[i__] = eltmp;
		    }
/* L40: */
		}
/* L50: */
	    }
	}
/* L60: */
    }

/*         Perform the Incomplete Cholesky decomposition by looping */
/*         over the rows. */
/*         Scale the first column.  Use the structure of A to pick out */
/*         the rows with something in column 1. */

    irbgn = ja[1] + 1;
    irend = ja[2] - 1;
    i__1 = irend;
    for (irr = irbgn; irr <= i__1; ++irr) {
	ir = ia[irr];
/*         Find the index into EL for EL(1,IR). */
/*         Hint: it's the second entry. */
	i__ = iel[ir] + 1;
	el[i__] /= d__[1];
/* L65: */
    }

    i__1 = *n;
    for (irow = 2; irow <= i__1; ++irow) {

/*         Update the IROW-th diagonal. */

	i__2 = irow - 1;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    r__[i__] = 0.f;
/* L66: */
	}
	ibgn = iel[irow] + 1;
	iend = iel[irow + 1] - 1;
	if (ibgn <= iend) {
/* LLL. OPTION ASSERT (NOHAZARD) */
/* DIR$ IVDEP */
/* VD$ NODEPCHK */
	    i__2 = iend;
	    for (i__ = ibgn; i__ <= i__2; ++i__) {
		r__[jel[i__]] = el[i__] * d__[jel[i__]];
		d__[irow] -= el[i__] * r__[jel[i__]];
/* L70: */
	    }

/*         Check to see if we have a problem with the diagonal. */

	    if (d__[irow] <= 0.f) {
		if (*iwarn == 0) {
		    *iwarn = irow;
		}
		d__[irow] = 1.f;
	    }
	}

/*         Update each EL(IROW+1:N,IROW), if there are any. */
/*         Use the structure of A to determine the Non-zero elements */
/*         of the IROW-th column of EL. */

	irbgn = ja[irow];
	irend = ja[irow + 1] - 1;
	i__2 = irend;
	for (irr = irbgn; irr <= i__2; ++irr) {
	    ir = ia[irr];
	    if (ir <= irow) {
		goto L100;
	    }
/*         Find the index into EL for EL(IR,IROW) */
	    ibgn = iel[ir] + 1;
	    iend = iel[ir + 1] - 1;
	    if (jel[ibgn] > irow) {
		goto L100;
	    }
	    i__3 = iend;
	    for (i__ = ibgn; i__ <= i__3; ++i__) {
		if (jel[i__] == irow) {
		    icend = iend;
L91:
		    if (jel[icend] >= irow) {
			--icend;
			goto L91;
		    }
/*         Sum up the EL(IR,1:IROW-1)*R(1:IROW-1) contributions. */
/* LLL. OPTION ASSERT (NOHAZARD) */
/* DIR$ IVDEP */
/* VD$ NODEPCHK */
		    i__4 = icend;
		    for (ic = ibgn; ic <= i__4; ++ic) {
			el[i__] -= el[ic] * r__[jel[ic]];
/* L80: */
		    }
		    el[i__] /= d__[irow];
		    goto L100;
		}
/* L90: */
	    }

/*         If we get here, we have real problems... */
	    fprintf(stderr, "A and EL data structure mismatch in row %i\n",
                    irow);
            fflush(stderr);
	    xermsg_("SLATEC", "SSICS", "see prev error",
                    &c__1, &c__2, 6, 5, 14);
L100:
	    ;
	}
/* L110: */
    }

/*         Replace diagonals by their inverses. */

/* VD$ CONCUR */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	d__[i__] = 1.f / d__[i__];
/* L120: */
    }
    return 0;
/* ------------- LAST LINE OF SSICS FOLLOWS ---------------------------- */
} /* ssics_ */

