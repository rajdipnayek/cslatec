/* dsilus.f -- translated by f2c (version 12.02.01).
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

/* DECK DSILUS */
/* Subroutine */ int dsilus_(integer *n, integer *nelt, integer *ia, integer *
	ja, doublereal *a, integer *isym, integer *nl, integer *il, integer *
	jl, doublereal *l, doublereal *dinv, integer *nu, integer *iu, 
	integer *ju, doublereal *u, integer *nrow, integer *ncol)
{
    /* System generated locals */
    integer i__1, i__2, i__3;

    /* Local variables */
    static integer i__, j, k, kc, kr, ibgn, iend, jbgn, jend, icol, indx;
    static doublereal temp;
    static integer irow, indx1, indx2, itemp, jtemp, indxc1, indxc2, indxr1, 
	    indxr2;

/* ***BEGIN PROLOGUE  DSILUS */
/* ***PURPOSE  Incomplete LU Decomposition Preconditioner SLAP Set Up. */
/*            Routine to generate the incomplete LDU decomposition of a */
/*            matrix.  The unit lower triangular factor L is stored by */
/*            rows and the unit upper triangular factor U is stored by */
/*            columns.  The inverse of the diagonal matrix D is stored. */
/*            No fill in is allowed. */
/* ***LIBRARY   SLATEC (SLAP) */
/* ***CATEGORY  D2E */
/* ***TYPE      DOUBLE PRECISION (SSILUS-S, DSILUS-D) */
/* ***KEYWORDS  INCOMPLETE LU FACTORIZATION, ITERATIVE PRECONDITION, */
/*             NON-SYMMETRIC LINEAR SYSTEM, SLAP, SPARSE */
/* ***AUTHOR  Greenbaum, Anne, (Courant Institute) */
/*           Seager, Mark K., (LLNL) */
/*             Lawrence Livermore National Laboratory */
/*             PO BOX 808, L-60 */
/*             Livermore, CA 94550 (510) 423-3141 */
/*             seager@llnl.gov */
/* ***DESCRIPTION */

/* *Usage: */
/*     INTEGER N, NELT, IA(NELT), JA(NELT), ISYM */
/*     INTEGER NL, IL(NL), JL(NL), NU, IU(NU), JU(NU) */
/*     INTEGER NROW(N), NCOL(N) */
/*     DOUBLE PRECISION A(NELT), L(NL), DINV(N), U(NU) */

/*     CALL DSILUS( N, NELT, IA, JA, A, ISYM, NL, IL, JL, L, */
/*    $    DINV, NU, IU, JU, U, NROW, NCOL ) */

/* *Arguments: */
/* N      :IN       Integer */
/*         Order of the Matrix. */
/* NELT   :IN       Integer. */
/*         Number of elements in arrays IA, JA, and A. */
/* IA     :IN       Integer IA(NELT). */
/* JA     :IN       Integer JA(NELT). */
/* A      :IN       Double Precision A(NELT). */
/*         These arrays should hold the matrix A in the SLAP Column */
/*         format.  See "Description", below. */
/* ISYM   :IN       Integer. */
/*         Flag to indicate symmetric storage format. */
/*         If ISYM=0, all non-zero entries of the matrix are stored. */
/*         If ISYM=1, the matrix is symmetric, and only the lower */
/*         triangle of the matrix is stored. */
/* NL     :OUT      Integer. */
/*         Number of non-zeros in the L array. */
/* IL     :OUT      Integer IL(NL). */
/* JL     :OUT      Integer JL(NL). */
/* L      :OUT      Double Precision L(NL). */
/*         IL, JL, L  contain the unit lower triangular factor of  the */
/*         incomplete decomposition  of some  matrix stored  in   SLAP */
/*         Row format.     The   Diagonal  of ones  *IS*  stored.  See */
/*         "DESCRIPTION", below for more details about the SLAP format. */
/* NU     :OUT      Integer. */
/*         Number of non-zeros in the U array. */
/* IU     :OUT      Integer IU(NU). */
/* JU     :OUT      Integer JU(NU). */
/* U      :OUT      Double Precision     U(NU). */
/*         IU, JU, U contain   the unit upper triangular factor of the */
/*         incomplete  decomposition    of some matrix  stored in SLAP */
/*         Column  format.   The Diagonal of ones   *IS*  stored.  See */
/*         "Description", below  for  more  details  about  the   SLAP */
/*         format. */
/* NROW   :WORK     Integer NROW(N). */
/*         NROW(I) is the number of non-zero elements in the I-th row */
/*         of L. */
/* NCOL   :WORK     Integer NCOL(N). */
/*         NCOL(I) is the number of non-zero elements in the I-th */
/*         column of U. */

/* *Description */
/*       IL, JL, L should contain the unit  lower triangular factor of */
/*       the incomplete decomposition of the A matrix  stored in SLAP */
/*       Row format.  IU, JU, U should contain  the unit upper factor */
/*       of the  incomplete decomposition of  the A matrix  stored in */
/*       SLAP Column format This ILU factorization can be computed by */
/*       the DSILUS routine. The diagonals (which are all one's) are */
/*       stored. */

/*       =================== S L A P Column format ================== */

/*       This routine  requires that  the matrix A  be stored in  the */
/*       SLAP Column format.  In this format the non-zeros are stored */
/*       counting down columns (except for  the diagonal entry, which */
/*       must appear first in each  "column")  and are stored  in the */
/*       double precision array A.   In other words,  for each column */
/*       in the matrix put the diagonal entry in  A.  Then put in the */
/*       other non-zero  elements going down  the column (except  the */
/*       diagonal) in order.   The  IA array holds the  row index for */
/*       each non-zero.  The JA array holds the offsets  into the IA, */
/*       A arrays  for  the  beginning  of each   column.   That  is, */
/*       IA(JA(ICOL)),  A(JA(ICOL)) points   to the beginning  of the */
/*       ICOL-th   column    in    IA and   A.      IA(JA(ICOL+1)-1), */
/*       A(JA(ICOL+1)-1) points to  the  end of the   ICOL-th column. */
/*       Note that we always have  JA(N+1) = NELT+1,  where N is  the */
/*       number of columns in  the matrix and NELT  is the number  of */
/*       non-zeros in the matrix. */

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

/* ***SEE ALSO  SILUR */
/* ***REFERENCES  1. Gene Golub and Charles Van Loan, Matrix Computations, */
/*                  Johns Hopkins University Press, Baltimore, Maryland, */
/*                  1983. */
/* ***ROUTINES CALLED  (NONE) */
/* ***REVISION HISTORY  (YYMMDD) */
/*   890404  DATE WRITTEN */
/*   890404  Previous REVISION DATE */
/*   890915  Made changes requested at July 1989 CML Meeting.  (MKS) */
/*   890922  Numerous changes to prologue to make closer to SLATEC */
/*           standard.  (FNF) */
/*   890929  Numerous changes to reduce SP/DP differences.  (FNF) */
/*   910411  Prologue converted to Version 4.0 format.  (BAB) */
/*   920511  Added complete declaration section.  (WRB) */
/*   920929  Corrected format of reference.  (FNF) */
/*   930701  Updated CATEGORY section.  (FNF, WRB) */
/* ***END PROLOGUE  DSILUS */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/* ***FIRST EXECUTABLE STATEMENT  DSILUS */

/*         Count number of elements in each row of the lower triangle. */

    /* Parameter adjustments */
    --ncol;
    --nrow;
    --dinv;
    --a;
    --ja;
    --ia;
    --l;
    --jl;
    --il;
    --u;
    --ju;
    --iu;

    /* Function Body */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	nrow[i__] = 0;
	ncol[i__] = 0;
/* L10: */
    }
/* VD$R NOCONCUR */
/* VD$R NOVECTOR */
    i__1 = *n;
    for (icol = 1; icol <= i__1; ++icol) {
	jbgn = ja[icol] + 1;
	jend = ja[icol + 1] - 1;
	if (jbgn <= jend) {
	    i__2 = jend;
	    for (j = jbgn; j <= i__2; ++j) {
		if (ia[j] < icol) {
		    ++ncol[icol];
		} else {
		    ++nrow[ia[j]];
		    if (*isym != 0) {
			++ncol[ia[j]];
		    }
		}
/* L20: */
	    }
	}
/* L30: */
    }
    ju[1] = 1;
    il[1] = 1;
    i__1 = *n;
    for (icol = 1; icol <= i__1; ++icol) {
	il[icol + 1] = il[icol] + nrow[icol];
	ju[icol + 1] = ju[icol] + ncol[icol];
	nrow[icol] = il[icol];
	ncol[icol] = ju[icol];
/* L40: */
    }

/*         Copy the matrix A into the L and U structures. */
    i__1 = *n;
    for (icol = 1; icol <= i__1; ++icol) {
	dinv[icol] = a[ja[icol]];
	jbgn = ja[icol] + 1;
	jend = ja[icol + 1] - 1;
	if (jbgn <= jend) {
	    i__2 = jend;
	    for (j = jbgn; j <= i__2; ++j) {
		irow = ia[j];
		if (irow < icol) {
/*         Part of the upper triangle. */
		    iu[ncol[icol]] = irow;
		    u[ncol[icol]] = a[j];
		    ++ncol[icol];
		} else {
/*         Part of the lower triangle (stored by row). */
		    jl[nrow[irow]] = icol;
		    l[nrow[irow]] = a[j];
		    ++nrow[irow];
		    if (*isym != 0) {
/*         Symmetric...Copy lower triangle into upper triangle as well. */
			iu[ncol[irow]] = icol;
			u[ncol[irow]] = a[j];
			++ncol[irow];
		    }
		}
/* L50: */
	    }
	}
/* L60: */
    }

/*         Sort the rows of L and the columns of U. */
    i__1 = *n;
    for (k = 2; k <= i__1; ++k) {
	jbgn = ju[k];
	jend = ju[k + 1] - 1;
	if (jbgn < jend) {
	    i__2 = jend - 1;
	    for (j = jbgn; j <= i__2; ++j) {
		i__3 = jend;
		for (i__ = j + 1; i__ <= i__3; ++i__) {
		    if (iu[j] > iu[i__]) {
			itemp = iu[j];
			iu[j] = iu[i__];
			iu[i__] = itemp;
			temp = u[j];
			u[j] = u[i__];
			u[i__] = temp;
		    }
/* L70: */
		}
/* L80: */
	    }
	}
	ibgn = il[k];
	iend = il[k + 1] - 1;
	if (ibgn < iend) {
	    i__2 = iend - 1;
	    for (i__ = ibgn; i__ <= i__2; ++i__) {
		i__3 = iend;
		for (j = i__ + 1; j <= i__3; ++j) {
		    if (jl[i__] > jl[j]) {
			jtemp = ju[i__];
			ju[i__] = ju[j];
			ju[j] = jtemp;
			temp = l[i__];
			l[i__] = l[j];
			l[j] = temp;
		    }
/* L90: */
		}
/* L100: */
	    }
	}
/* L110: */
    }

/*         Perform the incomplete LDU decomposition. */
    i__1 = *n;
    for (i__ = 2; i__ <= i__1; ++i__) {

/*           I-th row of L */
	indx1 = il[i__];
	indx2 = il[i__ + 1] - 1;
	if (indx1 > indx2) {
	    goto L200;
	}
	i__2 = indx2;
	for (indx = indx1; indx <= i__2; ++indx) {
	    if (indx == indx1) {
		goto L180;
	    }
	    indxr1 = indx1;
	    indxr2 = indx - 1;
	    indxc1 = ju[jl[indx]];
	    indxc2 = ju[jl[indx] + 1] - 1;
	    if (indxc1 > indxc2) {
		goto L180;
	    }
L160:
	    kr = jl[indxr1];
L170:
	    kc = iu[indxc1];
	    if (kr > kc) {
		++indxc1;
		if (indxc1 <= indxc2) {
		    goto L170;
		}
	    } else if (kr < kc) {
		++indxr1;
		if (indxr1 <= indxr2) {
		    goto L160;
		}
	    } else if (kr == kc) {
		l[indx] -= l[indxr1] * dinv[kc] * u[indxc1];
		++indxr1;
		++indxc1;
		if (indxr1 <= indxr2 && indxc1 <= indxc2) {
		    goto L160;
		}
	    }
L180:
	    l[indx] /= dinv[jl[indx]];
/* L190: */
	}

/*         I-th column of U */
L200:
	indx1 = ju[i__];
	indx2 = ju[i__ + 1] - 1;
	if (indx1 > indx2) {
	    goto L260;
	}
	i__2 = indx2;
	for (indx = indx1; indx <= i__2; ++indx) {
	    if (indx == indx1) {
		goto L240;
	    }
	    indxc1 = indx1;
	    indxc2 = indx - 1;
	    indxr1 = il[iu[indx]];
	    indxr2 = il[iu[indx] + 1] - 1;
	    if (indxr1 > indxr2) {
		goto L240;
	    }
L210:
	    kr = jl[indxr1];
L220:
	    kc = iu[indxc1];
	    if (kr > kc) {
		++indxc1;
		if (indxc1 <= indxc2) {
		    goto L220;
		}
	    } else if (kr < kc) {
		++indxr1;
		if (indxr1 <= indxr2) {
		    goto L210;
		}
	    } else if (kr == kc) {
		u[indx] -= l[indxr1] * dinv[kc] * u[indxc1];
		++indxr1;
		++indxc1;
		if (indxr1 <= indxr2 && indxc1 <= indxc2) {
		    goto L210;
		}
	    }
L240:
	    u[indx] /= dinv[iu[indx]];
/* L250: */
	}

/*         I-th diagonal element */
L260:
	indxr1 = il[i__];
	indxr2 = il[i__ + 1] - 1;
	if (indxr1 > indxr2) {
	    goto L300;
	}
	indxc1 = ju[i__];
	indxc2 = ju[i__ + 1] - 1;
	if (indxc1 > indxc2) {
	    goto L300;
	}
L270:
	kr = jl[indxr1];
L280:
	kc = iu[indxc1];
	if (kr > kc) {
	    ++indxc1;
	    if (indxc1 <= indxc2) {
		goto L280;
	    }
	} else if (kr < kc) {
	    ++indxr1;
	    if (indxr1 <= indxr2) {
		goto L270;
	    }
	} else if (kr == kc) {
	    dinv[i__] -= l[indxr1] * dinv[kc] * u[indxc1];
	    ++indxr1;
	    ++indxc1;
	    if (indxr1 <= indxr2 && indxc1 <= indxc2) {
		goto L270;
	    }
	}

L300:
	;
    }

/*         Replace diagonal elements by their inverses. */
/* VD$ VECTOR */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	dinv[i__] = 1. / dinv[i__];
/* L430: */
    }

    return 0;
/* ------------- LAST LINE OF DSILUS FOLLOWS ---------------------------- */
} /* dsilus_ */

