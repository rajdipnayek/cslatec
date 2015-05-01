/* scpplt.f -- translated by f2c (version 12.02.01).
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

/* DECK SCPPLT */
/* Subroutine */ int scpplt_(integer *n, integer *nelt, integer *ia, integer *
	ja, real *a, integer *isym, integer *iunit)
{
    /* Format strings */
    static char fmt_1000[] = "(/\002**** Picture of Column SLAP matrix follo"
	    "ws ****\002/\002 N, NELT and Density = \002,2i10,e16.7)";
    static char fmt_1010[] = "(4x,225(i1))";
    static char fmt_1020[] = "(1x,i3,a)";

    /* System generated locals */
    integer i__1, i__2;
    real r__1;

    /* Local variables */
    static integer i__, j, jbgn, jend, icol, nmax, irow;
    static char chmat[225*225];

    /* Fortran I/O blocks */
    static cilist io___9 = { 0, 0, 0, fmt_1000, 0 };
    static cilist io___10 = { 0, 0, 0, fmt_1010, 0 };
    static cilist io___11 = { 0, 0, 0, fmt_1020, 0 };


/* ***BEGIN PROLOGUE  SCPPLT */
/* ***PURPOSE  Printer Plot of SLAP Column Format Matrix. */
/*            Routine to print out a SLAP Column format matrix in a */
/*            "printer plot" graphical representation. */
/* ***LIBRARY   SLATEC (SLAP) */
/* ***CATEGORY  N1 */
/* ***TYPE      SINGLE PRECISION (SCPPLT-S, DCPPLT-D) */
/* ***KEYWORDS  DIAGNOSTICS, LINEAR SYSTEM, SLAP SPARSE */
/* ***AUTHOR  Seager, Mark K., (LLNL) */
/*             Lawrence Livermore National Laboratory */
/*             PO BOX 808, L-60 */
/*             Livermore, CA 94550 (510) 423-3141 */
/*             seager@llnl.gov */
/* ***DESCRIPTION */

/* *Usage: */
/*     INTEGER N, NELT, IA(NELT), JA(NELT), ISYM, IUNIT */
/*     REAL    A(NELT) */

/*     CALL SCPPLT( N, NELT, IA, JA, A, ISYM, IUNIT ) */

/* *Arguments: */
/* N      :IN       Integer */
/*         Order of the Matrix. */
/*         If N.gt.MAXORD, only the leading MAXORD x MAXORD */
/*         submatrix will be printed.  (Currently MAXORD = 225.) */
/* NELT   :IN       Integer. */
/*         Number of non-zeros stored in A. */
/* IA     :IN       Integer IA(NELT). */
/* JA     :IN       Integer JA(NELT). */
/* A      :IN       Real A(NELT). */
/*         These arrays should hold the matrix A in the SLAP */
/*         Column format.  See "Description", below. */
/* ISYM   :IN       Integer. */
/*         Flag to indicate symmetric storage format. */
/*         If ISYM=0, all non-zero entries of the matrix are stored. */
/*         If ISYM=1, the matrix is symmetric, and only the lower */
/*         triangle of the matrix is stored. */
/* IUNIT  :IN       Integer. */
/*         Fortran logical I/O device unit number to write the matrix */
/*         to.  This unit must be connected in a system dependent fashion */
/*         to a file or the console or you will get a nasty message */
/*         from the Fortran I/O libraries. */

/* *Description: */
/*       This routine prints out a SLAP  Column format matrix  to the */
/*       Fortran logical I/O unit   number  IUNIT.  The  numbers them */
/*       selves  are not printed  out, but   rather  a one  character */
/*       representation of the numbers.   Elements of the matrix that */
/*       are not represented in the (IA,JA,A)  arrays are  denoted by */
/*       ' ' character (a blank).  Elements of A that are *ZERO* (and */
/*       hence  should  really not be  stored) are  denoted  by a '0' */
/*       character.  Elements of A that are *POSITIVE* are denoted by */
/*       'D' if they are Diagonal elements  and '#' if  they are off */
/*       Diagonal  elements.  Elements of  A that are *NEGATIVE* are */
/*       denoted by 'N'  if they  are Diagonal  elements and  '*' if */
/*       they are off Diagonal elements. */

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

/* *Cautions: */
/*     This routine will attempt to write to the Fortran logical output */
/*     unit IUNIT, if IUNIT .ne. 0.  Thus, the user must make sure that */
/*     this logical unit is attached to a file or terminal before calling */
/*     this routine with a non-zero value for IUNIT.  This routine does */
/*     not check for the validity of a non-zero IUNIT unit number. */

/* *Portability: */
/*     This routine, as distributed, can generate lines up to 229 */
/*     characters long.  Some Fortran systems have more restricted */
/*     line lengths.  Change parameter MAXORD and the large number */
/*     in FORMAT 1010 to reduce this line length. */

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
/*   921007  Replaced hard-wired 225 with parameter MAXORD.  (FNF) */
/*   921021  Corrected syntax of CHARACTER declaration.  (FNF) */
/*   930701  Updated CATEGORY section.  (FNF, WRB) */
/* ***END PROLOGUE  SCPPLT */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Parameters .. */
/*     .. Local Scalars .. */
/*     .. Local Arrays .. */
/*     .. Intrinsic Functions .. */
/* ***FIRST EXECUTABLE STATEMENT  SCPPLT */

/*         Set up the character matrix... */

    /* Parameter adjustments */
    --a;
    --ja;
    --ia;

    /* Function Body */
    nmax = min(225,*n);
    i__1 = nmax;
    for (i__ = 1; i__ <= i__1; ++i__) {
	s_copy(chmat + (i__ - 1) * 225, " ", nmax, (ftnlen)1);
/* L10: */
    }
    i__1 = nmax;
    for (icol = 1; icol <= i__1; ++icol) {
	jbgn = ja[icol];
	jend = ja[icol + 1] - 1;
	i__2 = jend;
	for (j = jbgn; j <= i__2; ++j) {
	    irow = ia[j];
	    if (irow <= nmax) {
		if (*isym != 0) {
/*         Put in non-sym part as well... */
		    if (a[j] == 0.f) {
			*(unsigned char *)&chmat[(irow - 1) * 225 + (icol - 1)
				] = '0';
		    } else if (a[j] > 0.f) {
			*(unsigned char *)&chmat[(irow - 1) * 225 + (icol - 1)
				] = '#';
		    } else {
			*(unsigned char *)&chmat[(irow - 1) * 225 + (icol - 1)
				] = '*';
		    }
		}
		if (irow == icol) {
/*         Diagonal entry. */
		    if (a[j] == 0.f) {
			*(unsigned char *)&chmat[(irow - 1) * 225 + (icol - 1)
				] = '0';
		    } else if (a[j] > 0.f) {
			*(unsigned char *)&chmat[(irow - 1) * 225 + (icol - 1)
				] = 'D';
		    } else {
			*(unsigned char *)&chmat[(irow - 1) * 225 + (icol - 1)
				] = 'N';
		    }
		} else {
/*         Off-Diagonal entry */
		    if (a[j] == 0.f) {
			*(unsigned char *)&chmat[(irow - 1) * 225 + (icol - 1)
				] = '0';
		    } else if (a[j] > 0.f) {
			*(unsigned char *)&chmat[(irow - 1) * 225 + (icol - 1)
				] = '#';
		    } else {
			*(unsigned char *)&chmat[(irow - 1) * 225 + (icol - 1)
				] = '*';
		    }
		}
	    }
/* L20: */
	}
/* L30: */
    }

/*         Write out the heading. */
    io___9.ciunit = *iunit;
    s_wsfe(&io___9);
    do_fio(&c__1, (char *)&(*n), (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&(*nelt), (ftnlen)sizeof(integer));
    r__1 = (real) (*nelt) / (*n * *n);
    do_fio(&c__1, (char *)&r__1, (ftnlen)sizeof(real));
    e_wsfe();
    io___10.ciunit = *iunit;
    s_wsfe(&io___10);
    i__1 = nmax;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = i__ % 10;
	do_fio(&c__1, (char *)&i__2, (ftnlen)sizeof(integer));
    }
    e_wsfe();

/*         Write out the character representations matrix elements. */
    i__2 = nmax;
    for (irow = 1; irow <= i__2; ++irow) {
	io___11.ciunit = *iunit;
	s_wsfe(&io___11);
	do_fio(&c__1, (char *)&irow, (ftnlen)sizeof(integer));
	do_fio(&c__1, chmat + (irow - 1) * 225, nmax);
	e_wsfe();
/* L40: */
    }
    return 0;

/*      The following assumes MAXORD.le.225. */
/* ------------- LAST LINE OF SCPPLT FOLLOWS ---------------------------- */
} /* scpplt_ */

