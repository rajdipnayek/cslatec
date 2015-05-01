/* sbhin.f -- translated by f2c (version 12.02.01).
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

/* DECK SBHIN */
/* Subroutine */ int sbhin_(integer *n, integer *nelt, integer *ia, integer *
	ja, real *a, integer *isym, real *soln, real *rhs, integer *iunit, 
	integer *job)
{
    /* Format strings */
    static char fmt_9000[] = "(a80)";
    static char fmt_9010[] = "(5i14)";
    static char fmt_9020[] = "(a3,11x,4i14)";
    static char fmt_9030[] = "(2a16,2a20)";

    /* System generated locals */
    integer i__1, i__2, i__3;

    /* Local variables */
    static integer i__, j;
    static char code[3];
    static integer ibgn, iend, nele, icol, nind, ncol;
    static real temp;
    static integer npls, nrow, nline, itemp;
    static char title[80];
    static integer nrils, nnvls, jobret;
    static char rinfmt[16], rhsfmt[20], nvlfmt[20], pntfmt[16];
    static integer nrhsls;

    /* Fortran I/O blocks */
    static cilist io___1 = { 0, 0, 0, fmt_9000, 0 };
    static cilist io___3 = { 0, 0, 0, fmt_9010, 0 };
    static cilist io___9 = { 0, 0, 0, fmt_9020, 0 };
    static cilist io___15 = { 0, 0, 0, fmt_9030, 0 };
    static cilist io___21 = { 0, 0, 0, pntfmt, 0 };
    static cilist io___23 = { 0, 0, 0, rinfmt, 0 };
    static cilist io___24 = { 0, 0, 0, nvlfmt, 0 };
    static cilist io___25 = { 0, 5, 0, rhsfmt, 0 };


/* ***BEGIN PROLOGUE  SBHIN */
/* ***PURPOSE  Read a Sparse Linear System in the Boeing/Harwell Format. */
/*            The matrix is read in and if the right hand side is also */
/*            present in the input file then it too is read in.  The */
/*            matrix is then modified to be in the SLAP Column format. */
/* ***LIBRARY   SLATEC (SLAP) */
/* ***CATEGORY  N1 */
/* ***TYPE      SINGLE PRECISION (SBHIN-S, DBHIN-D) */
/* ***KEYWORDS  LINEAR SYSTEM, MATRIX READ, SLAP SPARSE */
/* ***AUTHOR  Seager, Mark K., (LLNL) */
/*             Lawrence Livermore National Laboratory */
/*             PO BOX 808, L-60 */
/*             Livermore, CA 94550 (510) 423-3141 */
/*             seager@llnl.gov */
/* ***DESCRIPTION */

/* *Usage: */
/*     INTEGER N, NELT, IA(NELT), JA(NELT), ISYM, IUNIT, JOB */
/*     REAL    A(NELT), SOLN(N), RHS(N) */

/*     CALL SBHIN( N, NELT, IA, JA, A, ISYM, SOLN, RHS, IUNIT, JOB ) */

/* *Arguments: */
/* N      :OUT      Integer */
/*         Order of the Matrix. */
/* NELT   :INOUT    Integer. */
/*         On input NELT is the maximum number of non-zeros that */
/*         can be stored in the IA, JA, A arrays. */
/*         On output NELT is the number of non-zeros stored in A. */
/* IA     :OUT      Integer IA(NELT). */
/* JA     :OUT      Integer JA(NELT). */
/* A      :OUT      Real A(NELT). */
/*         On output these arrays hold the matrix A in the SLAP */
/*         Triad format.  See "Description", below. */
/* ISYM   :OUT      Integer. */
/*         Flag to indicate symmetric storage format. */
/*         If ISYM=0, all non-zero entries of the matrix are stored. */
/*         If ISYM=1, the matrix is symmetric, and only the lower */
/*         triangle of the matrix is stored. */
/* SOLN   :OUT      Real SOLN(N). */
/*         The solution to the linear system, if present.  This array */
/*         is accessed if and only if JOB is set to read it in, see */
/*         below.  If the user requests that SOLN be read in, but it is */
/*         not in the file, then it is simply zeroed out. */
/* RHS    :OUT      Real RHS(N). */
/*         The right hand side vector.  This array is accessed if and */
/*         only if JOB is set to read it in, see below. */
/*         If the user requests that RHS be read in, but it is not in */
/*         the file, then it is simply zeroed out. */
/* IUNIT  :IN       Integer. */
/*         Fortran logical I/O device unit number to read the matrix */
/*         from.  This unit must be connected in a system dependent */
/*         fashion to a file, or you will get a nasty message */
/*         from the Fortran I/O libraries. */
/* JOB    :INOUT    Integer. */
/*         Flag indicating what I/O operations to perform. */
/*         On input JOB indicates what Input operations to try to */
/*         perform. */
/*         JOB = 0 => Read only the matrix. */
/*         JOB = 1 => Read matrix and RHS (if present). */
/*         JOB = 2 => Read matrix and SOLN (if present). */
/*         JOB = 3 => Read matrix, RHS and SOLN (if present). */
/*         On output JOB indicates what operations were actually */
/*         performed. */
/*         JOB = -3 => Unable to parse matrix "CODE" from input file */
/*                     to determine if only the lower triangle of matrix */
/*                     is stored. */
/*         JOB = -2 => Number of non-zeros (NELT) too large. */
/*         JOB = -1 => System size (N) too large. */
/*         JOB =  0 => Read in only the matrix. */
/*         JOB =  1 => Read in the matrix and RHS. */
/*         JOB =  2 => Read in the matrix and SOLN. */
/*         JOB =  3 => Read in the matrix, RHS and SOLN. */
/*         JOB = 10 => Read in only the matrix *STRUCTURE*, but no */
/*                     non-zero entries.  Hence, A(*) is not referenced */
/*                     and has the return values the same as the input. */
/*         JOB = 11 => Read in the matrix *STRUCTURE* and RHS. */
/*         JOB = 12 => Read in the matrix *STRUCTURE* and SOLN. */
/*         JOB = 13 => Read in the matrix *STRUCTURE*, RHS and SOLN. */

/* *Description: */
/*       The format for the input is as follows.  The first line contains */
/*       a title to identify the data file.  On the second line (5I4) are */
/*       counters: NLINE, NPLS, NRILS, NNVLS, NRHSLS. */
/*        NLINE  Number of data lines (after the header) in the file. */
/*        NPLS   Number of lines for the Column Pointer data in the file. */
/*        NRILS  Number of lines for the Row indices in the file. */
/*        NNVLS  Number of lines for the Matrix elements in the file. */
/*        NRHSLS Number of lines for the RHS in the file. */
/*       The third line (A3,11X,4I4) contains a symmetry code and some */
/*       additional counters: CODE, NROW, NCOL, NIND, NELE. */
/*       On the fourth line (2A16,2A20) are formats to be used to read */
/*       the following data: PNTFNT, RINFMT, NVLFMT, RHSFMT. */
/*       Following that are the blocks of data in the order indicated. */

/*       =================== S L A P Triad format =================== */
/*       This routine requires that the  matrix A be   stored in  the */
/*       SLAP  Triad format.  In  this format only the non-zeros  are */
/*       stored.  They may appear in  *ANY* order.  The user supplies */
/*       three arrays of  length NELT, where  NELT is  the number  of */
/*       non-zeros in the matrix: (IA(NELT), JA(NELT), A(NELT)).  For */
/*       each non-zero the user puts the row and column index of that */
/*       matrix element  in the IA and  JA arrays.  The  value of the */
/*       non-zero  matrix  element is  placed   in  the corresponding */
/*       location of the A array.   This is  an  extremely  easy data */
/*       structure to generate.  On  the  other hand it   is  not too */
/*       efficient on vector computers for  the iterative solution of */
/*       linear systems.  Hence,   SLAP changes   this  input    data */
/*       structure to the SLAP Column format  for  the iteration (but */
/*       does not change it back). */

/*       Here is an example of the  SLAP Triad   storage format for a */
/*       5x5 Matrix.  Recall that the entries may appear in any order. */

/*           5x5 Matrix      SLAP Triad format for 5x5 matrix on left. */
/*                              1  2  3  4  5  6  7  8  9 10 11 */
/*       |11 12  0  0 15|   A: 51 12 11 33 15 53 55 22 35 44 21 */
/*       |21 22  0  0  0|  IA:  5  1  1  3  1  5  5  2  3  4  2 */
/*       | 0  0 33  0 35|  JA:  1  2  1  3  5  3  5  2  5  4  1 */
/*       | 0  0  0 44  0| */
/*       |51  0 53  0 55| */

/* *Portability: */
/*         You must make sure that IUNIT is a valid Fortran logical */
/*         I/O device unit number and that the unit number has been */
/*         associated with a file or the console.  This is a system */
/*         dependent function. */

/* *Implementation note: */
/*         SOLN is not read by this version.  It will simply be */
/*         zeroed out if JOB = 2 or 3 and the returned value of */
/*         JOB will indicate SOLN has not been read. */
/* ***REFERENCES  (NONE) */
/* ***ROUTINES CALLED  (NONE) */
/* ***REVISION HISTORY  (YYMMDD) */
/*   881107  DATE WRITTEN */
/*   881213  Previous REVISION DATE */
/*   890915  Made changes requested at July 1989 CML Meeting.  (MKS) */
/*   890922  Numerous changes to prologue to make closer to SLATEC */
/*           standard.  (FNF) */
/*   890929  Numerous changes to reduce SP/DP differences.  (FNF) */
/*   910411  Prologue converted to Version 4.0 format.  (BAB) */
/*   911122  Added loop to zero out RHS if user wants to read RHS, but */
/*           it's not in the input file. (MKS) */
/*   911125  Minor improvements to prologue.  (FNF) */
/*   920511  Added complete declaration section.  (WRB) */
/*   921007  Corrected description of input format.  (FNF) */
/*   921208  Added Implementation Note and code to zero out SOLN.  (FNF) */
/*   930701  Updated CATEGORY section.  (FNF, WRB) */
/* ***END PROLOGUE  SBHIN */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. Intrinsic Functions .. */
/* ***FIRST EXECUTABLE STATEMENT  SBHIN */

/*         Read Matrices In BOEING-HARWELL format. */

/* TITLE  Header line to identify data file. */
/* NLINE  Number of data lines (after the header) in the file. */
/* NPLS   Number of lines for the Column Pointer data in the file. */
/* NRILS  Number of lines for the Row indices in the data file. */
/* NNVLS  Number of lines for the Matrix elements in the data file. */
/* NRHSLS Number of lines for the RHS in the data file. */
/* ---- Only those variables needed by SLAP are referenced. ---- */

    /* Parameter adjustments */
    --rhs;
    --soln;
    --a;
    --ja;
    --ia;

    /* Function Body */
    io___1.ciunit = *iunit;
    s_rsfe(&io___1);
    do_fio(&c__1, title, (ftnlen)80);
    e_rsfe();
    io___3.ciunit = *iunit;
    s_rsfe(&io___3);
    do_fio(&c__1, (char *)&nline, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&npls, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&nrils, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&nnvls, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&nrhsls, (ftnlen)sizeof(integer));
    e_rsfe();
    io___9.ciunit = *iunit;
    s_rsfe(&io___9);
    do_fio(&c__1, code, (ftnlen)3);
    do_fio(&c__1, (char *)&nrow, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&ncol, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&nind, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&nele, (ftnlen)sizeof(integer));
    e_rsfe();
    io___15.ciunit = *iunit;
    s_rsfe(&io___15);
    do_fio(&c__1, pntfmt, (ftnlen)16);
    do_fio(&c__1, rinfmt, (ftnlen)16);
    do_fio(&c__1, nvlfmt, (ftnlen)20);
    do_fio(&c__1, rhsfmt, (ftnlen)20);
    e_rsfe();

    if (nrow > *n) {
	*n = nrow;
	jobret = -1;
	goto L999;
    }
    if (nind > *nelt) {
	*nelt = nind;
	jobret = -2;
	goto L999;
    }

/*         Set the parameters. */

    *n = nrow;
    *nelt = nind;
    if (s_cmp(code, "RUA", (ftnlen)3, (ftnlen)3) == 0) {
	*isym = 0;
    } else if (s_cmp(code, "RSA", (ftnlen)3, (ftnlen)3) == 0) {
	*isym = 1;
    } else {
	jobret = -3;
	goto L999;
    }
    io___21.ciunit = *iunit;
    s_rsfe(&io___21);
    i__1 = *n + 1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	do_fio(&c__1, (char *)&ja[i__], (ftnlen)sizeof(integer));
    }
    e_rsfe();
    io___23.ciunit = *iunit;
    s_rsfe(&io___23);
    i__1 = *nelt;
    for (i__ = 1; i__ <= i__1; ++i__) {
	do_fio(&c__1, (char *)&ia[i__], (ftnlen)sizeof(integer));
    }
    e_rsfe();
    jobret = 10;
    if (nnvls > 0) {
	io___24.ciunit = *iunit;
	s_rsfe(&io___24);
	i__1 = *nelt;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    do_fio(&c__1, (char *)&a[i__], (ftnlen)sizeof(real));
	}
	e_rsfe();
	jobret = 0;
    }
    if (*job % 2 == 1) {

/*         User requests that the RHS be read in.  If it is in the input */
/*         file, read it in; otherwise just zero it out. */

	if (nrhsls > 0) {
	    s_rsfe(&io___25);
	    i__1 = *n;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		do_fio(&c__1, (char *)&rhs[i__], (ftnlen)sizeof(real));
	    }
	    e_rsfe();
	    ++jobret;
	} else {
	    i__1 = *n;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		rhs[i__] = 0.f;
/* L10: */
	    }
	}
    }
    if (*job == 2 || *job == 3) {

/*         User requests that the SOLN be read in. */
/*         Just zero out the array. */

	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    soln[i__] = 0.f;
/* L20: */
	}
    }

/*         Now loop through the IA array making sure that the diagonal */
/*         matrix element appears first in the column.  Then sort the */
/*         rest of the column in ascending order. */

/* VD$R NOCONCUR */
/* VD$R NOVECTOR */
    i__1 = *n;
    for (icol = 1; icol <= i__1; ++icol) {
	ibgn = ja[icol];
	iend = ja[icol + 1] - 1;
	i__2 = iend;
	for (i__ = ibgn; i__ <= i__2; ++i__) {
	    if (ia[i__] == icol) {

/*              Swap the diagonal element with the first element in the */
/*              column. */

		itemp = ia[i__];
		ia[i__] = ia[ibgn];
		ia[ibgn] = itemp;
		temp = a[i__];
		a[i__] = a[ibgn];
		a[ibgn] = temp;
		goto L40;
	    }
/* L30: */
	}
L40:
	++ibgn;
	if (ibgn < iend) {
	    i__2 = iend;
	    for (i__ = ibgn; i__ <= i__2; ++i__) {
		i__3 = iend;
		for (j = i__ + 1; j <= i__3; ++j) {
		    if (ia[i__] > ia[j]) {
			itemp = ia[i__];
			ia[i__] = ia[j];
			ia[j] = itemp;
			temp = a[i__];
			a[i__] = a[j];
			a[j] = temp;
		    }
/* L50: */
		}
/* L60: */
	    }
	}
/* L70: */
    }

/*         Set return flag. */
L999:
    *job = jobret;
    return 0;
/* ------------- LAST LINE OF SBHIN FOLLOWS ------------------------------ */
} /* sbhin_ */

