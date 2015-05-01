/* dtout.f -- translated by f2c (version 12.02.01).
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

/* DECK DTOUT */
/* Subroutine */ int dtout_(integer *n, integer *nelt, integer *ia, integer *
	ja, doublereal *a, integer *isym, doublereal *soln, doublereal *rhs, 
	integer *iunit, integer *job)
{
    /* Format strings */
    static char fmt_1000[] = "(5i10)";
    static char fmt_1010[] = "(1x,i5,1x,i5,1x,d16.7)";
    static char fmt_1020[] = "(1x,d16.7)";

    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__, irhs, isoln;

    /* Fortran I/O blocks */
    static cilist io___3 = { 0, 0, 0, fmt_1000, 0 };
    static cilist io___5 = { 0, 0, 0, fmt_1010, 0 };
    static cilist io___6 = { 0, 0, 0, fmt_1020, 0 };
    static cilist io___7 = { 0, 0, 0, fmt_1020, 0 };


/* ***BEGIN PROLOGUE  DTOUT */
/* ***PURPOSE  Write out SLAP Triad Format Linear System. */
/*            Routine to write out a SLAP Triad format matrix and right */
/*            hand side and solution to the system, if known. */
/* ***LIBRARY   SLATEC (SLAP) */
/* ***CATEGORY  N1 */
/* ***TYPE      DOUBLE PRECISION (STOUT-S, DTOUT-D) */
/* ***KEYWORDS  DIAGNOSTICS, LINEAR SYSTEM, SLAP SPARSE */
/* ***AUTHOR  Seager, Mark K., (LLNL) */
/*             Lawrence Livermore National Laboratory */
/*             PO BOX 808, L-60 */
/*             Livermore, CA 94550 (510) 423-3141 */
/*             seager@llnl.gov */
/* ***DESCRIPTION */

/* *Usage: */
/*     INTEGER N, NELT, IA(NELT), JA(NELT), ISYM, IUNIT, JOB */
/*     DOUBLE PRECISION A(NELT), SOLN(N), RHS(N) */

/*     CALL DTOUT( N, NELT, IA, JA, A, ISYM, SOLN, RHS, IUNIT, JOB ) */

/* *Arguments: */
/* N      :IN       Integer */
/*         Order of the Matrix. */
/* NELT   :IN       Integer. */
/*         Number of non-zeros stored in A. */
/* IA     :IN       Integer IA(NELT). */
/* JA     :IN       Integer JA(NELT). */
/* A      :IN       Double Precision A(NELT). */
/*         These arrays should hold the matrix A in the SLAP */
/*         Triad format.  See "Description", below. */
/* ISYM   :IN       Integer. */
/*         Flag to indicate symmetric storage format. */
/*         If ISYM=0, all non-zero entries of the matrix are stored. */
/*         If ISYM=1, the matrix is symmetric, and only the lower */
/*         triangle of the matrix is stored. */
/* SOLN   :IN       Double Precision SOLN(N). */
/*         The solution to the linear system, if known.  This array */
/*         is accessed if and only if JOB is set to print it out, */
/*         see below. */
/* RHS    :IN       Double Precision RHS(N). */
/*         The right hand side vector.  This array is accessed if and */
/*         only if JOB is set to print it out, see below. */
/* IUNIT  :IN       Integer. */
/*         Fortran logical I/O device unit number to write the matrix */
/*         to.  This unit must be connected in a system dependent fashion */
/*         to a file or the console or you will get a nasty message */
/*         from the Fortran I/O libraries. */
/* JOB    :IN       Integer. */
/*         Flag indicating what I/O operations to perform. */
/*         JOB = 0 => Print only the matrix. */
/*             = 1 => Print matrix and RHS. */
/*             = 2 => Print matrix and SOLN. */
/*             = 3 => Print matrix, RHS and SOLN. */

/* *Description: */
/*       The format for the output is as follows.  On  the first line */
/*       are counters and flags: N, NELT, ISYM, IRHS, ISOLN.  N, NELT */
/*       and ISYM are described above.  IRHS is  a flag indicating if */
/*       the RHS was  written out (1 is  yes, 0 is  no).  ISOLN  is a */
/*       flag indicating if the SOLN was written out  (1 is yes, 0 is */
/*       no).  The format for the fist line is: 5i10.  Then comes the */
/*       NELT Triad's IA(I), JA(I) and A(I), I = 1, NELT.  The format */
/*       for  these lines is   :  1X,I5,1X,I5,1X,D16.7.   Then  comes */
/*       RHS(I), I = 1, N, if IRHS = 1.  Then  comes SOLN(I), I  = 1, */
/*       N, if ISOLN = 1.  The format for these lines is: 1X,D16.7. */

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

/* *Cautions: */
/*     This routine will attempt to write to the Fortran logical output */
/*     unit IUNIT, if IUNIT .ne. 0.  Thus, the user must make sure that */
/*     this logical unit is attached to a file or terminal before calling */
/*     this routine with a non-zero value for IUNIT.  This routine does */
/*     not check for the validity of a non-zero IUNIT unit number. */
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
/*   921007  Changed E's to D's in formats.  (FNF) */
/*   930701  Updated CATEGORY section.  (FNF, WRB) */
/* ***END PROLOGUE  DTOUT */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/* ***FIRST EXECUTABLE STATEMENT  DTOUT */

/*         If RHS and SOLN are to be printed also. */
/*         Write out the information heading. */

    /* Parameter adjustments */
    --rhs;
    --soln;
    --a;
    --ja;
    --ia;

    /* Function Body */
    irhs = 0;
    isoln = 0;
    if (*job == 1 || *job == 3) {
	irhs = 1;
    }
    if (*job > 1) {
	isoln = 1;
    }
    io___3.ciunit = *iunit;
    s_wsfe(&io___3);
    do_fio(&c__1, (char *)&(*n), (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&(*nelt), (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&(*isym), (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&irhs, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&isoln, (ftnlen)sizeof(integer));
    e_wsfe();

/*         Write out the matrix non-zeros in Triad format. */
    i__1 = *nelt;
    for (i__ = 1; i__ <= i__1; ++i__) {
	io___5.ciunit = *iunit;
	s_wsfe(&io___5);
	do_fio(&c__1, (char *)&ia[i__], (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&ja[i__], (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&a[i__], (ftnlen)sizeof(doublereal));
	e_wsfe();
/* L10: */
    }

/*         If requested, write out the rhs. */
    if (irhs == 1) {
	io___6.ciunit = *iunit;
	s_wsfe(&io___6);
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    do_fio(&c__1, (char *)&rhs[i__], (ftnlen)sizeof(doublereal));
	}
	e_wsfe();
    }

/*         If requested, write out the solution. */
    if (isoln == 1) {
	io___7.ciunit = *iunit;
	s_wsfe(&io___7);
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    do_fio(&c__1, (char *)&soln[i__], (ftnlen)sizeof(doublereal));
	}
	e_wsfe();
    }
    return 0;
/* ------------- LAST LINE OF DTOUT FOLLOWS ---------------------------- */
} /* dtout_ */

