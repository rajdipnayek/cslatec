/* schkw.f -- translated by f2c (version 12.02.01).
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
static integer c__2 = 2;
static integer c__7 = 7;

/* DECK SCHKW */
/* Subroutine */ int schkw_(char *name__, integer *lociw, integer *leniw, 
	integer *locw, integer *lenw, integer *ierr, integer *iter, real *err,
	 ftnlen name_len)
{
    /* System generated locals */
    address a__1[7];
    integer i__1[7];
    char ch__1[89], ch__2[86];

    /* Local variables */
    static char xern1[8], xern2[8];
    extern doublereal r1mach_(integer *);
    static char xernam[8];
    extern /* Subroutine */ int xermsg_(char *, char *, char *, integer *, 
	    integer *, ftnlen, ftnlen, ftnlen);

    /* Fortran I/O blocks */
    static icilist io___3 = { 0, xern1, 0, "(I8)", 8, 1 };
    static icilist io___5 = { 0, xern2, 0, "(I8)", 8, 1 };
    static icilist io___6 = { 0, xern1, 0, "(I8)", 8, 1 };
    static icilist io___7 = { 0, xern2, 0, "(I8)", 8, 1 };


/* ***BEGIN PROLOGUE  SCHKW */
/* ***SUBSIDIARY */
/* ***PURPOSE  SLAP WORK/IWORK Array Bounds Checker. */
/*            This routine checks the work array lengths and interfaces */
/*            to the SLATEC error handler if a problem is found. */
/* ***LIBRARY   SLATEC (SLAP) */
/* ***CATEGORY  R2 */
/* ***TYPE      SINGLE PRECISION (SCHKW-S, DCHKW-D) */
/* ***KEYWORDS  ERROR CHECKING, SLAP, WORKSPACE CHECKING */
/* ***AUTHOR  Seager, Mark K., (LLNL) */
/*             Lawrence Livermore National Laboratory */
/*             PO BOX 808, L-60 */
/*             Livermore, CA 94550 (510) 423-3141 */
/*             seager@llnl.gov */
/* ***DESCRIPTION */

/* *Usage: */
/*     CHARACTER*(*) NAME */
/*     INTEGER LOCIW, LENIW, LOCW, LENW, IERR, ITER */
/*     REAL    ERR */

/*     CALL SCHKW( NAME, LOCIW, LENIW, LOCW, LENW, IERR, ITER, ERR ) */

/* *Arguments: */
/* NAME   :IN       Character*(*). */
/*         Name of the calling routine.  This is used in the output */
/*         message, if an error is detected. */
/* LOCIW  :IN       Integer. */
/*         Location of the first free element in the integer workspace */
/*         array. */
/* LENIW  :IN       Integer. */
/*         Length of the integer workspace array. */
/* LOCW   :IN       Integer. */
/*         Location of the first free element in the real workspace */
/*         array. */
/* LENRW  :IN       Integer. */
/*         Length of the real workspace array. */
/* IERR   :OUT      Integer. */
/*         Return error flag. */
/*               IERR = 0 => All went well. */
/*               IERR = 1 => Insufficient storage allocated for */
/*                           WORK or IWORK. */
/* ITER   :OUT      Integer. */
/*         Set to zero on return. */
/* ERR    :OUT      Real. */
/*         Set to the smallest positive magnitude if all went well. */
/*         Set to a very large number if an error is detected. */

/* ***REFERENCES  (NONE) */
/* ***ROUTINES CALLED  R1MACH, XERMSG */
/* ***REVISION HISTORY  (YYMMDD) */
/*   880225  DATE WRITTEN */
/*   881213  Previous REVISION DATE */
/*   890915  Made changes requested at July 1989 CML Meeting.  (MKS) */
/*   890922  Numerous changes to prologue to make closer to SLATEC */
/*           standard.  (FNF) */
/*   890929  Numerous changes to reduce SP/DP differences.  (FNF) */
/*   900805  Changed XERRWV calls to calls to XERMSG.  (RWC) */
/*   910411  Prologue converted to Version 4.0 format.  (BAB) */
/*   910502  Corrected XERMSG calls to satisfy Section 6.2.2 of ANSI */
/*           X3.9-1978.  (FNF) */
/*   910506  Made subsidiary.  (FNF) */
/*   920511  Added complete declaration section.  (WRB) */
/*   921015  Added code to initialize ITER and ERR when IERR=0.  (FNF) */
/* ***END PROLOGUE  SCHKW */
/*     .. Scalar Arguments .. */
/*     .. Local Scalars .. */
/*     .. External Functions .. */
/*     .. External Subroutines .. */
/* ***FIRST EXECUTABLE STATEMENT  SCHKW */

/*         Check the Integer workspace situation. */

    *ierr = 0;
    *iter = 0;
    *err = r1mach_(&c__1);
    if (*lociw > *leniw) {
	*ierr = 1;
	*err = r1mach_(&c__2);
	s_copy(xernam, name__, (ftnlen)8, name_len);
	s_wsfi(&io___3);
	do_fio(&c__1, (char *)&(*lociw), (ftnlen)sizeof(integer));
	e_wsfi();
	s_wsfi(&io___5);
	do_fio(&c__1, (char *)&(*leniw), (ftnlen)sizeof(integer));
	e_wsfi();
/* Writing concatenation */
	i__1[0] = 3, a__1[0] = "In ";
	i__1[1] = 8, a__1[1] = xernam;
	i__1[2] = 33, a__1[2] = ", INTEGER work array too short.  ";
	i__1[3] = 12, a__1[3] = "IWORK needs ";
	i__1[4] = 8, a__1[4] = xern1;
	i__1[5] = 17, a__1[5] = "; have allocated ";
	i__1[6] = 8, a__1[6] = xern2;
	s_cat(ch__1, a__1, i__1, &c__7, (ftnlen)89);
	xermsg_("SLATEC", "SCHKW", ch__1, &c__1, &c__1, (ftnlen)6, (ftnlen)5, 
		(ftnlen)89);
    }

/*         Check the Real workspace situation. */
    if (*locw > *lenw) {
	*ierr = 1;
	*err = r1mach_(&c__2);
	s_copy(xernam, name__, (ftnlen)8, name_len);
	s_wsfi(&io___6);
	do_fio(&c__1, (char *)&(*locw), (ftnlen)sizeof(integer));
	e_wsfi();
	s_wsfi(&io___7);
	do_fio(&c__1, (char *)&(*lenw), (ftnlen)sizeof(integer));
	e_wsfi();
/* Writing concatenation */
	i__1[0] = 3, a__1[0] = "In ";
	i__1[1] = 8, a__1[1] = xernam;
	i__1[2] = 30, a__1[2] = ", REAL work array too short.  ";
	i__1[3] = 12, a__1[3] = "RWORK needs ";
	i__1[4] = 8, a__1[4] = xern1;
	i__1[5] = 17, a__1[5] = "; have allocated ";
	i__1[6] = 8, a__1[6] = xern2;
	s_cat(ch__2, a__1, i__1, &c__7, (ftnlen)86);
	xermsg_("SLATEC", "SCHKW", ch__2, &c__1, &c__1, (ftnlen)6, (ftnlen)5, 
		(ftnlen)86);
    }
    return 0;
/* ------------- LAST LINE OF SCHKW FOLLOWS ---------------------------- */
} /* schkw_ */

