/* sclosm.f -- translated by f2c (version 12.02.01).
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
static integer c__100 = 100;
static integer c__2 = 2;

/* DECK SCLOSM */
/* Subroutine */ int sclosm_(integer *ipage)
{
    /* System generated locals */
    address a__1[2];
    integer i__1[2];
    char ch__1[40];
    cllist cl__1;

    /* Local variables */
    static integer ios;
    static char xern1[8];
    static integer ipagef;
    extern /* Subroutine */ int xermsg_(char *, char *, char *, integer *, 
	    integer *, ftnlen, ftnlen, ftnlen);

    /* Fortran I/O blocks */
    static icilist io___4 = { 0, xern1, 0, "(I8)", 8, 1 };


/* ***BEGIN PROLOGUE  SCLOSM */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to SPLP */
/* ***LIBRARY   SLATEC */
/* ***TYPE      ALL (SCLOSM-A) */
/* ***AUTHOR  (UNKNOWN) */
/* ***DESCRIPTION */

/*     1. UNLOAD, RELEASE, OR CLOSE UNIT NUMBER IPAGEF. */

/* ***SEE ALSO  SPLP */
/* ***ROUTINES CALLED  XERMSG */
/* ***REVISION HISTORY  (YYMMDD) */
/*   811215  DATE WRITTEN */
/*   890605  Corrected references to XERRWV.  (WRB) */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900402  Added TYPE section.  (WRB) */
/*   900510  Convert XERRWV calls to XERMSG calls.  (RWC) */
/* ***END PROLOGUE  SCLOSM */

/* ***FIRST EXECUTABLE STATEMENT  SCLOSM */
    ipagef = *ipage;
    cl__1.cerr = 1;
    cl__1.cunit = ipagef;
    cl__1.csta = "KEEP";
    ios = f_clos(&cl__1);
    if (ios != 0) {
	goto L100;
    }
    return 0;

L100:
    s_wsfi(&io___4);
    do_fio(&c__1, (char *)&ios, (ftnlen)sizeof(integer));
    e_wsfi();
/* Writing concatenation */
    i__1[0] = 32, a__1[0] = "IN SPLP, CLOSE HAS ERROR FLAG = ";
    i__1[1] = 8, a__1[1] = xern1;
    s_cat(ch__1, a__1, i__1, &c__2, (ftnlen)40);
    xermsg_("SLATEC", "SCLOSM", ch__1, &c__100, &c__1, (ftnlen)6, (ftnlen)6, (
	    ftnlen)40);
    return 0;
} /* sclosm_ */

