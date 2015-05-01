/* swritp.f -- translated by f2c (version 12.02.01).
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
static integer c__4 = 4;

/* DECK SWRITP */
/* Subroutine */ int swritp_(integer *ipage, integer *list, real *rlist, 
	integer *lpage, integer *irec)
{
    /* System generated locals */
    address a__1[4];
    integer i__1, i__2, i__3[4];
    char ch__1[40];

    /* Local variables */
    static integer i__, lpg;
    static char xern1[8], xern2[8];
    static integer irecn, ipagef;
    extern /* Subroutine */ int xermsg_(char *, char *, char *, integer *, 
	    integer *, ftnlen, ftnlen, ftnlen);

    /* Fortran I/O blocks */
    static cilist io___4 = { 1, 0, 0, 0, 0 };
    static cilist io___6 = { 1, 0, 0, 0, 0 };
    static icilist io___8 = { 0, xern1, 0, "(I8)", 8, 1 };
    static icilist io___10 = { 0, xern2, 0, "(I8)", 8, 1 };


/* ***BEGIN PROLOGUE  SWRITP */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to SPLP */
/* ***LIBRARY   SLATEC */
/* ***TYPE      SINGLE PRECISION (SWRITP-S, DWRITP-D) */
/* ***AUTHOR  (UNKNOWN) */
/* ***DESCRIPTION */

/*     WRITE RECORD NUMBER IRECN, OF LENGTH LPG, FROM STORAGE */
/*     ARRAY LIST(*) ONTO UNIT NUMBER IPAGEF. */
/*     WRITE RECORD NUMBER IRECN+1, OF LENGTH LPG, ONTO UNIT */
/*     NUMBER IPAGEF FROM THE STORAGE ARRAY RLIST(*). */

/*     TO CHANGE THIS PROGRAM UNIT TO DOUBLE PRECISION CHANGE */
/*     /REAL (12 BLANKS)/ TO /DOUBLE PRECISION/. */

/* ***SEE ALSO  SPLP */
/* ***ROUTINES CALLED  XERMSG */
/* ***REVISION HISTORY  (YYMMDD) */
/*   811215  DATE WRITTEN */
/*   890605  Corrected references to XERRWV.  (WRB) */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900328  Added TYPE section.  (WRB) */
/*   900510  Convert XERRWV calls to XERMSG calls.  (RWC) */
/* ***END PROLOGUE  SWRITP */
/* ***FIRST EXECUTABLE STATEMENT  SWRITP */
    /* Parameter adjustments */
    --rlist;
    --list;

    /* Function Body */
    ipagef = *ipage;
    lpg = *lpage;
    irecn = *irec;
    io___4.ciunit = ipagef;
    io___4.cirec = irecn;
    i__1 = s_wdue(&io___4);
    if (i__1 != 0) {
	goto L100;
    }
    i__2 = lpg;
    for (i__ = 1; i__ <= i__2; ++i__) {
	i__1 = do_uio(&c__1, (char *)&list[i__], (ftnlen)sizeof(integer));
	if (i__1 != 0) {
	    goto L100;
	}
    }
    i__1 = e_wdue();
    if (i__1 != 0) {
	goto L100;
    }
    io___6.ciunit = ipagef;
    io___6.cirec = irecn + 1;
    i__1 = s_wdue(&io___6);
    if (i__1 != 0) {
	goto L100;
    }
    i__2 = lpg;
    for (i__ = 1; i__ <= i__2; ++i__) {
	i__1 = do_uio(&c__1, (char *)&rlist[i__], (ftnlen)sizeof(real));
	if (i__1 != 0) {
	    goto L100;
	}
    }
    i__1 = e_wdue();
    if (i__1 != 0) {
	goto L100;
    }
    return 0;

L100:
    s_wsfi(&io___8);
    do_fio(&c__1, (char *)&lpg, (ftnlen)sizeof(integer));
    e_wsfi();
    s_wsfi(&io___10);
    do_fio(&c__1, (char *)&irecn, (ftnlen)sizeof(integer));
    e_wsfi();
/* Writing concatenation */
    i__3[0] = 15, a__1[0] = "IN SPLP, LGP = ";
    i__3[1] = 8, a__1[1] = xern1;
    i__3[2] = 9, a__1[2] = " IRECN = ";
    i__3[3] = 8, a__1[3] = xern2;
    s_cat(ch__1, a__1, i__3, &c__4, (ftnlen)40);
    xermsg_("SLATEC", "SWRITP", ch__1, &c__100, &c__1, (ftnlen)6, (ftnlen)6, (
	    ftnlen)40);
    return 0;
} /* swritp_ */

