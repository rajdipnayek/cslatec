/* svout.f -- translated by f2c (version 12.02.01).
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

/* DECK SVOUT */
/* Subroutine */ int svout_(integer *n, real *sx, char *ifmt, integer *idigit,
	 ftnlen ifmt_len)
{
    /* Format strings */
    static char fmt_1000[] = "(1x,i4,\002 - \002,i4,1p,10e12.3)";
    static char fmt_1001[] = "(1x,i4,\002 - \002,i4,1x,1p,8e14.5)";
    static char fmt_1002[] = "(1x,i4,\002 - \002,i4,1x,1p,6e18.9)";
    static char fmt_1003[] = "(1x,i4,\002 - \002,i4,1x,1p,5e24.13)";

    /* System generated locals */
    integer i__1, i__2, i__3;
    cilist ci__1;

    /* Local variables */
    static integer i__, j, k1, k2, lout;
    extern integer i1mach_(integer *);
    static integer ndigit;

    /* Fortran I/O blocks */
    static cilist io___6 = { 0, 0, 0, fmt_1000, 0 };
    static cilist io___8 = { 0, 0, 0, fmt_1001, 0 };
    static cilist io___9 = { 0, 0, 0, fmt_1002, 0 };
    static cilist io___10 = { 0, 0, 0, fmt_1003, 0 };
    static cilist io___11 = { 0, 0, 0, fmt_1000, 0 };
    static cilist io___12 = { 0, 0, 0, fmt_1001, 0 };
    static cilist io___13 = { 0, 0, 0, fmt_1002, 0 };
    static cilist io___14 = { 0, 0, 0, fmt_1003, 0 };


/* ***BEGIN PROLOGUE  SVOUT */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to SPLP */
/* ***LIBRARY   SLATEC */
/* ***TYPE      SINGLE PRECISION (SVOUT-S, DVOUT-D) */
/* ***AUTHOR  (UNKNOWN) */
/* ***DESCRIPTION */

/*     SINGLE PRECISION VECTOR OUTPUT ROUTINE. */

/*  INPUT.. */

/*  N,SX(*) PRINT THE SINGLE PRECISION ARRAY SX(I),I=1,...,N, ON */
/*          OUTPUT UNIT LOUT. THE HEADING IN THE FORTRAN FORMAT */
/*          STATEMENT IFMT(*), DESCRIBED BELOW, IS PRINTED AS A FIRST */
/*          STEP. THE COMPONENTS SX(I) ARE INDEXED, ON OUTPUT, */
/*          IN A PLEASANT FORMAT. */
/*  IFMT(*) A FORTRAN FORMAT STATEMENT. THIS IS PRINTED ON OUTPUT */
/*          UNIT LOUT WITH THE VARIABLE FORMAT FORTRAN STATEMENT */
/*                WRITE(LOUT,IFMT) */
/*  IDIGIT  PRINT AT LEAST ABS(IDIGIT) DECIMAL DIGITS PER NUMBER. */
/*          THE SUBPROGRAM WILL CHOOSE THAT INTEGER 4,6,10 OR 14 */
/*          WHICH WILL PRINT AT LEAST ABS(IDIGIT) NUMBER OF */
/*          PLACES.  IF IDIGIT.LT.0, 72 PRINTING COLUMNS ARE UTILIZED */
/*          TO WRITE EACH LINE OF OUTPUT OF THE ARRAY SX(*). (THIS */
/*          CAN BE USED ON MOST TIME-SHARING TERMINALS). IF */
/*          IDIGIT.GE.0, 133 PRINTING COLUMNS ARE UTILIZED. (THIS CAN */
/*          BE USED ON MOST LINE PRINTERS). */

/*  EXAMPLE.. */

/*  PRINT AN ARRAY CALLED (COSTS OF PURCHASES) OF LENGTH 100 SHOWING */
/*  6 DECIMAL DIGITS PER NUMBER. THE USER IS RUNNING ON A TIME-SHARING */
/*  SYSTEM WITH A 72 COLUMN OUTPUT DEVICE. */

/*     DIMENSION COSTS(100) */
/*     N = 100 */
/*     IDIGIT = -6 */
/*     CALL SVOUT(N,COSTS,'(''1COSTS OF PURCHASES'')',IDIGIT) */

/* ***SEE ALSO  SPLP */
/* ***ROUTINES CALLED  I1MACH */
/* ***REVISION HISTORY  (YYMMDD) */
/*   811215  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   891107  Added comma after 1P edit descriptor in FORMAT */
/*           statements.  (WRB) */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900328  Added TYPE section.  (WRB) */
/* ***END PROLOGUE  SVOUT */

/*     GET THE UNIT NUMBER WHERE OUTPUT WILL BE WRITTEN. */
/* ***FIRST EXECUTABLE STATEMENT  SVOUT */
    /* Parameter adjustments */
    --sx;

    /* Function Body */
    j = 2;
    lout = i1mach_(&j);
    ci__1.cierr = 0;
    ci__1.ciunit = lout;
    ci__1.cifmt = ifmt;
    s_wsfe(&ci__1);
    e_wsfe();
    if (*n <= 0) {
	return 0;
    }
    ndigit = *idigit;
    if (*idigit == 0) {
	ndigit = 4;
    }
    if (*idigit >= 0) {
	goto L80;
    }

    ndigit = -(*idigit);
    if (ndigit > 4) {
	goto L20;
    }

    i__1 = *n;
    for (k1 = 1; k1 <= i__1; k1 += 5) {
/* Computing MIN */
	i__2 = *n, i__3 = k1 + 4;
	k2 = min(i__2,i__3);
	io___6.ciunit = lout;
	s_wsfe(&io___6);
	do_fio(&c__1, (char *)&k1, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&k2, (ftnlen)sizeof(integer));
	i__2 = k2;
	for (i__ = k1; i__ <= i__2; ++i__) {
	    do_fio(&c__1, (char *)&sx[i__], (ftnlen)sizeof(real));
	}
	e_wsfe();
/* L10: */
    }
    return 0;

L20:
    if (ndigit > 6) {
	goto L40;
    }

    i__1 = *n;
    for (k1 = 1; k1 <= i__1; k1 += 4) {
/* Computing MIN */
	i__2 = *n, i__3 = k1 + 3;
	k2 = min(i__2,i__3);
	io___8.ciunit = lout;
	s_wsfe(&io___8);
	do_fio(&c__1, (char *)&k1, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&k2, (ftnlen)sizeof(integer));
	i__2 = k2;
	for (i__ = k1; i__ <= i__2; ++i__) {
	    do_fio(&c__1, (char *)&sx[i__], (ftnlen)sizeof(real));
	}
	e_wsfe();
/* L30: */
    }
    return 0;

L40:
    if (ndigit > 10) {
	goto L60;
    }

    i__1 = *n;
    for (k1 = 1; k1 <= i__1; k1 += 3) {
/* Computing MIN */
	i__2 = *n, i__3 = k1 + 2;
	k2 = min(i__2,i__3);
	io___9.ciunit = lout;
	s_wsfe(&io___9);
	do_fio(&c__1, (char *)&k1, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&k2, (ftnlen)sizeof(integer));
	i__2 = k2;
	for (i__ = k1; i__ <= i__2; ++i__) {
	    do_fio(&c__1, (char *)&sx[i__], (ftnlen)sizeof(real));
	}
	e_wsfe();
/* L50: */
    }
    return 0;

L60:
    i__1 = *n;
    for (k1 = 1; k1 <= i__1; k1 += 2) {
/* Computing MIN */
	i__2 = *n, i__3 = k1 + 1;
	k2 = min(i__2,i__3);
	io___10.ciunit = lout;
	s_wsfe(&io___10);
	do_fio(&c__1, (char *)&k1, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&k2, (ftnlen)sizeof(integer));
	i__2 = k2;
	for (i__ = k1; i__ <= i__2; ++i__) {
	    do_fio(&c__1, (char *)&sx[i__], (ftnlen)sizeof(real));
	}
	e_wsfe();
/* L70: */
    }
    return 0;

L80:
    if (ndigit > 4) {
	goto L100;
    }

    i__1 = *n;
    for (k1 = 1; k1 <= i__1; k1 += 10) {
/* Computing MIN */
	i__2 = *n, i__3 = k1 + 9;
	k2 = min(i__2,i__3);
	io___11.ciunit = lout;
	s_wsfe(&io___11);
	do_fio(&c__1, (char *)&k1, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&k2, (ftnlen)sizeof(integer));
	i__2 = k2;
	for (i__ = k1; i__ <= i__2; ++i__) {
	    do_fio(&c__1, (char *)&sx[i__], (ftnlen)sizeof(real));
	}
	e_wsfe();
/* L90: */
    }
    return 0;

L100:
    if (ndigit > 6) {
	goto L120;
    }

    i__1 = *n;
    for (k1 = 1; k1 <= i__1; k1 += 8) {
/* Computing MIN */
	i__2 = *n, i__3 = k1 + 7;
	k2 = min(i__2,i__3);
	io___12.ciunit = lout;
	s_wsfe(&io___12);
	do_fio(&c__1, (char *)&k1, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&k2, (ftnlen)sizeof(integer));
	i__2 = k2;
	for (i__ = k1; i__ <= i__2; ++i__) {
	    do_fio(&c__1, (char *)&sx[i__], (ftnlen)sizeof(real));
	}
	e_wsfe();
/* L110: */
    }
    return 0;

L120:
    if (ndigit > 10) {
	goto L140;
    }

    i__1 = *n;
    for (k1 = 1; k1 <= i__1; k1 += 6) {
/* Computing MIN */
	i__2 = *n, i__3 = k1 + 5;
	k2 = min(i__2,i__3);
	io___13.ciunit = lout;
	s_wsfe(&io___13);
	do_fio(&c__1, (char *)&k1, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&k2, (ftnlen)sizeof(integer));
	i__2 = k2;
	for (i__ = k1; i__ <= i__2; ++i__) {
	    do_fio(&c__1, (char *)&sx[i__], (ftnlen)sizeof(real));
	}
	e_wsfe();
/* L130: */
    }
    return 0;

L140:
    i__1 = *n;
    for (k1 = 1; k1 <= i__1; k1 += 5) {
/* Computing MIN */
	i__2 = *n, i__3 = k1 + 4;
	k2 = min(i__2,i__3);
	io___14.ciunit = lout;
	s_wsfe(&io___14);
	do_fio(&c__1, (char *)&k1, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&k2, (ftnlen)sizeof(integer));
	i__2 = k2;
	for (i__ = k1; i__ <= i__2; ++i__) {
	    do_fio(&c__1, (char *)&sx[i__], (ftnlen)sizeof(real));
	}
	e_wsfe();
/* L150: */
    }
    return 0;
} /* svout_ */

