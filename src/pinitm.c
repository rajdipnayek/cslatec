/* pinitm.f -- translated by f2c (version 12.02.01).
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

/* DECK PINITM */
/* Subroutine */ int pinitm_(integer *m, integer *n, real *sx, integer *ix, 
	integer *lmx, integer *ipagef)
{
    /* Initialized data */

    static real zero = 0.f;
    static real one = 1.f;

    static integer i__, lp4, n20012, n20008, nerr, iopt;
    extern /* Subroutine */ int xermsg_(char *, char *, char *, integer *, 
	    integer *, ftnlen, ftnlen, ftnlen);

/* ***BEGIN PROLOGUE  PINITM */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to SPLP */
/* ***LIBRARY   SLATEC */
/* ***TYPE      SINGLE PRECISION (PINITM-S, DPINTM-D) */
/* ***AUTHOR  Hanson, R. J., (SNLA) */
/*           Wisniewski, J. A., (SNLA) */
/* ***DESCRIPTION */

/*     PINITM LIMITS THE TYPE OF STORAGE TO A SEQUENTIAL SCHEME. */
/*     THE MATRIX IS STORED BY COLUMNS. */
/*     SPARSE MATRIX INITIALIZATION SUBROUTINE. */

/*            M=NUMBER OF ROWS OF THE MATRIX. */
/*            N=NUMBER OF COLUMNS OF THE MATRIX. */
/*  SX(*),IX(*)=THE WORK ARRAYS WHICH ARE USED TO STORE THE SPARSE */
/*              MATRIX.  THESE ARRAYS ARE AUTOMATICALLY MAINTAINED BY */
/*              THE PACKAGE FOR THE USER. */
/*          LMX=LENGTH OF THE WORK ARRAY SX(*). */
/*              LMX MUST BE AT LEAST N+7 WHERE */
/*              FOR GREATEST EFFICIENCY LMX SHOULD BE AT LEAST N+NZ+6 */
/*              WHERE NZ IS THE MAXIMUM NUMBER OF NONZEROES TO BE */
/*              STORED IN THE MATRIX.  VALUES OF LMX BETWEEN N+7 AND */
/*              N+NZ+6 WILL CAUSE DEMAND PAGING TO OCCUR. */
/*              THIS IS IMPLEMENTED BY THE PACKAGE. */
/*              IX(*) MUST BE DIMENSIONED AT LEAST LMX */
/*      IPAGEF=UNIT NUMBER WHERE DEMAND PAGES WILL BE STORED. */

/*     THIS SUBROUTINE IS A MODIFICATION OF THE SUBROUTINE LINITM, */
/*     SANDIA LABS. REPT. SAND78-0785. */
/*     MODIFICATIONS BY K.L. HIEBERT AND R.J. HANSON */
/*     REVISED 811130-1000 */
/*     REVISED YYMMDD-HHMM */

/* ***SEE ALSO  SPLP */
/* ***ROUTINES CALLED  XERMSG */
/* ***REVISION HISTORY  (YYMMDD) */
/*   811215  DATE WRITTEN */
/*   890831  Modified array declarations.  (WRB) */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ) */
/*   900328  Added TYPE section.  (WRB) */
/*   910403  Updated AUTHOR and DESCRIPTION sections.  (WRB) */
/* ***END PROLOGUE  PINITM */
    /* Parameter adjustments */
    --ix;
    --sx;

    /* Function Body */
/* ***FIRST EXECUTABLE STATEMENT  PINITM */
    iopt = 1;

/*     CHECK FOR INPUT ERRORS. */

    if (! (*m <= 0 || *n <= 0)) {
	goto L20002;
    }
    nerr = 55;
    xermsg_("SLATEC", "PINITM", "MATRIX DIMENSION M OR N .LE. 0.", &nerr, &
	    iopt, (ftnlen)6, (ftnlen)6, (ftnlen)31);

/*     VERIFY IF VALUE OF LMX IS LARGE ENOUGH. */

L20002:
    if (! (*lmx < *n + 7)) {
	goto L20005;
    }
    nerr = 55;
    xermsg_("SLATEC", "PINITM", "THE VALUE OF LMX IS TOO SMALL.", &nerr, &
	    iopt, (ftnlen)6, (ftnlen)6, (ftnlen)30);

/*     INITIALIZE DATA STRUCTURE INDEPENDENT VALUES. */

L20005:
    sx[1] = zero;
    sx[2] = zero;
    sx[3] = (real) (*ipagef);
    ix[1] = *lmx;
    ix[2] = *m;
    ix[3] = *n;
    ix[4] = 0;
    sx[*lmx - 1] = zero;
    sx[*lmx] = -one;
    ix[*lmx - 1] = -1;
    lp4 = *n + 4;

/*     INITIALIZE DATA STRUCTURE DEPENDENT VALUES. */

    i__ = 4;
    n20008 = lp4;
    goto L20009;
L20008:
    ++i__;
L20009:
    if (n20008 - i__ < 0) {
	goto L20010;
    }
    sx[i__] = zero;
    goto L20008;
L20010:
    i__ = 5;
    n20012 = lp4;
    goto L20013;
L20012:
    ++i__;
L20013:
    if (n20012 - i__ < 0) {
	goto L20014;
    }
    ix[i__] = lp4;
    goto L20012;
L20014:
    sx[*n + 5] = zero;
    ix[*n + 5] = 0;
    ix[*lmx] = 0;

/*     INITIALIZATION COMPLETE. */

    return 0;
} /* pinitm_ */

