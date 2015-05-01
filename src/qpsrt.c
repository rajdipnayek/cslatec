/* qpsrt.f -- translated by f2c (version 12.02.01).
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

/* DECK QPSRT */
/* Subroutine */ int qpsrt_(integer *limit, integer *last, integer *maxerr, 
	real *ermax, real *elist, integer *iord, integer *nrmax)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__, j, k, ido, ibeg, jbnd, isucc, jupbn;
    static real errmin, errmax;

/* ***BEGIN PROLOGUE  QPSRT */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to QAGE, QAGIE, QAGPE, QAGSE, QAWCE, QAWOE and */
/*            QAWSE */
/* ***LIBRARY   SLATEC */
/* ***TYPE      SINGLE PRECISION (QPSRT-S, DQPSRT-D) */
/* ***KEYWORDS  SEQUENTIAL SORTING */
/* ***AUTHOR  (UNKNOWN) */
/* ***DESCRIPTION */

/* 1.        QPSRT */
/*           Ordering Routine */
/*              Standard FORTRAN Subroutine */
/*              REAL Version */

/* 2.        PURPOSE */
/*              This routine maintains the descending ordering */
/*              in the list of the local error estimates resulting from */
/*              the interval subdivision process. At each call two error */
/*              estimates are inserted using the sequential search */
/*              method, top-down for the largest error estimate */
/*              and bottom-up for the smallest error estimate. */

/* 3.        CALLING SEQUENCE */
/*              CALL QPSRT(LIMIT,LAST,MAXERR,ERMAX,ELIST,IORD,NRMAX) */

/*           PARAMETERS (MEANING AT OUTPUT) */
/*              LIMIT  - INTEGER */
/*                       Maximum number of error estimates the list */
/*                       can contain */

/*              LAST   - INTEGER */
/*                       Number of error estimates currently */
/*                       in the list */

/*              MAXERR - INTEGER */
/*                       MAXERR points to the NRMAX-th largest error */
/*                       estimate currently in the list */

/*              ERMAX  - REAL */
/*                       NRMAX-th largest error estimate */
/*                       ERMAX = ELIST(MAXERR) */

/*              ELIST  - REAL */
/*                       Vector of dimension LAST containing */
/*                       the error estimates */

/*              IORD   - INTEGER */
/*                       Vector of dimension LAST, the first K */
/*                       elements of which contain pointers */
/*                       to the error estimates, such that */
/*                       ELIST(IORD(1)),... , ELIST(IORD(K)) */
/*                       form a decreasing sequence, with */
/*                       K = LAST if LAST.LE.(LIMIT/2+2), and */
/*                       K = LIMIT+1-LAST otherwise */

/*              NRMAX  - INTEGER */
/*                       MAXERR = IORD(NRMAX) */

/* ***SEE ALSO  QAGE, QAGIE, QAGPE, QAGSE, QAWCE, QAWOE, QAWSE */
/* ***ROUTINES CALLED  (NONE) */
/* ***REVISION HISTORY  (YYMMDD) */
/*   800101  DATE WRITTEN */
/*   890831  Modified array declarations.  (WRB) */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900328  Added TYPE section.  (WRB) */
/* ***END PROLOGUE  QPSRT */


/*           CHECK WHETHER THE LIST CONTAINS MORE THAN */
/*           TWO ERROR ESTIMATES. */

/* ***FIRST EXECUTABLE STATEMENT  QPSRT */
    /* Parameter adjustments */
    --iord;
    --elist;

    /* Function Body */
    if (*last > 2) {
	goto L10;
    }
    iord[1] = 1;
    iord[2] = 2;
    goto L90;

/*           THIS PART OF THE ROUTINE IS ONLY EXECUTED */
/*           IF, DUE TO A DIFFICULT INTEGRAND, SUBDIVISION */
/*           INCREASED THE ERROR ESTIMATE. IN THE NORMAL CASE */
/*           THE INSERT PROCEDURE SHOULD START AFTER THE */
/*           NRMAX-TH LARGEST ERROR ESTIMATE. */

L10:
    errmax = elist[*maxerr];
    if (*nrmax == 1) {
	goto L30;
    }
    ido = *nrmax - 1;
    i__1 = ido;
    for (i__ = 1; i__ <= i__1; ++i__) {
	isucc = iord[*nrmax - 1];
/* ***JUMP OUT OF DO-LOOP */
	if (errmax <= elist[isucc]) {
	    goto L30;
	}
	iord[*nrmax] = isucc;
	--(*nrmax);
/* L20: */
    }

/*           COMPUTE THE NUMBER OF ELEMENTS IN THE LIST TO */
/*           BE MAINTAINED IN DESCENDING ORDER. THIS NUMBER */
/*           DEPENDS ON THE NUMBER OF SUBDIVISIONS STILL */
/*           ALLOWED. */

L30:
    jupbn = *last;
    if (*last > *limit / 2 + 2) {
	jupbn = *limit + 3 - *last;
    }
    errmin = elist[*last];

/*           INSERT ERRMAX BY TRAVERSING THE LIST TOP-DOWN, */
/*           STARTING COMPARISON FROM THE ELEMENT ELIST(IORD(NRMAX+1)). */

    jbnd = jupbn - 1;
    ibeg = *nrmax + 1;
    if (ibeg > jbnd) {
	goto L50;
    }
    i__1 = jbnd;
    for (i__ = ibeg; i__ <= i__1; ++i__) {
	isucc = iord[i__];
/* ***JUMP OUT OF DO-LOOP */
	if (errmax >= elist[isucc]) {
	    goto L60;
	}
	iord[i__ - 1] = isucc;
/* L40: */
    }
L50:
    iord[jbnd] = *maxerr;
    iord[jupbn] = *last;
    goto L90;

/*           INSERT ERRMIN BY TRAVERSING THE LIST BOTTOM-UP. */

L60:
    iord[i__ - 1] = *maxerr;
    k = jbnd;
    i__1 = jbnd;
    for (j = i__; j <= i__1; ++j) {
	isucc = iord[k];
/* ***JUMP OUT OF DO-LOOP */
	if (errmin < elist[isucc]) {
	    goto L80;
	}
	iord[k + 1] = isucc;
	--k;
/* L70: */
    }
    iord[i__] = *last;
    goto L90;
L80:
    iord[k + 1] = *last;

/*           SET MAXERR AND ERMAX. */

L90:
    *maxerr = iord[*nrmax];
    *ermax = elist[*maxerr];
    return 0;
} /* qpsrt_ */

