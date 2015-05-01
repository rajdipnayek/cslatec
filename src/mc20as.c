/* mc20as.f -- translated by f2c (version 12.02.01).
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

/* DECK MC20AS */
/* Subroutine */ int mc20as_(integer *nc, integer *maxa, real *a, integer *
	inum, integer *jptr, integer *jnum, integer *jdisp)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer i__, j, k, ja, jb, kr;
    static real ace;
    static integer ice, jce, loc;
    static real acep;
    static integer icep, jcep, null;

/* ***BEGIN PROLOGUE  MC20AS */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to SPLP */
/* ***LIBRARY   SLATEC */
/* ***TYPE      SINGLE PRECISION (MC20AS-S, MC20AD-D) */
/* ***AUTHOR  (UNKNOWN) */
/* ***DESCRIPTION */

/*     THIS SUBPROGRAM IS A SLIGHT MODIFICATION OF A SUBPROGRAM */
/*     FROM THE C. 1979 AERE HARWELL LIBRARY.  THE NAME OF THE */
/*     CORRESPONDING HARWELL CODE CAN BE OBTAINED BY DELETING */
/*     THE FINAL LETTER =S= IN THE NAMES USED HERE. */
/*     REVISED SEP. 13, 1979. */

/*     ROYALTIES HAVE BEEN PAID TO AERE-UK FOR USE OF THEIR CODES */
/*     IN THE PACKAGE GIVEN HERE.  ANY PRIMARY USAGE OF THE HARWELL */
/*     SUBROUTINES REQUIRES A ROYALTY AGREEMENT AND PAYMENT BETWEEN */
/*     THE USER AND AERE-UK.  ANY USAGE OF THE SANDIA WRITTEN CODES */
/*     SPLP( ) (WHICH USES THE HARWELL SUBROUTINES) IS PERMITTED. */

/* ***SEE ALSO  SPLP */
/* ***ROUTINES CALLED  (NONE) */
/* ***REVISION HISTORY  (YYMMDD) */
/*   811215  DATE WRITTEN */
/*   890831  Modified array declarations.  (WRB) */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900402  Added TYPE section.  (WRB) */
/* ***END PROLOGUE  MC20AS */
/* ***FIRST EXECUTABLE STATEMENT  MC20AS */
    /* Parameter adjustments */
    --jptr;
    --a;
    --inum;
    --jnum;

    /* Function Body */
    null = -(*jdisp);
/* **      CLEAR JPTR */
    i__1 = *nc;
    for (j = 1; j <= i__1; ++j) {
	jptr[j] = 0;
/* L10: */
    }
/* **      COUNT THE NUMBER OF ELEMENTS IN EACH COLUMN. */
    i__1 = *maxa;
    for (k = 1; k <= i__1; ++k) {
	j = jnum[k] + *jdisp;
	++jptr[j];
/* L20: */
    }
/* **      SET THE JPTR ARRAY */
    k = 1;
    i__1 = *nc;
    for (j = 1; j <= i__1; ++j) {
	kr = k + jptr[j];
	jptr[j] = k;
	k = kr;
/* L30: */
    }

/* **      REORDER THE ELEMENTS INTO COLUMN ORDER.  THE ALGORITHM IS AN */
/*        IN-PLACE SORT AND IS OF ORDER MAXA. */
    i__1 = *maxa;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*        ESTABLISH THE CURRENT ENTRY. */
	jce = jnum[i__] + *jdisp;
	if (jce == 0) {
	    goto L50;
	}
	ace = a[i__];
	ice = inum[i__];
/*        CLEAR THE LOCATION VACATED. */
	jnum[i__] = null;
/*        CHAIN FROM CURRENT ENTRY TO STORE ITEMS. */
	i__2 = *maxa;
	for (j = 1; j <= i__2; ++j) {
/*        CURRENT ENTRY NOT IN CORRECT POSITION.  DETERMINE CORRECT */
/*        POSITION TO STORE ENTRY. */
	    loc = jptr[jce];
	    ++jptr[jce];
/*        SAVE CONTENTS OF THAT LOCATION. */
	    acep = a[loc];
	    icep = inum[loc];
	    jcep = jnum[loc];
/*        STORE CURRENT ENTRY. */
	    a[loc] = ace;
	    inum[loc] = ice;
	    jnum[loc] = null;
/*        CHECK IF NEXT CURRENT ENTRY NEEDS TO BE PROCESSED. */
	    if (jcep == null) {
		goto L50;
	    }
/*        IT DOES.  COPY INTO CURRENT ENTRY. */
	    ace = acep;
	    ice = icep;
	    jce = jcep + *jdisp;
/* L40: */
	}

L50:
	;
    }

/* **      RESET JPTR VECTOR. */
    ja = 1;
    i__1 = *nc;
    for (j = 1; j <= i__1; ++j) {
	jb = jptr[j];
	jptr[j] = ja;
	ja = jb;
/* L60: */
    }
    return 0;
} /* mc20as_ */

