/* idloc.f -- translated by f2c (version 12.02.01).
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

static integer c__55 = 55;
static integer c__1 = 1;

/* DECK IDLOC */
integer idloc_(integer *loc, doublereal *sx, integer *ix)
{
    /* System generated locals */
    integer ret_val, i__1;

    /* Local variables */
    static integer k, np, lpg, key, lmx, lmxm1, ipage, itemp;
    extern /* Subroutine */ int dprwpg_(integer *, integer *, integer *, 
	    doublereal *, integer *), xermsg_(char *, char *, char *, integer 
	    *, integer *, ftnlen, ftnlen, ftnlen);

/* ***BEGIN PROLOGUE  IDLOC */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to DSPLP */
/* ***LIBRARY   SLATEC */
/* ***TYPE      DOUBLE PRECISION (IPLOC-S, IDLOC-D) */
/* ***KEYWORDS  RELATIVE ADDRESS DETERMINATION FUNCTION, SLATEC */
/* ***AUTHOR  Boland, W. Robert, (LANL) */
/*           Nicol, Tom, (University of British Columbia) */
/* ***DESCRIPTION */

/*   Given a "virtual" location,  IDLOC returns the relative working */
/*   address of the vector component stored in SX, IX.  Any necessary */
/*   page swaps are performed automatically for the user in this */
/*   function subprogram. */

/*   LOC       is the "virtual" address of the data to be retrieved. */
/*   SX ,IX    represent the matrix where the data is stored. */

/* ***SEE ALSO  DSPLP */
/* ***ROUTINES CALLED  DPRWPG, XERMSG */
/* ***REVISION HISTORY  (YYMMDD) */
/*   890606  DATE WRITTEN */
/*   890606  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ) */
/*   910731  Added code to set IDLOC to 0 if LOC is non-positive.  (WRB) */
/* ***END PROLOGUE  IDLOC */
/* ***FIRST EXECUTABLE STATEMENT  IDLOC */
    /* Parameter adjustments */
    --ix;
    --sx;

    /* Function Body */
    if (*loc <= 0) {
	xermsg_("SLATEC", "IDLOC", "A value of LOC, the first argument, .LE."
		" 0 was encountered", &c__55, &c__1, (ftnlen)6, (ftnlen)5, (
		ftnlen)58);
	ret_val = 0;
	return ret_val;
    }

/*     Two cases exist:  (1.LE.LOC.LE.K) .OR. (LOC.GT.K). */

    k = ix[3] + 4;
    lmx = ix[1];
    lmxm1 = lmx - 1;
    if (*loc <= k) {
	ret_val = *loc;
	return ret_val;
    }

/*     Compute length of the page, starting address of the page, page */
/*     number and relative working address. */

    lpg = lmx - k;
    itemp = *loc - k - 1;
    ipage = itemp / lpg + 1;
    ret_val = itemp % lpg + k + 1;
    np = (i__1 = ix[lmxm1], abs(i__1));

/*     Determine if a page fault has occurred.  If so, write page NP */
/*     and read page IPAGE.  Write the page only if it has been */
/*     modified. */

    if (ipage != np) {
	if (sx[lmx] == 1.f) {
	    sx[lmx] = 0.f;
	    key = 2;
	    dprwpg_(&key, &np, &lpg, &sx[1], &ix[1]);
	}
	key = 1;
	dprwpg_(&key, &ipage, &lpg, &sx[1], &ix[1]);
    }
    return ret_val;
} /* idloc_ */

