/* dwnlt2.f -- translated by f2c (version 12.02.01).
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

/* DECK DWNLT2 */
logical dwnlt2_(integer *me, integer *mend, integer *ir, doublereal *factor, 
	doublereal *tau, doublereal *scale, doublereal *wic)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1;
    logical ret_val;

    /* Local variables */
    static integer j;
    static doublereal t, rn, sn;

/* ***BEGIN PROLOGUE  DWNLT2 */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to WNLIT */
/* ***LIBRARY   SLATEC */
/* ***TYPE      DOUBLE PRECISION (WNLT2-S, DWNLT2-D) */
/* ***AUTHOR  Hanson, R. J., (SNLA) */
/*           Haskell, K. H., (SNLA) */
/* ***DESCRIPTION */

/*     To test independence of incoming column. */

/*     Test the column IC to determine if it is linearly independent */
/*     of the columns already in the basis.  In the initial tri. step, */
/*     we usually want the heavy weight ALAMDA to be included in the */
/*     test for independence.  In this case, the value of FACTOR will */
/*     have been set to 1.E0 before this procedure is invoked. */
/*     In the potentially rank deficient problem, the value of FACTOR */
/*     will have been set to ALSQ=ALAMDA**2 to remove the effect of the */
/*     heavy weight from the test for independence. */

/*     Write new column as partitioned vector */
/*           (A1)  number of components in solution so far = NIV */
/*           (A2)  M-NIV components */
/*     And compute  SN = inverse weighted length of A1 */
/*                  RN = inverse weighted length of A2 */
/*     Call the column independent when RN .GT. TAU*SN */

/* ***SEE ALSO  DWNLIT */
/* ***ROUTINES CALLED  (NONE) */
/* ***REVISION HISTORY  (YYMMDD) */
/*   790701  DATE WRITTEN */
/*   890620  Code extracted from WNLIT and made a subroutine.  (RWC)) */
/*   900604  DP version created from SP version.  (RWC) */
/* ***END PROLOGUE  DWNLT2 */


/* ***FIRST EXECUTABLE STATEMENT  DWNLT2 */
    /* Parameter adjustments */
    --wic;
    --scale;

    /* Function Body */
    sn = 0.f;
    rn = 0.f;
    i__1 = *mend;
    for (j = 1; j <= i__1; ++j) {
	t = scale[j];
	if (j <= *me) {
	    t /= *factor;
	}
/* Computing 2nd power */
	d__1 = wic[j];
	t *= d__1 * d__1;

	if (j < *ir) {
	    sn += t;
	} else {
	    rn += t;
	}
/* L10: */
    }
/* Computing 2nd power */
    d__1 = *tau;
    ret_val = rn > sn * (d__1 * d__1);
    return ret_val;
} /* dwnlt2_ */

