/* hvnrm.f -- translated by f2c (version 12.02.01).
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

/* DECK HVNRM */
doublereal hvnrm_(real *v, integer *ncomp)
{
    /* System generated locals */
    integer i__1;
    real ret_val, r__1, r__2, r__3;

    /* Local variables */
    static integer k;

/* ***BEGIN PROLOGUE  HVNRM */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to DEABM, DEBDF and DERKF */
/* ***LIBRARY   SLATEC */
/* ***TYPE      SINGLE PRECISION (HVNRM-S, DHVNRM-D) */
/* ***AUTHOR  Watts, H. A., (SNLA) */
/* ***DESCRIPTION */

/*     Compute the maximum norm of the vector V(*) of length NCOMP and */
/*     return the result as HVNRM. */

/* ***SEE ALSO  DEABM, DEBDF, DERKF */
/* ***ROUTINES CALLED  (NONE) */
/* ***REVISION HISTORY  (YYMMDD) */
/*   800501  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890831  Modified array declarations.  (WRB) */
/*   891024  Changed routine name from VNORM to HVNRM.  (WRB) */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900328  Added TYPE section.  (WRB) */
/*   910722  Updated AUTHOR section.  (ALS) */
/* ***END PROLOGUE  HVNRM */
/* ***FIRST EXECUTABLE STATEMENT  HVNRM */
    /* Parameter adjustments */
    --v;

    /* Function Body */
    ret_val = 0.f;
    i__1 = *ncomp;
    for (k = 1; k <= i__1; ++k) {
/* L10: */
/* Computing MAX */
	r__2 = ret_val, r__3 = (r__1 = v[k], dabs(r__1));
	ret_val = dmax(r__2,r__3);
    }
    return ret_val;
} /* hvnrm_ */

