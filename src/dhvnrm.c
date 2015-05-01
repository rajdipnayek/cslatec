/* dhvnrm.f -- translated by f2c (version 12.02.01).
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

/* DECK DHVNRM */
doublereal dhvnrm_(doublereal *v, integer *ncomp)
{
    /* System generated locals */
    integer i__1;
    doublereal ret_val, d__1, d__2, d__3;

    /* Local variables */
    static integer k;

/* ***BEGIN PROLOGUE  DHVNRM */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to DDEABM, DDEBDF and DDERKF */
/* ***LIBRARY   SLATEC */
/* ***TYPE      DOUBLE PRECISION (HVNRM-S, DHVNRM-D) */
/* ***AUTHOR  Watts, H. A., (SNLA) */
/* ***DESCRIPTION */

/*     Compute the maximum norm of the vector V(*) of length NCOMP and */
/*     return the result as DHVNRM */

/* ***SEE ALSO  DDEABM, DDEBDF, DDERKF */
/* ***ROUTINES CALLED  (NONE) */
/* ***REVISION HISTORY  (YYMMDD) */
/*   820301  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890831  Modified array declarations.  (WRB) */
/*   891024  Changed references from DVNORM to DHVNRM.  (WRB) */
/*   891024  Changed routine name from DVNORM to DHVNRM.  (WRB) */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900328  Added TYPE section.  (WRB) */
/*   910722  Updated AUTHOR section.  (ALS) */
/* ***END PROLOGUE  DHVNRM */

/* ***FIRST EXECUTABLE STATEMENT  DHVNRM */
    /* Parameter adjustments */
    --v;

    /* Function Body */
    ret_val = 0.;
    i__1 = *ncomp;
    for (k = 1; k <= i__1; ++k) {
/* Computing MAX */
	d__2 = ret_val, d__3 = (d__1 = v[k], abs(d__1));
	ret_val = max(d__2,d__3);
/* L10: */
    }
    return ret_val;
} /* dhvnrm_ */

