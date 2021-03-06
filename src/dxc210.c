/* dxc210.f -- translated by f2c (version 12.02.01).
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

/* Common Block Declarations */

struct {
    integer nlg102, mlg102, lg102[21];
} dxblk3_;

#define dxblk3_1 dxblk3_

/* Table of constant values */

static doublereal c_b9 = 10.;
static integer c__208 = 208;
static integer c__1 = 1;

/* DECK DXC210 */
/* Subroutine */ int dxc210_(integer *k, doublereal *z__, integer *j, integer 
	*ierror)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1;

    /* Local variables */
    static integer i__, m, ja, ka, ic, id, ii, it, ka1, ka2, nm1, np1;
    extern /* Subroutine */ int xermsg_(char *, char *, char *, integer *, 
	    integer *, ftnlen, ftnlen, ftnlen);

/* ***BEGIN PROLOGUE  DXC210 */
/* ***PURPOSE  To provide double-precision floating-point arithmetic */
/*            with an extended exponent range. */
/* ***LIBRARY   SLATEC */
/* ***CATEGORY  A3D */
/* ***TYPE      DOUBLE PRECISION (XC210-S, DXC210-D) */
/* ***KEYWORDS  EXTENDED-RANGE DOUBLE-PRECISION ARITHMETIC */
/* ***AUTHOR  Lozier, Daniel W., (National Bureau of Standards) */
/*           Smith, John M., (NBS and George Mason University) */
/* ***DESCRIPTION */
/*     INTEGER K, J */
/*     DOUBLE PRECISION Z */

/*                  GIVEN K THIS SUBROUTINE COMPUTES J AND Z */
/*                  SUCH THAT  RADIX**K = Z*10**J, WHERE Z IS IN */
/*                  THE RANGE 1/10 .LE. Z .LT. 1. */
/*                  THE VALUE OF Z WILL BE ACCURATE TO FULL */
/*                  DOUBLE-PRECISION PROVIDED THE NUMBER */
/*                  OF DECIMAL PLACES IN THE LARGEST */
/*                  INTEGER PLUS THE NUMBER OF DECIMAL */
/*                  PLACES CARRIED IN DOUBLE-PRECISION DOES NOT */
/*                  EXCEED 60. DXC210 IS CALLED BY SUBROUTINE */
/*                  DXCON WHEN NECESSARY. THE USER SHOULD */
/*                  NEVER NEED TO CALL DXC210 DIRECTLY. */

/* ***SEE ALSO  DXSET */
/* ***REFERENCES  (NONE) */
/* ***ROUTINES CALLED  XERMSG */
/* ***COMMON BLOCKS    DXBLK3 */
/* ***REVISION HISTORY  (YYMMDD) */
/*   820712  DATE WRITTEN */
/*   890126  Revised to meet SLATEC CML recommendations.  (DWL and JMS) */
/*   901019  Revisions to prologue.  (DWL and WRB) */
/*   901106  Changed all specific intrinsics to generic.  (WRB) */
/*           Corrected order of sections in prologue and added TYPE */
/*           section.  (WRB) */
/*           CALLs to XERROR changed to CALLs to XERMSG.  (WRB) */
/*   920127  Revised PURPOSE section of prologue.  (DWL) */
/* ***END PROLOGUE  DXC210 */

/*   THE CONDITIONS IMPOSED ON NLG102, MLG102, AND LG102 BY */
/* THIS SUBROUTINE ARE */

/*     (1) NLG102 .GE. 2 */

/*     (2) MLG102 .GE. 1 */

/*     (3) 2*MLG102*(MLG102 - 1) .LE. 2**NBITS - 1 */

/* THESE CONDITIONS MUST BE MET BY APPROPRIATE CODING */
/* IN SUBROUTINE DXSET. */

/* ***FIRST EXECUTABLE STATEMENT  DXC210 */
    *ierror = 0;
    if (*k == 0) {
	goto L70;
    }
    m = dxblk3_1.mlg102;
    ka = abs(*k);
    ka1 = ka / m;
    ka2 = ka % m;
    if (ka1 >= m) {
	goto L60;
    }
    nm1 = dxblk3_1.nlg102 - 1;
    np1 = dxblk3_1.nlg102 + 1;
    it = ka2 * dxblk3_1.lg102[np1 - 1];
    ic = it / m;
    id = it % m;
    *z__ = (doublereal) id;
    if (ka1 > 0) {
	goto L20;
    }
    i__1 = nm1;
    for (ii = 1; ii <= i__1; ++ii) {
	i__ = np1 - ii;
	it = ka2 * dxblk3_1.lg102[i__ - 1] + ic;
	ic = it / m;
	id = it % m;
	*z__ = *z__ / m + id;
/* L10: */
    }
    ja = ka * dxblk3_1.lg102[0] + ic;
    goto L40;
L20:
    i__1 = nm1;
    for (ii = 1; ii <= i__1; ++ii) {
	i__ = np1 - ii;
	it = ka2 * dxblk3_1.lg102[i__ - 1] + ka1 * dxblk3_1.lg102[i__] + ic;
	ic = it / m;
	id = it % m;
	*z__ = *z__ / m + id;
/* L30: */
    }
    ja = ka * dxblk3_1.lg102[0] + ka1 * dxblk3_1.lg102[1] + ic;
L40:
    *z__ /= m;
    if (*k > 0) {
	goto L50;
    }
    *j = -ja;
    d__1 = -(*z__);
    *z__ = pow_dd(&c_b9, &d__1);
    goto L80;
L50:
    *j = ja + 1;
    d__1 = *z__ - 1.;
    *z__ = pow_dd(&c_b9, &d__1);
    goto L80;
L60:
/*   THIS ERROR OCCURS IF K EXCEEDS  MLG102**2 - 1  IN MAGNITUDE. */

    xermsg_("SLATEC", "DXC210", "K too large", &c__208, &c__1, (ftnlen)6, (
	    ftnlen)6, (ftnlen)11);
    *ierror = 208;
    return 0;
L70:
    *j = 0;
    *z__ = 1.;
L80:
    return 0;
} /* dxc210_ */

