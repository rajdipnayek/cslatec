/* indxc.f -- translated by f2c (version 12.02.01).
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
    integer npp, k;
    real eps, cnv;
    integer nm, ncmplx, ik;
} cblkt_;

#define cblkt_1 cblkt_

/* Table of constant values */

static integer c__2 = 2;

/* DECK INDXC */
/* Subroutine */ int indxc_(integer *i__, integer *ir, integer *idxc, integer 
	*nc)
{
/* ***BEGIN PROLOGUE  INDXC */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to BLKTRI */
/* ***LIBRARY   SLATEC */
/* ***TYPE      INTEGER (INDXC-I) */
/* ***AUTHOR  (UNKNOWN) */
/* ***SEE ALSO  BLKTRI */
/* ***ROUTINES CALLED  (NONE) */
/* ***COMMON BLOCKS    CBLKT */
/* ***REVISION HISTORY  (YYMMDD) */
/*   801001  DATE WRITTEN */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900402  Added TYPE section.  (WRB) */
/* ***END PROLOGUE  INDXC */
/* ***FIRST EXECUTABLE STATEMENT  INDXC */
    *nc = pow_ii(&c__2, ir);
    *idxc = *i__;
    if (*idxc + *nc - 1 - cblkt_1.nm <= 0) {
	goto L102;
    } else {
	goto L101;
    }
L101:
    *nc = 0;
L102:
    return 0;
} /* indxc_ */

