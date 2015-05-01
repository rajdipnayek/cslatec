/* indxb.f -- translated by f2c (version 12.02.01).
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

/* DECK INDXB */
/* Subroutine */ int indxb_(integer *i__, integer *ir, integer *idx, integer *
	idp)
{
    /* Local variables */
    static integer id, ipl, izh;

/* ***BEGIN PROLOGUE  INDXB */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to BLKTRI */
/* ***LIBRARY   SLATEC */
/* ***TYPE      INTEGER (INDXB-I) */
/* ***AUTHOR  (UNKNOWN) */
/* ***SEE ALSO  BLKTRI */
/* ***ROUTINES CALLED  (NONE) */
/* ***COMMON BLOCKS    CBLKT */
/* ***REVISION HISTORY  (YYMMDD) */
/*   801001  DATE WRITTEN */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900402  Added TYPE section.  (WRB) */
/*   920422  Added statement so IDX would always be defined.  (WRB) */
/* ***END PROLOGUE  INDXB */

/* ***FIRST EXECUTABLE STATEMENT  INDXB */
    *idx = *i__;
    *idp = 0;
    if (*ir < 0) {
	goto L107;
    } else if (*ir == 0) {
	goto L101;
    } else {
	goto L103;
    }
L101:
    if (*i__ - cblkt_1.nm <= 0) {
	goto L102;
    } else {
	goto L107;
    }
L102:
    *idx = *i__;
    *idp = 1;
    return 0;
L103:
    izh = pow_ii(&c__2, ir);
    id = *i__ - izh - izh;
    *idx = id + id + (*ir - 1) * cblkt_1.ik + *ir + (cblkt_1.ik - *i__) / izh 
	    + 4;
    ipl = izh - 1;
    *idp = izh + izh - 1;
    if (*i__ - ipl - cblkt_1.nm <= 0) {
	goto L105;
    } else {
	goto L104;
    }
L104:
    *idp = 0;
    return 0;
L105:
    if (*i__ + ipl - cblkt_1.nm <= 0) {
	goto L107;
    } else {
	goto L106;
    }
L106:
    *idp = cblkt_1.nm + ipl - *i__ + 1;
L107:
    return 0;
} /* indxb_ */

