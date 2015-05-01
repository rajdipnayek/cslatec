/* dfulmt.f -- translated by f2c (version 12.02.01).
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

/* DECK DFULMT */
/* Subroutine */ int dfulmt_(integer *i__, integer *j, doublereal *aij, 
	integer *indcat, doublereal *prgopt, doublereal *dattrv, integer *
	iflag)
{
    static integer lp, key, nerr, next;
    static doublereal zero;
    static integer level;
    extern /* Subroutine */ int xermsg_(char *, char *, char *, integer *, 
	    integer *, ftnlen, ftnlen, ftnlen);

/* ***BEGIN PROLOGUE  DFULMT */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to DSPLP */
/* ***LIBRARY   SLATEC */
/* ***TYPE      DOUBLE PRECISION (FULMAT-S, DFULMT-D) */
/* ***AUTHOR  (UNKNOWN) */
/* ***DESCRIPTION */

/*     DECODES A STANDARD TWO-DIMENSIONAL FORTRAN ARRAY PASSED */
/*     IN THE ARRAY DATTRV(IA,*).  THE ROW DIMENSION IA AND THE */
/*     MATRIX DIMENSIONS MRELAS AND NVARS MUST SIMULTANEOUSLY BE */
/*     PASSED USING THE OPTION ARRAY, PRGOPT(*).  IT IS AN ERROR */
/*     IF THIS DATA IS NOT PASSED TO DFULMT( ). */
/*     EXAMPLE-- (FOR USE TOGETHER WITH DSPLP().) */
/*      EXTERNAL DUSRMT */
/*      DIMENSION DATTRV(IA,*) */
/*      PRGOPT(01)=7 */
/*      PRGOPT(02)=68 */
/*      PRGOPT(03)=1 */
/*      PRGOPT(04)=IA */
/*      PRGOPT(05)=MRELAS */
/*      PRGOPT(06)=NVARS */
/*      PRGOPT(07)=1 */
/*     CALL DSPLP(  ... DFULMT INSTEAD OF DUSRMT...) */

/* ***SEE ALSO  DSPLP */
/* ***ROUTINES CALLED  XERMSG */
/* ***REVISION HISTORY  (YYMMDD) */
/*   811215  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ) */
/*   900328  Added TYPE section.  (WRB) */
/* ***END PROLOGUE  DFULMT */
/* ***FIRST EXECUTABLE STATEMENT  DFULMT */
    /* Parameter adjustments */
    --iflag;
    --dattrv;
    --prgopt;

    /* Function Body */
    if (! (iflag[1] == 1)) {
	goto L50;
    }
/*     INITIALIZE POINTERS TO PROCESS FULL TWO-DIMENSIONAL FORTRAN */
/*     ARRAYS. */
    zero = 0.;
    lp = 1;
L10:
    next = (integer) prgopt[lp];
    if (! (next <= 1)) {
	goto L20;
    }
    nerr = 29;
    level = 1;
    xermsg_("SLATEC", "DFULMT", "IN DSPLP, ROW DIM., MRELAS, NVARS ARE MISSI"
	    "NG FROM PRGOPT.", &nerr, &level, (ftnlen)6, (ftnlen)6, (ftnlen)58)
	    ;
    iflag[1] = 3;
    goto L110;
L20:
    key = (integer) prgopt[lp + 1];
    if (! (key != 68)) {
	goto L30;
    }
    lp = next;
    goto L10;
L30:
    if (! (prgopt[lp + 2] == zero)) {
	goto L40;
    }
    lp = next;
    goto L10;
L40:
    iflag[2] = 1;
    iflag[3] = 1;
    iflag[4] = (integer) prgopt[lp + 3];
    iflag[5] = (integer) prgopt[lp + 4];
    iflag[6] = (integer) prgopt[lp + 5];
    goto L110;
L50:
    if (! (iflag[1] == 2)) {
	goto L100;
    }
L60:
    *i__ = iflag[2];
    *j = iflag[3];
    if (! (*j > iflag[6])) {
	goto L70;
    }
    iflag[1] = 3;
    goto L110;
L70:
    if (! (*i__ > iflag[5])) {
	goto L80;
    }
    iflag[2] = 1;
    iflag[3] = *j + 1;
    goto L60;
L80:
    *aij = dattrv[iflag[4] * (*j - 1) + *i__];
    iflag[2] = *i__ + 1;
    if (! (*aij == zero)) {
	goto L90;
    }
    goto L60;
L90:
    *indcat = 0;
    goto L110;
L100:
L110:
    return 0;
} /* dfulmt_ */

