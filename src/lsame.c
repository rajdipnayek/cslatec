/* lsame.f -- translated by f2c (version 12.02.01).
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

/* DECK LSAME */
logical lsame_(char *ca, char *cb, ftnlen ca_len, ftnlen cb_len)
{
    /* Initialized data */

    static logical first = TRUE_;

    /* System generated locals */
    logical ret_val;

    /* Local variables */
    static integer ioff;

/* ***BEGIN PROLOGUE  LSAME */
/* ***SUBSIDIARY */
/* ***PURPOSE  Test two characters to determine if they are the same */
/*            letter, except for case. */
/* ***LIBRARY   SLATEC */
/* ***CATEGORY  R, N3 */
/* ***TYPE      LOGICAL (LSAME-L) */
/* ***KEYWORDS  CHARACTER COMPARISON, LEVEL 2 BLAS, LEVEL 3 BLAS */
/* ***AUTHOR  Hanson, R., (SNLA) */
/*           Du Croz, J., (NAG) */
/* ***DESCRIPTION */

/*  LSAME  tests if CA is the same letter as CB regardless of case. */
/*  CB is assumed to be an upper case letter. LSAME returns .TRUE. if */
/*  CA is either the same as CB or the equivalent lower case letter. */

/*  N.B. This version of the code is correct for both ASCII and EBCDIC */
/*       systems.  Installers must modify the routine for other */
/*       character-codes. */

/*       For CDC systems using 6-12 bit representations, the system- */
/*       specific code in comments must be activated. */

/*  Parameters */
/*  ========== */

/*  CA     - CHARACTER*1 */
/*  CB     - CHARACTER*1 */
/*           On entry, CA and CB specify characters to be compared. */
/*           Unchanged on exit. */

/* ***REFERENCES  (NONE) */
/* ***ROUTINES CALLED  (NONE) */
/* ***REVISION HISTORY  (YYMMDD) */
/*   860720  DATE WRITTEN */
/*   910606  Modified to meet SLATEC prologue standards.  Only comment */
/*           lines were modified.  (BKS) */
/*   910607  Modified to handle ASCII and EBCDIC codes.  (WRB) */
/*   930201  Tests for equality and equivalence combined.  (RWC and WRB) */
/* ***END PROLOGUE  LSAME */
/*     .. Scalar Arguments .. */
/*     .. Local Scalars .. */
/*     .. Intrinsic Functions .. */
/*     .. Save statement .. */
/*     .. Data statements .. */
/* ***FIRST EXECUTABLE STATEMENT  LSAME */
    if (first) {
	ioff = 'a' - 'A';
    }

    first = FALSE_;

/*     Test if the characters are equal or equivalent. */

    ret_val = *(unsigned char *)ca == *(unsigned char *)cb || *(unsigned char 
	    *)ca - ioff == *(unsigned char *)cb;

    return ret_val;

/*  The following comments contain code for CDC systems using 6-12 bit */
/*  representations. */

/*     .. Parameters .. */
/*     INTEGER                ICIRFX */
/*     PARAMETER            ( ICIRFX=62 ) */
/*     .. Scalar Arguments .. */
/*     CHARACTER*1            CB */
/*     .. Array Arguments .. */
/*     CHARACTER*1            CA(*) */
/*     .. Local Scalars .. */
/*     INTEGER                IVAL */
/*     .. Intrinsic Functions .. */
/*     INTRINSIC              ICHAR, CHAR */
/*     .. Executable Statements .. */
/*     INTRINSIC              ICHAR, CHAR */

/*     See if the first character in string CA equals string CB. */

/*     LSAME = CA(1) .EQ. CB .AND. CA(1) .NE. CHAR(ICIRFX) */

/*     IF (LSAME) RETURN */

/*     The characters are not identical. Now check them for equivalence. */
/*     Look for the 'escape' character, circumflex, followed by the */
/*     letter. */

/*     IVAL = ICHAR(CA(2)) */
/*     IF (IVAL.GE.ICHAR('A') .AND. IVAL.LE.ICHAR('Z')) THEN */
/*        LSAME = CA(1) .EQ. CHAR(ICIRFX) .AND. CA(2) .EQ. CB */
/*     ENDIF */

/*     RETURN */

/*     End of LSAME. */

} /* lsame_ */

