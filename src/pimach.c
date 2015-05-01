/* pimach.f -- translated by f2c (version 12.02.01).
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

/* DECK PIMACH */
doublereal pimach_(real *dum)
{
    /* System generated locals */
    real ret_val;

/* ***BEGIN PROLOGUE  PIMACH */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to HSTCSP, HSTSSP and HWSCSP */
/* ***LIBRARY   SLATEC */
/* ***TYPE      SINGLE PRECISION (PIMACH-S) */
/* ***AUTHOR  (UNKNOWN) */
/* ***DESCRIPTION */

/*     This subprogram supplies the value of the constant PI correct to */
/*     machine precision where */

/*     PI=3.1415926535897932384626433832795028841971693993751058209749446 */

/* ***SEE ALSO  HSTCSP, HSTSSP, HWSCSP */
/* ***ROUTINES CALLED  (NONE) */
/* ***REVISION HISTORY  (YYMMDD) */
/*   801001  DATE WRITTEN */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900402  Added TYPE section.  (WRB) */
/* ***END PROLOGUE  PIMACH */

/* ***FIRST EXECUTABLE STATEMENT  PIMACH */
    ret_val = 3.14159265358979f;
    return ret_val;
} /* pimach_ */

