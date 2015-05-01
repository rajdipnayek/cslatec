/* mpblas.f -- translated by f2c (version 12.02.01).
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
    integer mpb, mpt, mpm, mplun, mpmxr, mpr[30];
} mpcom_;

#define mpcom_1 mpcom_

/* Table of constant values */

static integer c__8 = 8;
static integer c__2 = 2;
static integer c__4 = 4;
static integer c__14 = 14;
static integer c__1 = 1;
static integer c__9 = 9;

/* DECK MPBLAS */
/* Subroutine */ int mpblas_(integer *i1)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    extern integer i1mach_(integer *);
    static integer mpbexp;
    extern /* Subroutine */ int xermsg_(char *, char *, char *, integer *, 
	    integer *, ftnlen, ftnlen, ftnlen);

/* ***BEGIN PROLOGUE  MPBLAS */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to DQDOTA and DQDOTI */
/* ***LIBRARY   SLATEC */
/* ***TYPE      ALL (MPBLAS-A) */
/* ***AUTHOR  (UNKNOWN) */
/* ***DESCRIPTION */

/*     This subroutine is called to set up Brent's 'mp' package */
/*     for use by the extended precision inner products from the BLAS. */

/*     In the SLATEC library we require the Extended Precision MP number */
/*     to have a mantissa twice as long as Double Precision numbers. */
/*     The calculation of MPT (and MPMXR which is the actual array size) */
/*     in this routine will give 2x (or slightly more) on the machine */
/*     that we are running on.  The INTEGER array size of 30 was chosen */
/*     to be slightly longer than the longest INTEGER array needed on */
/*     any machine that we are currently aware of. */

/* ***SEE ALSO  DQDOTA, DQDOTI */
/* ***REFERENCES  R. P. Brent, A Fortran multiple-precision arithmetic */
/*                 package, ACM Transactions on Mathematical Software 4, */
/*                 1 (March 1978), pp. 57-70. */
/*               R. P. Brent, MP, a Fortran multiple-precision arithmetic */
/*                 package, Algorithm 524, ACM Transactions on Mathema- */
/*                 tical Software 4, 1 (March 1978), pp. 71-81. */
/* ***ROUTINES CALLED  I1MACH, XERMSG */
/* ***COMMON BLOCKS    MPCOM */
/* ***REVISION HISTORY  (YYMMDD) */
/*   791001  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890531  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900402  Added TYPE section.  (WRB) */
/*   920501  Reformatted the REFERENCES section.  (WRB) */
/*   930124  Increased Array size in MPCON for SUN -r8, and calculate */
/*               size for Quad Precision for 2x DP.  (RWC) */
/* ***END PROLOGUE  MPBLAS */
/* ***FIRST EXECUTABLE STATEMENT  MPBLAS */
    *i1 = 1;

/*     For full extended precision accuracy, MPB should be as large as */
/*     possible, subject to the restrictions in Brent's paper. */

/*     Statements below are for an integer wordlength of  48, 36, 32, */
/*     24, 18, and 16.  Pick one, or generate a new one. */
/*       48     MPB = 4194304 */
/*       36     MPB =   65536 */
/*       32     MPB =   16384 */
/*       24     MPB =    1024 */
/*       18     MPB =     128 */
/*       16     MPB =      64 */

    mpbexp = i1mach_(&c__8) / 2 - 2;
    mpcom_1.mpb = pow_ii(&c__2, &mpbexp);

/*     Set up remaining parameters */
/*                  UNIT FOR ERROR MESSAGES */
    mpcom_1.mplun = i1mach_(&c__4);
/*                  NUMBER OF MP DIGITS */
    mpcom_1.mpt = ((i1mach_(&c__14) << 1) + mpbexp - 1) / mpbexp;
/*                  DIMENSION OF R */
    mpcom_1.mpmxr = mpcom_1.mpt + 4;

    if (mpcom_1.mpmxr > 30) {
	xermsg_("SLATEC", "MPBLAS", "Array space not sufficient for Quad Pre"
		"cision 2x Double Precision, Proceeding.", &c__1, &c__1, (
		ftnlen)6, (ftnlen)6, (ftnlen)78);
	mpcom_1.mpt = 26;
	mpcom_1.mpmxr = 30;
    }
/*                  EXPONENT RANGE */
/* Computing MIN */
    i__1 = 32767, i__2 = i1mach_(&c__9) / 4 - 1;
    mpcom_1.mpm = min(i__1,i__2);
    return 0;
} /* mpblas_ */

