/* dxset.f -- translated by f2c (version 12.02.01).
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
    integer nbitsf;
} dxblk1_;

#define dxblk1_1 dxblk1_

struct {
    doublereal radix, radixl, rad2l, dlg10r;
    integer l, l2, kmax;
} dxblk2_;

#define dxblk2_1 dxblk2_

struct {
    integer nlg102, mlg102, lg102[21];
} dxblk3_;

#define dxblk3_1 dxblk3_

/* Table of constant values */

static integer c__10 = 10;
static integer c__14 = 14;
static integer c__15 = 15;
static integer c__16 = 16;
static integer c__8 = 8;
static integer c__201 = 201;
static integer c__1 = 1;
static integer c__202 = 202;
static integer c__203 = 203;
static integer c__2 = 2;
static integer c__204 = 204;
static integer c__205 = 205;
static integer c__206 = 206;

/* DECK DXSET */
/* Subroutine */ int dxset_(integer *irad, integer *nradpl, doublereal *dzero,
	 integer *nbits, integer *ierror)
{
    /* Initialized data */

    static integer log102[20] = { 301,29,995,663,981,195,213,738,894,724,493,
	    26,768,189,881,462,108,541,310,428 };
    static integer iflag = 0;

    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1;

    /* Local variables */
    static integer i__, j, k, ic, nb, ii, kk, it, lx, np1, lg102x, log2r, 
	    iradx;
    extern integer i1mach_(integer *);
    static integer nrdplc, lgtemp[20], iminex, imaxex;
    extern /* Subroutine */ int xermsg_(char *, char *, char *, integer *, 
	    integer *, ftnlen, ftnlen, ftnlen);
    static integer nbitsx;
    static doublereal dzerox;

/* ***BEGIN PROLOGUE  DXSET */
/* ***PURPOSE  To provide double-precision floating-point arithmetic */
/*            with an extended exponent range. */
/* ***LIBRARY   SLATEC */
/* ***CATEGORY  A3D */
/* ***TYPE      DOUBLE PRECISION (XSET-S, DXSET-D) */
/* ***KEYWORDS  EXTENDED-RANGE DOUBLE-PRECISION ARITHMETIC */
/* ***AUTHOR  Lozier, Daniel W., (National Bureau of Standards) */
/*           Smith, John M., (NBS and George Mason University) */
/* ***DESCRIPTION */

/*   SUBROUTINE  DXSET  MUST BE CALLED PRIOR TO CALLING ANY OTHER */
/* EXTENDED-RANGE SUBROUTINE. IT CALCULATES AND STORES SEVERAL */
/* MACHINE-DEPENDENT CONSTANTS IN COMMON BLOCKS. THE USER MUST */
/* SUPPLY FOUR CONSTANTS THAT PERTAIN TO HIS PARTICULAR COMPUTER. */
/* THE CONSTANTS ARE */

/*          IRAD = THE INTERNAL BASE OF DOUBLE-PRECISION */
/*                 ARITHMETIC IN THE COMPUTER. */
/*        NRADPL = THE NUMBER OF RADIX PLACES CARRIED IN */
/*                 THE DOUBLE-PRECISION REPRESENTATION. */
/*         DZERO = THE SMALLEST OF 1/DMIN, DMAX, DMAXLN WHERE */
/*                 DMIN = THE SMALLEST POSITIVE DOUBLE-PRECISION */
/*                 NUMBER OR AN UPPER BOUND TO THIS NUMBER, */
/*                 DMAX = THE LARGEST DOUBLE-PRECISION NUMBER */
/*                 OR A LOWER BOUND TO THIS NUMBER, */
/*                 DMAXLN = THE LARGEST DOUBLE-PRECISION NUMBER */
/*                 SUCH THAT LOG10(DMAXLN) CAN BE COMPUTED BY THE */
/*                 FORTRAN SYSTEM (ON MOST SYSTEMS DMAXLN = DMAX). */
/*         NBITS = THE NUMBER OF BITS (EXCLUSIVE OF SIGN) IN */
/*                 AN INTEGER COMPUTER WORD. */

/* ALTERNATIVELY, ANY OR ALL OF THE CONSTANTS CAN BE GIVEN */
/* THE VALUE 0 (0.0D0 FOR DZERO). IF A CONSTANT IS ZERO, DXSET TRIES */
/* TO ASSIGN AN APPROPRIATE VALUE BY CALLING I1MACH */
/* (SEE P.A.FOX, A.D.HALL, N.L.SCHRYER, ALGORITHM 528 FRAMEWORK */
/* FOR A PORTABLE LIBRARY, ACM TRANSACTIONS ON MATH SOFTWARE, */
/* V.4, NO.2, JUNE 1978, 177-188). */

/*   THIS IS THE SETTING-UP SUBROUTINE FOR A PACKAGE OF SUBROUTINES */
/* THAT FACILITATE THE USE OF EXTENDED-RANGE ARITHMETIC. EXTENDED-RANGE */
/* ARITHMETIC ON A PARTICULAR COMPUTER IS DEFINED ON THE SET OF NUMBERS */
/* OF THE FORM */

/*               (X,IX) = X*RADIX**IX */

/* WHERE X IS A DOUBLE-PRECISION NUMBER CALLED THE PRINCIPAL PART, */
/* IX IS AN INTEGER CALLED THE AUXILIARY INDEX, AND RADIX IS THE */
/* INTERNAL BASE OF THE DOUBLE-PRECISION ARITHMETIC.  OBVIOUSLY, */
/* EACH REAL NUMBER IS REPRESENTABLE WITHOUT ERROR BY MORE THAN ONE */
/* EXTENDED-RANGE FORM.  CONVERSIONS BETWEEN  DIFFERENT FORMS ARE */
/* ESSENTIAL IN CARRYING OUT ARITHMETIC OPERATIONS.  WITH THE CHOICE */
/* OF RADIX WE HAVE MADE, AND THE SUBROUTINES WE HAVE WRITTEN, THESE */
/* CONVERSIONS ARE PERFORMED WITHOUT ERROR (AT LEAST ON MOST COMPUTERS). */
/* (SEE SMITH, J.M., OLVER, F.W.J., AND LOZIER, D.W., EXTENDED-RANGE */
/* ARITHMETIC AND NORMALIZED LEGENDRE POLYNOMIALS, ACM TRANSACTIONS ON */
/* MATHEMATICAL SOFTWARE, MARCH 1981). */

/*   AN EXTENDED-RANGE NUMBER  (X,IX)  IS SAID TO BE IN ADJUSTED FORM IF */
/* X AND IX ARE ZERO OR */

/*           RADIX**(-L) .LE. ABS(X) .LT. RADIX**L */

/* IS SATISFIED, WHERE L IS A COMPUTER-DEPENDENT INTEGER DEFINED IN THIS */
/* SUBROUTINE. TWO EXTENDED-RANGE NUMBERS IN ADJUSTED FORM CAN BE ADDED, */
/* SUBTRACTED, MULTIPLIED OR DIVIDED (IF THE DIVISOR IS NONZERO) WITHOUT */
/* CAUSING OVERFLOW OR UNDERFLOW IN THE PRINCIPAL PART OF THE RESULT. */
/* WITH PROPER USE OF THE EXTENDED-RANGE SUBROUTINES, THE ONLY OVERFLOW */
/* THAT CAN OCCUR IS INTEGER OVERFLOW IN THE AUXILIARY INDEX. IF THIS */
/* IS DETECTED, THE SOFTWARE CALLS XERROR (A GENERAL ERROR-HANDLING */
/* FORTRAN SUBROUTINE PACKAGE). */

/*   MULTIPLICATION AND DIVISION IS PERFORMED BY SETTING */

/*                 (X,IX)*(Y,IY) = (X*Y,IX+IY) */
/* OR */
/*                 (X,IX)/(Y,IY) = (X/Y,IX-IY). */

/* PRE-ADJUSTMENT OF THE OPERANDS IS ESSENTIAL TO AVOID */
/* OVERFLOW OR  UNDERFLOW OF THE PRINCIPAL PART. SUBROUTINE */
/* DXADJ (SEE BELOW) MAY BE CALLED TO TRANSFORM ANY EXTENDED- */
/* RANGE NUMBER INTO ADJUSTED FORM. */

/*   ADDITION AND SUBTRACTION REQUIRE THE USE OF SUBROUTINE DXADD */
/* (SEE BELOW).  THE INPUT OPERANDS NEED NOT BE IN ADJUSTED FORM. */
/* HOWEVER, THE RESULT OF ADDITION OR SUBTRACTION IS RETURNED */
/* IN ADJUSTED FORM.  THUS, FOR EXAMPLE, IF (X,IX),(Y,IY), */
/* (U,IU),  AND (V,IV) ARE IN ADJUSTED FORM, THEN */

/*                 (X,IX)*(Y,IY) + (U,IU)*(V,IV) */

/* CAN BE COMPUTED AND STORED IN ADJUSTED FORM WITH NO EXPLICIT */
/* CALLS TO DXADJ. */

/*   WHEN AN EXTENDED-RANGE NUMBER IS TO BE PRINTED, IT MUST BE */
/* CONVERTED TO AN EXTENDED-RANGE FORM WITH DECIMAL RADIX.  SUBROUTINE */
/* DXCON IS PROVIDED FOR THIS PURPOSE. */

/*   THE SUBROUTINES CONTAINED IN THIS PACKAGE ARE */

/*     SUBROUTINE DXADD */
/* USAGE */
/*                  CALL DXADD(X,IX,Y,IY,Z,IZ,IERROR) */
/*                  IF (IERROR.NE.0) RETURN */
/* DESCRIPTION */
/*                  FORMS THE EXTENDED-RANGE SUM  (Z,IZ) = */
/*                  (X,IX) + (Y,IY).  (Z,IZ) IS ADJUSTED */
/*                  BEFORE RETURNING. THE INPUT OPERANDS */
/*                  NEED NOT BE IN ADJUSTED FORM, BUT THEIR */
/*                  PRINCIPAL PARTS MUST SATISFY */
/*                  RADIX**(-2L).LE.ABS(X).LE.RADIX**(2L), */
/*                  RADIX**(-2L).LE.ABS(Y).LE.RADIX**(2L). */

/*     SUBROUTINE DXADJ */
/* USAGE */
/*                  CALL DXADJ(X,IX,IERROR) */
/*                  IF (IERROR.NE.0) RETURN */
/* DESCRIPTION */
/*                  TRANSFORMS (X,IX) SO THAT */
/*                  RADIX**(-L) .LE. ABS(X) .LT. RADIX**L. */
/*                  ON MOST COMPUTERS THIS TRANSFORMATION DOES */
/*                  NOT CHANGE THE MANTISSA OF X PROVIDED RADIX IS */
/*                  THE NUMBER BASE OF DOUBLE-PRECISION ARITHMETIC. */

/*     SUBROUTINE DXC210 */
/* USAGE */
/*                  CALL DXC210(K,Z,J,IERROR) */
/*                  IF (IERROR.NE.0) RETURN */
/* DESCRIPTION */
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

/*     SUBROUTINE DXCON */
/* USAGE */
/*                  CALL DXCON(X,IX,IERROR) */
/*                  IF (IERROR.NE.0) RETURN */
/* DESCRIPTION */
/*                  CONVERTS (X,IX) = X*RADIX**IX */
/*                  TO DECIMAL FORM IN PREPARATION FOR */
/*                  PRINTING, SO THAT (X,IX) = X*10**IX */
/*                  WHERE 1/10 .LE. ABS(X) .LT. 1 */
/*                  IS RETURNED, EXCEPT THAT IF */
/*                  (ABS(X),IX) IS BETWEEN RADIX**(-2L) */
/*                  AND RADIX**(2L) THEN THE REDUCED */
/*                  FORM WITH IX = 0 IS RETURNED. */

/*     SUBROUTINE DXRED */
/* USAGE */
/*                  CALL DXRED(X,IX,IERROR) */
/*                  IF (IERROR.NE.0) RETURN */
/* DESCRIPTION */
/*                  IF */
/*                  RADIX**(-2L) .LE. (ABS(X),IX) .LE. RADIX**(2L) */
/*                  THEN DXRED TRANSFORMS (X,IX) SO THAT IX=0. */
/*                  IF (X,IX) IS OUTSIDE THE ABOVE RANGE, */
/*                  THEN DXRED TAKES NO ACTION. */
/*                  THIS SUBROUTINE IS USEFUL IF THE */
/*                  RESULTS OF EXTENDED-RANGE CALCULATIONS */
/*                  ARE TO BE USED IN SUBSEQUENT ORDINARY */
/*                  DOUBLE-PRECISION CALCULATIONS. */

/* ***REFERENCES  Smith, Olver and Lozier, Extended-Range Arithmetic and */
/*                 Normalized Legendre Polynomials, ACM Trans on Math */
/*                 Softw, v 7, n 1, March 1981, pp 93--105. */
/* ***ROUTINES CALLED  I1MACH, XERMSG */
/* ***COMMON BLOCKS    DXBLK1, DXBLK2, DXBLK3 */
/* ***REVISION HISTORY  (YYMMDD) */
/*   820712  DATE WRITTEN */
/*   881020  Revised to meet SLATEC CML recommendations.  (DWL and JMS) */
/*   901019  Revisions to prologue.  (DWL and WRB) */
/*   901106  Changed all specific intrinsics to generic.  (WRB) */
/*           Corrected order of sections in prologue and added TYPE */
/*           section.  (WRB) */
/*           CALLs to XERROR changed to CALLs to XERMSG.  (WRB) */
/*   920127  Revised PURPOSE section of prologue.  (DWL) */
/* ***END PROLOGUE  DXSET */


/*   LOG102 CONTAINS THE FIRST 60 DIGITS OF LOG10(2) FOR USE IN */
/* CONVERSION OF EXTENDED-RANGE NUMBERS TO BASE 10 . */

/* FOLLOWING CODING PREVENTS DXSET FROM BEING EXECUTED MORE THAN ONCE. */
/* THIS IS IMPORTANT BECAUSE SOME SUBROUTINES (SUCH AS DXNRMP AND */
/* DXLEGF) CALL DXSET TO MAKE SURE EXTENDED-RANGE ARITHMETIC HAS */
/* BEEN INITIALIZED. THE USER MAY WANT TO PRE-EMPT THIS CALL, FOR */
/* EXAMPLE WHEN I1MACH IS NOT AVAILABLE. SEE CODING BELOW. */
/* ***FIRST EXECUTABLE STATEMENT  DXSET */
    *ierror = 0;
    if (iflag != 0) {
	return 0;
    }
    iradx = *irad;
    nrdplc = *nradpl;
    dzerox = *dzero;
    iminex = 0;
    imaxex = 0;
    nbitsx = *nbits;
/* FOLLOWING 5 STATEMENTS SHOULD BE DELETED IF I1MACH IS */
/* NOT AVAILABLE OR NOT CONFIGURED TO RETURN THE CORRECT */
/* MACHINE-DEPENDENT VALUES. */
    if (iradx == 0) {
	iradx = i1mach_(&c__10);
    }
    if (nrdplc == 0) {
	nrdplc = i1mach_(&c__14);
    }
    if (dzerox == 0.) {
	iminex = i1mach_(&c__15);
    }
    if (dzerox == 0.) {
	imaxex = i1mach_(&c__16);
    }
    if (nbitsx == 0) {
	nbitsx = i1mach_(&c__8);
    }
    if (iradx == 2) {
	goto L10;
    }
    if (iradx == 4) {
	goto L10;
    }
    if (iradx == 8) {
	goto L10;
    }
    if (iradx == 16) {
	goto L10;
    }
    xermsg_("SLATEC", "DXSET", "IMPROPER VALUE OF IRAD", &c__201, &c__1, (
	    ftnlen)6, (ftnlen)5, (ftnlen)22);
    *ierror = 201;
    return 0;
L10:
    log2r = 0;
    if (iradx == 2) {
	log2r = 1;
    }
    if (iradx == 4) {
	log2r = 2;
    }
    if (iradx == 8) {
	log2r = 3;
    }
    if (iradx == 16) {
	log2r = 4;
    }
    dxblk1_1.nbitsf = log2r * nrdplc;
    dxblk2_1.radix = (doublereal) iradx;
    dxblk2_1.dlg10r = d_lg10(&dxblk2_1.radix);
    if (dzerox != 0.) {
	goto L14;
    }
/* Computing MIN */
    i__1 = (1 - iminex) / 2, i__2 = (imaxex - 1) / 2;
    lx = min(i__1,i__2);
    goto L16;
L14:
    lx = (integer) (d_lg10(&dzerox) * .5 / dxblk2_1.dlg10r);
/* RADIX**(2*L) SHOULD NOT OVERFLOW, BUT REDUCE L BY 1 FOR FURTHER */
/* PROTECTION. */
    --lx;
L16:
    dxblk2_1.l2 = lx << 1;
    if (lx >= 4) {
	goto L20;
    }
    xermsg_("SLATEC", "DXSET", "IMPROPER VALUE OF DZERO", &c__202, &c__1, (
	    ftnlen)6, (ftnlen)5, (ftnlen)23);
    *ierror = 202;
    return 0;
L20:
    dxblk2_1.l = lx;
    dxblk2_1.radixl = pow_di(&dxblk2_1.radix, &dxblk2_1.l);
/* Computing 2nd power */
    d__1 = dxblk2_1.radixl;
    dxblk2_1.rad2l = d__1 * d__1;
/*    IT IS NECESSARY TO RESTRICT NBITS (OR NBITSX) TO BE LESS THAN SOME */
/* UPPER LIMIT BECAUSE OF BINARY-TO-DECIMAL CONVERSION. SUCH CONVERSION */
/* IS DONE BY DXC210 AND REQUIRES A CONSTANT THAT IS STORED TO SOME FIXED */
/* PRECISION. THE STORED CONSTANT (LOG102 IN THIS ROUTINE) PROVIDES */
/* FOR CONVERSIONS ACCURATE TO THE LAST DECIMAL DIGIT WHEN THE INTEGER */
/* WORD LENGTH DOES NOT EXCEED 63. A LOWER LIMIT OF 15 BITS IS IMPOSED */
/* BECAUSE THE SOFTWARE IS DESIGNED TO RUN ON COMPUTERS WITH INTEGER WORD */
/* LENGTH OF AT LEAST 16 BITS. */
    if (15 <= nbitsx && nbitsx <= 63) {
	goto L30;
    }
    xermsg_("SLATEC", "DXSET", "IMPROPER VALUE OF NBITS", &c__203, &c__1, (
	    ftnlen)6, (ftnlen)5, (ftnlen)23);
    *ierror = 203;
    return 0;
L30:
    i__1 = nbitsx - 1;
    dxblk2_1.kmax = pow_ii(&c__2, &i__1) - dxblk2_1.l2;
    nb = (nbitsx - 1) / 2;
    dxblk3_1.mlg102 = pow_ii(&c__2, &nb);
    if (1 <= nrdplc * log2r && nrdplc * log2r <= 120) {
	goto L40;
    }
    xermsg_("SLATEC", "DXSET", "IMPROPER VALUE OF NRADPL", &c__204, &c__1, (
	    ftnlen)6, (ftnlen)5, (ftnlen)24);
    *ierror = 204;
    return 0;
L40:
    dxblk3_1.nlg102 = nrdplc * log2r / nb + 3;
    np1 = dxblk3_1.nlg102 + 1;

/*   AFTER COMPLETION OF THE FOLLOWING LOOP, IC CONTAINS */
/* THE INTEGER PART AND LGTEMP CONTAINS THE FRACTIONAL PART */
/* OF LOG10(IRADX) IN RADIX 1000. */
    ic = 0;
    for (ii = 1; ii <= 20; ++ii) {
	i__ = 21 - ii;
	it = log2r * log102[i__ - 1] + ic;
	ic = it / 1000;
	lgtemp[i__ - 1] = it % 1000;
/* L50: */
    }

/*   AFTER COMPLETION OF THE FOLLOWING LOOP, LG102 CONTAINS */
/* LOG10(IRADX) IN RADIX MLG102. THE RADIX POINT IS */
/* BETWEEN LG102(1) AND LG102(2). */
    dxblk3_1.lg102[0] = ic;
    i__1 = np1;
    for (i__ = 2; i__ <= i__1; ++i__) {
	lg102x = 0;
	i__2 = nb;
	for (j = 1; j <= i__2; ++j) {
	    ic = 0;
	    for (kk = 1; kk <= 20; ++kk) {
		k = 21 - kk;
		it = (lgtemp[k - 1] << 1) + ic;
		ic = it / 1000;
		lgtemp[k - 1] = it % 1000;
/* L60: */
	    }
	    lg102x = (lg102x << 1) + ic;
/* L70: */
	}
	dxblk3_1.lg102[i__ - 1] = lg102x;
/* L80: */
    }

/* CHECK SPECIAL CONDITIONS REQUIRED BY SUBROUTINES... */
    if (nrdplc < dxblk2_1.l) {
	goto L90;
    }
    xermsg_("SLATEC", "DXSET", "NRADPL .GE. L", &c__205, &c__1, (ftnlen)6, (
	    ftnlen)5, (ftnlen)13);
    *ierror = 205;
    return 0;
L90:
    if (dxblk2_1.l * 6 <= dxblk2_1.kmax) {
	goto L100;
    }
    xermsg_("SLATEC", "DXSET", "6*L .GT. KMAX", &c__206, &c__1, (ftnlen)6, (
	    ftnlen)5, (ftnlen)13);
    *ierror = 206;
    return 0;
L100:
    iflag = 1;
    return 0;
} /* dxset_ */

