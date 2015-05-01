/* la05ed.f -- translated by f2c (version 12.02.01).
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
    doublereal small;
    integer lp, lenl, lenu, ncp, lrow, lcol;
} la05dd_;

#define la05dd_1 la05dd_

/* DECK LA05ED */
/* Subroutine */ int la05ed_(doublereal *a, integer *irn, integer *ip, 
	integer *n, integer *iw, integer *ia, logical *reals)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer j, k, kl, kn, nz, ipi;

/* ***BEGIN PROLOGUE  LA05ED */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to DSPLP */
/* ***LIBRARY   SLATEC */
/* ***TYPE      DOUBLE PRECISION (LA05ES-S, LA05ED-D) */
/* ***AUTHOR  (UNKNOWN) */
/* ***DESCRIPTION */

/*     THIS SUBPROGRAM IS A SLIGHT MODIFICATION OF A SUBPROGRAM */
/*     FROM THE C. 1979 AERE HARWELL LIBRARY.  THE NAME OF THE */
/*     CORRESPONDING HARWELL CODE CAN BE OBTAINED BY DELETING */
/*     THE FINAL LETTER =D= IN THE NAMES USED HERE. */
/*     REVISED SEP. 13, 1979. */

/*     ROYALTIES HAVE BEEN PAID TO AERE-UK FOR USE OF THEIR CODES */
/*     IN THE PACKAGE GIVEN HERE.  ANY PRIMARY USAGE OF THE HARWELL */
/*     SUBROUTINES REQUIRES A ROYALTY AGREEMENT AND PAYMENT BETWEEN */
/*     THE USER AND AERE-UK.  ANY USAGE OF THE SANDIA WRITTEN CODES */
/*     DSPLP( ) (WHICH USES THE HARWELL SUBROUTINES) IS PERMITTED. */

/* ***SEE ALSO  DSPLP */
/* ***ROUTINES CALLED  (NONE) */
/* ***COMMON BLOCKS    LA05DD */
/* ***REVISION HISTORY  (YYMMDD) */
/*   811215  DATE WRITTEN */
/*   890831  Modified array declarations.  (WRB) */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900402  Added TYPE section.  (WRB) */
/* ***END PROLOGUE  LA05ED */
/* ***FIRST EXECUTABLE STATEMENT  LA05ED */
    /* Parameter adjustments */
    --iw;
    --ip;
    --irn;
    --a;

    /* Function Body */
    ++la05dd_1.ncp;
/*     COMPRESS FILE OF POSITIVE INTEGERS. ENTRY J STARTS AT IRN(IP(J)) */
/*  AND CONTAINS IW(J) INTEGERS,J=1,N. OTHER COMPONENTS OF IRN ARE ZERO. */
/*  LENGTH OF COMPRESSED FILE PLACED IN LROW IF REALS IS .TRUE. OR LCOL */
/*  OTHERWISE. */
/*  IF REALS IS .TRUE. ARRAY A CONTAINS A FILE ASSOCIATED WITH IRN */
/*  AND THIS IS COMPRESSED TOO. */
/*  A,IRN,IP,IW,IA ARE INPUT/OUTPUT VARIABLES. */
/*  N,REALS ARE INPUT/UNCHANGED VARIABLES. */

    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
/* STORE THE LAST ELEMENT OF ENTRY J IN IW(J) THEN OVERWRITE IT BY -J. */
	nz = iw[j];
	if (nz <= 0) {
	    goto L10;
	}
	k = ip[j] + nz - 1;
	iw[j] = irn[k];
	irn[k] = -j;
L10:
	;
    }
/* KN IS THE POSITION OF NEXT ENTRY IN COMPRESSED FILE. */
    kn = 0;
    ipi = 0;
    kl = la05dd_1.lcol;
    if (*reals) {
	kl = la05dd_1.lrow;
    }
/* LOOP THROUGH THE OLD FILE SKIPPING ZERO (DUMMY) ELEMENTS AND */
/*     MOVING GENUINE ELEMENTS FORWARD. THE ENTRY NUMBER BECOMES */
/*     KNOWN ONLY WHEN ITS END IS DETECTED BY THE PRESENCE OF A NEGATIVE */
/*     INTEGER. */
    i__1 = kl;
    for (k = 1; k <= i__1; ++k) {
	if (irn[k] == 0) {
	    goto L30;
	}
	++kn;
	if (*reals) {
	    a[kn] = a[k];
	}
	if (irn[k] >= 0) {
	    goto L20;
	}
/* END OF ENTRY. RESTORE IRN(K), SET POINTER TO START OF ENTRY AND */
/*     STORE CURRENT KN IN IPI READY FOR USE WHEN NEXT LAST ENTRY */
/*     IS DETECTED. */
	j = -irn[k];
	irn[k] = iw[j];
	ip[j] = ipi + 1;
	iw[j] = kn - ipi;
	ipi = kn;
L20:
	irn[kn] = irn[k];
L30:
	;
    }
    if (*reals) {
	la05dd_1.lrow = kn;
    }
    if (! (*reals)) {
	la05dd_1.lcol = kn;
    }
    return 0;
} /* la05ed_ */

