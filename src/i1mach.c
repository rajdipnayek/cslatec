/* i1mach.f -- translated by f2c (version 12.02.01).
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

/* Table of constant values */

static integer c__32767 = 32767;
static integer c_b8 = 16777215;
static integer c__16405 = 16405;
static integer c_b12 = 9876536;
static integer c__0 = 0;
static integer c_b18 = 4194303;
static integer c__9 = 9;
static integer c__1 = 1;
static integer c__3 = 3;

/* DECK I1MACH */
integer i1mach_(integer *i__)
{
    /* Initialized data */

    static integer t3e[3] = { 9777664,5323660,46980 };
    static integer sc = 0;

    /* Format strings */
    static char fmt_9010[] = "(/\002 Adjust autodoubled I1MACH by uncommenti"
	    "ng data\002/\002 statements appropriate for your machine and set"
	    "ting\002/\002 IMACH(I) = IMACH(I+3) for I = 11, 12, and 13.\002)";
    static char fmt_9020[] = "(/\002 Adjust I1MACH by uncommenting data stat"
	    "ements\002/\002 appropriate for your machine.\002)";

    /* System generated locals */
    integer ret_val;
    static integer equiv_0[16];
    static real equiv_1[2];

    /* Local variables */
    static integer j, k, i3;
#define imach (equiv_0)
#define rmach (equiv_1)
    extern /* Subroutine */ int i1mcr1_(integer *, integer *, integer *, 
	    integer *, integer *);
#define small ((integer *)equiv_1)
#define output (equiv_0 + 3)

    /* Fortran I/O blocks */
    static cilist io___7 = { 0, 6, 0, fmt_9010, 0 };
    static cilist io___11 = { 0, 6, 0, fmt_9020, 0 };
    static cilist io___12 = { 0, 6, 0, 0, 0 };


/* ***BEGIN PROLOGUE  I1MACH */
/* ***PURPOSE  Return integer machine dependent constants. */
/* ***LIBRARY   SLATEC */
/* ***CATEGORY  R1 */
/* ***TYPE      INTEGER (I1MACH-I) */
/* ***KEYWORDS  MACHINE CONSTANTS */
/* ***AUTHOR  Fox, P. A., (Bell Labs) */
/*           Hall, A. D., (Bell Labs) */
/*           Schryer, N. L., (Bell Labs) */
/* ***DESCRIPTION */

/*   I1MACH can be used to obtain machine-dependent parameters for the */
/*   local machine environment.  It is a function subprogram with one */
/*   (input) argument and can be referenced as follows: */

/*        K = I1MACH(I) */

/*   where I=1,...,16.  The (output) value of K above is determined by */
/*   the (input) value of I.  The results for various values of I are */
/*   discussed below. */

/*   I/O unit numbers: */
/*     I1MACH( 1) = the standard input unit. */
/*     I1MACH( 2) = the standard output unit. */
/*     I1MACH( 3) = the standard punch unit. */
/*     I1MACH( 4) = the standard error message unit. */

/*   Words: */
/*     I1MACH( 5) = the number of bits per integer storage unit. */
/*     I1MACH( 6) = the number of characters per integer storage unit. */

/*   Integers: */
/*     assume integers are represented in the S-digit, base-A form */

/*                sign ( X(S-1)*A**(S-1) + ... + X(1)*A + X(0) ) */

/*                where 0 .LE. X(I) .LT. A for I=0,...,S-1. */
/*     I1MACH( 7) = A, the base. */
/*     I1MACH( 8) = S, the number of base-A digits. */
/*     I1MACH( 9) = A**S - 1, the largest magnitude. */

/*   Floating-Point Numbers: */
/*     Assume floating-point numbers are represented in the T-digit, */
/*     base-B form */
/*                sign (B**E)*( (X(1)/B) + ... + (X(T)/B**T) ) */

/*                where 0 .LE. X(I) .LT. B for I=1,...,T, */
/*                0 .LT. X(1), and EMIN .LE. E .LE. EMAX. */
/*     I1MACH(10) = B, the base. */

/*   Single-Precision: */
/*     I1MACH(11) = T, the number of base-B digits. */
/*     I1MACH(12) = EMIN, the smallest exponent E. */
/*     I1MACH(13) = EMAX, the largest exponent E. */

/*   Double-Precision: */
/*     I1MACH(14) = T, the number of base-B digits. */
/*     I1MACH(15) = EMIN, the smallest exponent E. */
/*     I1MACH(16) = EMAX, the largest exponent E. */

/*   To alter this function for a particular environment, the desired */
/*   set of DATA statements should be activated by removing the C from */
/*   column 1.  Also, the values of I1MACH(1) - I1MACH(4) should be */
/*   checked for consistency with the local operating system. */

/* ***REFERENCES  P. A. Fox, A. D. Hall and N. L. Schryer, Framework for */
/*                 a portable library, ACM Transactions on Mathematical */
/*                 Software 4, 2 (June 1978), pp. 177-188. */
/* ***ROUTINES CALLED  (NONE) */
/* ***REVISION HISTORY  (YYMMDD) */
/*   750101  DATE WRITTEN */
/*   891012  Added VAX G-floating constants.  (WRB) */
/*   891012  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900618  Added DEC RISC constants.  (WRB) */
/*   900723  Added IBM RS 6000 constants.  (WRB) */
/*   901009  Correct I1MACH(7) for IBM Mainframes. Should be 2 not 16. */
/*           (RWC) */
/*   910710  Added HP 730 constants.  (SMR) */
/*   911114  Added Convex IEEE constants.  (WRB) */
/*   920121  Added SUN -r8 compiler option constants.  (WRB) */
/*   920229  Added Touchstone Delta i860 constants.  (WRB) */
/*   920501  Reformatted the REFERENCES section.  (WRB) */
/*   920625  Added Convex -p8 and -pd8 compiler option constants. */
/*           (BKS, WRB) */
/*   930201  Added DEC Alpha and SGI constants.  (RWC and WRB) */
/*   930618  Corrected I1MACH(5) for Convex -p8 and -pd8 compiler */
/*           options.  (DWL, RWC and WRB). */
/*   150501  Replaced with alternative implementation from BLAS. */
/* ***END PROLOGUE  I1MACH */

/*  THIS VERSION ADAPTS AUTOMATICALLY TO MOST CURRENT MACHINES, */
/*  INCLUDING AUTO-DOUBLE COMPILERS. */
/*  TO COMPILE ON OLDER MACHINES, ADD A C IN COLUMN 1 */
/*  ON THE NEXT LINE */
/*  AND REMOVE THE C FROM COLUMN 1 IN ONE OF THE SECTIONS BELOW. */
/*  CONSTANTS FOR EVEN OLDER MACHINES CAN BE OBTAINED BY */
/*          mail netlib@research.bell-labs.com */
/*          send old1mach from blas */
/*  PLEASE SEND CORRECTIONS TO dmg OR ehg@bell-labs.com. */

/*     MACHINE CONSTANTS FOR THE HONEYWELL DPS 8/70 SERIES. */

/*      DATA IMACH( 1) /    5 / */
/*      DATA IMACH( 2) /    6 / */
/*      DATA IMACH( 3) /   43 / */
/*      DATA IMACH( 4) /    6 / */
/*      DATA IMACH( 5) /   36 / */
/*      DATA IMACH( 6) /    4 / */
/*      DATA IMACH( 7) /    2 / */
/*      DATA IMACH( 8) /   35 / */
/*      DATA IMACH( 9) / O377777777777 / */
/*      DATA IMACH(10) /    2 / */
/*      DATA IMACH(11) /   27 / */
/*      DATA IMACH(12) / -127 / */
/*      DATA IMACH(13) /  127 / */
/*      DATA IMACH(14) /   63 / */
/*      DATA IMACH(15) / -127 / */
/*      DATA IMACH(16) /  127 /, SC/987/ */

/*     MACHINE CONSTANTS FOR PDP-11 FORTRANS SUPPORTING */
/*     32-BIT INTEGER ARITHMETIC. */

/*      DATA IMACH( 1) /    5 / */
/*      DATA IMACH( 2) /    6 / */
/*      DATA IMACH( 3) /    7 / */
/*      DATA IMACH( 4) /    6 / */
/*      DATA IMACH( 5) /   32 / */
/*      DATA IMACH( 6) /    4 / */
/*      DATA IMACH( 7) /    2 / */
/*      DATA IMACH( 8) /   31 / */
/*      DATA IMACH( 9) / 2147483647 / */
/*      DATA IMACH(10) /    2 / */
/*      DATA IMACH(11) /   24 / */
/*      DATA IMACH(12) / -127 / */
/*      DATA IMACH(13) /  127 / */
/*      DATA IMACH(14) /   56 / */
/*      DATA IMACH(15) / -127 / */
/*      DATA IMACH(16) /  127 /, SC/987/ */

/*     MACHINE CONSTANTS FOR THE UNIVAC 1100 SERIES. */

/*     NOTE THAT THE PUNCH UNIT, I1MACH(3), HAS BEEN SET TO 7 */
/*     WHICH IS APPROPRIATE FOR THE UNIVAC-FOR SYSTEM. */
/*     IF YOU HAVE THE UNIVAC-FTN SYSTEM, SET IT TO 1. */

/*      DATA IMACH( 1) /    5 / */
/*      DATA IMACH( 2) /    6 / */
/*      DATA IMACH( 3) /    7 / */
/*      DATA IMACH( 4) /    6 / */
/*      DATA IMACH( 5) /   36 / */
/*      DATA IMACH( 6) /    6 / */
/*      DATA IMACH( 7) /    2 / */
/*      DATA IMACH( 8) /   35 / */
/*      DATA IMACH( 9) / O377777777777 / */
/*      DATA IMACH(10) /    2 / */
/*      DATA IMACH(11) /   27 / */
/*      DATA IMACH(12) / -128 / */
/*      DATA IMACH(13) /  127 / */
/*      DATA IMACH(14) /   60 / */
/*      DATA IMACH(15) /-1024 / */
/*      DATA IMACH(16) / 1023 /, SC/987/ */

    if (sc != 987) {
/*        *** CHECK FOR AUTODOUBLE *** */
	small[1] = 0;
	*rmach = 1e13f;
	if (small[1] != 0) {
/*           *** AUTODOUBLED *** */
	    if (small[0] == 1117925532 && small[1] == -448790528 || small[1] 
		    == 1117925532 && small[0] == -448790528) {
/*               *** IEEE *** */
		imach[9] = 2;
		imach[13] = 53;
		imach[14] = -1021;
		imach[15] = 1024;
	    } else if (small[0] == -2065213935 && small[1] == 10752) {
/*               *** VAX WITH D_FLOATING *** */
		imach[9] = 2;
		imach[13] = 56;
		imach[14] = -127;
		imach[15] = 127;
	    } else if (small[0] == 1267827943 && small[1] == 704643072) {
/*               *** IBM MAINFRAME *** */
		imach[9] = 16;
		imach[13] = 14;
		imach[14] = -64;
		imach[15] = 63;
	    } else {
		s_wsfe(&io___7);
		e_wsfe();
		s_stop("777", (ftnlen)3);
	    }
	    imach[10] = imach[13];
	    imach[11] = imach[14];
	    imach[12] = imach[15];
	} else {
	    *rmach = 1234567.f;
	    if (small[0] == 1234613304) {
/*               *** IEEE *** */
		imach[9] = 2;
		imach[10] = 24;
		imach[11] = -125;
		imach[12] = 128;
		imach[13] = 53;
		imach[14] = -1021;
		imach[15] = 1024;
		sc = 987;
	    } else if (small[0] == -1271379306) {
/*               *** VAX *** */
		imach[9] = 2;
		imach[10] = 24;
		imach[11] = -127;
		imach[12] = 127;
		imach[13] = 56;
		imach[14] = -127;
		imach[15] = 127;
		sc = 987;
	    } else if (small[0] == 1175639687) {
/*               *** IBM MAINFRAME *** */
		imach[9] = 16;
		imach[10] = 6;
		imach[11] = -64;
		imach[12] = 63;
		imach[13] = 14;
		imach[14] = -64;
		imach[15] = 63;
		sc = 987;
	    } else if (small[0] == 1251390520) {
/*              *** CONVEX C-1 *** */
		imach[9] = 2;
		imach[10] = 24;
		imach[11] = -128;
		imach[12] = 127;
		imach[13] = 53;
		imach[14] = -1024;
		imach[15] = 1023;
	    } else {
		for (i3 = 1; i3 <= 3; ++i3) {
		    j = small[0] / 10000000;
		    k = small[0] - j * 10000000;
		    if (k != t3e[i3 - 1]) {
			goto L20;
		    }
		    small[0] = j;
/* L10: */
		}
/*              *** CRAY T3E *** */
		imach[0] = 5;
		imach[1] = 6;
		imach[2] = 0;
		imach[3] = 0;
		imach[4] = 64;
		imach[5] = 8;
		imach[6] = 2;
		imach[7] = 63;
		i1mcr1_(&imach[8], &k, &c__32767, &c_b8, &c_b8);
		imach[9] = 2;
		imach[10] = 53;
		imach[11] = -1021;
		imach[12] = 1024;
		imach[13] = 53;
		imach[14] = -1021;
		imach[15] = 1024;
		goto L35;
L20:
		i1mcr1_(&j, &k, &c__16405, &c_b12, &c__0);
		if (small[0] != j) {
		    s_wsfe(&io___11);
		    e_wsfe();
		    s_stop("777", (ftnlen)3);
		}
/*              *** CRAY 1, XMP, 2, AND 3 *** */
		imach[0] = 5;
		imach[1] = 6;
		imach[2] = 102;
		imach[3] = 6;
		imach[4] = 46;
		imach[5] = 8;
		imach[6] = 2;
		imach[7] = 45;
		i1mcr1_(&imach[8], &k, &c__0, &c_b18, &c_b8);
		imach[9] = 2;
		imach[10] = 47;
		imach[11] = -8188;
		imach[12] = 8189;
		imach[13] = 94;
		imach[14] = -8141;
		imach[15] = 8189;
		goto L35;
	    }
	}
	imach[0] = 5;
	imach[1] = 6;
	imach[2] = 7;
	imach[3] = 6;
	imach[4] = 32;
	imach[5] = 4;
	imach[6] = 2;
	imach[7] = 31;
	imach[8] = 2147483647;
L35:
	sc = 987;
    }
    if (*i__ < 1 || *i__ > 16) {
	goto L40;
    }
    ret_val = imach[*i__ - 1];
    return ret_val;
L40:
    s_wsle(&io___12);
    do_lio(&c__9, &c__1, "I1MACH(I): I =", (ftnlen)14);
    do_lio(&c__3, &c__1, (char *)&(*i__), (ftnlen)sizeof(integer));
    do_lio(&c__9, &c__1, " is out of bounds.", (ftnlen)18);
    e_wsle();
    s_stop("", (ftnlen)0);
    return ret_val;
} /* i1mach_ */

#undef output
#undef small
#undef rmach
#undef imach


/* Subroutine */ int i1mcr1_(integer *a, integer *a1, integer *b, integer *
	c__, integer *d__)
{
/* *** SPECIAL COMPUTATION FOR OLD CRAY MACHINES **** */
    *a1 = (*b << 24) + *c__;
    *a = (*a1 << 24) + *d__;
    return 0;
} /* i1mcr1_ */

