/* r1mach.f -- translated by f2c (version 12.02.01).
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

static integer c__16 = 16;
static integer c__0 = 0;
static integer c__32751 = 32751;
static integer c_b11 = 16777215;
static integer c__15520 = 15520;
static integer c__15536 = 15536;
static integer c__16339 = 16339;
static integer c_b20 = 4461392;
static integer c_b21 = 10451455;
static integer c__16405 = 16405;
static integer c_b24 = 9876536;
static integer c__8195 = 8195;
static integer c_b30 = 8388608;
static integer c__1 = 1;
static integer c__24574 = 24574;
static integer c_b34 = 16777214;
static integer c__16338 = 16338;
static integer c__16383 = 16383;
static integer c_b42 = 10100890;
static integer c_b43 = 8715216;
static integer c__9 = 9;
static integer c__3 = 3;

/* DECK R1MACH */
doublereal r1mach_(integer *i__)
{
    /* Initialized data */

    static integer t3e[3] = { 9777664,5323660,46980 };
    static integer sc = 0;

    /* Format strings */
    static char fmt_9010[] = "(/\002 Adjust autodoubled R1MACH by getting da"
	    "ta\002/\002 appropriate for your machine from D1MACH.\002)";
    static char fmt_9020[] = "(/\002 Adjust R1MACH by uncommenting data stat"
	    "ements\002/\002 appropriate for your machine.\002)";

    /* System generated locals */
    real ret_val;
    static real equiv_4[6];

    /* Local variables */
    static integer j, k, l;
#define log10 ((integer *)equiv_4 + 4)
#define large ((integer *)equiv_4 + 1)
#define rmach (equiv_4)
#define small ((integer *)equiv_4)
#define diver ((integer *)equiv_4 + 3)
#define right ((integer *)equiv_4 + 2)
    extern /* Subroutine */ int i1mcra_(integer *, integer *, integer *, 
	    integer *, integer *);

    /* Fortran I/O blocks */
    static cilist io___9 = { 0, 6, 0, fmt_9010, 0 };
    static cilist io___13 = { 0, 6, 0, fmt_9020, 0 };
    static cilist io___14 = { 0, 6, 0, 0, 0 };


/* ***BEGIN PROLOGUE  R1MACH */
/* ***PURPOSE  Return floating point machine dependent constants. */
/* ***LIBRARY   SLATEC */
/* ***CATEGORY  R1 */
/* ***TYPE      SINGLE PRECISION (R1MACH-S, D1MACH-D) */
/* ***KEYWORDS  MACHINE CONSTANTS */
/* ***AUTHOR  Fox, P. A., (Bell Labs) */
/*           Hall, A. D., (Bell Labs) */
/*           Schryer, N. L., (Bell Labs) */
/* ***DESCRIPTION */

/*   R1MACH can be used to obtain machine-dependent parameters for the */
/*   local machine environment.  It is a function subprogram with one */
/*   (input) argument, and can be referenced as follows: */

/*        A = R1MACH(I) */

/*   where I=1,...,5.  The (output) value of A above is determined by */
/*   the (input) value of I.  The results for various values of I are */
/*   discussed below. */

/*   R1MACH(1) = B**(EMIN-1), the smallest positive magnitude. */
/*   R1MACH(2) = B**EMAX*(1 - B**(-T)), the largest magnitude. */
/*   R1MACH(3) = B**(-T), the smallest relative spacing. */
/*   R1MACH(4) = B**(1-T), the largest relative spacing. */
/*   R1MACH(5) = LOG10(B) */

/*   Assume single precision numbers are represented in the T-digit, */
/*   base-B form */

/*              sign (B**E)*( (X(1)/B) + ... + (X(T)/B**T) ) */

/*   where 0 .LE. X(I) .LT. B for I=1,...,T, 0 .LT. X(1), and */
/*   EMIN .LE. E .LE. EMAX. */

/*   The values of B, T, EMIN and EMAX are provided in I1MACH as */
/*   follows: */
/*   I1MACH(10) = B, the base. */
/*   I1MACH(11) = T, the number of base-B digits. */
/*   I1MACH(12) = EMIN, the smallest exponent E. */
/*   I1MACH(13) = EMAX, the largest exponent E. */

/*   To alter this function for a particular environment, the desired */
/*   set of DATA statements should be activated by removing the C from */
/*   column 1.  Also, the values of R1MACH(1) - R1MACH(4) should be */
/*   checked for consistency with the local operating system. */

/* ***REFERENCES  P. A. Fox, A. D. Hall and N. L. Schryer, Framework for */
/*                 a portable library, ACM Transactions on Mathematical */
/*                 Software 4, 2 (June 1978), pp. 177-188. */
/* ***ROUTINES CALLED  XERMSG */
/* ***REVISION HISTORY  (YYMMDD) */
/*   790101  DATE WRITTEN */
/*   890213  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ) */
/*   900618  Added DEC RISC constants.  (WRB) */
/*   900723  Added IBM RS 6000 constants.  (WRB) */
/*   910710  Added HP 730 constants.  (SMR) */
/*   911114  Added Convex IEEE constants.  (WRB) */
/*   920121  Added SUN -r8 compiler option constants.  (WRB) */
/*   920229  Added Touchstone Delta i860 constants.  (WRB) */
/*   920501  Reformatted the REFERENCES section.  (WRB) */
/*   920625  Added CONVEX -p8 and -pd8 compiler option constants. */
/*           (BKS, WRB) */
/*   930201  Added DEC Alpha and SGI constants.  (RWC and WRB) */
/*   150501  Replaced with alternative implementation from BLAS. */
/* ***END PROLOGUE  R1MACH */

/*     needs to be (2) for AUTODOUBLE, HARRIS SLASH 6, ... */
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
/*      DATA RMACH(1) / O402400000000 / */
/*      DATA RMACH(2) / O376777777777 / */
/*      DATA RMACH(3) / O714400000000 / */
/*      DATA RMACH(4) / O716400000000 / */
/*      DATA RMACH(5) / O776464202324 /, SC/987/ */

/*     MACHINE CONSTANTS FOR PDP-11 FORTRANS SUPPORTING */
/*     32-BIT INTEGERS (EXPRESSED IN INTEGER AND OCTAL). */
/*      DATA SMALL(1) /    8388608 / */
/*      DATA LARGE(1) / 2147483647 / */
/*      DATA RIGHT(1) /  880803840 / */
/*      DATA DIVER(1) /  889192448 / */
/*      DATA LOG10(1) / 1067065499 /, SC/987/ */
/*      DATA RMACH(1) / O00040000000 / */
/*      DATA RMACH(2) / O17777777777 / */
/*      DATA RMACH(3) / O06440000000 / */
/*      DATA RMACH(4) / O06500000000 / */
/*      DATA RMACH(5) / O07746420233 /, SC/987/ */

/*     MACHINE CONSTANTS FOR THE UNIVAC 1100 SERIES. */
/*      DATA RMACH(1) / O000400000000 / */
/*      DATA RMACH(2) / O377777777777 / */
/*      DATA RMACH(3) / O146400000000 / */
/*      DATA RMACH(4) / O147400000000 / */
/*      DATA RMACH(5) / O177464202324 /, SC/987/ */

    if (sc != 987) {
/*        *** CHECK FOR AUTODOUBLE *** */
	small[1] = 0;
	rmach[0] = 1e13f;
	if (small[1] != 0) {
/*           *** AUTODOUBLED *** */
	    if (small[0] == 1117925532 && small[1] == -448790528) {
/*              *** IEEE BIG ENDIAN *** */
		small[0] = 1048576;
		small[1] = 0;
		large[0] = 2146435071;
		large[1] = -1;
		right[0] = 1017118720;
		right[1] = 0;
		diver[0] = 1018167296;
		diver[1] = 0;
		log10[0] = 1070810131;
		log10[1] = 1352628735;
	    } else if (small[1] == 1117925532 && small[0] == -448790528) {
/*              *** IEEE LITTLE ENDIAN *** */
		small[1] = 1048576;
		small[0] = 0;
		large[1] = 2146435071;
		large[0] = -1;
		right[1] = 1017118720;
		right[0] = 0;
		diver[1] = 1018167296;
		diver[0] = 0;
		log10[1] = 1070810131;
		log10[0] = 1352628735;
	    } else if (small[0] == -2065213935 && small[1] == 10752) {
/*              *** VAX WITH D_FLOATING *** */
		small[0] = 128;
		small[1] = 0;
		large[0] = -32769;
		large[1] = -1;
		right[0] = 9344;
		right[1] = 0;
		diver[0] = 9472;
		diver[1] = 0;
		log10[0] = 546979738;
		log10[1] = -805796613;
	    } else if (small[0] == 1267827943 && small[1] == 704643072) {
/*              *** IBM MAINFRAME *** */
		small[0] = 1048576;
		small[1] = 0;
		large[0] = 2147483647;
		large[1] = -1;
		right[0] = 856686592;
		right[1] = 0;
		diver[0] = 873463808;
		diver[1] = 0;
		log10[0] = 1091781651;
		log10[1] = 1352628735;
	    } else {
		s_wsfe(&io___9);
		e_wsfe();
		s_stop("777", (ftnlen)3);
	    }
	} else {
	    rmach[0] = 1234567.f;
	    if (small[0] == 1234613304) {
/*              *** IEEE *** */
		small[0] = 8388608;
		large[0] = 2139095039;
		right[0] = 864026624;
		diver[0] = 872415232;
		log10[0] = 1050288283;
	    } else if (small[0] == -1271379306) {
/*              *** VAX *** */
		small[0] = 128;
		large[0] = -32769;
		right[0] = 13440;
		diver[0] = 13568;
		log10[0] = 547045274;
	    } else if (small[0] == 1175639687) {
/*              *** IBM MAINFRAME *** */
		small[0] = 1048576;
		large[0] = 2147483647;
		right[0] = 990904320;
		diver[0] = 1007681536;
		log10[0] = 1091781651;
	    } else if (small[0] == 1251390520) {
/*              *** CONVEX C-1 *** */
		small[0] = 8388608;
		large[0] = 2147483647;
		right[0] = 880803840;
		diver[0] = 889192448;
		log10[0] = 1067065499;
	    } else {
		for (l = 1; l <= 3; ++l) {
		    j = small[0] / 10000000;
		    k = small[0] - j * 10000000;
		    if (k != t3e[l - 1]) {
			goto L20;
		    }
		    small[0] = j;
/* L10: */
		}
/*              *** CRAY T3E *** */
		i1mcra_(small, &k, &c__16, &c__0, &c__0);
		i1mcra_(large, &k, &c__32751, &c_b11, &c_b11);
		i1mcra_(right, &k, &c__15520, &c__0, &c__0);
		i1mcra_(diver, &k, &c__15536, &c__0, &c__0);
		i1mcra_(log10, &k, &c__16339, &c_b20, &c_b21);
		goto L30;
L20:
		i1mcra_(&j, &k, &c__16405, &c_b24, &c__0);
		if (small[0] != j) {
		    s_wsfe(&io___13);
		    e_wsfe();
		    s_stop("777", (ftnlen)3);
		}
/*              *** CRAY 1, XMP, 2, AND 3 *** */
		i1mcra_(small, &k, &c__8195, &c_b30, &c__1);
		i1mcra_(large, &k, &c__24574, &c_b11, &c_b34);
		i1mcra_(right, &k, &c__16338, &c_b30, &c__0);
		i1mcra_(diver, &k, &c__16339, &c_b30, &c__0);
		i1mcra_(log10, &k, &c__16383, &c_b42, &c_b43);
	    }
	}
L30:
	sc = 987;
    }
/*     SANITY CHECK */
    if (rmach[3] >= 1.f) {
	s_stop("776", (ftnlen)3);
    }
    if (*i__ < 1 || *i__ > 5) {
	s_wsle(&io___14);
	do_lio(&c__9, &c__1, "R1MACH(I): I =", (ftnlen)14);
	do_lio(&c__3, &c__1, (char *)&(*i__), (ftnlen)sizeof(integer));
	do_lio(&c__9, &c__1, " is out of bounds.", (ftnlen)18);
	e_wsle();
	s_stop("", (ftnlen)0);
    }
    ret_val = rmach[*i__ - 1];
    return ret_val;
} /* r1mach_ */

#undef right
#undef diver
#undef small
#undef rmach
#undef large
#undef log10


/* Subroutine */ int i1mcra_(integer *a, integer *a1, integer *b, integer *
	c__, integer *d__)
{
/* *** SPECIAL COMPUTATION FOR CRAY MACHINES **** */
    *a1 = (*b << 24) + *c__;
    *a = (*a1 << 24) + *d__;
    return 0;
} /* i1mcra_ */

