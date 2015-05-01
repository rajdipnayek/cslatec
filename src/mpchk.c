/* mpchk.f -- translated by f2c (version 12.02.01).
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
    integer b, t, m, lun, mxr, r__[30];
} mpcom_;

#define mpcom_1 mpcom_

/* Table of constant values */

static integer c__4 = 4;
static integer c__1 = 1;

/* DECK MPCHK */
/* Subroutine */ int mpchk_(integer *i__, integer *j)
{
    /* Format strings */
    static char fmt_30[] = "(\002 *** B =\002,i10,\002 ILLEGAL IN CALL TO MP"
	    "CHK,\002/\002 PERHAPS NOT SET BEFORE CALL TO AN MP ROUTINE **"
	    "*\002)";
    static char fmt_50[] = "(\002 *** T =\002,i10,\002 ILLEGAL IN CALL TO MP"
	    "CHK,\002/\002 PERHAPS NOT SET BEFORE CALL TO AN MP ROUTINE **"
	    "*\002)";
    static char fmt_70[] = "(\002 *** M .LE. T IN CALL TO MPCHK,\002/\002 PE"
	    "RHAPS NOT SET BEFORE CALL TO AN MP ROUTINE ***\002)";
    static char fmt_90[] = "(\002 *** B TOO LARGE IN CALL TO MPCHK ***\002)";
    static char fmt_110[] = "(\002 *** MXR TOO SMALL OR NOT SET TO DIM(R) BE"
	    "FORE CALL\002,\002 TO AN MP ROUTINE *** \002/\002 *** MXR SHOULD"
	    " BE AT LEAST\002,i3,\002*T +\002,i4,\002 =\002,i6,\002  ***\002"
	    "/\002 *** ACTUALLY MXR =\002,i10,\002, AND T =\002,i10,\002  **"
	    "*\002)";

    /* Local variables */
    static integer ib, mx;
    extern /* Subroutine */ int mperr_(void);
    extern integer i1mach_(integer *);

    /* Fortran I/O blocks */
    static cilist io___1 = { 0, 0, 0, fmt_30, 0 };
    static cilist io___2 = { 0, 0, 0, fmt_50, 0 };
    static cilist io___3 = { 0, 0, 0, fmt_70, 0 };
    static cilist io___5 = { 0, 0, 0, fmt_90, 0 };
    static cilist io___7 = { 0, 0, 0, fmt_110, 0 };


/* ***BEGIN PROLOGUE  MPCHK */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to DQDOTA and DQDOTI */
/* ***LIBRARY   SLATEC */
/* ***TYPE      ALL (MPCHK-A) */
/* ***AUTHOR  (UNKNOWN) */
/* ***DESCRIPTION */

/*  Checks legality of B, T, M, MXR and LUN which should be set */
/*  in COMMON. The condition on MXR (the dimension of the EP arrays) */
/*  is that  MXR .GE. (I*T + J) */

/* ***SEE ALSO  DQDOTA, DQDOTI, MPBLAS */
/* ***ROUTINES CALLED  I1MACH, MPERR */
/* ***COMMON BLOCKS    MPCOM */
/* ***REVISION HISTORY  (YYMMDD) */
/*   791001  DATE WRITTEN */
/*   ??????  Modified for use with BLAS.  Blank COMMON changed to named */
/*           COMMON.  R given dimension 12. */
/*   891009  Removed unreferenced statement label.  (WRB) */
/*   891009  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900402  Added TYPE section.  (WRB) */
/*   930124  Increased Array size in MPCON for SUN -r8.  (RWC) */
/* ***END PROLOGUE  MPCHK */
/* ***FIRST EXECUTABLE STATEMENT  MPCHK */
    mpcom_1.lun = i1mach_(&c__4);
/* NOW CHECK LEGALITY OF B, T AND M */
    if (mpcom_1.b > 1) {
	goto L40;
    }
    io___1.ciunit = mpcom_1.lun;
    s_wsfe(&io___1);
    do_fio(&c__1, (char *)&mpcom_1.b, (ftnlen)sizeof(integer));
    e_wsfe();
    mperr_();
L40:
    if (mpcom_1.t > 1) {
	goto L60;
    }
    io___2.ciunit = mpcom_1.lun;
    s_wsfe(&io___2);
    do_fio(&c__1, (char *)&mpcom_1.t, (ftnlen)sizeof(integer));
    e_wsfe();
    mperr_();
L60:
    if (mpcom_1.m > mpcom_1.t) {
	goto L80;
    }
    io___3.ciunit = mpcom_1.lun;
    s_wsfe(&io___3);
    e_wsfe();
    mperr_();
/* 8*B*B-1 SHOULD BE REPRESENTABLE, IF NOT WILL OVERFLOW */
/* AND MAY BECOME NEGATIVE, SO CHECK FOR THIS */
L80:
    ib = (mpcom_1.b << 2) * mpcom_1.b - 1;
    if (ib > 0 && (ib << 1) + 1 > 0) {
	goto L100;
    }
    io___5.ciunit = mpcom_1.lun;
    s_wsfe(&io___5);
    e_wsfe();
    mperr_();
/* CHECK THAT SPACE IN COMMON IS SUFFICIENT */
L100:
    mx = *i__ * mpcom_1.t + *j;
    if (mpcom_1.mxr >= mx) {
	return 0;
    }
/* HERE COMMON IS TOO SMALL, SO GIVE ERROR MESSAGE. */
    io___7.ciunit = mpcom_1.lun;
    s_wsfe(&io___7);
    do_fio(&c__1, (char *)&(*i__), (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&(*j), (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&mx, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&mpcom_1.mxr, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&mpcom_1.t, (ftnlen)sizeof(integer));
    e_wsfe();
    mperr_();
    return 0;
} /* mpchk_ */

