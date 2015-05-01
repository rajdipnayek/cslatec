/* splpup.f -- translated by f2c (version 12.02.01).
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

static integer c__1 = 1;
static integer c__10 = 10;
static integer c__3 = 3;
static integer c__11 = 11;
static integer c__7 = 7;
static integer c__12 = 12;
static integer c__13 = 13;
static integer c__8 = 8;
static integer c__5 = 5;
static integer c__9 = 9;
static integer c__22 = 22;

/* DECK SPLPUP */
/* Subroutine */ int splpup_(S_fp usrmat, integer *mrelas, integer *nvars__, 
	real *prgopt, real *dattrv, real *bl, real *bu, integer *ind, integer 
	*info, real *amat, integer *imat, logical *sizeup, real *asmall, real 
	*abig)
{
    /* System generated locals */
    address a__1[3], a__2[7], a__3[5];
    integer i__1, i__2[3], i__3[7], i__4, i__5[5];
    char ch__1[56], ch__2[130], ch__3[54], ch__4[128], ch__5[73], ch__6[74];

    /* Local variables */
    static integer i__, j;
    static real aij, amn, amx, xval, zero;
    static char xern1[8], xern2[8], xern3[16], xern4[16];
    static integer iflag[10], index, itcnt, itmax;
    static logical first;
    static integer iplace, indcat;
    extern /* Subroutine */ int pchngs_(integer *, real *, integer *, real *, 
	    integer *, integer *), xermsg_(char *, char *, char *, integer *, 
	    integer *, ftnlen, ftnlen, ftnlen), pnnzrs_(integer *, real *, 
	    integer *, real *, integer *, integer *);

    /* Fortran I/O blocks */
    static icilist io___4 = { 0, xern1, 0, "(I8)", 8, 1 };
    static icilist io___5 = { 0, xern1, 0, "(I8)", 8, 1 };
    static icilist io___7 = { 0, xern3, 0, "(1PE15.6)", 16, 1 };
    static icilist io___9 = { 0, xern4, 0, "(1PE15.6)", 16, 1 };
    static icilist io___11 = { 0, xern1, 0, "(I8)", 8, 1 };
    static icilist io___12 = { 0, xern1, 0, "(I8)", 8, 1 };
    static icilist io___13 = { 0, xern3, 0, "(1PE15.6)", 16, 1 };
    static icilist io___14 = { 0, xern4, 0, "(1PE15.6)", 16, 1 };
    static icilist io___23 = { 0, xern1, 0, "(I8)", 8, 1 };
    static icilist io___25 = { 0, xern2, 0, "(I8)", 8, 1 };
    static icilist io___29 = { 0, xern1, 0, "(I8)", 8, 1 };


/* ***BEGIN PROLOGUE  SPLPUP */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to SPLP */
/* ***LIBRARY   SLATEC */
/* ***TYPE      SINGLE PRECISION (SPLPUP-S, DPLPUP-D) */
/* ***AUTHOR  (UNKNOWN) */
/* ***DESCRIPTION */

/*     THE EDITING REQUIRED TO CONVERT THIS SUBROUTINE FROM SINGLE TO */
/*     DOUBLE PRECISION INVOLVES THE FOLLOWING CHARACTER STRING CHANGES. */

/*     USE AN EDITING COMMAND (CHANGE) /STRING-1/(TO)STRING-2/. */
/*     /REAL (12 BLANKS)/DOUBLE PRECISION/. */

/*     REVISED 810613-1130 */
/*     REVISED YYMMDD-HHMM */

/*     THIS SUBROUTINE COLLECTS INFORMATION ABOUT THE BOUNDS AND MATRIX */
/*     FROM THE USER.  IT IS PART OF THE SPLP( ) PACKAGE. */

/* ***SEE ALSO  SPLP */
/* ***ROUTINES CALLED  PCHNGS, PNNZRS, XERMSG */
/* ***REVISION HISTORY  (YYMMDD) */
/*   811215  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890605  Corrected references to XERRWV.  (WRB) */
/*   890605  Removed unreferenced labels.  (WRB) */
/*   891009  Removed unreferenced variables.  (WRB) */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ) */
/*   900328  Added TYPE section.  (WRB) */
/*   900510  Convert XERRWV calls to XERMSG calls, changed do-it-yourself */
/*           DO loops to DO loops.  (RWC) */
/*   900602  Get rid of ASSIGNed GOTOs.  (RWC) */
/* ***END PROLOGUE  SPLPUP */

/* ***FIRST EXECUTABLE STATEMENT  SPLPUP */
    /* Parameter adjustments */
    --imat;
    --amat;
    --ind;
    --bu;
    --bl;
    --dattrv;
    --prgopt;

    /* Function Body */
    zero = 0.f;

/*     CHECK USER-SUPPLIED BOUNDS */

/*     CHECK THAT IND(*) VALUES ARE 1,2,3 OR 4. */
/*     ALSO CHECK CONSISTENCY OF UPPER AND LOWER BOUNDS. */

    i__1 = *nvars__;
    for (j = 1; j <= i__1; ++j) {
	if (ind[j] < 1 || ind[j] > 4) {
	    s_wsfi(&io___4);
	    do_fio(&c__1, (char *)&j, (ftnlen)sizeof(integer));
	    e_wsfi();
/* Writing concatenation */
	    i__2[0] = 32, a__1[0] = "IN SPLP, INDEPENDENT VARIABLE = ";
	    i__2[1] = 8, a__1[1] = xern1;
	    i__2[2] = 16, a__1[2] = " IS NOT DEFINED.";
	    s_cat(ch__1, a__1, i__2, &c__3, (ftnlen)56);
	    xermsg_("SLATEC", "SPLPUP", ch__1, &c__10, &c__1, (ftnlen)6, (
		    ftnlen)6, (ftnlen)56);
	    *info = -10;
	    return 0;
	}

	if (ind[j] == 3) {
	    if (bl[j] > bu[j]) {
		s_wsfi(&io___5);
		do_fio(&c__1, (char *)&j, (ftnlen)sizeof(integer));
		e_wsfi();
		s_wsfi(&io___7);
		do_fio(&c__1, (char *)&bl[j], (ftnlen)sizeof(real));
		e_wsfi();
		s_wsfi(&io___9);
		do_fio(&c__1, (char *)&bu[j], (ftnlen)sizeof(real));
		e_wsfi();
/* Writing concatenation */
		i__3[0] = 23, a__2[0] = "IN SPLP, LOWER BOUND = ";
		i__3[1] = 16, a__2[1] = xern3;
		i__3[2] = 19, a__2[2] = " AND UPPER BOUND = ";
		i__3[3] = 16, a__2[3] = xern4;
		i__3[4] = 28, a__2[4] = " FOR INDEPENDENT VARIABLE = ";
		i__3[5] = 8, a__2[5] = xern1;
		i__3[6] = 20, a__2[6] = " ARE NOT CONSISTENT.";
		s_cat(ch__2, a__2, i__3, &c__7, (ftnlen)130);
		xermsg_("SLATEC", "SPLPUP", ch__2, &c__11, &c__1, (ftnlen)6, (
			ftnlen)6, (ftnlen)130);
		return 0;
	    }
	}
/* L10: */
    }

    i__1 = *nvars__ + *mrelas;
    for (i__ = *nvars__ + 1; i__ <= i__1; ++i__) {
	if (ind[i__] < 1 || ind[i__] > 4) {
	    s_wsfi(&io___11);
	    i__4 = i__ - *nvars__;
	    do_fio(&c__1, (char *)&i__4, (ftnlen)sizeof(integer));
	    e_wsfi();
/* Writing concatenation */
	    i__2[0] = 30, a__1[0] = "IN SPLP, DEPENDENT VARIABLE = ";
	    i__2[1] = 8, a__1[1] = xern1;
	    i__2[2] = 16, a__1[2] = " IS NOT DEFINED.";
	    s_cat(ch__3, a__1, i__2, &c__3, (ftnlen)54);
	    xermsg_("SLATEC", "SPLPUP", ch__3, &c__12, &c__1, (ftnlen)6, (
		    ftnlen)6, (ftnlen)54);
	    *info = -12;
	    return 0;
	}

	if (ind[i__] == 3) {
	    if (bl[i__] > bu[i__]) {
		s_wsfi(&io___12);
		do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
		e_wsfi();
		s_wsfi(&io___13);
		do_fio(&c__1, (char *)&bl[i__], (ftnlen)sizeof(real));
		e_wsfi();
		s_wsfi(&io___14);
		do_fio(&c__1, (char *)&bu[i__], (ftnlen)sizeof(real));
		e_wsfi();
/* Writing concatenation */
		i__3[0] = 23, a__2[0] = "IN SPLP, LOWER BOUND = ";
		i__3[1] = 16, a__2[1] = xern3;
		i__3[2] = 19, a__2[2] = " AND UPPER BOUND = ";
		i__3[3] = 16, a__2[3] = xern4;
		i__3[4] = 26, a__2[4] = " FOR DEPENDANT VARIABLE = ";
		i__3[5] = 8, a__2[5] = xern1;
		i__3[6] = 20, a__2[6] = " ARE NOT CONSISTENT.";
		s_cat(ch__4, a__2, i__3, &c__7, (ftnlen)128);
		xermsg_("SLATEC", "SPLPUP", ch__4, &c__13, &c__1, (ftnlen)6, (
			ftnlen)6, (ftnlen)128);
		*info = -13;
		return 0;
	    }
	}
/* L20: */
    }

/*     GET UPDATES OR DATA FOR MATRIX FROM THE USER */

/*     GET THE ELEMENTS OF THE MATRIX FROM THE USER.  IT WILL BE STORED */
/*     BY COLUMNS USING THE SPARSE STORAGE CODES OF RJ HANSON AND */
/*     JA WISNIEWSKI. */

    iflag[0] = 1;

/*     KEEP ACCEPTING ELEMENTS UNTIL THE USER IS FINISHED GIVING THEM. */
/*     LIMIT THIS LOOP TO 2*NVARS*MRELAS ITERATIONS. */

    itmax = (*nvars__ << 1) * *mrelas + 1;
    itcnt = 0;
    first = TRUE_;

/*     CHECK ON THE ITERATION COUNT. */

L30:
    ++itcnt;
    if (itcnt > itmax) {
	xermsg_("SLATEC", "SPLPUP", "IN SPLP, MORE THAN 2*NVARS*MRELAS ITERA"
		"TIONS DEFINING OR UPDATING MATRIX DATA.", &c__7, &c__1, (
		ftnlen)6, (ftnlen)6, (ftnlen)78);
	*info = -7;
	return 0;
    }

    aij = zero;
    (*usrmat)(&i__, &j, &aij, &indcat, &prgopt[1], &dattrv[1], iflag);
    if (iflag[0] == 1) {
	iflag[0] = 2;
	goto L30;
    }

/*     CHECK TO SEE THAT THE SUBSCRIPTS I AND J ARE VALID. */

    if (i__ < 1 || i__ > *mrelas || j < 1 || j > *nvars__) {

/*        CHECK ON SIZE OF MATRIX DATA */
/*        RECORD THE LARGEST AND SMALLEST(IN MAGNITUDE) NONZERO ELEMENTS. */

	if (iflag[0] == 3) {
	    if (*sizeup && dabs(aij) != zero) {
		if (first) {
		    amx = dabs(aij);
		    amn = dabs(aij);
		    first = FALSE_;
		} else if (dabs(aij) > amx) {
		    amx = dabs(aij);
		} else if (dabs(aij) < amn) {
		    amn = dabs(aij);
		}
	    }
	    goto L40;
	}

	s_wsfi(&io___23);
	do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
	e_wsfi();
	s_wsfi(&io___25);
	do_fio(&c__1, (char *)&j, (ftnlen)sizeof(integer));
	e_wsfi();
/* Writing concatenation */
	i__5[0] = 21, a__3[0] = "IN SPLP, ROW INDEX = ";
	i__5[1] = 8, a__3[1] = xern1;
	i__5[2] = 19, a__3[2] = " OR COLUMN INDEX = ";
	i__5[3] = 8, a__3[3] = xern2;
	i__5[4] = 17, a__3[4] = " IS OUT OF RANGE.";
	s_cat(ch__5, a__3, i__5, &c__5, (ftnlen)73);
	xermsg_("SLATEC", "SPLPUP", ch__5, &c__8, &c__1, (ftnlen)6, (ftnlen)6,
		 (ftnlen)73);
	*info = -8;
	return 0;
    }

/*     IF INDCAT=0 THEN SET A(I,J)=AIJ. */
/*     IF INDCAT=1 THEN ACCUMULATE ELEMENT, A(I,J)=A(I,J)+AIJ. */

    if (indcat == 0) {
	pchngs_(&i__, &aij, &iplace, &amat[1], &imat[1], &j);
    } else if (indcat == 1) {
	index = -(i__ - 1);
	pnnzrs_(&index, &xval, &iplace, &amat[1], &imat[1], &j);
	if (index == i__) {
	    aij += xval;
	}
	pchngs_(&i__, &aij, &iplace, &amat[1], &imat[1], &j);
    } else {
	s_wsfi(&io___29);
	do_fio(&c__1, (char *)&indcat, (ftnlen)sizeof(integer));
	e_wsfi();
/* Writing concatenation */
	i__2[0] = 27, a__1[0] = "IN SPLP, INDICATION FLAG = ";
	i__2[1] = 8, a__1[1] = xern1;
	i__2[2] = 39, a__1[2] = " FOR MATRIX DATA MUST BE EITHER 0 OR 1.";
	s_cat(ch__6, a__1, i__2, &c__3, (ftnlen)74);
	xermsg_("SLATEC", "SPLPUP", ch__6, &c__9, &c__1, (ftnlen)6, (ftnlen)6,
		 (ftnlen)74);
	*info = -9;
	return 0;
    }

/*     CHECK ON SIZE OF MATRIX DATA */
/*     RECORD THE LARGEST AND SMALLEST(IN MAGNITUDE) NONZERO ELEMENTS. */

    if (*sizeup && dabs(aij) != zero) {
	if (first) {
	    amx = dabs(aij);
	    amn = dabs(aij);
	    first = FALSE_;
	} else if (dabs(aij) > amx) {
	    amx = dabs(aij);
	} else if (dabs(aij) < amn) {
	    amn = dabs(aij);
	}
    }
    if (iflag[0] != 3) {
	goto L30;
    }

L40:
    if (*sizeup && ! first) {
	if (amn < *asmall || amx > *abig) {
	    xermsg_("SLATEC", "SPLPUP", "IN SPLP, A MATRIX ELEMENT'S SIZE IS"
		    " OUT OF THE SPECIFIED RANGE.", &c__22, &c__1, (ftnlen)6, (
		    ftnlen)6, (ftnlen)63);
	    *info = -22;
	    return 0;
	}
    }
    return 0;
} /* splpup_ */

