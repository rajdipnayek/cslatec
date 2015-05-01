/* xsetua.f -- translated by f2c (version 12.02.01).
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
static integer c__2 = 2;
static logical c_true = TRUE_;
static integer c__5 = 5;

/* DECK XSETUA */
/* Subroutine */ int xsetua_(integer *iunita, integer *n)
{
    /* System generated locals */
    address a__1[2];
    integer i__1[2], i__2;
    char ch__1[37];

    /* Local variables */
    static integer i__, junk;
    static char xern1[8];
    static integer index;
    extern integer j4save_(integer *, integer *, logical *);
    extern /* Subroutine */ int xermsg_(char *, char *, char *, integer *, 
	    integer *, ftnlen, ftnlen, ftnlen);

    /* Fortran I/O blocks */
    static icilist io___2 = { 0, xern1, 0, "(I8)", 8, 1 };


/* ***BEGIN PROLOGUE  XSETUA */
/* ***PURPOSE  Set logical unit numbers (up to 5) to which error */
/*            messages are to be sent. */
/* ***LIBRARY   SLATEC (XERROR) */
/* ***CATEGORY  R3B */
/* ***TYPE      ALL (XSETUA-A) */
/* ***KEYWORDS  ERROR, XERROR */
/* ***AUTHOR  Jones, R. E., (SNLA) */
/* ***DESCRIPTION */

/*     Abstract */
/*        XSETUA may be called to declare a list of up to five */
/*        logical units, each of which is to receive a copy of */
/*        each error message processed by this package. */
/*        The purpose of XSETUA is to allow simultaneous printing */
/*        of each error message on, say, a main output file, */
/*        an interactive terminal, and other files such as graphics */
/*        communication files. */

/*     Description of Parameters */
/*      --Input-- */
/*        IUNIT - an array of up to five unit numbers. */
/*                Normally these numbers should all be different */
/*                (but duplicates are not prohibited.) */
/*        N     - the number of unit numbers provided in IUNIT */
/*                must have 1 .LE. N .LE. 5. */

/* ***REFERENCES  R. E. Jones and D. K. Kahaner, XERROR, the SLATEC */
/*                 Error-handling Package, SAND82-0800, Sandia */
/*                 Laboratories, 1982. */
/* ***ROUTINES CALLED  J4SAVE, XERMSG */
/* ***REVISION HISTORY  (YYMMDD) */
/*   790801  DATE WRITTEN */
/*   861211  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900510  Change call to XERRWV to XERMSG.  (RWC) */
/*   920501  Reformatted the REFERENCES section.  (WRB) */
/* ***END PROLOGUE  XSETUA */
/* ***FIRST EXECUTABLE STATEMENT  XSETUA */

    /* Parameter adjustments */
    --iunita;

    /* Function Body */
    if (*n < 1 || *n > 5) {
	s_wsfi(&io___2);
	do_fio(&c__1, (char *)&(*n), (ftnlen)sizeof(integer));
	e_wsfi();
/* Writing concatenation */
	i__1[0] = 29, a__1[0] = "INVALID NUMBER OF UNITS, N = ";
	i__1[1] = 8, a__1[1] = xern1;
	s_cat(ch__1, a__1, i__1, &c__2, (ftnlen)37);
	xermsg_("SLATEC", "XSETUA", ch__1, &c__1, &c__2, (ftnlen)6, (ftnlen)6,
		 (ftnlen)37);
	return 0;
    }

    i__2 = *n;
    for (i__ = 1; i__ <= i__2; ++i__) {
	index = i__ + 4;
	if (i__ == 1) {
	    index = 3;
	}
	junk = j4save_(&index, &iunita[i__], &c_true);
/* L10: */
    }
    junk = j4save_(&c__5, n, &c_true);
    return 0;
} /* xsetua_ */

