/* usrmat.f -- translated by f2c (version 12.02.01).
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

/* DECK USRMAT */
/* Subroutine */ int usrmat_(integer *i__, integer *j, real *aij, integer *
	indcat, real *prgopt, real *dattrv, integer *iflag)
{
    static integer l;

/* ***BEGIN PROLOGUE  USRMAT */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to SPLP */
/* ***LIBRARY   SLATEC */
/* ***TYPE      SINGLE PRECISION (USRMAT-S, DUSRMT-D) */
/* ***AUTHOR  (UNKNOWN) */
/* ***DESCRIPTION */

/*   The user may supply this code */

/* ***SEE ALSO  SPLP */
/* ***ROUTINES CALLED  (NONE) */
/* ***REVISION HISTORY  (YYMMDD) */
/*   811215  DATE WRITTEN */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900328  Added TYPE section.  (WRB) */
/* ***END PROLOGUE  USRMAT */

/* ***FIRST EXECUTABLE STATEMENT  USRMAT */
    /* Parameter adjustments */
    --iflag;
    --dattrv;
    --prgopt;

    /* Function Body */
    if (iflag[1] == 1) {

/*     THIS IS THE INITIALIZATION STEP.  THE VALUES OF IFLAG(K),K=2,3,4, */
/*     ARE RESPECTIVELY THE COLUMN INDEX, THE ROW INDEX (OR THE NEXT COL. */
/*     INDEX), AND THE POINTER TO THE MATRIX ENTRY'S VALUE WITHIN */
/*     DATTRV(*).  ALSO CHECK (DATTRV(1)=0.) SIGNIFYING NO DATA. */
	if (dattrv[1] == 0.f) {
	    *i__ = 0;
	    *j = 0;
	    iflag[1] = 3;
	} else {
	    iflag[2] = -dattrv[1];
	    iflag[3] = dattrv[2];
	    iflag[4] = 3;
	}

	return 0;
    } else {
	*j = iflag[2];
	*i__ = iflag[3];
	l = iflag[4];
	if (*i__ == 0) {

/*     SIGNAL THAT ALL OF THE NONZERO ENTRIES HAVE BEEN DEFINED. */
	    iflag[1] = 3;
	    return 0;
	} else if (*i__ < 0) {

/*     SIGNAL THAT A SWITCH IS MADE TO A NEW COLUMN. */
	    *j = -(*i__);
	    *i__ = dattrv[l];
	    ++l;
	}

	*aij = dattrv[l];

/*     UPDATE THE INDICES AND POINTERS FOR THE NEXT ENTRY. */
	iflag[2] = *j;
	iflag[3] = dattrv[l + 1];
	iflag[4] = l + 2;

/*     INDCAT=0 DENOTES THAT ENTRIES OF THE MATRIX ARE ASSIGNED THE */
/*     VALUES FROM DATTRV(*).  NO ACCUMULATION IS PERFORMED. */
	*indcat = 0;
	return 0;
    }
    return 0;
} /* usrmat_ */

