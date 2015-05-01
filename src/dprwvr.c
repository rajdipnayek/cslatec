/* dprwvr.f -- translated by f2c (version 12.02.01).
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

/* DECK DPRWVR */
/* Subroutine */ int dprwvr_(integer *key, integer *ipage, integer *lpg, 
	doublereal *sx, integer *ix)
{
    /* Initialized data */

    static doublereal zero = 0.;
    static doublereal one = 1.;

    static integer iaddr;
    static logical first;
    static integer ipagef;
    extern /* Subroutine */ int dreadp_(integer *, integer *, doublereal *, 
	    integer *, integer *), sopenm_(integer *, integer *);
    static integer istart;
    extern /* Subroutine */ int dwritp_(integer *, integer *, doublereal *, 
	    integer *, integer *);

/* ***BEGIN PROLOGUE  DPRWVR */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to DSPLP */
/* ***LIBRARY   SLATEC */
/* ***TYPE      DOUBLE PRECISION (PRWVIR-S, DPRWVR-D) */
/* ***AUTHOR  Hanson, R. J., (SNLA) */
/*           Wisniewski, J. A., (SNLA) */
/* ***DESCRIPTION */

/*     DPRWVR LIMITS THE TYPE OF STORAGE TO A SEQUENTIAL SPARSE MATRIX */
/*     STORAGE SCHEME.  THE PAGE STORAGE IS ON RANDOM ACCESS DISK. */
/*     DPRWVR IS PART OF THE SPARSE LP PACKAGE, DSPLP. */

/*     KEY       IS A FLAG WHICH INDICATES WHETHER A READ OR WRITE */
/*               OPERATION IS TO BE PERFORMED. A VALUE OF KEY=1 INDICATES */
/*               A READ. A VALUE OF KEY=2 INDICATES A WRITE. */
/*     IPAGE     IS THE PAGE OF MATRIX MN WE ARE ACCESSING. */
/*     LPG       IS THE LENGTH OF THE PAGE. */
/*   SX(*),IX(*) IS THE MATRIX DATA. */

/*     THIS SUBROUTINE IS A MODIFICATION OF THE SUBROUTINE LRWVIR, */
/*     SANDIA LABS. REPT. SAND78-0785. */
/*     MODIFICATIONS BY K.L. HIEBERT AND R.J. HANSON */

/* ***SEE ALSO  DSPLP */
/* ***ROUTINES CALLED  DREADP, DWRITP, SOPENM */
/* ***REVISION HISTORY  (YYMMDD) */
/*   811215  DATE WRITTEN */
/*   891009  Removed unreferenced variables.  (WRB) */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900328  Added TYPE section.  (WRB) */
/*   910403  Updated AUTHOR and DESCRIPTION sections.  (WRB) */
/* ***END PROLOGUE  DPRWVR */
    /* Parameter adjustments */
    --ix;
    --sx;

    /* Function Body */
/* ***FIRST EXECUTABLE STATEMENT  DPRWVR */

/*     COMPUTE STARTING ADDRESS OF PAGE. */

    ipagef = (integer) sx[3];
    istart = ix[3] + 5;

/*     OPEN RANDOM ACCESS FILE NUMBER IPAGEF, IF FIRST PAGE WRITE. */

    first = sx[4] == zero;
    if (! first) {
	goto L20002;
    }
    sopenm_(&ipagef, lpg);
    sx[4] = one;

/*     PERFORM EITHER A READ OR A WRITE. */

L20002:
    iaddr = (*ipage << 1) - 1;
    if (! (*key == 1)) {
	goto L20005;
    }
    dreadp_(&ipagef, &ix[istart], &sx[istart], lpg, &iaddr);
    goto L20006;
L20005:
    if (! (*key == 2)) {
	goto L10001;
    }
    dwritp_(&ipagef, &ix[istart], &sx[istart], lpg, &iaddr);
L10001:
L20006:
    return 0;
} /* dprwvr_ */

